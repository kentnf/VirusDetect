#!/usr/bin/perl

package align;

use strict;
use warnings;
use Exporter;
use FindBin;
use lib "$FindBin::RealBin";
use Util;


# our @ISA = qw( Exporter );

# our @EXPORT_OK = qw( renameFasta bwa_remove removeRedundancy filter_SAM files_combine2 );

=head2

 align_to_reference -- align read to reference, output is sam
 
=cut
sub align_to_reference
{
	my ($align_program, $input_file, $reference, $output_file, $parameters, $temp_folder, $debug) = @_;

	# align for bwa
	if ($align_program =~ m/bwa/)
	{
		my $sai = $temp_folder."/bwa.sai";
		my $log = $temp_folder."/bwa.log";
		Util::process_cmd("$align_program index -p $reference -a bwtsw $reference 2> $log", $debug) unless -s "$reference.amb";
		Util::process_cmd("$align_program aln $parameters $reference $input_file 1> $sai 2>> $log", $debug);
		Util::process_cmd("$align_program samse -n 10000 -s $reference $sai $input_file 1> $output_file 2>> $log", $debug);
	}

	# align for bowtie2
	if ($align_program =~ m/bowtie2/) 
	{
		my $build_program = $align_program."-build";
		Util::process_cmd("$build_program $reference $reference", $debug) unless -s "$reference.1.bt2";
		Util::process_cmd("$align_program $parameters -x $reference -U $input_file -S $output_file", $debug);		
	}

	return 1;
}

=head

=cut
sub generate_unmapped_reads
{
	my ($input_SAM, $output_reads) = @_;

	my %mapped_reads;
	my ($num_unmapped_reads, $ratio) = (0, 0);

	my $in = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$output_reads) || die $!;
	while(<$in>)
	{
		chomp;
		next if $_ =~ m/^@/;
		my @a = split(/\t/, $_);

		if ( $a[1] == 4 ) {
			print $out "\@$a[0]\n$a[9]\n+\n$a[10]\n";
			$num_unmapped_reads++;
		} else {
			$mapped_reads{$a[0]} = 1;
		}
	}
	$in->close;
	$out->close;

	$ratio = ( $num_unmapped_reads / ($num_unmapped_reads + scalar(keys %mapped_reads))) * 100;
	$ratio = sprintf("%.2f", $ratio);
	$ratio = $ratio."%";
	
	return ($num_unmapped_reads, $ratio);
}

=head2

  filter_SAM -- filter sam file ,unmapped reads, 2nd hits

=cut
sub filter_SAM
{
        my $input_SAM = shift;

        my $temp_SAM = $input_SAM.".temp";
        my ($total_count, $filtered_count) = (0, 0);
        my ($query_col, $opt_col) = (0, 11);    # query and option column number for sam
        my $max_distance = 2;                   # set $max_distance for all selected hits
        my $bestEditDist = -1;                  # set best edit distance
        my @alignment = ();                     # alignment to array
        my $pre_query_name = '';                # previous query name

        my $in  = IO::File->new($input_SAM) || die $!;
        my $out = IO::File->new(">".$temp_SAM) || die $!;
        while(<$in>)
        {
                chomp;
                if ($_ =~ m/^@/) { print $out $_."\n"; next; }
                my @a = split(/\t/, $_);
                if ( $a[1] == 4 ) { $filtered_count++; }
                else {  print $out $_."\n"; }
                $total_count++;
        }
        #my ($total_count, $kept_align) = (0,0);
	my $kept_align = 0;
        $in->close;
        $out->close;
        Util::process_cmd("mv $temp_SAM $input_SAM");
        #print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as unmapped reads, only for BWA\n";



        $in  = IO::File->new($input_SAM) || die $!;
        $out = IO::File->new(">".$temp_SAM) || die $!;
        while(<$in>)
        {
                chomp;
                if ($_ =~ m/^@/) { 
			next; #print $out $_."\n"; next; 
		}
                my @a = split(/\t/, $_);

                my $query_name = $a[$query_col];

                if ($query_name ne $pre_query_name)
                {
                        # parse the pre results
                        foreach my $align (@alignment)
                        {
                                my $editDistance;
                                if ($align =~ m/\tNM:i:(\d+)/) { $editDistance = $1; }
                                else { die "Error, this alignment info do not have edit distance : $align\n"; }
                                if ($editDistance == $bestEditDist) { print $out $align."\n"; $kept_align++; }
                        }

                        # init vars;
                        @alignment = ();
                        $bestEditDist = -1;
                        $pre_query_name = $query_name;
                }

                my $distance;
                if ($_ =~ /\tNM:i:(\d+)/) { $distance = $1; }
                else { die "Error, this alignment info do not have edit distance : $_\n"; }
                next if $distance >= $max_distance;
                if ($bestEditDist == -1) { $bestEditDist = $distance; }
                if ($distance < $bestEditDist) { $bestEditDist = $distance; }
                push (@alignment, $_);

                $total_count++;
        }
        $in->close;

        # parse final query recoed
        if (scalar(@alignment) > 0)
        {
                foreach my $align (@alignment)
                {
                        my $editDistance;
                        if ($align =~ m/\tNM:i:(\d+)/) { $editDistance = $1; }
                        else { die "Error, this alignment info do not have edit distance : $align\n"; }
                        if ($editDistance == $bestEditDist) { print $out $align."\n"; $kept_align++; }
                }
        }
        $out->close;
        $filtered_count = $total_count - $kept_align;

	Util::process_cmd("mv $temp_SAM $input_SAM");	

        #print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as 2ndhits reads, only for BWA\n";
}

=head2
	pileup_filter -- filter the pileup file by coverage (sRNA aligned)
	20140712: the input and output file is same
=cut
sub pileup_filter
{
	my ($input_pileup, $virus_seq_info, $coverage, $output_pileup) = @_;

	# put plant virus sequence length to hash
	my %seq_len;
	my $fh1 = IO::File->new($virus_seq_info) || die $!;
	while(<$fh1>) {
		chomp;
		next if $_ =~ m/^#/;
		# ID Len Desc
		# AB000048 2007 Parvovirus Feline panleukopenia virus gene ...... 1 Vertebrata
		my @a = split(/\t/, $_);
		$seq_len{$a[0]} = $a[1];
	}
	$fh1->close;

	# filter the pileup file by coverage length
	my $pre_id;
	my @pileup_info;
	my $out = IO::File->new(">".$output_pileup) || die $!;
	my $in = IO::File->new($input_pileup) || die $!;
	while(<$in>)
	{
		chomp;
		# ref pos base depth match Qual?
		# AB000282 216 T 1 ^!, ~
		my @a = split(/\t/, $_);
		die "[EROR]Pileup File, Line $_\n" if scalar @a < 5;
		my ($id, $pos, $base, $depth) = ($a[0], $a[1], $a[2], $a[3]);

		if ( defined $pre_id && $id ne $pre_id)
		{

			# parse previous ref pileup info
			die "[ERROR]undef ref id $id for lengh\n" unless defined $seq_len{$id};
			my $len_cutoff = $seq_len{$pre_id} * $coverage;
			my $len = scalar(@pileup_info);
			if ($len > $len_cutoff) {
				print $out join("\n", @pileup_info),"\n";
			}

			@pileup_info = ();
			$pre_id = $id;
		}

		push(@pileup_info, $_);
		$pre_id = $id unless defined $pre_id;
	}
	$in->close;

	# parse info for last reference
	my $len_cutoff = $seq_len{$pre_id} * $coverage;
	my $len = scalar(@pileup_info);
	if ($len > $len_cutoff) {
		print $out join("\n", @pileup_info),"\n";
	}	

	$out->close;
}

=head2
pileup_to_contig -- convert pileup file to contig sequences
=cut
sub pileup_to_contig
{
	my ($input_pileup, $output_seq, $len_cutoff, $depth_cutoff, $prefix) = @_;

	my %contig;
	my $contig_num = 0;
	$prefix = 'ALIGNED';

	my $fh = IO::File->new($input_pileup) || die $!;
	my $out = IO::File->new(">".$output_seq) || die $!;

	my $line = <$fh>; chomp($line); die "[ERROR]Pileup File, 1st Line $line\n" if $line =~ m/^#/;
	my @t = split(/\t/, $line); die "[EROR]Pileup File, Line $line\n" if scalar @t < 5;
	my ($pre_id, $pre_pos, $base, $dep, $align) = ($t[0], $t[1], $t[2], $t[3], $t[4]);
        my $match_depth = $align =~ tr/.,/.,/; my $match_max = $match_depth;
        my $match_baseA = $align =~ tr/Aa/Aa/;
        my $match_baseT = $align =~ tr/Tt/Tt/;
        my $match_baseC = $align =~ tr/Cc/Cc/;
        my $match_baseG = $align =~ tr/Gg/Gg/;
        if ($match_baseA >= $match_max) { $match_max = $match_baseA; $base = 'A'; }
        if ($match_baseT >= $match_max) { $match_max = $match_baseT; $base = 'T'; }
        if ($match_baseC >= $match_max) { $match_max = $match_baseC; $base = 'C'; }
        if ($match_baseG >= $match_max) { $match_max = $match_baseG; $base = 'G'; }
	my $seq = $base;
	my $len = 1;
	my ($id, $pos, $depth, $qual);

	while(<$fh>)
	{
		chomp;
		# ref pos base depth match Qual?
		# AB000282 216 T 1 ^!, ~
		my @a = split(/\t/, $_);
		die "[EROR]Pileup File, Line $_\n" if scalar @a < 5;
		($id, $pos, $base, $depth, $align, $qual) = ($a[0], $a[1], $a[2], $a[3], $a[4], $a[5]);

		# get best base by depth and quality
		$align =~ s/\^\S//g;
		$align =~ s/\d+\S//g;
		$align =~ s/\$//g;
		$align =~ s/\+//g;
		$align =~ s/-//g;
		die "[ERROR]Align Len $_\n$align\n" if length($align) != $depth;
		die "[ERROR]Qua Len" if length($qual) != $depth;
		$align =~ tr/atcg/ATCG/;

		my @m = split(//, $align);
		my @n = split(//, $qual);
		my %binfo = ();
		my $maxb = 1;
		for(my $i=0; $i<$depth; $i++)
		{
			my $b = $m[$i];
			my $q = ord($n[$i]);
			$b = $base if ($b eq "," || $b eq ".");

			if (defined $binfo{$b}{'freq'}) {
				$binfo{$b}{'freq'}++;
				if ($binfo{$b}{'freq'} > $maxb) { $maxb = $binfo{$b}{'freq'}; }
			} else {
				$binfo{$b}{'freq'} = 1;
			}

			if (defined $binfo{$b}{'qual'}) {
				$binfo{$b}{'qual'}+=$q;
			} else {
				$binfo{$b}{'qual'} = $q;
			}
		}

		my @maxb;
		foreach my $b (sort keys %binfo) {
			if ($binfo{$b}{'freq'} == $maxb) {
				push(@maxb, $b)
			}
		}

		my $maxq = 0;
		foreach my $b (@maxb) {
			if ($binfo{$b}{'qual'} > $maxq) {
				$maxq = $binfo{$b}{'qual'};
				$base = $b;
			}	
		}

		#gene the contig sequence
		if ( $id eq $pre_id ) # same reference
		{
			if ($pos == $pre_pos + 1) # extend contig
			{
				$len++;
				$dep+=$depth;
				$seq.=$base;
				$pre_pos = $pos;
				$pre_id = $id;
			}
			else # end pre contig, start new one
			{

				if ( $len > $len_cutoff && ($dep/$len > $depth_cutoff)) {
					die "[ERROR]Already defined contig $pre_id $pre_pos $seq\n" if defined $contig{$pre_id."#".$pre_pos};
					$contig{$pre_id."#".$pre_pos} = $seq;
					$contig_num++;
					print $out ">".$prefix.$contig_num." ".$pre_id." ".$pre_pos."\n".$seq."\n";
				}

				$len = 1;
				$dep = $depth;
				$seq = $base;
				$pre_pos = $pos;
				$pre_id = $id;
			}
		}
		else
		{
			# end pre contig, start new one
			if ( $len > $len_cutoff && ($dep/$len > $depth_cutoff)) {
				die "[ERROR]Already defined contig $pre_id $pre_pos $seq\n" if defined $contig{$pre_id."#".$pre_pos};
				$contig{$pre_id."#".$pre_pos} = $seq;
				$contig_num++;
				print $out ">".$prefix.$contig_num." ".$pre_id." ".$pre_pos."\n".$seq."\n";
			}

			$len = 1;
			$dep = $depth;
			$seq = $base;
			$pre_pos = $pos;
			$pre_id = $id;
		}
	}
	$fh->close;

	# parse the last one
	if ( $len > $len_cutoff && ( $dep/$len > $depth_cutoff )) {
		die "[ERROR]Already defined contig $pre_id $pre_pos $seq\n" if defined $contig{$pre_id."#".$pre_pos};
		$contig{$pre_id."#".$pre_pos} = $seq;
		$contig_num++;
		print $out ">".$prefix.$contig_num." ".$pre_id." ".$pre_pos."\n".$seq."\n";
	}

	$out->close;
}

=head
# input reads;
my $sample = $ARGV[0];
my $output_contig = $ARGV[1];
die "No sample input\n" unless -e $sample;
#die "No sample input\n" unless -s $output_contig;	

# then range for k-mer, and coverage
my $kmer_start = 9;
my $kmer_end   = 19;
my $coverage_start = 5;
my $coverage_end   = 15;

# objective type (n50, maxLen, avgLen) for optimization
my $objective_type='maxLen';

# set folder and file path
my $WORKING_DIR = cwd();		# current folder
my $bin_dir  = ${FindBin::RealBin};	# bin folder
my $temp_dir = $WORKING_DIR."/temp";	# temp folder
mkdir $temp_dir unless -s $temp_dir;
=cut

=head2
 velvet_optimiser -- combine velvet optimiser and run velvet script
=cut
sub velvet_optimiser_combine
{
	my ($sample, $output_contig, $kmer_start, $kmer_end, $coverage_start, $coverage_end, $objective_type, $bin_dir, $temp_dir, $debug) = @_;

	my $current_folder;
   	my $statfile;
    	my $objective;						# the objective for each run
    	my $max_objective = 0;					# the maximum objective
    	my $opt_kmer = $kmer_start; 				# the optimization kmer_length when objective is max
    	my $opt_coverage = $coverage_start;			# the optimization coverage when objective is max
    	my $opt_avgLen = 0;                			# the avg Length when objective is max
		
	# optimize k-mer length using fixed coverage
	for(my $i=$kmer_start; $i<= $kmer_end; $i=$i+2) 
	{
		runVelvet($sample, $i, $coverage_start, $bin_dir, $temp_dir, $debug);
		$current_folder = $sample."_".$i."_".$coverage_start;
		$statfile = $current_folder."/contigs.fa";
		my $aa = contigStats($statfile);
		if ( defined $aa->{$objective_type} ) 
		{
			$objective = $aa->{$objective_type};
			
			# output best hash length, coverage, objective, avgLength, contigs num
			if ( $objective > $max_objective ) 
			{
				$max_objective = $objective;
				$opt_kmer = $i;
				$opt_avgLen = $aa->{avgLen};
			}
		}
		Util::process_cmd("rm $current_folder -r", $debug);		
	} 

	# optimize coverage using fixed k-mer length
	# coverage start from 7 (coverage_start + 2), because 5 has been computed, how about 6? 
	for(my $j=$coverage_start+2; $j<=$coverage_end; $j=$j+1) 
	{  
		runVelvet($sample, $opt_kmer, $j, $bin_dir, $temp_dir, $debug);
		$current_folder = $sample."_".$opt_kmer."_".$j;
		$statfile = $current_folder. "/contigs.fa";
		my $aa = contigStats($statfile);

		if ( defined $aa->{$objective_type} ) 
		{
			$objective=$aa->{$objective_type};

			# output the best hast length, coverage, objective, avglength, contigs num
			if($objective>$max_objective)
			{
				# print OUT1 "yes"."\n";# print yes of the optimization value improved
				$max_objective = $objective;
				$opt_coverage = $j;
				$opt_avgLen = $aa->{avgLen};
			}
		}
		Util::process_cmd("rm $current_folder -r");
	}       

	print "Best Parameters for $sample:\nKmer: $opt_kmer\nCoverage: $opt_coverage\nMax $objective_type : $max_objective\nBest avgLen: $opt_avgLen\n" if $debug;

	# run velvet again using best parameters
	runVelvet($sample, $opt_kmer, $opt_coverage, $bin_dir, $temp_dir, $debug);
	$current_folder = $sample."_".$opt_kmer."_".$opt_coverage;
	Util::process_cmd("mv $current_folder/contigs.fa $output_contig", $debug);
}

=head2
 runVelvet -- run Velvet with differnt parameters
=cut
sub runVelvet 
{
	my ($sample, $kmer_length, $cov_cutoff, $bin_dir, $temp_dir, $debug) = @_;
	my $outputDir = $sample."_".$kmer_length."_".$cov_cutoff;
	my $file_type = Util::detect_FileType($sample);
	Util::process_cmd($bin_dir."/velveth $outputDir $kmer_length -$file_type $sample >> $temp_dir/velvet.log", $debug);
	Util::process_cmd($bin_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 >> $temp_dir/velvet.log", $debug);	
	return 1;
}

=head2
 contigStats -- get contig stats info
=cut
sub contigStats 
{
	
	my $file = shift;
	#my $minsize = shift;
	
	#print "In contigStats with $file, $minsize\n" if $interested;
	
	my $numseq=0;
	my $avglen=0;
	my $minlen=1E9;
	my $maxlen=0;
	my @len;
	my $toosmall=0;
	my $nn=0;
	
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	while(my $seq = $in->next_seq()){
		my $L = $seq->length;
		#if($L < $minsize){
			#$toosmall ++;
			#next;
		#}
		#count Ns
		my $s = $seq->seq;
		my $n = $s =~ s/N/N/gi;
		$n ||= 0;
		$nn += $n;
		#count seqs and other stats
		$numseq ++;
		$avglen += $L;
		$maxlen = $L if $L > $maxlen;
		$minlen = $L if $L < $minlen;
		push @len, $L;
	}
	@len = sort { $a <=> $b } @len;
	my $cum = 0;
	my $n50 = 0;
	for my $i (0 .. $#len){
		$cum += $len[$i];
		if($cum >= $avglen/2) {
			$n50 = $len[$i];
			last;
		}
	}
	
	my %out;
	if($numseq > 0)
	{
		$out{numSeqs} = $numseq;
		$out{numBases} = $avglen;
		$out{numOK} = ($avglen - $nn);
		$out{numNs} = $nn;
		$out{minLen} = $minlen;
		$out{avgLen} = $avglen/$numseq;
		$out{maxLen} = $maxlen;
		$out{n50} = $n50;
		#$out{minsize} = $minsize;
		$out{numTooSmall} = $toosmall;
	}
	else 
	{
		$out{$numseq} = 0;
	}
	
	#print "Leaving contigstats!\n" if $interested;
	return (\%out);
}

=head2
 
 remove_redundancy -- remove redundancy contig 
=cut
sub remove_redundancy
{
	my ($contig_file, $read_file, $parameters, $max_end_clip, $min_overlap, $bin_dir, $temp_dir, $debug) = @_;

	my $min_identity = '';	# pass this value to other subroutine for furture functions

	# step 1. remove low complexity sequence using dust, 2 remov XN reads
	dust_parse($contig_file, $bin_dir, $debug);
	trim_XNseq($contig_file, 0.8, 40);

	# step 2. get number of contigs 
	my ($before_contig_num, $after_contig_num);
	$before_contig_num = count_seq($contig_file);	# get the seq number before remove redundancy
	$after_contig_num  = $before_contig_num - 1;	# this is seq number after remove redundancy	
	
	# if the new contig number != old contig number, continue remove redundancy
	while( $after_contig_num != $before_contig_num )
	{ 
		remove_redundancy_main($contig_file, $parameters, $max_end_clip, $min_overlap, $min_identity, $bin_dir);
		$before_contig_num = $after_contig_num; 	# renew contig_num1
		$after_contig_num  = count_seq($contig_file); 	# get seq number of new contig
	}

	if ($after_contig_num > 1) {
		Util::print_user_submessage("$after_contig_num of uniq contigs were generated");
	} elsif ($after_contig_num == 1) {
		Util::print_user_submessage("$after_contig_num of uniq contig was generated");
	} elsif ($after_contig_num == 0) {
		Util::print_user_submessage("None of uniq contig was generated");
	}

	# finish remove redundancy, next for base correction
	my $cpu_num = 20;
	if ($parameters =~ m/-a (\d+)/) { $cpu_num = $1; }

	base_correction($read_file, $contig_file, $cpu_num, $bin_dir, $temp_dir, $debug);
}

=head2 
 base_correction -- correct contig base using reads aligned to contigs
 # aligment -> sam -> bam -> sorted bam -> pileup -> Consensus sequence
=cut

sub base_correction
{
	my ($read_file, $contig_file, $cpu_num, $bin_dir, $temp_dir, $debug) = @_;

	my $file_type = Util::detect_FileType($read_file);
	my $format = '';
	$format = '-f' if $file_type eq 'fasta';
	die "[ERR]Undef file type for $read_file\n" if ($file_type ne "fasta" && $file_type ne "fastq");
	#aligment -> sam -> bam -> sorted bam -> pileup
        Util::process_cmd("$bin_dir/bowtie-build --quiet -f $contig_file $contig_file", $debug);
        Util::process_cmd("$bin_dir/samtools faidx $contig_file 2> $temp_dir/samtools.log", $debug);
        Util::process_cmd("$bin_dir/bowtie --quiet -v 1 -p $cpu_num $format -S -a --best $contig_file $read_file $read_file.sam", $debug);
        Util::process_cmd("$bin_dir/samtools view -bS $read_file.sam > $read_file.bam 2> $temp_dir/samtools.log");
        Util::process_cmd("$bin_dir/samtools sort $read_file.bam $read_file.sorted 2> $temp_dir/samtools.log");
        Util::process_cmd("$bin_dir/samtools mpileup -f $contig_file $read_file.sorted.bam > $read_file.pileup 2> $temp_dir/samtools.log");

        my $file_size = -s "$read_file.pileup";
        if( $file_size == 0 ){                     # if file size = 0, create blank file, and exit the loop
                Util::process_cmd("touch $contig_file.blank");
                next;
        }

	pileup_to_contig("$read_file.pileup", $contig_file, 40, 1, 'UNIQUE');
        #Util::process_cmd("java -cp $BIN_DIR extractConsensus $sample 1 40 $i");
        #renameFasta("$sample.contigs$i.fa", "$sample.$input_suffix", $contig_prefix);

	# remove temp files
        system("rm $read_file.sam");
        system("rm $read_file.bam");
        system("rm $read_file.sorted.bam");
        system("rm $read_file.pileup"); 	# must delete this file for next cycle remove redundancy
        system("rm $contig_file.fai");
        system("rm $temp_dir/*.ebwt");
}
=cut

=head2
 count_seq
=cut
sub count_seq
{
	my $input = shift;
	my $seq_num = 0;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input);
	while(my $input = $in->next_seq)
	{
		$seq_num++;
	}
	return $seq_num;
}

=head2
 dust_parse -- using dust to remove low complexity contigs 
=cut
sub dust_parse
{
	my ($input_contig, $bin_dir, $debug) = @_;
	Util::process_cmd("$bin_dir/dust $input_contig 1> $input_contig.masked 2> $input_contig.dust.log", $debug);
	Util::process_cmd("mv $input_contig.masked $input_contig", $debug);
}

=head
 trim_XNseq -- trim X N base in contigs
=cut
sub trim_XNseq
{
	my ($input_contigs, $max_Xratio, $min_BaseNum) = @_;

	$max_Xratio = 0.8 unless defined $max_Xratio; 
	$min_BaseNum = 40 unless defined $min_BaseNum; 

	my $output_seq = '';
	my ($seq_length, $seq_id, $sequence);
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_contigs);
	while(my $inseq = $in->next_seq)
	{
		$seq_length= $inseq->length;
		$seq_id    = $inseq->id;
		$sequence  = $inseq->seq;
		my $xn_num = $sequence =~ tr/XxNn/XxNn/;
		
		if (($xn_num/$seq_length) >= $max_Xratio || ($seq_length - $xn_num) < $min_BaseNum )
		{
			#$bad_seq_num++;
		}
		else
		{
			$output_seq.=">$seq_id\n$sequence\n";
		}
	}

	my $fh = IO::File->new(">".$input_contigs) || die $!;
	print $fh $output_seq;
	$fh->close;
}

=head2

=cut			 

sub remove_redundancy_main
{
	my ($input_contig, $parameters, $max_end_clip, $min_overlap, $min_identity, $bin_dir) = @_;

	# step 1. sort seq by length
	# create array for store sequence info
	# [    Seq1     ,      Seq2     ,      Seq3    ]
	# [ID, Len, Seq], [ID, Len, Seq], [ID, Len, Seq]
	my @all_data; 
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_contig);
	while(my $inseq = $in->next_seq)
	{
		push(@all_data, [$inseq->id, $inseq->length, $inseq->seq]);
	}
	@all_data = sort { -1*($a->[1] <=> $b->[1]) } @all_data;

	# step 2. remove redundancy
	# put un-redundancy sequence to hash : inset
	# key: seqID
	# value: sequence
	#
	# put redundancy sequnece to vars : restset
	# fasta format
	# >seqID \n sequence \n
	# 
	# contig_conunt : inset set number
	my %inset;
	my $restset = ''; 
	my $contig_count=1;

	for my $tr (@all_data) 
	{
		if (scalar(keys %inset)  == 0) 
		{
			# put the longest sequence to hash
			$inset{$tr->[0]} = $tr->[2]; 
		}
		else
		{
			my $query = ">$tr->[0]\n$tr->[2]\n";

			# gene redundancy stat between query and %inset sequences  
			# return_string is tab delimit file
			# 1. seqID
			# 2. r or n, then combine the sequence according to r and n
			my $return_string = ifRedundant(\%inset, \$query, $input_contig, $parameters, $max_end_clip, $min_overlap, $min_identity, $bin_dir);
			my @return_col = split(/\t/, $return_string);

			if ($return_col[1] eq "r") 
			{
				# remove redundancy sequence
				$restset.=$query;
			}
			elsif ( $return_col[1] eq "n")
			{
				# add un-redundancy sequence to hash
				$contig_count++;
				$inset{$tr->[0]} = $tr->[2]; 
			}
			else
			{
				# partily redundancy , combine then replace
				# the format return_string is different
				# 1. hit_name:query_name 
				my @names = split(/\:/, $return_col[0]);#µÚÒ»ÁÐÊÇ(hit_name:query_name)		
				$inset{$names[0]} = $return_col[1];; #Ô­À´hit_name¶ÔÓ¦µÄ¼ÇÂ¼±»ÐÂÐòÁÐ¸²¸Ç
				$restset .= $query;
			}
		}
	}

	# step 3 output contigs after remove redundancy
	my $out1 = IO::File->new(">".$input_contig) || die $!; 
	while (my ($k,$v) = each %inset) { print $out1 ">$k\n$v\n";  }
	$out1->close; 

	# print "[DEBUG]\t".$input."\t".$contig_count."\n";
	
	return 1;
}

=head2

 ifRedundant: check if the query sequence is redundancy after compared with inset sequences 

=cut
sub ifRedundant 
{
	my ($inset, $query, $input_contig, $blast_parameters, $max_end_clip, $min_overlap, $min_identity, $bin_dir) = @_; 
	
	# save query and hit seqeunces to files
	my $query_seq_file = $input_contig."_query";
	my $hit_seq_file   = $input_contig."_tem";
	my $blast_output   = $input_contig."_tem.paired";

	my $fh1 = IO::File->new(">".$query_seq_file) || die $!;
	print $fh1 $$query;
	$fh1->close;
	
	my $fh2 = IO::File->new(">".$hit_seq_file) || die $!;
	while (my ($k,$v) = each %$inset) { print $fh2 ">$k\n$v\n"; }
	$fh2->close;
	
    	# perform blast. 
	# using process_cmd() could debug ouput result
	system("$bin_dir/formatdb -i $hit_seq_file -p F");
	my $blast_program = $bin_dir."/megablast";

	my $blast_param = "-i $query_seq_file -d $hit_seq_file -o $blast_output $blast_parameters";
	system($blast_program." ".$blast_param) && die "Error at blast command: $blast_param\n";
	
	# get redundancy info from blast result
	my $result = findRedundancy($inset, $query, $blast_output, $max_end_clip, $min_overlap, $min_identity);

	# for debug
	# if($$query =~ />NOVEL1\n/){
	# 	print STDERR $result."good\n";
	#	die "this is >NOVEL1";
	# }
	
	unlink ($query_seq_file, $hit_seq_file, $blast_output, "$hit_seq_file.nhr", "$hit_seq_file.nin", "$hit_seq_file.nsq");
	return $result;
}

=head2

 findRedundancy : get all HSP for one query against multiply hits

=cut
sub findRedundancy
{
	my ($inset, $query, $blast_output, $max_end_clip, $min_overlap, $min_identity) = @_;

	my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end);

	my %hsp = ();  # key: query, value: hsp for one query
	my $hsp_count=0; 
	my $query_sequenc; 
	my $subject_sequenc; 	
	my $is_hsp = 1;

	my $bfh = IO::File->new($blast_output) || "Can not open blast result file: $blast_output $!\n";
	while(<$bfh>)
	{
		chomp;
		if (/Query=\s(\S+)/ || eof) # new Query or End of file, output result for previous hsp
		{
			#########################################
			#  previous hsp info to hash		#
			#########################################
			if ($hsp_count > 0)#Èç¹û²»ÊÇµÚÒ»´Î³
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#########################################
			# parse all hsp info for prev query-hit	#
			#########################################
			if (scalar(keys(%hsp)) > 0)
			{
				foreach my $hsp_num (sort {$a<=>$b} keys %hsp)
				{
					my @one_hit = split(/\t/, $hsp{$hsp_num});
					unless(scalar(@one_hit) == 13) { die "Error! the parsed blast may lost some info:\n $hsp{$hsp_num} \n"; }
					
					my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand,
					    $query_start, $query_end, $hit_start, $hit_end) =
					   ($one_hit[0], $one_hit[1], $one_hit[2], $one_hit[3], $one_hit[4], $one_hit[5], $one_hit[6], $one_hit[7],
					    $one_hit[8], $one_hit[9], $one_hit[10],$one_hit[11],$one_hit[12]);

					my $query_to_end = $query_length - $query_end;
					my $hit_to_end = $hit_length - $hit_end;
					my $hit_to_start = $hit_length - $hit_start;

					# set min identity for different hsp length
					$identity =~ s/%//;
					if($hsp_length <= 50) {$min_identity = 95;}
					elsif($hsp_length > 50 && $hsp_length <= 100) {$min_identity = 96;}
					else{$min_identity = 97;}
                    
					# the if-eles order is importent (below)
					# The query-subject are not redundancy if the identitiy is low
					next if ($identity < $min_identity); 

					# the query is included in subject
					if ($query_start -1 <= $max_end_clip  && $query_to_end <= $max_end_clip)  
					{
						#my $hit_seq = $inset->{$hit_name}; 
					    	#print "type1\t".$hit_name."\t".$query_name."\t".$hit_seq."\n";#Êä³öºÏ²¢ÐÅÏ¢¹©ÈË¹¤Ð£¶Ô£¬µ÷ÊÔÓÃ
						return $hit_name."\tr";
					}

					# the query-subject are partially overlapped, need to combined, 
					if($hsp_length >= $min_overlap)
					{
						my $combined_seq;
					    	(my $query_seq = $$query) =~ s/^>[^\n]+\n//;# remove query name
					        $query_seq =~ s/\s//g;
						my $hit_seq = $inset->{$hit_name}; 
						if ($strand==1)
						{
						    	# Query on left
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_start <= $max_end_clip)
							{ 
						      		my $query_string = substr($query_seq, 0, $query_start); 							   
								my $hit_string = substr($hit_seq, $hit_start, $hit_to_start);
								#print "type2\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#Êä³öºÏ²¢ÐÅÏ¢¹©ÈË¹¤Ð£¶Ô£¬µ÷ÊÔÓÃ
								$combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
								return $combined_seq;
							} 
							# hit on the left
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
								my $hit_string = substr($hit_seq, 0, $hit_end); 
								my $query_string = substr($query_seq, $query_end, $query_to_end); 
								#print "type3\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#Êä³öºÏ²¢ÐÅÏ¢¹©ÈË¹¤Ð£¶Ô£¬µ÷ÊÔÓÃ
								$combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
								return $combined_seq;
							}
						}
						if ($strand==-1)
						{
							# Query one the left
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
								my $query_string = substr($query_seq, 0, $query_start); 							   
								my $hit_string = substr($hit_seq, 0, $hit_end-1);
								rcSeq(\$hit_string, 'rc'); 
								#print "type4\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#Êä³öºÏ²¢ÐÅÏ¢¹©ÈË¹¤Ð£¶Ô£¬µ÷ÊÔÓÃ
								$combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
								return $combined_seq;
							} 
							# hit one the left
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_start-1 <= $max_end_clip)
							{ 
								my $hit_string = substr($hit_seq, $hit_start, $hit_to_start); 
								rcSeq(\$hit_string, 'rc');
								my $query_string = substr($query_seq, $query_end, $query_to_end); 
								#print "type5\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#Êä³öºÏ²¢ÐÅÏ¢¹©ÈË¹¤Ð£¶Ô£¬µ÷ÊÔÓÃ
								$combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
								return $combined_seq;
							} 
						}
					}
				}
				return $hit_name."\tn";
			}
			
			#####################################
			#  start a new query hsp	    #
			#####################################
			%hsp = ();$hsp_count = 0;
			$query_name = $1; $query_length = ""; $hit_name = ""; $hit_length = "";
		}
		elsif (/\s+\((\S+)\sletters\)/) # get Query Length
		{
			$query_length = $1;
			$query_length =~ s/,//g;
		}
		
		elsif (/>(\S+)/)# Hit Name
		{
			#########################################
			#  put previous hsp info to hash        #
			#########################################
			if ($hsp_count > 0 || eof)
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;	
				$is_hsp = 0;
			}
			#################################
			#  se4t a new hit	        #
			#################################
		    $hit_name = $1; $hit_length = "";
		}
		elsif (/\s+Length = (\d+)/)
		{
				$hit_length = $1;
				$hit_length =~ s/,//g;
		}

		elsif (/Score =\s+(\S+) bits.+Expect(\(\d+\))? = (\S+)/)# hsp info
		{
			if ($hsp_count > 0 && $is_hsp == 1)
			{	
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
								$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
								$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#################################
			#  start a new hsp		#
			#################################
			$is_hsp = 1;
			$hsp_count++; 
			$score = $1; $evalue = $3;
			$evalue = "1$evalue" if ($evalue =~ m/^e/);
			$query_start = 0; $query_end = 0; $hit_start = 0; $hit_end = 0;

		}
		elsif (/\s+Identities = (\d+)\/(\d+)\s+\((\S+)\)/ && $hsp_count >= 1)
		{
			$identity = $1/$2*100;
			$identity = sprintf("%.".(2)."f", $identity);
			if ( $1 > $2 ) { $hsp_length = $1; } else { $hsp_length = $2; }
		}

		elsif (/\s+Strand = (\S+) \/ (\S+)/ && $hsp_count >= 1) 
		{
			if ( $2 eq "Plus" ) { $strand = 1;} else { $strand = -1;}
		}

		elsif (/Query\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1)
		{
			if ($query_start == 0) { $query_start = $1; $query_start =~ s/,//g;}
			$query_end = $3;
			$query_end =~ s/,//g;
		}

		elsif (/Sbjct\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1) 
		{
			if ( $strand == -1 )
			{
				if ($hit_end == 0) { $hit_end = $1; $hit_end =~ s/,//g;  };
				$hit_start = $3;
				$hit_start =~ s/,//g;
			}
			else
			{
				if ($hit_start == 0) { $hit_start = $1; $hit_start =~ s/,//g; };
				$hit_end = $3;
				$hit_end =~ s/,//g;
			}
		}
		else
		{
			next;
		}
	}
	$bfh->close;
	return "null\tn";
}

=head2
 reverse comp seq
=cut

sub rcSeq 
{
        my ($seq_r, $tag) = @_;
        defined $tag or $tag = 'rc';

        my ($Is_r, $Is_c) = (0,0);
        $tag =~ /r/i and $Is_r = 1;
        $tag =~ /c/i and $Is_c = 1;

        !$Is_r and !$Is_c and die "Wrong Input for function rcSeq! $!\n";
        $Is_r and $$seq_r = reverse ($$seq_r);
        $Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDH/tgcaayrmkvbhdTGCAAYRMKVBHD/; # edit on 2010-11-14;

        return 0;
}

1;
