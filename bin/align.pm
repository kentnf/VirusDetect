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

  filter_SAM -- filter sam file

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

1;
