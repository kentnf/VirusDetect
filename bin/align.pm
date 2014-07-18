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
		Util::process_cmd("$align_program aln $parameters $reference $input_file 1> $sai 2>> $log", $debug) unless -s $sai;
		Util::process_cmd("$align_program samse -n 10000 -s $reference $sai $input_file 1> $output_file 2>> $log", $debug) unless -s $output_file;
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

=head
sub renameFasta
{
	my ($input_fasta_file, $output_fasta_file, $prefix) = @_;
    
	my $seq_num = 0;
	my $out = IO::File->new(">".$output_fasta_file) || die $!;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_fasta_file);
    while(my $inseq = $in->next_seq)
	{
		$seq_num++;
		print $out ">".$prefix.$seq_num."\n".$inseq->seq."\n";
	}
	$out->close;
}

sub bwa_remove{
    my ($file_list, $reference, $host_reference, $parameters_bwa_align) = @_;
    my $ref = (split(/\//, $reference))[-1];
    Util::process_cmd("$BIN_DIR/bowtie2-build -f $ref $reference");
    #Util::process_cmd("$BIN_DIR/bowtie2-build -f $ref $reference") unless (-s $ref.".1.bt2");
    my $sample;
    my $i=0;
    open(IN, "$file_list") || die $!;
    while (<IN>) {
        
        chomp;
        $i=$i+1;
        $sample=$_;
        print "#processing sample $i by $0: $sample\n";
        
        # command lines
        Util::process_cmd("$BIN_DIR/bowtie2 --quiet -N $max_dist -p $thread_num -L $len_seed --sensitive -q -x $DATABASE_DIR"."/vrl_plant -U $sample -S $sample.sam", 1);
        #Util::process_cmd("$BIN_DIR/bowtie2 --quiet -N $max_dist -p $thread_num -L $len_seed --sensitive -q -x $DATABASE_DIR"."/vrl_plant -U $sample -S $sample.sam") unless (-s "$sample.sam");
        # generate unmapped reads
        #my ($num_unmapped_reads, $ratio) = generate_unmapped_reads("$sample.sam", "$sample.unmapped");
        #my ($input_SAM, $output_reads) = @_;
        
        my %mapped_reads;
        my ($num_unmapped_reads, $ratio) = (0, 0);
        
        my $in  = IO::File->new("$sample.sam") || die $!;
        my $out = IO::File->new(">"."$sample.unmapped") || die $!;
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
        ($num_unmapped_reads, $ratio) = ("$sample.sam", "$sample.unmapped");
        Util::print_user_submessage("$num_unmapped_reads ($ratio) of reads can not be aligned to host reference");
        
        # remove temp file for read alignment
        #system("rm $sample.sai");
        #system("rm $sample.sam");
    }
    close(IN);
}

sub removeRedundancy{
    my ($file_list, $file_type, $input_suffix, $contig_prefix, $parameters_remove_redundancy) = @_;
    open(IN1,$file_list) || die "Can't open the file $file_list\n";
    my ($j, $sample, $contig_file);
    $j=0;
    while(<IN1>){
        $j=$j+1;
        chomp;
        $sample = $_;
        #print "#processing sample $j by $0: $sample\n";
        $contig_file = $sample.".".$input_suffix;
        
        # get aligned files size, do not remove redundancy if file size is 0
        my $file_size = -s "$contig_file";
        next if $file_size == 0;
#        if ($input_suffix =~ "combined") {
#            my $in  = IO::File->new("$sample.$input_suffix") || die $!;
#            my $out = IO::File->new(">"."$sample.$input_suffix.fasta") || die $!;
#            while(<$in>)
#            {
#                chomp;
#                next if $_ =~ m/^@/;
#                my @a = split(/\t/, $_);
#                
#                
#                if ( scalar(@a > 8) ) {
#                    print $out ">"."\@$a[0]\n";
#                    print $out "\@$a[9]\n";
#                }
#            }
#            $in->close;
#            $out->close;
#            Util::process_cmd("mv $sample.$input_suffix.fasta $sample.$input_suffix");
#        }
        # if file has sequence, move it to temp folder to remove redundancy
        # 1. remove simple repeate sequence using mask
        Util::process_cmd("$BIN_DIR/dust $sample.$input_suffix 1> $sample.masked 2> $tf/dust.log");
        Util::process_cmd("$BIN_DIR/trim_XNseq1.pl $sample.masked $sample.$input_suffix 0.8 40 > $sample.$input_suffix.1");
        Util::process_cmd("rm $sample.masked");
        
        my ($before_contig_num, $after_contig_num, $i);
        $i = 1;								# get the number of contig file, default is 1
        $before_contig_num = `grep -c \'>\' $sample.$input_suffix.$i`;	# get the seq number before remove redundancy
        $after_contig_num  = 0;						# this is seq number after remove redundancy
        
        if( $before_contig_num == 0 ){                          # if file size = 0, create blank file, and exit the loop
            Util::process_cmd("touch $sample.$input_suffix");
            next;
        }
        
        # if the new contig number != old contig number, continue remove redundancy
        while( $after_contig_num != $before_contig_num )
        {
            # the default output is the input_inset;
            Util::process_cmd ("$BIN_DIR/removeRedundancy.pl --input $sample.$input_suffix.$i --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
            Util::process_cmd("rm $sample.$input_suffix.$i");# rm old file
            my $remove_redundancy_result = "$sample.$input_suffix.$i"."_inset";
            
            $i++;
            $before_contig_num = $after_contig_num; # renew contig_num1
            # renew contig_num2
            $after_contig_num =  `grep -c \'>\' $remove_redundancy_result`; # get seq number of new contig
            chomp($after_contig_num);
            Util::process_cmd("mv $remove_redundancy_result $sample.$input_suffix.$i");
        }
        
        if ($after_contig_num > 1) {
            Util::print_user_submessage("$after_contig_num of unique contigs were generated");
        } elsif ($after_contig_num == 1) {
            Util::print_user_submessage("$after_contig_num of unique contig was generated");
        } elsif ($after_contig_num == 0) {
            Util::print_user_submessage("None of unique contig was generated, should be error!");
        }
        
        # finish remove redundancy, next for base correction
        my $sample_reference = "$sample.$input_suffix.$i";	# sample_reference file after remove Redundancy
        my $sample_reads     = $sample;				# read file, need to re-aligned to sample_reference file
        
        #aligment -> sam -> bam -> sorted bam -> pileup
        my $format = "-q"; if ($file_type eq "fasta") {$format="-f"};
        Util::process_cmd("$BIN_DIR/bowtie-build --quiet -f $sample_reference $sample") unless (-e "$sample.1.amb");
        Util::process_cmd("$BIN_DIR/samtools faidx $sample_reference 2> $tf/samtools.log") unless (-e "$sample_reference.fai");
        Util::process_cmd("$BIN_DIR/bowtie --quiet $sample -v 1 -p $cpu_num $format $sample -S -a --best $sample.sam") unless (-s "$sample.sam");
        Util::process_cmd("$BIN_DIR/samtools view -bt $sample_reference.fai $sample.sam > $sample.bam 2> $tf/samtools.log") unless (-s "$sample.bam");
        Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted 2> $tf/samtools.log") unless (-s "$sample.sorted.bam");
        Util::process_cmd("$BIN_DIR/samtools mpileup -f $sample_reference $sample.sorted.bam > $sample.pileup 2> $tf/samtools.log") unless (-s "$sample.pileup");
        
        $file_size = -s "$sample.pileup";		# get file size
        if( $file_size == 0 ){				# if file size = 0, create blank file, and exit the loop
            Util::process_cmd("touch $sample.$input_suffix");
            next;
        }
        
        $i++;
        Util::process_cmd("java -cp $BIN_DIR extractConsensus $sample 1 40 $i");
        #renameFasta("$sample.contigs$i.fa", "$sample.$input_suffix", $contig_prefix);
        #my ($input_fasta_file, $output_fasta_file, $prefix) = @_;
        renameFasta("$sample.contigs$i.fa", "$sample.$input_suffix",$contig_prefix);
        
        # remove temp files
        system("rm $sample.sam");
        system("rm $sample.bam");
        system("rm $sample.sorted.bam");
        system("rm $sample.pileup"); # must delete this file for next cycle remove redundancy
        system("rm $sample_reference");
        system("rm $sample_reference.fai");
        system("rm $tf/*.ebwt");
        system("rm $sample.contigs$i.fa");
    }
    close(IN1);
}

=head
sub Velvet_Optimiser_combined{
    my ($parameters, $file_list, $input_suffix, $file_type, $bjective_type, $hash_end, $coverage_end, $output_suffix);
    my $sample;
    my $sampleNum=0;
    my $current_folder;
    my $statfile;
    my $objective;						# the objective for each run
    my $max_objective;						# the maximum objective
    my $opt_hash_length=$hash_start; 				# the optimization hash_length when objective is max
    my $opt_coverage=$coverage_start;				# the optimization coverage when objective is max
    my $opt_avgLen=0;                				# the avg Length when objective is max
    open(IN, "$file_list");
    open(OUT1, ">$tf/optimization.log") || die $!; 		# save optimization information
    open(OUT2, ">$tf/optimization.result") || die $!;		# save final optimization result
    while (<IN>) {
		$sampleNum = $sampleNum+1;
		chomp;
		$sample = $_;
		#print "#processing sample $sampleNum: $sample\n";
		
		# reset parameters
		$max_objective = 0;
		$opt_hash_length = $hash_start;
		$opt_coverage = $coverage_start;
		
		# optimize k-mer length using fixed coverage
		for(my $i=$hash_start; $i<= $hash_end; $i=$i+2) {
			runVelvet($sample, $i, $coverage_start);
			$current_folder = $sample."_".$i."_".$coverage_start;
			$statfile = $current_folder."/contigs.fa";
			my $aa = contigStats($statfile);        # return hash reference
			if ( defined $aa->{$objective_type} ) {
				$objective = $aa->{$objective_type};
				print OUT1 $i."\t".$coverage_start."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";
                
				# output best hash length, coverage, objective, avgLength, contigs num
				if ( $objective > $max_objective ) {
					print OUT1 "yes"."\n"; #print yes if the optimization value improved(higher than before)
					$max_objective = $objective;
					$opt_hash_length = $i;
					$opt_avgLen=$aa->{avgLen};
				}
			}
			Util::process_cmd("rm $current_folder -r");
		}
		# optimize coverage using fixed k-mer length
		for(my $j=$coverage_start+2; $j<=$coverage_end; $j=$j+1) {  # start from 7, 注意这里从7开始，因为5上面都算过了
			runVelvet($sample,$opt_hash_length,$j);
			$current_folder=$sample."_".$opt_hash_length."_".$j;
			$statfile=$current_folder. "/contigs.fa";
			my $aa=contigStats($statfile);
			if ( defined $aa->{$objective_type} ) {
				$objective=$aa->{$objective_type};
				print OUT1 $opt_hash_length."\t".$j."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";
                
				# output the best hast length, coverage, objective, avglength, contigs num
				if($objective>$max_objective){
					print OUT1 "yes"."\n";# print yes of the optimization value improved
					$max_objective=$objective;
					$opt_coverage=$j;
					$opt_avgLen=$aa->{avgLen};
				}
			}
			Util::process_cmd("rm $current_folder -r");
		}
		print OUT2 $sample."\t".$opt_hash_length."\t".$opt_coverage."\t".$max_objective."\t".$opt_avgLen."\n";
	}
	close(IN);
	close(OUT1);
	close(OUT2);
    my $i=0;
	my $resultDir;
    
	open(IN, "$parameters");
	while (<IN>) {
		$i=$i+1;
		chomp;
		my @a = split(/\t/, $_);				# array: sample, Hash_len, Coverage_cutoff
		runVelvet1($a[0],$a[1],$a[2]);				# assemblied reads into contigs using Velvet
        
		# move the assemblied contigs to ouput file
		$resultDir = $a[0]."_".$a[1]."_".$a[2];
		my $contigName = $a[0].".".$output_suffix;		# output file name
		Util::process_cmd("mv $resultDir/contigs.fa $contigName");# move assemblied contigs to output file
		#system("rm -r $resultDir");				# remove assemblied output folder
        
		my $count = `grep -c \'>\' $contigName `;		#  stat the number of contigs
		chomp($count);
		Util::print_user_submessage("Assemblied Contig No: ".$count);
	}
	close(IN);
}
sub runVelvet {
	($sample1, $hash_length, $cov_cutoff) = @_;
	my $outputDir=$sample1."_".$hash_length."_".$cov_cutoff;
    
	my $file;
	if ($input_suffix)	{ $file = "$sample1.$input_suffix"; }
	else 			{ $file = $sample1; }
    
	Util::process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $file >> $tf/velvet.log");
	Util::process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 >> $tf/velvet.log");
}
sub runVelvet1 {
    ($sample, $hash_length, $cov_cutoff) = @_;
	my $outputDir = $sample."_".$hash_length."_".$cov_cutoff;
	my $file;
	if ($input_suffix)	{ $file = "$sample.$input_suffix"; }
	else			{ $file = $sample; }
	Util::process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $file > $tf/velvet.log");
	Util::process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 > $tf/velvet.log");
}

sub contigStats {
	
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

=cut

1;
