#!/usr/bin/perl
package align;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use Bio::SeqIO;
use lib "$FindBin::RealBin/bin/PerlLib";
use Util;
use IO::File;
use Exporter;


our @ISA = qw( Exporter );

our @EXPORT_OK = qw( renameFasta bwa_remove removeRedundancy filter_SAM files_combine2 );
################################
our $WORKING_DIR   = cwd();				# set current folder as working folder
our $DATABASE_DIR  = $WORKING_DIR."/databases";	# set database folder
our $BIN_DIR	  = ${FindBin::RealBin}."/bin";		# set script folder
our $TEMP_DIR	  = $WORKING_DIR."/temp";		# set temp folder
our $tf = $TEMP_DIR;
our $output_suffix;	# output suffix
# basic options
my $file_type= "fastq";				# [autodetect for it] input file type, fasta or fastq
my $reference= "vrl_plant";			# virus sequence
my $host_reference;       			# host reference
my $thread_num = 8; 				# thread number
my $file_list;

# paras for BWA
my $max_dist = 1;  				# max edit distance
my $max_open = 1;  				# max gap opening
my $max_extension = 1; 				# max gap extension (gap length)
my $len_seed = 15; 				# bwa seed length
my $dist_seed = 1; 				# bwa seed max edit distance

# paras for megablast detection (remove redundancy )
my $strand_specific;  				# switch for strand specific transcriptome data?
my $min_overlap = 30; 				# minimum overlap for hsp combine
my $max_end_clip = 6; 				# max end clip for hsp combine
my $mis_penalty = -1;     			# megablast mismatch penlty, minus integer
my $gap_cost = 2;         			# megablast gap open cost, plus integer
my $gap_extension = 1;    			# megablast gap extension cost, plus integer
my $cpu_num = 8;

# paras for blast && identification
my $word_size = 11;
my $exp_value = 1e-5;				#
my $percent_identity = 25;			# tblastx “‘µ∞∞◊÷ –Ú¡–¿¥±»∂‘ ±hspµƒ◊Ó–°Õ¨“ª–‘
my $mis_penalty_b = -1;				# megablast mismatch penlty, minus integer
my $gap_cost_b = 2;				# megablast gap open cost, plus integer
my $gap_extension_b = 1;			# megablast gap extension cost, plus integer

my $filter_query = "F";				# megablast switch for remove simple sequence
my $hits_return = 500;				# megablast number of hit returns

# paras for result filter
my $hsp_cover = 0.75;
my $coverage_cutoff = 0.1;			# coverage cutoff for final result
my $depth_cutoff = 5;				# depth cutoff for final result

# disabled parameters or used as fixed value
my $input_suffix='clean'; 			# input_suffix, disabled
my $coverage=0.3;  				# √øÃı≤Œøº–Ú¡–»Áπ˚±ªreads∏≤∏«µƒ≤ø∑÷’º»´≥§±»¿˝µƒ„–÷µ
my $objective_type='maxLen';			# objective type for Velvet assembler: n50°¢maxLen, avgLen
my $diff_ratio= 0.25;
my $diff_contig_cover = 0.5;
my $diff_contig_length= 100;

# get input paras #
GetOptions(
'file_type=s'	=> 	\$file_type,
'reference=s'	=> 	\$reference,
'host_reference=s' => 	\$host_reference,
'thread_num=i' => 	\$thread_num,
'file_list=s' => \$file_list,

'max_dist=i' => 	\$max_dist,
'max_open=i' => 	\$max_open,
'max_extension=i' => 	\$max_extension,
'len_seed=i' => 	\$len_seed,
'dist_seed=i' => 	\$dist_seed,


'strand_specific!' => 	\$strand_specific,
'min_overlap=i' => 	\$min_overlap,
'max_end_clip=i' => 	\$max_end_clip,
'cpu_num=i' 		=> \$cpu_num,
'mis_penalty=i' => 	\$mis_penalty,
'gap_cost=i' => 	\$gap_cost,
'gap_extension=i' => 	\$gap_extension,

'word_size=i' =>  	\$word_size,
'exp_value-s' =>  	\$exp_value,
'percent_identity=s' => 	\$percent_identity,	# tblastx“‘µ∞∞◊÷ –Ú¡–¿¥±»∂‘ ±hspµƒ◊Ó–°Õ¨“ª–‘
'mis_penalty_b=i' => 	\$mis_penalty_b,
'gap_cost_b=i' => 	\$gap_cost_b,
'gap_extension_b=i' => 	\$gap_extension_b,

'hsp_cover=s' =>	\$hsp_cover,
'diff_ratio=s' => 	\$diff_ratio,
'diff_contig_cover=s' =>\$diff_contig_cover,
'diff_contig_length=s'=>\$diff_contig_length,

'coverage_cutoff=f' =>	\$coverage_cutoff,
'depth_cutoff=f' =>	\$depth_cutoff,
'output_suffix=s'	=> \$output_suffix,

'velvet_dir=s'		=> \$velvet_dir,
'parameters=s' 		=> \$parameters,
'objective_type=s'	=> \$objective_type,
'hash_end=i'		=> \$hash_end,
'coverage_end=i'	=> \$coverage_end
);

sub bowtie2_align
{
	my ($input, $output, $index, $parameters) = @_;

	
}

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

sub filter_SAM
{
	my $input_SAM = shift;
	my $temp_SAM = $input_SAM.".temp";
	my ($total_count, $filtered_count) = (0, 0);
    my ($query_col, $opt_col) = (0, 11); 	# query and option column number for sam
	my $max_distance = 2;			# set $max_distance for all selected hits
	my $bestEditDist = -1;			# set best edit distance
	my @alignment = ();			# alignment to array
	my $pre_query_name = '';		# previous query name
	
    
	my $in  = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$temp_SAM) || die $!;
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^@/) { print $out $_."\n"; next; }
		my @a = split(/\t/, $_);
		if ( $a[1] == 4 ) { $filtered_count++; }
		else {	print $out $_."\n"; }
		$total_count++;
	}
    my ($total_count, $kept_align) = (0,0);
	$in->close;
	$out->close;
	Util::process_cmd("mv $temp_SAM $input_SAM");
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as unmapped reads, only for BWA\n";
    
	my $in  = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$temp_SAM) || die $!;
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^@/) { next; #print $out $_."\n"; next;
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
	my $filtered_count = $total_count - $kept_align;
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as 2ndhits reads, only for BWA\n";
}
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

1;