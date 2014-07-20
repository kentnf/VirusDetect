#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use IO::File;
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/bin";
use Util;
use align;

my $usage = <<_EOUSAGE_;
########################################################################
# virus_detect.pl --reference [FILE] [option] input
#  
# Basic Options:
#  --reference		The name of a fasta file containing the virus 
#                       reference sequences  [vrl_plant]
#  --host_reference	Name of host reference file for subtraction [Null]
#  --thread_num		Number of threads (multi-threading mode) [8] 
# 
# BWA-related options (align sRNA to reference virus database or host 
#  sequences):
#  --max_dist		Maximum edit distance [1]  
#  --max_open		Maximum number of gap opens [1]  
#  --max_extension	Maximum number of gap extensions [1]  
#  --len_seed		Take the first INT subsequence as seed [15] 
#  --dist_seed		Maximum edit distance in the seed [1]  
# 
# Megablast-related options (remove redundancy within virus contigs):
#  --strand_specific	Only for sequences assembled from strand-
#                        specific RNA-seq  [Not selected]
#  --min_overlap	The minimum overlap length between two 
#                        contigs to be combined [30]
#  --max_end_clip	The maximum length of end clips [6]
#  --mis_penalty	Penalty for a nucleotide mismatch [-1]
#  --gap_cost		Cost to open a gap [2] 
#  --gap_extension	Cost to extend a gap [1] 
#
# Megablast-related options (align virus contigs to reference virus 
#  database for virus identification):
#  --word_size	     	[11] 
#  --exp_value	     	[1e-5]
#  --percent_identity 	[80] 
#  --mis_penalty_b   	Penalty for a nucleotide mismatch [-1]
#  --gap_cost_b      	Cost to open a gap [2] 
#  --gap_extension_b 	Cost to extend a gap [1]
#
# Result filter options:
#  --hsp_cover		Coverage cutoff of a reported virus contig by
#                        reference virus sequences [0.75]
#  --coverage_cutoff	Coverage cutoff of a reported virus reference 
#                        sequences by assembled virus contigs [0.1] 
#  --depth_cutoff	Depth cutoff of a reported virus reference [5]
########################################################################
_EOUSAGE_
;

my $file_list;

# basic options
my $reference= "vrl_plant";			# virus sequence
my $host_reference;       			# host reference
my $thread_num = 8; 				# thread number

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
	'r|reference=s'	=> 		\$reference,
	'h|host_reference=s' => 	\$host_reference,
	't|thread_num=i' => 		\$thread_num,

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
	'percent_identity=s' => 	\$percent_identity,
	'mis_penalty_b=i' => 	\$mis_penalty_b,
	'gap_cost_b=i' => 	\$gap_cost_b,
	'gap_extension_b=i' => 	\$gap_extension_b,

	'hsp_cover=s' =>	\$hsp_cover,
	'diff_ratio=s' => 	\$diff_ratio,
	'diff_contig_cover=s' =>\$diff_contig_cover,
	'diff_contig_length=s'=>\$diff_contig_length,

	'coverage_cutoff=f' =>	\$coverage_cutoff,
	'depth_cutoff=f' =>	\$depth_cutoff
);

my $debug = 1;

# check input file
die "Please input one file:\n\n$usage\n" unless scalar(@ARGV) == 1;
die "Input file is not exist, please check it.\n\n$usage\n" unless -s $ARGV[0];

# check filetype

foreach my $sample (@ARGV) 
{
	my $file_type = Util::detect_FileType($sample);
	my $sample_base = basename($sample);

	# set path and folder folder
	my $WORKING_DIR   = cwd();					# set current folder as working folder
	my $DATABASE_DIR  = ${FindBin::RealBin}."/databases";		# set database folder
	my $BIN_DIR       = ${FindBin::RealBin}."/bin";			# set script folder 
	my $TEMP_DIR      = $WORKING_DIR."/".$sample_base."_temp";	# set temp folder
	$reference	  = $DATABASE_DIR."/".$reference;		# set reference
	my $seq_info	  = $DATABASE_DIR."/vrl_genbank.info";		# set vrl info
	print "Working: $WORKING_DIR\nDatabase: $DATABASE_DIR\nBin: $BIN_DIR\nTemp: $TEMP_DIR\n" if $debug;

	# create temp folder and create link for sample
	Util::process_cmd("mkdir $TEMP_DIR", $debug) unless -e $TEMP_DIR;
	Util::process_cmd("ln -s $WORKING_DIR/$sample $TEMP_DIR/$sample_base", $debug) unless -e "$TEMP_DIR/$sample_base";
	$sample = "$TEMP_DIR/$sample_base";				# change the sample name to the linke from this step

	# set parameters for remove reduncancy (rr)
	my $rr_blast_word_size = int($min_overlap/3);
	my $rr_hits_returns = 10;
	my $rr_blast_parameters = "-F F -a $thread_num -W $rr_blast_word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $rr_hits_returns";
	if ($strand_specific) { $rr_blast_parameters .=" -S 1"; }

	# part A: 1. align reads to plant virus; 2. extract aligned seqs; 3. remove redundancy contigs
	my $align_parameters = "-n $max_dist -o $max_open -e $max_extension -i 0 -l $len_seed -k $dist_seed -t $thread_num";
	my $align_program    = "$BIN_DIR/bwa";

	print "Align Program: $align_program\nAlign Parameters: $align_parameters\nAlign Input File: $sample\nAlign Output File: $sample.sam\n" if $debug;

	print ("####################################################################\nprocess sample $sample_base\n");
	Util::print_user_message("Align reads to reference virus sequence database");
	align::align_to_reference($align_program, $sample, $reference, "$sample.sam", $align_parameters, $TEMP_DIR, $debug);
	align::filter_SAM($sample.".sam");	# filter out unmapped, 2nd hits, only keep the best hit
	Util::process_cmd("$BIN_DIR/samtools view -bt $reference.fai $sample.sam > $sample.bam 2> $TEMP_DIR/samtools.log", $debug) unless (-s "$sample.bam");
	Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted 2> $TEMP_DIR/samtools.log", $debug) unless (-s "$sample.sorted.bam");
	Util::process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pileup 2> $TEMP_DIR/samtools.log", $debug) unless (-s "$sample.pre.pileup");
	align::pileup_filter("$sample.pre.pileup", "$seq_info", "$coverage", "$sample.pileup", $debug) unless (-s "$sample.pileup");	# filter pileup file 
	align::pileup_to_contig("$sample.pileup", "$sample.aligned", 40, 1, 'ALIGNED') unless -s "$sample.aligned"; # input, output, min_len, min_depth, prefix
	align::removeRedundancy("$sample.aligned", $sample, $rr_blast_parameters, $max_end_clip, $min_overlap, $BIN_DIR, $TEMP_DIR, $debug);
	
	

	# part B: 1. remove host related reads  2. de novo assembly 3. remove redundancy contigs
	# parameter for velvet: $sample, $output_contig, $kmer_start, $kmer_end, $coverage_start, $coverage_end, $objective_type, $bin_dir, $temp_dir, $debug
	if( $host_reference ){
		Util::print_user_message("Align reads to host reference sequences");
		align::align_to_reference($align_program, $sample, $host_reference, "$sample.sam", $align_parameters, $TEMP_DIR, $debug);
		align::generate_unmapped_reads("$sample.sam", "$sample.unmapped");
		Util::print_user_message("De novo assembly");
		align::velvet_optimiser_combine("$sample.unmapped", "$sample.assembled", 9, 19, 5, 25, $objective_type, $BIN_DIR, $TEMP_DIR, $debug);
	}	
	else
	{
		Util::print_user_message("De novo assembly");
		align::velvet_optimiser_combine($sample, "$sample.assembled", 9, 19, 5, 25, $objective_type, $BIN_DIR, $TEMP_DIR, $debug);
	}

    	# removeRedundancy($file_list, $file_type, "assemblied", "NOVEL", $parameters_remove_redundancy);
	
	# combine the known and unknown virus, remove redundancy of combined results, it must be using strand_specific parameter
	Util::print_user_message("Remove redundancies in virus contigs");
	Util::process_cmd("cat $sample.aligned $sample.assembled > $sample.combined", $debug);

	# $cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix combined ".
    	#   				"--contig_prefix CONTIG --strand_specific $parameters_remove_redundancy";
	# Util::process_cmd($cmd_removeRedundancy);
    	# removeRedundancy($file_list, $file_type, "combined", "CONTIG", $parameters_remove_redundancy);
	# identify the virus
    
	Util::print_user_message("Virus identification");
	my $cmd_identify = "$BIN_DIR/virus_identify.pl ";
	$cmd_identify .= "--word_size $word_size --exp_value $exp_value --identity_percen $percent_identity ";
	$cmd_identify .= "--cpu_num $thread_num --mis_penalty $mis_penalty_b --gap_cost $gap_cost_b --gap_extension $gap_extension_b ";
	$cmd_identify .= "--hsp_cover $hsp_cover --diff_ratio $diff_ratio --diff_contig_cover $diff_contig_cover --diff_contig_length $diff_contig_length ";
	$cmd_identify .= "--coverage_cutoff $coverage_cutoff --depth_cutoff $depth_cutoff $sample $sample.combined";
	Util::process_cmd($cmd_identify, $debug);

	# delete temp files and log files 
	# system("rm -r *.log");
	# system("rm -r temp");

	Util::print_user_message("Finished");
	print ("####################################################################\n\n");

}
