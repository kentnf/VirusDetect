#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use Bio::SeqIO;
use lib "$FindBin::RealBin/bin/PerlLib";
use Util;
use IO::File;
use align qw( removeRedundancy bwa_remove );

my $usage = <<_EOUSAGE_;
# 12345
########################################################################
# virus_detect.pl --file_type [String] --reference [FILE] [option] input
#  
# Basic Options:
#  --file_type		Format of input file (fastq or fasta)  [fastq]
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

################################
##  set file folder path      ##
################################
my $WORKING_DIR   = cwd();				# set current folder as working folder
my $DATABASE_DIR  = $WORKING_DIR."/databases";	# set database folder
my $BIN_DIR	  = ${FindBin::RealBin}."/bin";		# set script folder 
my $TEMP_DIR	  = $WORKING_DIR."/temp";		# set temp folder
my $tf = $TEMP_DIR;

# basic options
our $file_type= "fastq";				# [autodetect for it] input file type, fasta or fastq
our $reference= "vrl_plant";			# virus sequence
our $host_reference;       			# host reference
our $thread_num = 8; 				# thread number
our $file_list;
our $output_suffix;

# paras for BWA
our $max_dist = 1;  				# max edit distance
our $max_open = 1;  				# max gap opening
our $max_extension = 1; 				# max gap extension (gap length)
our $len_seed = 15; 				# bwa seed length
our $dist_seed = 1; 				# bwa seed max edit distance

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
	'depth_cutoff=f' =>	\$depth_cutoff
);

# check input file
die "Please input one file:\n\n$usage\n" unless scalar(@ARGV) == 1;
die "Input file is not exist, please check it.\n\n$usage\n" unless -s $ARGV[0];

# check filetype 
my $file_type_check = Util::detect_FileType($ARGV[0]);
if ($file_type_check ne $file_type){
	print "Warning, file type is not correct, and has been changed to $file_type_check\n";
	$file_type = $file_type_check;
}

# set paramsters
my $parameters_remove_redundancy = "--min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $thread_num ".
				   "--mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension";
if ($strand_specific) { $parameters_remove_redundancy.=" --strand_specific"; }

my $parameters_bwa_align = "--max_dist $max_dist --max_open $max_open --max_extension $max_extension ".
			   "--len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num";


#Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
#reference_in            comma-separated list of files with ref sequences
#bt2_index_base          write bt2 data to files with this dir/basename
sub files_combine2
{
	my $file_list = shift;
	my $fh = IO::File->new($file_list) || die "Can not open innput file list $file_list $!\n";
	while(<$fh>)
	{
		chomp;
		my $file = $_;
		my $file_aligned	= $file.".aligned";
		my $file_assemblied	= $file.".assemblied";
		my $file_combined	= $file.".combined";
		Util::process_cmd("cat $file_aligned $file_assemblied > $file_combined");
	}
	$fh->close;
}

# main
main: {
	# create temp folder
	# create softlink for input file
	# create temp file list for input file
	# *the file in file list has full path
	Util::process_cmd("mkdir $TEMP_DIR") unless -e $TEMP_DIR;
	my $file_list = "$TEMP_DIR/temp_file_list";
	my $fh_list = IO::File->new(">".$file_list) || die "Can not open temp file list $file_list $!\n";
	print $fh_list $TEMP_DIR."/".$ARGV[0]."\n"; 
	$fh_list->close;
	Util::process_cmd("ln -s $WORKING_DIR/$ARGV[0] $TEMP_DIR/$ARGV[0]") unless -e "$TEMP_DIR/$ARGV[0]";
    my $ref_list = "$TEMP_DIR/temp_ref_list";
	my $refopen = IO::File->new(">".$ref_list) || die "Can not open reference file list $file_list $!\n";
    print $refopen $DATABASE_DIR."/".$reference."\n";
    $refopen->close;
	# detect known virus. Including:
	# 1. align small RNA to plant virus
	# 2. remove redundancy for aligned contigs 
			  print ("####################################################################\nprocess sample $ARGV[0]\n");
	Util::print_user_message("Align reads to reference virus sequence database");
    
    my $cmd_align = "$BIN_DIR/alignAndCorrect.pl --file_list $file_list --reference $DATABASE_DIR/$reference --coverage $coverage --output_suffix aligned";
	#my $cmd_align = "$BIN_DIR/alignAndCorrect.pl --file_list $file_list --reference $ref_list --coverage $coverage --output_suffix aligned";
	Util::process_cmd($cmd_align, 1);
    removeRedundancy($file_list, $file_type, "aligned", "KNOWN", $parameters_remove_redundancy);
	#my $cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix aligned ".
	#			   "--contig_prefix KNOWN $parameters_remove_redundancy";
    #Util::process_cmd($cmd_removeRedundancy);
    #removeRedundancy($file_list, $file_type, "aligned", "KNOWN", $parameters_remove_redundancy);

	# detect unknown virus, then assembly them, it including
	# 1. remove host related reads  
	# 2. de novo assembly
	# 3. remove redundancy contigs after assembly
	
	if( $host_reference ){
		Util::print_user_message("Align reads to host reference sequences");
		#Util::process_cmd("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$host_reference $parameters_bwa_align");
		bwa_remove($file_list, $reference, $host_reference, $parameters_bwa_align);

		# the input suffix of unmapped reads is 'unmapped'
		# the seq from sam file is fastq format, no matter the format of input file
		Util::print_user_message("De novo assembly");
		Util::process_cmd("$BIN_DIR/Velvet_Optimiser_combined.pl --parameters $TEMP_DIR/optimization.result --file_list $file_list --input_suffix unmapped --file_type fastq --objective_type $objective_type --hash_end 19 --coverage_end 25 --output_suffix assemblied");
	}	
	else
	{
		Util::print_user_message("De novo assembly");
		Util::process_cmd("$BIN_DIR/Velvet_Optimiser_combined.pl --parameters $TEMP_DIR/optimization.result --file_list $file_list --file_type $file_type --objective_type $objective_type --hash_end 19 --coverage_end 25 --output_suffix assemblied");
	}
    #the script is messing up here in this removeRedundancy line:
    removeRedundancy($file_list, $file_type, "assemblied", "NOVEL", $parameters_remove_redundancy);
	
	# combine the known and unknown virus, remove redundancy of combined results, it must be using strand_specific parameter
	Util::print_user_message("Remove redundancies in virus contigs");
	files_combine2($file_list);
	#$cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix combined ".
    #   				"--contig_prefix CONTIG --strand_specific $parameters_remove_redundancy";
	#Util::process_cmd($cmd_removeRedundancy);
    removeRedundancy($file_list, $file_type, "combined", "CONTIG", $parameters_remove_redundancy);
	# identify the virus
    
	Util::print_user_message("Virus identification");
	my $cmd_identify = "$BIN_DIR/virus_identify.pl ";
       	$cmd_identify .= "--file_list $file_list --file_type $file_type --reference $reference --contig_type combined ";
	$cmd_identify .= "--word_size $word_size --exp_value $exp_value --percent_identity $percent_identity ";
	$cmd_identify .= "--cpu_num $thread_num --mis_penalty $mis_penalty_b --gap_cost $gap_cost_b --gap_extension $gap_extension_b ";
	$cmd_identify .= "--hsp_cover $hsp_cover --diff_ratio $diff_ratio --diff_contig_cover $diff_contig_cover --diff_contig_length $diff_contig_length ";
	$cmd_identify .= "--coverage_cutoff $coverage_cutoff --depth_cutoff $depth_cutoff ";
	Util::process_cmd($cmd_identify);

	# delete temp files and log files 
	# system("rm -r *.log");
	# system("rm -r temp");

	Util::print_user_message("Finished");
	print ("####################################################################\n\n");
}
