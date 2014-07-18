#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::RealBin/PerlLib";
use Util;
use Bio::SeqIO;
use Cwd;

my $usage = <<_EOUSAGE_;

#########################################################################################################
# removeRedundancy_batch.pl --file_list <FILE> --file_type <String> --input_suffix <String> --contig_type <String> --contig_prefix <String> 
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#
# Required(4):
#  --file_list       A txt file containing a list of input file names without any suffix
#  --file_type       The format of input files(fastq or fasta)
#  --input_suffix    A file name suffix of input data 
#  --contig_type     The contig type for redundancy removal (aligned or assembled)
#  --contig_prefix   A prefix of each contig name in the data (fasta format)
#
# Megablast-related options(7):
#  --strand_specific Only for sequences assembled from strand specific RNA-seq [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [4]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]  
#########################################################################################################

_EOUSAGE_
	;
	
#################
##   global    ##
#################
# parameters for input files
our $file_list;		# list of input file name for virus detection
our $file_type;		# fastq or fasta
our $input_suffix;	# suffix file name : contigs1.fa
our $contig_type; 	# contig type : aligned还是assembled
our $contig_prefix;	# fasta 文件改名时，每条记录用到的前缀

# parameters for megablast (remove redundancy)
our $strand_specific;	#
our $min_overlap = 30;	# hsp combine
our $max_end_clip = 6; 	# hsp combine
our $cpu_num = 8;	# megablast: thread
our $mis_penalty = -1;	# megablast: penalty for mismatch
our $gap_cost = 2;	# megablast: penalty for gap open
our $gap_extension = 1;	# megablast: penalty for gap extension

################################
# set path for file and folders#
################################
our $WORKING_DIR	= cwd();				# current folder
our $DATABASE_DIR	= ${FindBin::RealBin}."/../databases";	# database folder
our $BIN_DIR		= ${FindBin::RealBin};			# program folder
our $TEMP_DIR 		= $WORKING_DIR."/temp";			# temp folder
my $tf = $TEMP_DIR;						# short name of temp folder

#################
# parameters	#
#################
GetOptions( 
	'file_list=s'		=> \$file_list,
	'file_type=s' 		=> \$file_type,
	'input_suffix=s' 	=> \$input_suffix,
	'contig_type=s' 	=> \$contig_type,
	'contig_prefix=s' 	=> \$contig_prefix,

	'strand_specific!' 	=> \$strand_specific,
	'min_overlap=i' 	=> \$min_overlap,
	'max_end_clip=i' 	=> \$max_end_clip,
	'cpu_num=i' 		=> \$cpu_num,
	'mis_penalty=i' 	=> \$mis_penalty,
	'gap_cost=i' 		=> \$gap_cost,
	'gap_extension=i' 	=> \$gap_extension			 			 
);
			 
die $usage unless ($file_list && $input_suffix && $contig_prefix);	# required parameters

#################
# main 		#
#################
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

	# if file has sequence, move it to temp folder to remove redundancy
	# 1. remove simple repeate sequence using mask
	Util::process_cmd("$BIN_DIR/dust $sample.$input_suffix 1> $sample.masked 2> $tf/dust.log");
	Util::process_cmd("$BIN_DIR/trim_XNseq1.pl $sample.masked $sample.$input_suffix 0.8 40 > $sample.$input_suffix.1");
	Util::process_cmd("rm $sample.masked");

	my ($before_contig_num, $after_contig_num, $i);
	$i = 1;								# get the number of contig file, default is 1
	$before_contig_num = `grep -c \'>\' $sample.$input_suffix.$i`;	# get the seq number before remove redundancy
	$after_contig_num  = 0;						# this is seq number after remove redundancy	
	
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
		Util::print_user_submessage("$after_contig_num of uniq contigs were generated");
	} elsif ($after_contig_num == 1) {
		Util::print_user_submessage("$after_contig_num of uniq contig was generated");
	} elsif ($after_contig_num == 0) {
		Util::print_user_submessage("None of uniq contig was generated");
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
	renameFasta("$sample.contigs$i.fa", "$sample.$input_suffix", $contig_prefix);

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

#system("rm *.log");

#################
# subroutine	#
#################
=head2

 renameFasta: rename the fasta file with spefic prefix

=cut
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

=head
our $max_Xratio = 80/100; 
our $min_BaseNum = 40; 
#将fasta文件中的所有序列满足，[XxNn]碱基比例小于一个ratio，或者非[XxNn]碱基总数大于一个绝对值的序列输出到output
 
if (@ARGV < 2)
{
  print "usage: trm_XNseq1.pl input1 input2 max_Xratio min_BaseNum > output\n";
  exit(0);
}

our $input1 = $ARGV[0];
our $input2 = $ARGV[1];
$max_Xratio = $ARGV[2];
$min_BaseNum = $ARGV[3];

my %seq_required;
my ($pre_name, $seq) = ('', ''); 
open(IN1, "$input1");
while (<IN1>) {
	if (/^>(\S+)/) {
		my $current_name = $1; 
		if ($pre_name ne '') {
			my $is_out = 1; #表示该序列应该输出，这个变量调试用
			$seq =~ s/\s//g; 
			my $xnnum = ($seq =~ tr/XxNn/XxNn/); #得到非4种常见碱基[XxNn]的碱基总数
			my $seqlen = length($seq); 
			if ($xnnum/$seqlen >= $max_Xratio) {#[XxNn]碱基比例大于一个ratio
				$is_out = 0; #表示该序列不应该输出
			}elsif ($seqlen-$xnnum < $min_BaseNum) {#或者[ATCG]碱基总数不够一定数量
				$is_out = 0; #表示该序列不应该输出
			}else{
				#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
				#print ">$pre_name\n$seq\n"; #都不满足的就输出到标准输出
				defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#建立query和query_length之间的映射
			}
			#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; #调试用
		}
		$pre_name=$current_name; $seq = ''; 
	}else{
		$seq .= $_; 
	}
}
#不要忘记处理剩下的
if ($pre_name ne '') {
	my $is_out = 1; 
	$seq =~ s/\s//g; 
	my $xnnum = ($seq =~ tr/XxNn/XxNn/);
	my $seqlen = length($seq);
	if ($xnnum/$seqlen >= $max_Xratio) {
		$is_out = 0;
	}elsif ($seqlen-$xnnum < $min_BaseNum) {
		$is_out = 0;
	}else{
		#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
		#print ">$pre_name\n$seq\n"; #都不满足的就输出到标准输出
		defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#建立query和query_length之间的映射
	}
	#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; 
}
close(IN1);

#从input2中把在%seq_required中的序列提取出来
open(IN2, $input2);

my $flag = "off";
while(<IN2>) {
	if($_ =~ m/^>/) {
		my $head = $_;
		chomp($head);
		$head=~s/>//;

		if(defined $seq_required{$head}) {#如果包括这个name
			print $_;#就输出
			$flag = "on";#同时改变标志，表示后面的序列需要继续向OUT1输出
		}
		else {#如果不包括这个name
			#print OUT $_;#就输出到OUT
			$flag = "off";#同时改变标志，表示后面的序列需要继续向OUT2输出
		}
	}
	else {
		if($flag eq "on") {#表示为"on"
			print $_;#后面的序列需要继续输出
		}
	}
}
close(IN2);
=cut


