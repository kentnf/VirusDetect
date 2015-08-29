#!/usr/bin/env perl

=head

 update information
 
 1. input fasta/fastq file, do not need file-type, contig-type, file-list parameters. 
    add parameters for identify novel virus?
   
 2. delete the reference parameters, because it should have sequence information, do not need changed 
    but user can input their own reference information ? 
    or generated reference form the sequence ?

 3. add debug and force [debug] parameters, 
    force, it will using knwon contig to perform novel identification 

=cut

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use IO::File;
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin";
use Util;

my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_itentify.pl [options] --reference reference input_read contig
#  		
#  		the contig was assembled from input_read
#
#                   --diff_ratio --word_size [INT] --exp_value <Float> --identity_percen <Float>
#                   --cpu_num  [INT] --mis_penalty [INT] --gap_cost[INT] --gap_extension [INT]
#
# Options(3):
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta] 
#  --diff_ratio The hits with distance less than 0.25 will be combined into one  [0.25] 
# 
# blast-related options(7):
#  --word_size [11] 
#  --exp_value [1e-5]
#  --identity_percen [25] 
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]
#
# New options(5):
#  --hsp_cover		The coverage of hsp should be more than this cutoff for query or hit [0.75]
#  --diff_contig_cover	The coverage for different contigs [0.5]
#  --diff_contig_length	The length of different contigs [100 (bp)]
#  --identity_drop_off	The percentage identity drop off compared with the highest identity [5 (%)]
#
#  --coverage_cutoff	The coverage cutoff for final ouput [0.1] 
#  --depth_cutoff	The depth cutoff for final ouput    [5]
#  --depth_norm		Normalized depth by library size [off]
#  --novel_len_cutoff	Force unconsidered contig longer than this cutoff as novel contigs [1000]
###########################################################################################

_EOUSAGE_

;

################################
## set folder and file path   ##
################################

my $WORKING_DIR = cwd();				# current folder : working folder
my $BIN_DIR = ${FindBin::RealBin};			# bin folder
my $result_dir = $WORKING_DIR."/result";		# result folder
my $tf = $WORKING_DIR."/temp";				# temp folder

my $DATABASE_DIR = ${FindBin::RealBin}."/../databases";	# database folder
my $seq_info  = $DATABASE_DIR."/vrl_genbank.info.gz";	# virus sequence info
my $reference; # = $DATABASE_DIR."/vrl_plant";       	# virus sequence
my $prot_tab  = $DATABASE_DIR."/vrl_idmapping.gz";	# virus protein table
Util::process_cmd("$BIN_DIR/formatdb -i $reference -p F") unless (-e "$reference.nhr");

# format of seq_info
# AB000048 \t 2007 \t Parvovirus \t Feline panleukopenia virus gene for nonstructural protein 1, complete cds, isolate: 483. \t 1 \t Vertebrata \n

###############################
##  global vars		     ##
###############################
my $hsp_cover = 0.75;		# for blast filter
my $drop_off = 5;		# for blast filter
my $diff_ratio= 0.25;		# ratio for number of diff contig
my $diff_contig_cover = 0.5;	# for hit filter
my $diff_contig_length = 100;	# for hit filter
my $coverage_cutoff = 0.1;	# coverage cutoff for final result
my $depth_cutoff = 5;		# depth cutoff for final result
my $depth_norm = 0;		# normalization depth by library size
my $novel_len_cutoff = 1000;	# 

my $word_size = 11;
my $cpu_num = 8;		# megablast: thread number
my $mis_penalty = -1;		# megablast: penalty for mismatch
my $gap_cost = 2;		# megablast: penalty for gap open
my $gap_extension = 1;		# megablast: penalty for gap extension
my $exp_value = 1e-5;		#
my $identity_percen = 25;	# tblastx: hsp_identity cutoff for protein blast

my $filter_query = "F";		# megablast: F - disable remove simple sequence
my $hits_return = 500;		# megablast: hit number
my $input_suffix='';

my ($debug, $debug_force, $novel_check);

########################
# get input parameters #
########################
GetOptions( 
	'reference=s' 		=> \$reference,
	'diff_ratio=f' 		=> \$diff_ratio,
	'hsp_cover=f'		=> \$hsp_cover,
	'diff_contig_cover=f'	=> \$diff_contig_cover,
	'diff_contig_length=i'	=> \$diff_contig_length,
	'word_size=i' 		=> \$word_size,
	'exp_value=f' 		=> \$exp_value,
	'identity_percen=f' 	=> \$identity_percen,
	'cpu_num=i' 		=> \$cpu_num,
	'mis_penalty=i' 	=> \$mis_penalty,
	'gap_cost=i' 		=> \$gap_cost,
	'gap_extension=i' 	=> \$gap_extension,
	'coverage_cutoff=f'	=> \$coverage_cutoff,
	'depth_cutoff=f'	=> \$depth_cutoff,
	'depth_norm'		=> \$depth_norm,
	'novel_len_cutoff=i'    => \$novel_len_cutoff,
	'd|debug'		=> \$debug,
	'f|force'		=> \$debug_force,
	'n|novel-check'		=> \$novel_check,
);

if (scalar(@ARGV) != 2) {
	print $usage;
	exit;
}

main: {
	my ($sample, $contig) = @ARGV;
	my $file_type = Util::detect_FileType($sample);
	die "[ERROR]undef input reads: $sample\n$usage" unless -s $sample;
	die "[ERROR]undef input contig: $contig\n$usage" unless -s $contig;

	# create result folder according to sample name, copy contig to result folder
	my $sample_base = basename($sample);
	my $sample_dir = $result_dir."_".$sample_base;
	Util::process_cmd("mkdir $sample_dir", $debug) unless -e $sample_dir;
	Util::process_cmd("cp $contig $sample_dir/contig_sequences.fa", $debug);

	# load contig info to hash; key: seqID, seq, length; value: seq, length
	my %contig_info = Util::load_seq($contig);
	my %reference_info = Util::load_seq($reference);
	my %reference_prot_info = Util::load_seq($reference."_prot");

	# load virus seqinfo to hash
	# key: seqID, Length, type, desc, version, host_type
	my %virus_info = Util::load_virus_info($seq_info, $prot_tab);

	# Part A
	# 1. blast contigs against reference, parse and filter blast results to table format (NOT USE SEARCH::IO)
	#    please check util::filter_blast_table for filter rules ($hsp_cover: 0.75, $drop_off: 5)
	#    hsp_cover/query_len >= 0.75 || hsp_cover/hit_len >= 0.75
	#    for each contig, if the best identity hit to viral seq is 90%, the min identify for other hits viral seq
	#    should be 85% (90% - 5%). that is 5 dropoff works
	#
	#    output blast table format:
	#    query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t
	#    query_start \t query_end \t hit_start \t hit_end \t identity2 \t aligned_query \t aligned_hit \t aligned_string\n
	# 
	# 2. then find known contig according to blast results, using identify and ratio as cutoff
	#    key: contig_id; value: 1
	#    1. get the all the hsps with > 60 identity for one query
	#    2. combine overlapped hsp for each query, then get total coverage from non-overlapped hsp of each query
	#       this total cov means number of base in query contig could aligned to viral seqs
	#    3. get ratio of hsp = total cov / query length * 100, if the ratio > 50%, the contig was as known
	#    4. get know blast result with known contigs
	#
	# 3. exit if known contig identified
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_output  = "$contig.blastn.paired";
	my $blast_param   = "-F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";
	Util::process_cmd("$blast_program -i $contig -d $reference -o $blast_output $blast_param", $debug) unless -s $blast_output;
	my $blast_table  = Util::parse_blast_to_table($blast_output, $blast_program);

	Util::save_file($blast_table, "$sample.blastn.table1");	# checking files

	my $nn1 = Util::line_num($blast_table);
	   $blast_table  = Util::filter_blast_table($blast_table, $hsp_cover, $drop_off, $blast_program);

	Util::save_file($blast_table, "$sample.blastn.table");	# checking files

	my ($known_contig, $known_blast_table) = Util::find_known_contig($blast_table, 60, 50);

	unlink($blast_output) unless $debug;
	if (scalar(keys(%$known_contig)) == 0) { system("touch $sample_dir/no_virus_detected"); exit; } 

	# report information for error checking:
	if ($debug) {
		print "=== For all contigs ===\nBlast Table: $nn1\nFiltered Blast Table: ".Util::line_num($blast_table);
		print "Contig Num:".scalar(keys(%contig_info))."\nKnown Contigs:".scalar(keys(%$known_contig))."\n";
		print "Blast Table for known contig: ".Util::line_num($known_blast_table)."\n";
	}

	# 4. get hit coverage information, then remove redundancy hit
	#    this function looks redundance, it has save filter compared to step 3, but the filtered
	#    result only works on known coverage and block computation
	my ($known_coverage, $known_block) = Util::get_hit_coverage($known_blast_table, 60, 0.5, 1);

	Util::save_file($known_coverage, "$sample.known.cov");  # checking files

	# 5. remove redundancy hit.
	#    after remove redundancy, the similar viral seq will be removed
	#    Diff_ratio : 
	#    Diff_contig_cover : 
	#    Diff_contig_length : 
	#    the output file is known identified
	my $known_identified = Util::remove_redundancy_hit($known_coverage, $known_blast_table, $diff_ratio, $diff_contig_cover, $diff_contig_length);

	Util::save_file($known_identified, "$sample.known.identified");

	# 6. align read to contigs get depth, then compute reference depth using contig depth
	# %contig_depth -> key: contig ID, depth, cover; value: depth, coverage
	# %known_depth  -> key1: ref_ID; 
	# 		   key2; depth, norm
	# 		   value: raw depth, normalized depth
	# filter the depth and coverage by depth cutoff and coverage cutoff, save the file
	my %contig_depth = get_contig_mapped_depth($contig, $sample, $cpu_num, $file_type);
	my %known_depth  = correct_depth($known_identified, \%contig_depth, \%contig_info, $sample, $depth_norm); # input sample will provide lib_size for depth normalization

	# 7 the known identified table was filter again using depth and coverage
	#   define coverage: query contig cov / total viral seq, 
	#   define depth:
	#   if the coverage meet cutoff, and raw & normalized depth meet the cutoff, identified virus will reported (final table)
	my ($final_known_ctgs, $removed_ctgs);
	($known_identified, $final_known_ctgs, $removed_ctgs) = filter_by_coverage_depth($known_identified, \%known_depth, $coverage_cutoff, $depth_cutoff);
	Util::save_file($known_identified, "$sample.known.identified_with_depth");

	my ($known_contig_table, $known_contig_blast_table, $known_reference) =  combine_table1($known_identified, $known_blast_table, \%contig_info, \%virus_info, \%reference_info);
	my $known_contig_blast_sam = Util::blast_table_to_sam($known_contig_blast_table);
	Util::save_file($known_contig_table, "$sample_dir/$sample_base.known.xls");
	Util::save_file($known_reference, "$sample_dir/known.reference.fa");
	Util::save_file($known_contig_blast_sam, "$sample_dir/$sample_base.known.sam");
	
	my $known_num = 0; ($known_num, $known_identified) = arrange_col2($known_identified, \%virus_info);

	if ( length($known_identified) > 1 ) { 
		Util::plot_result($known_identified, $known_contig_blast_table, $sample_dir, 'known'); 
	} else {
		unlink "$sample_dir/$sample_base.known.xls" if -e "$sample_dir/$sample_base.known.xls";
	}
	
	Util::print_user_submessage("$known_num of known virus has been identified");

	# output final-known/half-knwon/novel contigs;
	#    final-known contigs: contigs in final result known.xls
	#    half-known contigs : contigs not in final result, but meet the cutoff of 60% identity to virus, and 50% coverage of itself
	#			  --- it include contig also pass remove redundancy, but failed on coverage and depth filter
	#    novel contigs : need re-define 
	my $known_contig_f = "$sample.known.contigs";
	my $novel_contig = "$sample.novel.contigs";
	my $half_known_contig_f = "$sample.unconsidered.contigs";
	my $fh1 = IO::File->new(">".$known_contig_f) || die $!;
	my $fh2 = IO::File->new(">".$novel_contig) || die $!;
	my $fh3 = IO::File->new(">".$half_known_contig_f) || die $!;
                
	foreach my $cid (sort keys %contig_info) {
		if (defined $$final_known_ctgs{$cid}) {
			print $fh1 ">$cid\n$contig_info{$cid}{'seq'}\n";
		} elsif (defined $$known_contig{$cid} ) {
			if (defined $$removed_ctgs{$cid}) {
				print $fh2 ">$cid\n$contig_info{$cid}{'seq'}\n";
			} elsif ($contig_info{$cid}{'length'} >= $novel_len_cutoff) { 
				print $fh2 ">$cid\n$contig_info{$cid}{'seq'}\n";
			} else {
				print $fh3 ">$cid\n$contig_info{$cid}{'seq'}\n";
			}
		} else {
			print $fh2 ">$cid\n$contig_info{$cid}{'seq'}\n";
		}
	}
	$fh1->close;
	$fh2->close;
	$fh3->close;

	Util::process_cmd("cp $known_contig_f $sample_dir/contig_sequences.known.fa");
	Util::process_cmd("cp $novel_contig   $sample_dir/contig_sequences.novel.fa");
	Util::process_cmd("cp $half_known_contig_f $sample_dir/contig_sequences.unconsidered.fa");

	if ( $novel_check && -s $novel_contig )
	{
		# compare noval contigs against virus database using tblastx
		$blast_output = "$sample.novel.paired";
		$blast_program = $BIN_DIR."/blastall -p blastx";
		$blast_param = "-F $filter_query -a $cpu_num -e $exp_value";				
		my $reference_prot = $reference."_prot"; 
		Util::process_cmd("$blast_program -i $novel_contig -d $reference_prot -o $blast_output $blast_param", $debug) unless -s $blast_output;
		my $blast_novel_table = Util::parse_blast_to_table($blast_output, $blast_program);
		my $mm1 = Util::line_num($blast_novel_table);
		   $blast_novel_table = Util::filter_blast_table($blast_novel_table, 0, $drop_off, $blast_program);
		unlink($blast_output) unless $debug;

		# report information
		print "=== For novel contig ===\n Blast Table: $mm1\nFiltered Blast Table: ".Util::line_num($blast_novel_table)."\n" if $debug;

		# exit if no blast table was generated
		unless (length $blast_novel_table > 1 )
		{
			Util::process_cmd("touch $sample_dir/no_novel_virus_detected", $debug);
			Util::process_cmd("cp $novel_contig $sample_dir/unknown.contigs.fa", $debug);
			exit;
		}
	
		# 4. get coverage remove redundancy 
		my ($novel_coverage, $novel_block) = Util::get_hit_coverage($blast_novel_table, $identity_percen, 0, 1);
		my $novel_identified = Util::remove_redundancy_hit($novel_coverage, $blast_novel_table, $diff_ratio, $diff_contig_cover, $diff_contig_length);		

		# 5. get depth 
		my %novel_depth  = correct_depth($novel_identified, \%contig_depth, \%contig_info, $sample, $depth_norm);
		($novel_identified, $final_known_ctgs, $removed_ctgs) = filter_by_coverage_depth($novel_identified, \%novel_depth, $coverage_cutoff, $depth_cutoff);	
		
		# combine
		my ($novel_contig_table, $novel_contig_blast_table, $novel_reference) =  combine_table1($novel_identified, $blast_novel_table, \%contig_info, \%virus_info, \%reference_prot_info);
		my $novel_contig_blast_sam = Util::blast_table_to_sam($novel_contig_blast_table);
		Util::save_file($novel_contig_table, "$sample_dir/$sample_base.novel.xls");
		Util::save_file($novel_reference, "$sample_dir/novel.reference.fa");
        	Util::save_file($novel_contig_blast_sam, "$sample_dir/$sample_base.novel.sam");

		my $novel_num = 0; ($novel_num, $novel_identified) = arrange_col2($novel_identified, \%virus_info);

		if ( length($novel_identified) > 1 ) { 
			Util::plot_result($novel_identified, $novel_contig_blast_table, $sample_dir, 'novel'); 
		} else {
                	unlink "$sample_dir/$sample_base.novel.xls" if -e "$sample_dir/$sample_base.novel.xls";
                }

                Util::print_user_submessage("$novel_num of novel virus has been identified");
	}
}

# put folder new folder

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 arrange_col1 : discard in the new
 
 just re-order the contig dataset, 
 then sort by some cols
=cut
sub arrange_col1 
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
		# line format
		# 0 - Contig_ID
		# 1 - Contig_length
		# 2 - Hit_ID
		# 3 - Hit_length
		# 4 - strand
		# 5 - Contig_start
		# 6 - Contig_end
		# 7 - Hit_start 
		# 8 - Hit_end
		# 9 - identity
		#10 - genus
		#11 - description
		#12 - contig_sequence
		print "STOPPoint";
		print $_."\n"; exit;
		my @a = split(/\t/, $_);
		push(@all_data, [@a[0,19,1,2,3,17,18,9,10,11,12,13,6,8]]);	# re-order the data, then put them to array 
	}
	$in->close;
	
	@all_data = sort { ($a->[5] cmp $b->[5]) || ($a->[3] cmp $b->[3])} @all_data; # sort according to Genus and hit

	my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";	
	print $out "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";
	foreach my $data ( @all_data ) {
		print $out join("\t", @$data)."\n";
	}
	$out->close;
}
=head2

 arrange_col2 : arrange col of file identified1 to identified
 input:  sample.known.identified1 or sample.novel.identified1
 output: sample.known.identified  or sample.novel.identified  

 Filter the result using coverage and depth setted by input parameter
 sort result by genus, then by read_cov(bp) 

=cut
sub arrange_col2
{
	my ($known_identified, $virus_info, $sample, $coverage_norm) = @_;

	my @all_data;

	chomp($known_identified);
	my @a = split(/\n/, $known_identified);

	foreach my $line (@a)
	{
		# 0 [col1] - virus seq ID
		# 1 [col2] - virus seq length
		# 2 [col3] - virus seq cover length
		# 3 [col4] - coverage (cover_length/seq_length: col3/col2)
		# 4 [col5] - contigs
		# 5 [col6] - contig number
		# 6 [col7] - raw depth (read depth, from pileup file)
		# 7 [col8] - normalized depth

		my @ta = split(/\t/, $line);
		#my $coverage= 1.0 * $a[7] / $a[1];					# get coverage, shan's method
		my $coverage = $ta[2] / $ta[1];						# changed by kentnf
		my $depth = $ta[6];
		my $norm_depth = $ta[7];
		my $genus = $$virus_info{$ta[0]}{'genus'};				# col 9
		my $desc  = $$virus_info{$ta[0]}{'desc'};				# col 10
		push(@all_data, [@ta[0,1,2],$coverage,@ta[4,5,6,7],$genus,$desc]); 	# re-order the data, then put them to array 
	}
	
	@all_data = sort { ($a->[8] cmp $b->[8]) || ($b->[2] cmp $a->[2])} @all_data; # sort according to Genus and read_cov(bp)

	my $output_identified = '';

	foreach my $data ( @all_data ) {
		$output_identified.= join("\t", @$data)."\n";
	}

	my $num = scalar(@all_data);
	return ($num, $output_identified);
}


=head2

 get_contig_mapped_depth -- align reads to contig, then get mean mapped depth, total mapped length, coverage

=cut
sub get_contig_mapped_depth
{
	my ($contig, $sample, $cpu_num, $file_type) = @_;


	# using bw
	my $sai = $sample."_bwa.sai";
	my $log = $sample."_bwa.log";
	my $parameters = "-n 1 -o 1 -e 1 -i 0 -l 15 -k 1 -t $cpu_num";
	my $bwa_mhit_param = "-n 10000";
	Util::process_cmd("$BIN_DIR/bwa index -p $contig -a bwtsw $contig 1> $log 2>> $log", $debug);
	Util::process_cmd("$BIN_DIR/bwa aln $parameters $contig $sample 1> $sai 2>> $log", $debug);
	Util::process_cmd("$BIN_DIR/bwa samse $bwa_mhit_param $contig $sai $sample 1> $sample.sam 2>> $log", $debug);
	Util::xa2multi("$sample.sam");
	
	# using bowtie
	# my $format = ''; if( $file_type eq "fasta" ){ $format = "-f" };
	# Util::process_cmd("$BIN_DIR/bowtie-build --quiet $contig $contig", $debug);
	# Util::process_cmd("$BIN_DIR/bowtie --quiet $contig -v 1 -p $cpu_num -a --best --strata $format $sample -S --sam-nohead $sample.sam", $debug);

	Util::process_cmd("$BIN_DIR/samtools faidx $contig");
	Util::process_cmd("$BIN_DIR/samtools view -bt $contig.fai $sample.sam > $sample.bam 2>$sample.samtools.log");
	Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted 2>$sample.samtools.log");
	Util::process_cmd("$BIN_DIR/samtools mpileup -f $contig $sample.sorted.bam > $sample.pileup 2>$sample.samtools.log"); 

	# parse pileup file to save depth info to hash
	# key: ref_id, mean, total. cover
	# value: mean, total, cover
	# cover: number base covered for each contig
	# total: total depth for each contig
	# mean: mean depth for each contig (total/cover)
	my %depth = Util::pileup_depth("$sample.pileup");

	# unlink temp file
	unlink("$sample.pileup", "$sample.sorted.bam", "$sample.bam", "$contig.fai", "$sample.sam") unless $debug;

	return %depth;
}

=head2

 correct_depth -- using contig depth to compute reference depth

=cut
sub correct_depth
{	
	my ($known_identified, $contig_depth, $contig_info, $sample, $ref_depth_norm) = @_;

	# get total read num
	my $lib_size = 0;
	my $fh = IO::File->new($sample) || die $!;
	while(<$fh>) {
		my $id = $_;
		if ($id =~ m/^>/) { <$fh>; }
		if ($id =~ m/^@/) { <$fh>; <$fh>; <$fh>; }
		$lib_size++;
        }
	$fh->close;

	# convert contig depth to reference depth
	my %known_depth;
	my %known_depth_norm;

	chomp($known_identified);
	my @l = split(/\n/, $known_identified);
	foreach my $l (@l)
	{
		my @a = split(/\t/, $l);
		my @contigs = split(/,/, $a[4]);

		my $contig_len = 0;
		my $total_len = 0;

		foreach my $cid (@contigs)
		{
			die "Error, undef contig length $cid\n" unless defined $$contig_info{$cid}{'length'};
			die "Error, undef contig total  $cid\n" unless defined $$contig_depth{$cid}{'total'};
			$contig_len = $contig_len + $$contig_info{$cid}{'length'};	# combine all contig length for each reference
			$total_len  = $total_len  + $$contig_depth{$cid}{'total'};	# conbine all contig depth for each reference
		}

		die "Error in contig length: $contig_len\n" unless $contig_len > 0;
		die "Error in total length: $total_len\n" unless $total_len > 0;

		# normalization depth using RPM method
		# depth = (all_contig_depth / all_contig_length)
		# norm_depth = (all_contig_depth * 1e+6) / (all_contig_length * lib_size)
		my $ref_depth = $total_len / $contig_len;
		my $norm_depth = (1e+6 * $total_len) / ($contig_len * $lib_size);

		# select depth by user provided parameters	
		$known_depth{$a[0]}{'norm'} = $norm_depth;
		$known_depth{$a[0]}{'depth'} = $ref_depth;
	}
	
	return %known_depth;
}

=head2

 filter_by_coverage_depth -- using the coverage and depth to filter ref

=cut
sub filter_by_coverage_depth
{
	my ($known_identified, $known_depth, $coverage_cutoff, $depth_cutoff) = @_;

	my $output_identified = '';	# output identified content
	my %known_contigs;		# known contigs after filter
	my %removed_contigs;		# removed contigs for novel

	chomp($known_identified);
        my @a = split(/\n/, $known_identified);
        foreach my $line (@a)
        {
                my @ta = split(/\t/, $line);
		my @tb = split(/,/, $ta[4]);
		die "[ERR]no support contig: $line\n" if scalar @tb < 1;
	
                if ($ta[3] > $coverage_cutoff) {
			die "[ERR]Undef depth for $ta[0]\n" unless defined $$known_depth{$ta[0]}{'norm'};
                        if ( $$known_depth{$ta[0]}{'norm'} > $depth_cutoff || $$known_depth{$ta[0]}{'depth'} > $depth_cutoff ) {
                               $output_identified.=$line."\t".$$known_depth{$ta[0]}{'depth'}."\t".$$known_depth{$ta[0]}{'norm'}."\n";
				foreach my $ctg (@tb) { $known_contigs{$ctg} = 1; }
                        } else {
				foreach my $ctg (@tb) { $removed_contigs{$ctg} = 1; }
			}
                } else {
			foreach my $ctg (@tb) { $removed_contigs{$ctg} = 1; }
		}
        }

	my $org_removed_num = scalar(keys(%removed_contigs));
	my $org_known_num   = scalar(keys(%known_contigs));
	my $ccn = 0;	# correct contig num
	# check if removed contigs has knwon one, then delete it from removed
	foreach my $ctg (sort keys %removed_contigs) {
		delete $removed_contigs{$ctg} and $ccn++ if defined $known_contigs{$ctg};
			
	}
	#print "removed:$org_removed_num\tknown:$org_known_num\tcorrect:$ccn\n";
	return ($output_identified, \%known_contigs, \%removed_contigs);
}
=head2 

 combine_table1 -- combine the result to generate final table

 source: query_filter2

# annotation of soure query_filter2 (by kentnf)
# 1. the source query_fitler2 is not modified
# 2. put all the identified ref ID to hash %hit_index
# 3. find the best result of each query according to order, the 1st one should be best one, output it 
# 4. if the best one is identified in %hit_index, find the other best with same identify and evalue, output it
# 5. if the best one is not identified in %hit_index, find next best identified in %hit_index, output it, then 
#    find the other hit meet with identified in %hit_index, with same or high evalue or identified
    
# annotation of modified (by kentnf)
# 1. find the best result of each query according to order, the 1st one should be best one, output it
# 2. output all the query if it's hit identified in %hit_index

# 2nd modification 
# 1. only output all the query if it's hit identified in %hit_index

my ($known_contig_table, $knwon_contig_blast_table, $known_reference) =  combine_table1($known_identified, $known_blast_table, \%contig_info, \%virus_info, \%reference_info);


=cut 
sub combine_table1
{
	my ($known_identified, $known_blast_table, $contig_info, $virus_info, $reference_info) = @_;
	
	# put the refID to hash
	# key: refID (The reference is virus reference sequence)
	# value: 1
	my $reference_seq = '';
	my %hit_index;
	chomp($known_identified);
	my @a = split(/\n/, $known_identified);

	foreach my $line (@a)
        {
		next if $line =~ m/^#/;
                # HQ593108        7458    96      0.0113971574148565      CONTIG282       1
                my @ta = split(/\t/, $line);
		die "[ERR]Undef Ref Seq: $ta[0]\n" unless defined $$reference_info{$ta[0]}{'seq'};
		$reference_seq.= ">$ta[0]\n$$reference_info{$ta[0]}{'seq'}\n";
                $hit_index{$ta[0]} = 1;
        }

        # main
        my $output_table = "#Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";
	my $output_blast_table = '';
        chomp($known_blast_table);
        my @b = split(/\n/, $known_blast_table);
        foreach my $line (@b)
        {
		next if $line =~ m/^#/;
		# query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t
		# query_start \t query_end \t hit_start \t hit_end \t identity2 \t aligned_query \t aligned_hit \t aligned_string\n
		my @ta = split(/\t/, $line);

		if ( defined $hit_index{$ta[2]} ) 
		{
			die "[ERR]Undef seq for $ta[0]\n"   unless defined $$contig_info{$ta[0]}{'seq'};
			die "[ERR]Undef Genus for $ta[2]\n" unless defined $$virus_info{$ta[2]}{'genus'};
			die "[ERR]Undef Desc for $ta[2]\n"  unless defined $$virus_info{$ta[2]}{'desc'};
		
			$output_table.="$ta[0]\t$$contig_info{$ta[0]}{'seq'}\t$ta[1]\t$ta[2]\t$ta[3]\t$$virus_info{$ta[2]}{'genus'}\t$$virus_info{$ta[2]}{'desc'}\t";
			$output_table.="$ta[9]\t$ta[10]\t$ta[11]\t$ta[12]\t$ta[13]\t$ta[6]\t$ta[8]\n";
			$output_blast_table.=$line."\n";
		}
	}

	return ($output_table, $output_blast_table, $reference_seq);
}


