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
use FindBin;
use lib "$FindBin::RealBin/PerlLib";
use Util;
use File::Basename;

my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_itentify.pl [options] input_read contig
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
#  --identity_percen [80] 
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
#
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
my $seq_info = $DATABASE_DIR."/vrl_genbank.info";	# virus sequence info
my $reference = $DATABASE_DIR."/vrl_plant";       	# virus sequence
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
	#'reference=s' 		=> \$reference,
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
	my $sample_dir = $result_dir."_".$sample_base."_new";
	Util::process_cmd("mkdir $sample_dir", $debug) unless -e $sample_dir;
	Util::process_cmd("cp $contig $sample_dir/contig_sequences.fa", $debug);

	# load contig info to hash; key: seqID, seq, length; value: seq, length
	my %contig_info = Util::load_seq($contig);
	my %reference_info = Util::load_seq($reference);

	# load virus seqinfo to hash
	# key: seqID, Length, type, desc, version, host_type
	my %virus_info = Util::load_virus_info($seq_info);

	# Part A
	# 1. blast contigs against reference, parse and filter blast results to table format (NOT USE SEARCH::IO)
	#    please check util::filter_blast_table for filter rules
	#
	#    output blast table format:
	#    query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t
	#    query_start \t query_end \t hit_start \t hit_end \t identity2 \t aligned_query \t aligned_hit \t aligned_string\n
	# 
	# 2. then find known contig according to blast results
	#    key: contig_id; value: 1
	#    1. get the all the hsps with > 60 identity for one query
	#    2. combine all the hsps for one query
	#    3. get ratio of combined hsp and query
	#    4. get know if ratio, other is novel
	#
	# 3. exit if known contig identified
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_output  = "$contig.blastn.paired";
	my $blast_param   = "-F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";
	Util::process_cmd("$blast_program -i $contig -d $reference -o $blast_output $blast_param", $debug) unless -s $blast_output;
	my $blast_table  = Util::parse_blast_to_table($blast_output, $blast_program);
	my $nn1 = Util::line_num($blast_table);
	   $blast_table  = Util::filter_blast_table($blast_table, $hsp_cover, $drop_off, $blast_program);
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
	my ($known_coverage, $known_block) = Util::get_hit_coverage($known_blast_table, 60, 0.5, 1);
	my $known_identified = Util::remove_redundancy_hit($known_coverage, $known_blast_table, $diff_ratio, $diff_contig_cover, $diff_contig_length);

	# 5. aligne read to contigs get depth, then compute reference depth using contig depth
	# key: contig ID, depth, coverage; value: depth, coverage
	my %contig_depth = get_contig_mapped_depth($contig, $sample, $cpu_num, $file_type);
	my %known_depth  = correct_depth($known_identified, \%contig_depth, \%contig_info);

	$known_identified = filter_by_coverage_depth($known_identified, \%known_depth, $coverage_cutoff, $depth_cutoff);

	my ($known_contig_table, $known_contig_blast_table, $known_reference) =  combine_table1($known_identified, $known_blast_table, \%contig_info, \%virus_info, \%reference_info);
	my $known_contig_blast_sam = Util::blast_table_to_sam($known_contig_blast_table);
	Util::save_file($known_contig_table, "$sample_dir/$sample_base.new.known.xls");
	Util::save_file($known_reference, "$sample_dir/new.known.reference.fa");
	Util::save_file($known_contig_blast_sam, "$sample_dir/$sample_base.new.known.sam");
	
	my $known_num = 0; ($known_num, $known_identified) = arrange_col2($known_identified, \%virus_info);
	
	if ( length($known_identified) > 1 ) { 
		Util::plot_result($known_identified, $known_contig_blast_table, $sample_dir, 'known'); 
	} else {
		unlink "$sample_dir/$sample_base.known.xls" if -e "$sample_dir/$sample_base.known.xls";
	}
	
	Util::print_user_submessage("$known_num of known virus has been identified");

	if ( $novel_check )
	{
		# compare noval contigs against virus database using tblastx
		my $novel_contig = "$sample.novel.contigs";
		my $fh = IO::File->new(">".$novel_contig) || die $!;
		foreach my $cid (sort keys %contig_info) {
			if ( $debug_force ) {
				print $fh ">$cid\n$contig_info{$cid}{'seq'}\n" if defined $$known_contig{$cid};
			} else {
				print $fh ">$cid\n$contig_info{$cid}{'seq'}\n" unless defined $$known_contig{$cid};
			}
		}
		$fh->close;
	
		$blast_output = "$sample.novel.paired";
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-F $filter_query -a $cpu_num -e $exp_value";				
		Util::process_cmd("$blast_program -i $novel_contig -d $reference -o $blast_output $blast_param", $debug) unless -s $blast_output;
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
			next;
		}
	
		# 4. get coverage remove redundancy 
		my ($novel_coverage, $novel_block) = Util::get_hit_coverage($blast_novel_table, $identity_percen, 0, 1);
		my $novel_identified = Util::remove_redundancy_hit($novel_coverage, $blast_novel_table, $diff_ratio, $diff_contig_cover, $diff_contig_length);		

		# 5. get depth 
		my %novel_depth  = correct_depth($novel_identified, \%contig_depth, \%contig_info);
		$novel_identified = filter_by_coverage_depth($novel_identified, \%novel_depth, $coverage_cutoff, $depth_cutoff);	
		
		# combine
		my ($novel_contig_table, $novel_contig_blast_table, $novel_reference) =  combine_table1($novel_identified, $blast_novel_table, \%contig_info, \%virus_info, \%reference_info);
		my $novel_contig_blast_sam = Util::blast_table_to_sam($novel_contig_blast_table);
		Util::save_file($novel_contig_table, "$sample_dir/$sample_base.new.novel.xls");
		Util::save_file($novel_reference, "$sample_dir/new.novel.reference.fa");
        	Util::save_file($novel_contig_blast_sam, "$sample_dir/$sample_base.new.novel.sam");

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
 arrange_col1 : 

 # input format
 # Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

 just re-order the contig dataset, 
 the sort by some cols

=cut
sub arrange_col1 
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
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
 sort result by some col 

=cut
sub arrange_col2
{
	my ($known_identified, $virus_info) = @_;

	my @all_data;

	chomp($known_identified);
	my @a = split(/\n/, $known_identified);
	foreach my $line (@a)
	{
		my @ta = split(/\t/, $line);
		#my $coverage= 1.0 * $a[7] / $a[1];			# get coverage, shan's method
		my $coverage = $ta[2] / $ta[1];				# changed by kentnf
		my $depth = $ta[6];
		my $genus = $$virus_info{$ta[0]}{'genus'};
		my $desc  = $$virus_info{$ta[0]}{'desc'};
		push(@all_data, [@ta[0,1,2],$coverage,@ta[4,5,6],$genus,$desc]);  # re-order the data, then put them to array 
	}
	
	@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] cmp $a->[2])} @all_data; # sort according to Genus and hit_covered(bp)

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

	my $format = ''; if( $file_type eq "fasta" ){ $format = "-f" };
	Util::process_cmd("$BIN_DIR/bowtie-build --quiet $contig $contig", $debug);
	Util::process_cmd("$BIN_DIR/bowtie --quiet $contig -v 1 -p $cpu_num -a --best --strata $format $sample -S --sam-nohead $sample.sam", $debug);
	Util::process_cmd("$BIN_DIR/samtools faidx $contig");
	Util::process_cmd("$BIN_DIR/samtools view -bt $contig.fai $sample.sam > $sample.bam");
	Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	Util::process_cmd("$BIN_DIR/samtools mpileup -f $contig $sample.sorted.bam > $sample.pileup"); 
	my %depth = Util::pileup_depth("$sample.pileup");

	# unlink temp file
	unlink("$sample.pileup", "$sample.sorted.bam", "$sample.bam", "$contig.fai", "$sample.sam") unless $debug;
	unlink("$contig.1.ebwt", "$contig.2.ebwt", "$contig.3.ebwt", "$contig.4.ebwt", "$contig.rev.1.ebwt", "$contig.rev.2.ebwt") unless $debug;

	return %depth;
}

=head2

 correct_depth -- using contig depth to compute reference depth

=cut
sub correct_depth
{	
	my ($known_identified, $contig_depth, $contig_info) = @_;

	my %known_depth;

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
			$contig_len = $contig_len + $$contig_info{$cid}{'length'};
			$total_len  = $total_len  + $$contig_depth{$cid}{'total'};
		}

		die "Error in contig length: $contig_len\n" unless $contig_len > 0;
		die "Error in total length: $total_len\n" unless $total_len > 0;

		my $ref_depth = $total_len/$contig_len;
		$known_depth{$a[0]} = $ref_depth;
	}
	
	return %known_depth;
}

=head2

 filter_by_coverage_depth -- using the coverage and depth to filter ref

=cut
sub filter_by_coverage_depth
{
	my ($known_identified, $known_depth, $coverage_cutoff, $depth_cutoff) = @_;

	my $output_identified = '';	

	chomp($known_identified);
        my @a = split(/\n/, $known_identified);
        foreach my $line (@a)
        {
                my @ta = split(/\t/, $line);
	
                if ($ta[3] > $coverage_cutoff) {
			die "[ERR]Undef depth for $ta[0]\n" unless defined $$known_depth{$ta[0]};
                        if ( $$known_depth{$ta[0]} > $depth_cutoff) 
			{
                               $output_identified.=$line."\t".$$known_depth{$ta[0]}."\n";
                        }
                }
        }

	return $output_identified;
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


