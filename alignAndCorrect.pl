#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::RealBin/PerlLib";
use Util;
use File::Basename;
use IO::File;
use Bio::SeqIO;

my $usage = <<_EOUSAGE_;

#########################################################################################
# alignAndCorrect.pl --file_list <FILE> --reference <FILE> --coverage <Float>
#                 		 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#
# Required(3):
#  --file_list A txt file containing a list of input file names without any suffix
#  --reference A fasta file containing all the reference sequences
#  --coverage  A reference sequence covered by reads will be output as results
#
# BWA-related options(6): #maybe try to figure out how to get the parallels for these options in bowtie2
#  --max_dist      Maximum edit distance [1]
#  --max_open      Maximum number of gap opens [1]
#  --max_extension Maximum number of gap extensions [1]
#  --len_seed      Take the first INT subsequence as seed [15]
#  --dist_seed     Maximum edit distance in the seed [1]
#  --thread_num    Number of threads (multi-threading mode) [8]
###########################################################################################

_EOUSAGE_

;
#################
#  global vars  #
#################
our $file_list;		# list of sample name
our $reference;		# reference sequence, fasta format, full path
our $coverage;  	# cutoff of mapping reads / reference
our $max_dist = 1; 	# BWA allowed max distance
our $max_open = 1; 	# BWA allowed max gap opens
our $max_extension = 1; # BWA allowed max gap extension
our $len_seed = 15; 	# length of seed
our $dist_seed = 1; 	# BWA allowed seed max distance
our $thread_num = 8; 	# thread number
our $output_suffix;	# output suffix

################################
# set folder and file path     #
################################
our $WORKING_DIR	= cwd()."/";				# current folder
our $TEMP_DIR		= $WORKING_DIR."temp";			# temp foldr
our $DATABASE_DIR	= $WORKING_DIR."databases";	# database folder
our $BIN_DIR		= ${FindBin::RealBin};			# programs folder
our $seq_info		= $DATABASE_DIR."/vrl_genbank.info";	# genbank info

my $tf = $TEMP_DIR; # short name of temp folder for easy to do

###################
# get input paras #
###################
GetOptions(
'file_list=s'		=> \$file_list,
'reference=s' 		=> \$reference,
'coverage=f' 		=> \$coverage,
'max_dist=i' 		=> \$max_dist,
'max_open=i' 		=> \$max_open,
'max_extension=i' 	=> \$max_extension,
'len_seed=i' 		=> \$len_seed,
'dist_seed=i' 		=> \$dist_seed,
'thread_num=i' 		=> \$thread_num,
'output_suffix=s'	=> \$output_suffix
);

die $usage unless ($file_list && $reference && $coverage && $output_suffix);	# required parameters

#################
# main          #
#################
sub filter_SAM
{
	my $input_SAM = shift;
	my $temp_SAM = $input_SAM.".temp";
    my ($query_col, $opt_col) = (0, 11); 	# query and option column number for sam
	my $max_distance = 2;			# set $max_distance for all selected hits
	my $bestEditDist = -1;			# set best edit distance
	my @alignment = ();			# alignment to array
    my $pre_query_name="";
	my ($total_count, $filtered_count, $kept_align) = (0, 0, 0);
    
	my $in  = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$temp_SAM) || die $!;
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^@/) { print $out $_."\n"; next; }
		my @a = split(/\t/, $_);
        my $query_name = $a[$query_col];
		if ($query_name ne $pre_query_name)
        {   # parse the pre results
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
		if ( $a[1] == 4 ) { $filtered_count++; }
		else {	print $out $_."\n"; }
		$total_count++;
	}
	# parse final query recoed
	if (scalar(@alignment) > 0)
	{
 		foreach my $align (@alignment)
		{
			my $editDistance1;
			if ($align =~ m/\tNM:i:(\d+)/) { $editDistance1 = $1; }
			else { die "Error, this alignment info do not have edit distance : $align\n"; }
			if ($editDistance1 == $bestEditDist) { print $out $align."\n"; $kept_align++; }
		}
	}
    
    $in->close;
	my $filtered_count = $total_count - $kept_align;
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as 2ndhits reads, only for BWA\n";
	$out->close;
	Util::process_cmd("mv $temp_SAM $input_SAM");
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as unmapped reads, only for BWA\n";
}
main: {
    
    # create bwa and fasta index file for reference file (plant virus)
    my $ref;
    my $file_size;
    $ref = $DATABASE_DIR."/".(split(/\//, $reference))[-1];
    Util::process_cmd("$BIN_DIR/bowtie2-build -f $ref --quiet $reference", 1) unless (-s $ref.".1.bt2");
    #print $reference;
    Util::process_cmd("$BIN_DIR/samtools faidx $ref");
    my $ref_dir=$WORKING_DIR.$ref;
    # parse samples in file list
    
    my ($i, $sample, $output_file);
    $i=0;
    open(IN, "$file_list") || die $!;
    while (<IN>) {
		$i=$i+1;
		chomp;
		$sample = $_;
		$output_file = $sample.".".$output_suffix;
        #$output_file = $sample.".sam";
        #print "OUTPUTFILE".$output_file;
		die "Error, file $sample does not exist\n" unless -s $sample;
		#print "# processing the $i sample -- $sample using $0\n";
		#aligment -> sam -> bam -> sorted bam -> pileup
		Util::process_cmd("$BIN_DIR/bowtie2 --quiet -N $max_dist -p $thread_num -L $len_seed --sensitive -q -x $DATABASE_DIR"."/vrl_plant -U $sample -S $output_file", 1);
		# filter out unmapped reads
		# filter out 2nd hits, and only keep the best hits of reads alignment to reference
		filter_SAM($output_file);
        
		# sort sam to bam and generate pileup file
        Util::process_cmd("$BIN_DIR/samtools view -bt $reference.fai $output_file > $sample.bam 2> $tf/samtools.log");
		#Util::process_cmd("$BIN_DIR/samtools view -bt $reference.fai $output_file > $sample.bam 2> $tf/samtools.log") unless (-s "$sample.bam");
		Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted 2> $tf/samtools.log") unless (-s "$sample.sorted");
		Util::process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pre.pileup 2> $tf/samtools.log") unless (-s "$sample.pre.pileup");
        
		# if coverage of any reference is lower than ratio cutoff (default : 0.3 ), remove the alignment for this reference
		my $depth_cutoff;
		$depth_cutoff = pileup_filter("$sample.pre.pileup", "$seq_info", "$coverage", "$sample.pileup");
        Util::process_cmd("mv $sample.pre.pileup $sample.pileup", 1);
		# *** notice, could rewrite this function ***
		# get the size of pileup file
        $file_size = -s "$sample.pileup";
        print "FILE SIZE".$file_size;
		if($file_size != 0)
		{
			# get continous fragment for depth>=1 and length>=40, output file name has '1'
			Util::process_cmd("java -cp $BIN_DIR extractConsensus $sample 0 40 1");
			my $result_contigs = $sample.".contigs1.fa";			# this is a corrected sequences
			my $count = `grep \'>\' $result_contigs | wc -l`; chomp($count);# get seq number
			Util::print_user_submessage("$count of contigs were retrieved from aligned reads");
            print "COUNT".$count;
			if($count!=0){ 							# if seq number > 0, move seq file to align folder
				#system("$BIN_DIR/renameFasta.pl --inputfile $result_contigs --outputfile $sample.contigs2.fa --prefix ALIGNED");# the seq name should be formatted
				renameFasta($result_contigs, "$sample.contigs2.fa", "ALIGNED");
				Util::process_cmd("mv $sample.contigs2.fa $output_file");# move the seq file to folder
                print "OUTPUT FILE".$output_file;
				Util::process_cmd("rm $result_contigs");
			}
			else
			{
				Util::process_cmd("mv $result_contigs $output_file");	# if seq number = 0, move file to align folder
			}
		}
		else{
			Util::print_user_submessage("None of contigs were retrieved from aligned reads");# file size is 0 or sample.pileup is not exist
			Util::process_cmd("touch $output_file");
		}
		
		#system("rm $sample.sai");
		#system("rm $sample.sam");
		#system("rm $sample.bam");
		#system("rm $sample.sorted.bam");
		#system("rm $sample.pre.pileup");
		#system("rm $sample.pileup");
		#system("rm $tf/bowtie.log");
	}
	close(IN);
}

#################
# subroutine    #
#################

sub pileup_filter
{
	my ($input_pileup, $virus_seq_info, $coverage, $output_pileup) = @_;
    
	# put plant virus sequence length to hash
	my %seq_len;
	my $fh1 = IO::File->new($virus_seq_info) || die $!;
	while(<$fh1>) {
		chomp;
		next if $_ =~ m/^#/;
            # ID        Len               Desc
            # AB000048  2007  Parvovirus  Feline panleukopenia virus gene ......  1  Vertebrata
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
		# ref             pos     base    depth   match   Qual?
		# AB000282        216     T       1       ^!,     ~
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



