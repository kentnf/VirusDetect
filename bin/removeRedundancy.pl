#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::RealBin";
use Util;
use Bio::SeqIO;
use Cwd;
	
################################
# set path for file and folders#
################################
my $bin_dir	= ${FindBin::RealBin};		# program folder

# check input parameters
if (scalar(@ARGV) != 2) {
	print "USAGE $0 contig_file read_file\n";
	exit;
}
my $contig_file = $ARGV[0];
my $read_file = $ARGV[1];

my $debug = 1;

# parameters for megablast (remove redundancy)
my $min_overlap = 30;   # hsp combine
my $max_end_clip = 6;   # hsp combine

my $strand_specific = 1;
my $cpu_num = 8;        # megablast: thread
my $mis_penalty = -1;   # megablast: penalty for mismatch
my $gap_cost = 2;       # megablast: penalty for gap open
my $gap_extension = 1;  # megablast: penalty for gap extension

my $parameters = "-F F -W 10 -b 10 -a $cpu_num -q $mis_penalty -G $gap_cost -E $gap_extension";
if ($strand_specific) { $parameters .= " -S 1"; }

remove_redundancy($contig_file, $read_file, $parameters, $max_end_clip, $min_overlap, $bin_dir, $debug);

=head2
 
 remove_redundancy -- remove redundancy contig 

=cut
sub remove_redundancy
{
	my ($contig_file, $read_file, $parameters, $max_end_clip, $min_overlap, $bin_dir, $debug) = @_;

	my $min_identity = '';	# pass this value to other subroutine for furture functions

	# step 1. remove low complexity sequence using dust, 2 remov XN reads
	dust_parse($contig_file, $bin_dir);
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
	base_correction($read_file, $contig_file);
}


=head2 
 base_correction -- correct contig base using reads aligned to contigs
 # aligment -> sam -> bam -> sorted bam -> pileup -> Consensus sequence
=cut

sub base_correction
{
	my ($read_file, $contig_file, $bin_dir, $debug);


	my $cpu_num = 8;

	#aligment -> sam -> bam -> sorted bam -> pileup
        Util::process_cmd("$bin_dir/bowtie-build --quiet -f $contig_file $contig_file") unless (-e "$sample.1.ebwt");
        #Util::process_cmd("$BIN_DIR/samtools faidx $sample_reference 2> $tf/samtools.log") unless (-e "$sample_reference.fai");
        Util::process_cmd("$bin_dir/bowtie --quiet $sample -v 1 -p $cpu_num $format $sample -S -a --best $sample.sam", $debug);
        Util::process_cmd("$bin_dir/samtools view -bS $read_file.sam > $read_file.bam 2> $tf/samtools.log");
        Util::process_cmd("$bin_dir/samtools sort $read_file.bam $read_file.sorted 2> $tf/samtools.log");
        Util::process_cmd("$bin_dir/samtools mpileup -f $contig_file $read_file.sorted.bam > $read_file.pileup 2> $tf/samtools.log");

        $file_size = -s "$sample.pileup";               # get file size
        if( $file_size == 0 ){                          # if file size = 0, create blank file, and exit the loop
                Util::process_cmd("touch $sample.$input_suffix");
                next;
        }

        Util::process_cmd("java -cp $BIN_DIR extractConsensus $sample 1 40 $i");
        renameFasta("$sample.contigs$i.fa", "$sample.$input_suffix", $contig_prefix);

	# remove temp files
        system("rm $read_file.sam");
        system("rm $read_file.bam");
        system("rm $read_file.sorted.bam");
        system("rm $read_file.pileup"); 	# must delete this file for next cycle remove redundancy
        system("rm $contig_file");
        system("rm $contig_file.fai");
        system("rm $temp_dir/*.ebwt");
        system("rm $sample.contigs$i.fa");
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
	my ($input_contig, $bin_dir) = @_;
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
	my $out1 = IO::File->new(">".$contig_file) || die $!; 
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

sub process_cmd 
{
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
