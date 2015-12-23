#!/usr/bin/perl

=head1
 sRNAtool -- tools for sRNA data preparation
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:g:i:j:k:l:m:n:o:p:q:r:s:t:v:w:x:y:z:fuh', \%options);
rmadp(\%options, \@ARGV);

#########################
# kentnf: subroutine	#
#########################

=head2
 rmadp -- remove adapter 
=cut
sub rmadp
{
	my ($options, $files) = @_;
	
	my $subUsage = qq'
USAGE: $0 -s adapter_sequence -l minimum_length sRNA1 sRNA2 ... sRNAn

';

	# checking parameters
	print $subUsage and exit unless defined $$files[0];
	foreach my $f ( @$files ) { print "[ERR]no file $f\n" and exit unless -s $$files[0]; }

	my ($adp_3p, $adp_3p_len, $adp_3p_sub, $adp_5p, $adp_5p_len, $adp_5p_sub, $adp_5p_len_long, $adp_5p_sub_long,);
	$adp_3p_len = 11;	# default 3p adp length for trimming
	$adp_5p_len = 5;	# default 5p adp length for trimming, match and 5'end match and trim
	$adp_5p_len_long = 9;	# default 5p adp length for trimming, match and trim
	if (defined $$options{'s'} || defined $$options{'p'}) { } else { print $subUsage and exit; } 
	if (defined $$options{'s'}) {
		$adp_3p = $$options{'s'};
		print "[ERR]short adapter 3p\n" and exit if length($adp_3p) < 9;
		#$adp_3p_len = $$options{'l'} if (defined $$options{'l'} && $$options{'l'} > 5);
		$adp_3p_sub = substr($adp_3p, 0, $adp_3p_len);
	} 
	if (defined $$options{'p'} ) {
		$adp_5p = $$options{'p'};
		print "[ERR]short adapter 5p\n" and exit if length($adp_5p) < 9;
		$adp_5p_sub = substr($adp_5p, -$adp_5p_len);
		$adp_5p_sub_long = substr($adp_5p, -$adp_5p_len_long);
	} 
	my $long_match_end_enable = 1;
	
	my $distance = 1;
	$distance = int($$options{'d'}) if defined $$options{'d'} && $$options{'d'} >= 0;

	my $cutoff_min_len = 15;
	$cutoff_min_len =  $$options{'l'} if (defined $$options{'l'} );

	# set hash for sRNA size distribution
	my %sRNA_len_dist;	# key1: sample; key2: size; value: count
	my %sRNA_size;		# key: size; value: 1


	my $report_file = 'report_sRNA_trim.txt';
	print "[ERR]report file exist\n" if -s $report_file;

	my $report_info = '';
	$report_info.= "#sRNA\ttotal\t3Punmatch\t3Pnull\t3Pmatch\tbaseN\tshort\tcleaned\n";

	foreach my $inFile ( @$files ) 
	{
		my $prefix = $inFile;
		$prefix =~ s/\.gz$//; $prefix =~ s/\.fastq$//; $prefix =~ s/\.fq$//; $prefix =~ s/\.fasta$//; $prefix =~ s/\.fa$//;

		my $outFile1 = $prefix.".clean.fq";
		#my $outFile2 = $prefix.".clean.report.txt";
		#print "[ERR]out file exist\n" and exit if -s $outFile1;
		my $out1 = IO::File->new(">".$outFile1) || die $!;
		#my $out2 = IO::File->new(">".$outFile2) || die $!;
		#print $out2 "#ReadID\tlength\tstart\tend\t3p_edit_distacne\tlabel\n";

		my $fh;
		if ($inFile =~ m/\.gz$/) { 
			$fh = IO::File->new("gunzip -c $inFile | ") || die $!;
		} else { 
			$fh = IO::File->new($inFile) || die $!;
		}
	
		my ($format, $id1, $seq, $id2, $qul, $read_len);
		my ($total_num,$unmatch_5p_num,$null_5p_num,$match_5p_num,
			$unmatch_3p_num,$null_3p_num,$match_3p_num,
			$baseN_num,$short_num,
			$clean_num,$adp_clean,$adp5p_clean,$adp3p_clean)=(0,0,0,0,0,0,0,0,0,0,0,0,0);
		my ($mode_5p_a, $mode_5p_b, $mode_5p_a_match) = (0, 0, 0);
		while(<$fh>)
		{
			$id1 = $_;	chomp($id1);
			$seq = <$fh>;	chomp($seq);	$seq = uc($seq);
			$read_len = length($seq);
	
			if      ($id1 =~ m/^>/) { $format = 'fasta'; $id1 =~ s/^>//; }
			elsif   ($id1 =~ m/^@/) { $format = 'fastq'; $id1 =~ s/^@//; }
			else    { die "[ERR]seq format $id1\n"; }

			# match 3' adapter to reads,
			# this method will find the best adapter
			my ($pos_3p, $match_ed);
			if (defined $adp_3p) {
				for (my $d=0; $d <=$distance; $d++) 
				{
					for (my $i=0; $i<(length($seq)-$adp_3p_len+1); $i++)
					{
						my $read_substr = substr($seq, $i, $adp_3p_len);
                        			my $edit_distance = hamming($read_substr,$adp_3p_sub);
                        			if ( $edit_distance <= $d ){
                                			$pos_3p = $i;
							last;
                        			}
					}
					$match_ed = $d;
					last if defined $pos_3p;
				}
			}
			$pos_3p = length($seq) unless defined $pos_3p;

			# match 5' adapter to reads
			my $pos_5p = 0;
			if (defined $adp_5p)
			{
				# find the chimera locate in 5p using longer adp seed.
				# be careful, must analyze them before trimming
				my $long_pre_match_end = 0;
				my $long_match_end = 0;

				if ($long_match_end_enable == 1) {
					while($seq =~ /\Q$adp_5p_sub_long\E/g) {
						$long_match_end = pos($seq) - 1;
						$long_pre_match_end = $long_match_end + 1;
					}
			
					# report the pattern of matched reads		
					if ($long_pre_match_end > 0) {
                                        	my $sub1 = substr($seq, 0, $long_pre_match_end);
						my $sub2 = substr($adp_5p, -$long_pre_match_end);
                                        	my $mm = 'unmatch';
						my $ed = 0;
                                        	if ($sub1 eq $sub2) { 
							$mode_5p_a_match++; 
							$mm = "match"; 
						} else {
							$ed = hamming($sub1, $sub2);
						}
						# my $perc = sprintf("%.2f", ($ed/length($sub1)) );					
	
						# rules to perform trim
						# 1 trim match and unmatch if length <=20 or >=25 
						# 2 trim if edit distance <= 3 (21-24nt)
						if (length($seq)>=21 && length($seq)<=24 && $ed > 3) {
							$long_pre_match_end = 0;
						}
                                        	# print "$seq, $long_pre_match_end, $sub1, $sub2, $mm, $ed, $perc\n";
					}
				}

				# find the position locate with 5 base seed
				my ($match_start, $match_end, $match_len, $match_seq, $match_adp, $match_5p_ed, $pre_match_end);
				$match_start = 0;
				$match_end = 0;
				$pre_match_end = 0;
				while($seq =~ m/\Q$adp_5p_sub\E/g) {
					$match_end = pos($seq) - 1;
					#print "$seq, $adp_5p_sub\n$match_end\n";
					$match_len = $match_end - $match_start + 1;
					$match_seq = substr($seq, $match_start, $match_len);
					$match_adp = substr($adp_5p, -length($match_seq));
					#print "$match_seq, $match_adp\n";
					$match_5p_ed = hamming($match_seq, $match_adp);
					#print "$match_5p_ed\n";
					if ($match_5p_ed > 0) { last; }
					$match_start = $match_end + 1;
					$pre_match_end = $match_start;
				}

				# locate the best match 
				if ($long_pre_match_end > $pre_match_end) {
					$pos_5p = $long_pre_match_end;
					$mode_5p_a++;
				} else {
					$pos_5p = $pre_match_end;
					$mode_5p_b++;
				}
			}

			# read id and qul from fastq file
			if ( $format eq 'fastq' ) {
				$id2 = <$fh>;   chomp($id2);	
				$qul = <$fh>;   chomp($qul);
			}

			# output result 
			my $label = '';
			if (defined $adp_5p) {
				if ($pos_5p == 0 ) {
					$label.=",5p_unmatch";
					$unmatch_5p_num++;
				} elsif ($pos_5p == length($seq)) {
					$label.=",5p_null";
					$null_5p_num++;
				} else {
					$label.=",5p_match"; 
					$match_5p_num++;
				}
			}

			if (defined $adp_3p) {
				if ( $pos_3p == 0 ) {
					$label.=",3p_null";
					$null_3p_num++;
				} elsif ($pos_3p == length($seq)) {
					$label.=",3p_unmatch";
					$match_ed = "NA";
					$unmatch_3p_num++;
				} else {
					$label.=",3p_match";
					$match_3p_num++;
				}
			}

			#if ($label =~ m/3p_unmatch/) {
			#	print $out1 "@".$id1."\n".$seq."\n".$id2."\n".$qul."\n";
			#}

			if ( ($pos_3p > $pos_5p) && ($label !~ m/3p_null/) && ($label !~ m/3p_unmatch/) && ($label !~ m/5p_null/) )
			{
				my $trimmed_len = $pos_3p - $pos_5p;
				my $trimmed_seq = substr($seq, $pos_5p, $trimmed_len);
				my $trimmed_qul = substr($qul, $pos_5p, $trimmed_len) if $format eq 'fastq';
				my $baseN = $trimmed_seq =~ tr/N/N/;
				if ($baseN > 0) {
					$label.=",baseN";
					$baseN_num++;
				} elsif (length($trimmed_seq) < $cutoff_min_len ) {
					$label.=",short";
					$short_num++;
				} else {
					$clean_num++;
					if (defined $adp_5p && defined $adp_3p && $pos_5p > 0 && $pos_3p > 0) {
						$adp_clean++;
					} elsif (defined $adp_5p && $pos_5p > 0) {
						$adp5p_clean++;
					} elsif (defined $adp_3p && $pos_3p > 0) {
						$adp3p_clean++;
					}

					if ( $format eq 'fastq' ) {
						print $out1 "@".$id1."\n".$trimmed_seq."\n".$id2."\n".$trimmed_qul."\n";
					} else {
						print $out1 ">".$id1."\n".$trimmed_seq."\n";
					}

					# save hash for size distribution
					$sRNA_size{$trimmed_len} = 1 unless defined $sRNA_size{$trimmed_len};
					defined $sRNA_len_dist{$prefix}{$trimmed_len} and $sRNA_len_dist{$prefix}{$trimmed_len}++ or $sRNA_len_dist{$prefix}{$trimmed_len} = 1;
				}
			} 

			$label =~ s/^,//;
			$total_num++;
			$match_ed = 'NA' unless defined $match_ed;
			#print $out2 "$id1\t$read_len\t$pos_5p\t$pos_3p\t$match_ed\t$label\n";
 
		}
		$fh->close;
		$out1->close;
		#$out2->close;
		$report_info.="$inFile\t$total_num\t".
			"$unmatch_3p_num\t$null_3p_num\t$match_3p_num\t".
			"$baseN_num\t$short_num\t$clean_num\n";

		# debug 5p info
		#print "$inFile\t$mode_5p_a\t$mode_5p_a_match\t$mode_5p_b\n";
	}

	# report sRNA trim information
	my $outr = IO::File->new(">$report_file") || die $!;
	print $outr $report_info;
	$outr->close;

	# report sRNA length distribution
	my $out_len = IO::File->new(">sRNA_length.txt") || die $!;
	print $out_len "#sample";
	foreach my $len (sort {$a<=>$b} keys %sRNA_size ) {
		print $out_len "\t".$len;
	}
	print $out_len "\n";
	foreach my $sample (sort keys %sRNA_len_dist) {
		print $out_len $sample;
		foreach my $len (sort {$a<=>$b} keys %sRNA_size ) {
			my $num = 0;
			$num = $sRNA_len_dist{$sample}{$len} if defined $sRNA_len_dist{$sample}{$len};
			print $out_len "\t".$num;
		}	
		print $out_len "\n";	
	}
	$out_len->close;

}

sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] )  }

