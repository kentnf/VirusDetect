package Util;

use strict;
use warnings;
use Bio::SeqIO; 
use Bio::Graphics;
use Bio::SeqFeature::Generic;

sub line_num {
	my $input = shift;
	my $num = $input =~ tr/\n/\n/;
	return $num;
}

# check if the program avaiable 
sub check_program {
	my ($bin, $program) = @_;
}

# load fasta sequence to hash
sub load_seq {
	my $input_seq = shift;
	my %seq_info;

	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_seq);
	while(my $inseq = $in->next_seq)
	{
		$seq_info{$inseq->id}{'seq'} = $inseq->seq;
		$seq_info{$inseq->id}{'length'} = $inseq->length;
	}
	return %seq_info;
}

# remove duplicated reads
sub remove_dup {
	my $input_seq = shift;
	my $uniq_file = $input_seq.'.uniq';

	# hash for store unique seq
	my %uniq_seq;

	open(OUT, ">".$uniq_file) || die $!;
	open(IN, $input_seq) || die $!;
	while(<IN>) {
		my $file_type;
		my $sid = $_; chomp($sid);
		if		($sid =~ m/^>/) { $file_type = 'fasta'; }
		elsif	($sid =~ m/^@/) { $file_type = 'fastq'; }
		else	{ die "[ERR], can not detect the file type for file: $input_seq\n"; }
		my $seq = <IN>; chomp($seq);

		my $out_seq = "$sid\n$seq\n";
		if ($file_type eq 'fastq') { $out_seq.=<IN>; $out_seq.=<IN>; }

		unless (defined $uniq_seq{$seq}) {
			print OUT $out_seq;	
			$uniq_seq{$seq} = 1;
		}
	}
	close(IN);
	close(OUT);
	#process_cmd("mv $uniq_file $input_seq");
}

sub process_cmd {
	my ($cmd, $debug) = @_;	
	if (defined $debug && $debug > 0) { print "[CMD FOR DEBUG]: $cmd\n"; }
	my $ret = system($cmd);	
	if ($ret) {
		die "[ERR IN CMD]: \n$cmd\n Died with ret: $ret\n";
	}
	return($ret);
}

sub detect_FileType {
	my $file = shift;
	die "[ERR]Undef $file\n" unless -s $file; 
	open(FH, $file) || die $!;
	my $line = <FH>;
	my $file_type;
	if	($line =~ m/^>/) { $file_type = 'fasta'; }
	elsif	($line =~ m/^@/) { $file_type = 'fastq'; }
	else	{ die "[ERR], can not detect the file type for file: $file\n"; }
	return $file_type;
}

# determine the dataType (sRNA, mRNA) by sequence length
sub detect_DataType {
	my $file = shift;
	die "[ERR]Undef $file\n" unless -s $file;
	my $n = 0;
	my $total = 0;
	open(FH, $file) || die $!;
	while(<FH>) {
		my $file_type;
		my $sid = $_;	chomp($sid); 
		if  	($sid =~ m/^>/) { $file_type = 'fasta'; }
		elsif   ($sid =~ m/^@/) { $file_type = 'fastq'; }
		else    { die "[ERR], can not detect the file type for file: $file\n"; }
		my $seq = <FH>;	chomp($seq);
		$total += length($seq);
		if ($file_type eq 'fastq') { <FH>; <FH>; }
		$n++;
		last if $n > 1000;
	}

	my $data_type = 'sRNA';
	my $mean = $total / $n;
	if ($mean >= 36) {
		$data_type = 'mRNA';
	}

	return $data_type;
}

sub detect_seqNum {
	my $file = shift;
	die "[ERR]Undef $file\n" unless -s $file;
	my $seq_num = 0;
	my $file_type;
	open(FH, $file) || die $!;
	while(<FH>)
	{
		my $id = $_;
		if      ($id =~ m/^>/) { $file_type = 'fasta'; }
		elsif   ($id =~ m/^@/) { $file_type = 'fastq'; }
		else    { die "[ERR], can not detect the file type for file: $file\n"; }
		<FH>;
		$seq_num++;
		if ($file_type eq 'fastq') { <FH>; <FH>; }
	}
	close(FH);
	return $seq_num;
}

sub print_user_message {
	my @message = @_;
	#print "\n";
	my $time = get_time();
	foreach my $line (@message) { 
		print $time." ".$line."\n";
	}
	#print "\n";
}

sub print_user_submessage{
	my @message = @_;
	foreach my $line (@message) { 
		print "   ".$line."\n";
	}
}

sub get_time {
	my $time = `date +"%D %T"`;
	chomp $time;
	$time = "[$time]";
	return $time
}

sub save_file {
	my ($content, $file) = @_;
	if (defined $content && length($content) > 1) {
		open(FH, ">".$file) || die $!;
		print FH $content;
		close(FH);
	}
}

=head2
 host_subtraction -- subtract host derived contigs
=cut
sub host_subtraction
{
	my ($contigs, $input_blast_table, $identity, $min_cov) = @_;
      
    # put blast info to hash
    # %blk ->
    # key: query ID
    # value: [Query start, Query End], [Query start, Query End] ....
    # %query_length -> keu: query ID, value: query length
    # %host -> key: contig ID, value: 1
    my %blk;
    my %query_len;
    my %host;
    
	chomp($input_blast_table);
	my @a = split(/\n/, $input_blast_table);
	foreach my $line (@a)
	{
		# $query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end, $identity2, $aligned_query, $aligned_hit, $aligned_string#
		chomp($line); 
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);

		if ($ta[5]>=$identity) {   
			push(@{$blk{$ta[0]}}, [@ta[9,10]]); # create blk hash for query ID and [query_start,query_end], query is contig
			defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];
		}
	}

	foreach my $tk (sort keys %blk)
	{   
		my @o; # put non-overlap blocks to array @o for one query
		for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}})
		{   
			if (scalar(@o) == 0)
			{   
				push(@o, [@$ar]);
			} 
			elsif ($o[-1][1] >= $ar->[0]-1) # combine overlapped block
			{   
				$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1];
			}
			else
			{   
				push(@o, [@$ar]);
			}
		}

		my $total_cov = 0; # compute total coverage according to non-overlap blocks for one query 
		for my $ar ( @o )
		{     
			# print information for error checking
			# #query_name query_len block_start block_end block_length
			# print join("\t", $tk, $query_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1)."\n"; 
			$total_cov += ($ar->[1]-$ar->[0]+1);
		}
        
		my $ratio=int($total_cov/$query_len{$tk}*10000+0.5)/100;
		#print "XXX:$tk\t$ratio\n";
        if($ratio >= $min_cov) {   
			defined $host{$tk} or $host{$tk} = 1;
		}
    }

    # filter assembled contigs to subtract host 
	my $host_rm_seq = '';
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$contigs);
	while(my $inseq = $in->next_seq) {
		my $id = $inseq->id;
		my $sq = $inseq->seq;
		#print ">$id\n$sq\n";
		next if defined $host{$id};
		$host_rm_seq .= ">$id\n$sq\n";
	}

	my $fh = IO::File->new(">".$contigs) || die $!;
	print $fh $host_rm_seq;
	$fh->close;

	return %host;
}

=head2
 xa2multi -- convert xa tag to multiply hit
=cut
sub xa2multi
{
	my $input_sam = shift;
	my $tem_sam = $input_sam.".temp";

	my $out = IO::File->new(">".$tem_sam) || die $!;
	my $in = IO::File->new($input_sam) || die $!;
	while (<$in>) 
	{
		if (/\tXA:Z:(\S+)/)
		{
			my $l = $1;
			print $out $_;
			my @t = split("\t", $_);

			while ($l =~ /([^,;]+),([-+]\d+),([^,]+),(\d+);/g)
			{
				my $mchr = ($t[6] eq $1)? '=' : $t[6]; # FIXME: TLEN/ISIZE is not calculated!
				my $seq = $t[9];
				my $phred = $t[10];
				# if alternative alignment has other orientation than primary, 
				# then print the reverse (complement) of sequence and phred string
				if ((($t[1]&0x10)>0) xor ($2<0)) {
					$seq = reverse $seq;
					$seq =~ tr/ACGTacgt/TGCAtgca/;
					$phred = reverse $phred;
				}
				print $out (join("\t", $t[0], ($t[1]&0x6e9)|($2<0?0x10:0), $1, abs($2), 0, $3, @t[6..7], 0, $seq, $phred, "NM:i:$4"), "\n");
			}
		}
		else
		{
			print $out $_;
		}
	}
	$in->close;
	$out->close;
	process_cmd("mv $tem_sam $input_sam");
}

=head2
 load_virus_info -- load virus information from vrl_genebank file 
 discard for high memory usage

 fetch_virus_info -- fetch virus information from vrl_genebank file according to virus id list
 this will only use less memory for fetch identified virus info

=cut
sub fetch_virus_info
{
	my ($file, $prot_tab, $iden_vid) = @_;

	# load protein tab info to hash
	my %vid_pid; # key: virus ID, value: pid1 # pid2 # ... # pidN
	if (defined $prot_tab && -s $prot_tab) {
		my $fh1;
		if ($prot_tab =~ m/\.gz$/) {
			open($fh1, "-|", "gzip -cd $prot_tab") || die $!
		} else {
			open($fh1, $prot_tab) || die $!;
		}
		while(<$fh1>) {
			chomp;
			next if $_ =~ m/^#/;
			my @a = split(/\t/, $_);
			my $nid = $a[0];
			my $pid = $a[1];
			if (defined $$iden_vid{$pid}) {
				if ( defined $vid_pid{$nid} ) {
					$vid_pid{$nid}.= "#".$pid;
				} else {
					$vid_pid{$nid} = $pid;
				}
			}
		}
		$fh1->close;
	}

	# load virus info to hash
	my %virus_info;
	my $fh;	
	if  ($file =~ m/\.gz$/) {
		open($fh, "-|", "gzip -cd $file") || die $!;
	} else {
		open($fh, $file) || die $!;
	}

	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		# AB000048	2007	Parvovirus	Feline panleukopenia virus gene for nonstructural protein 1, complete cds, isolate: 483.	1	Vertebrata
		if (defined $$iden_vid{$a[0]}) {
			$virus_info{$a[0]}{'length'}	= $a[1];
			$virus_info{$a[0]}{'genus'} 	= $a[2];
			$virus_info{$a[0]}{'desc'} 		= $a[3];
			$virus_info{$a[0]}{'verison'} 	= $a[4];
			$virus_info{$a[0]}{'host_type'} = $a[5];
		}

		if (defined $vid_pid{$a[0]}) {
			my @b = split(/#/, $vid_pid{$a[0]});
			foreach my $pid (@b) {
				$virus_info{$pid}{'length'}    = $a[1];
				$virus_info{$pid}{'genus'}     = $a[2];
				$virus_info{$pid}{'desc'}      = $a[3];
				$virus_info{$pid}{'verison'}   = $a[4];
				$virus_info{$pid}{'host_type'} = $a[5];
			}
		}
	}
	close($fh);
	return %virus_info;
}

=head2 

 parse_blast_to_table -- this script is used for parse blastn result with fast speed than Bio::SearchIO

 input file is blast output, output is tabular information
 bioperl: $hsp->start('hit') always smaller than $hsp->end('hit'), this is different with blast table output
 This script use the bioperl format as output, it will have addtional strand info

 modified by kentnf
 20140715: 1. combine blast_parse_table22.pl and blast_parse_table24.pl
		   2. make it to subroutine parse_blast_to_table
=cut

sub parse_blast_to_table 
{
	my ($input_blast, $program) = @_;

	my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end);
	my ($identity2, $aligned_query, $aligned_hit, $aligned_string); # will output four additional col
	my ($query_strand, $hit_strand);				# for tblastx

	my (%hsp, $hsp_count); 
	%hsp = (); $hsp_count=0; 

	my $query_sequenc; 
	my $subject_sequenc; 	
	my $is_hsp = 1;

	# for different program
	my $blast = '';
	if    ($program =~ m/megablast/ || $program =~ m/blastn/) { $blast = 'blastn'; }
	elsif ($program =~ m/tblastx/) { $blast = 'tblastx'; }
	elsif ($program =~ m/blastx/)  { $blast = 'blastx'; }
	else  {die "[ERR]Undef blast program $program\n"; }

	# output table title
	my $output_table = "#query_name\tquery_length\thit_name\thit_length\thsp_length\tidentity\tevalue\tscore\tstrand".
			   "\tquery_start\tquery_end\thit_start\thit_end\tidentity2\taligned_query\taligned_hit\taligned_string\n";

	open(IN, "$input_blast") || die $!;
	while (<IN>) 
	{
		chomp;
		if (/Query=\s(\S+)/ || eof) # output all hsp info if the line has new Query name, or End of file 
		{
			# previous hsp information
			if ($hsp_count > 0) 
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
				$hsp{$hsp_count} = $hsp_info;
			}

			#########################################
			#  output all previous query hit hsp    #
			#########################################
			if (scalar(keys(%hsp)) > 0)
			{
				my $aa = scalar(keys(%hsp));
				foreach my $hsp_num (sort {$a<=>$b} keys %hsp)
				{
					$output_table.= $hsp{$hsp_num}."\n";
				}
			}

			# start new query record
			%hsp = ();$hsp_count = 0;

			$query_name = $1; $query_length = ""; $hit_name = ""; $hit_length = "";
		}
		elsif (/\s+\((\S+)\sletters\)/) # get Query Length
		{
			$query_length = $1;
			$query_length =~ s/,//g;
		}
		elsif (/>(\S+)/)# Hit Name, parse hit info
		{
			# previous hsp information
			if ($hsp_count > 0 || eof)
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
				$hsp{$hsp_count} = $hsp_info;
				$is_hsp = 0;			
			}

			# start new hit
			$hit_name = $1; $hit_length = "";

		}
		elsif (/\s+Length = (\d+)/)
		{
			$hit_length = $1;
			$hit_length =~ s/,//g;
		}
		elsif (/Score =\s+(\S+) bits.+Expect(\(\d+\))? = (\S+)/) # parse hsp info
		{
			if ($hsp_count > 0 && $is_hsp == 1)
			{	
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
			        	        $hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
				                $query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
				$hsp{$hsp_count} = $hsp_info;
			}

			#  new hsp record
			$is_hsp = 1;
			$hsp_count++; 
			$score = $1; $evalue = $3;
			$evalue = "1$evalue" if ($evalue =~ m/^e/);
			$query_start = 0; $query_end = 0; $hit_start = 0; $hit_end = 0;
			$aligned_query = ""; $aligned_hit= ""; $aligned_string= "";

		}
	
		elsif (/\s+Identities = (\d+)\/(\d+)\s+\((\S+)\)/ && $hsp_count >= 1) 	# identity info, hsp length
		{
			$identity = $1/$2*100;
			$identity = sprintf("%.".(2)."f", $identity);
			if ( $1 > $2 ) { $hsp_length = $1; $hsp_length =~ s/,//g; } else { $hsp_length = $2; $hsp_length =~ s/,//g; }
			$identity2 = $1."/".$2."(".$3.")";
		}
		elsif (/\s+Strand = (\S+) \/ (\S+)/ && $hsp_count >= 1 ) 		# strand info
		{
			die "[ERR]Unmatch blast:$blast\t$_\n" if $blast ne 'blastn';
			if ( $2 eq "Plus" ) { $strand = 1;} else { $strand = -1;}
		}
				
		elsif (/\s+Frame = ([\+\-]\d) \/ ([\+\-]\d)/ && $hsp_count >= 1)
		{
			die "[ERR]Unmatch blast:$blast\t$_\n" if $blast ne 'tblastx';
                	my ($a, $b)=($1,$2);
                	$strand = "[".$a."/".$b."]";
                	if ( $a =~ /\+/ ) { $query_strand = 1;} else { $query_strand = -1;}
                	if ( $b =~ /\+/ ) { $hit_strand = 1;} else { $hit_strand = -1;}
        	}
		
		elsif (/\s+Frame = (\S+)/ && $hsp_count >= 1)
		{
			die "[ERR]Unmatch blast:$blast\t$_\n" if $blast ne 'blastx';
			my $a = $1;
			if ( $a =~ /\+/ ) { $strand = 1;} else { $strand = -1;}
		}

		elsif (/Query\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1) 		# query sequence, && aligned string
		{
			if ($query_start == 0) { $query_start = $1; $query_start =~ s/,//g;}
			$query_end = $3;
			$query_end =~ s/,//g;

			if (defined $query_strand && $query_strand == -1 && $blast eq 'tblastx')
			{
				if ($query_end == 0) { $query_end = $1; $query_end =~ s/,//g;};
				$query_start = $3;
				$query_start =~ s/,//g;
			}

			chomp($_);
			$aligned_query .= $_.",";
			my $next_line=<IN>;
			chomp($next_line); 
			$aligned_string .= $next_line.",";
		}
		elsif (/Sbjct\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1)		# subject sequence
		{
			if (defined $hit_strand && $blast eq 'tblastx') 
			{
				$strand = $hit_strand;
			}

			if ( $strand == -1 )
			{
				if ($hit_end == 0) { $hit_end = $1; $hit_end =~ s/,//g; };
				$hit_start = $3;
				$hit_start =~ s/,//g;
			}
			else
			{
				if ($hit_start == 0) { $hit_start = $1; $hit_start =~ s/,//g;};
				$hit_end = $3;
				$hit_end =~ s/,//g;
			}
			chomp($_);
			$aligned_hit .= $_.",";
		}
		else
		{
			next;
		}
	}
	close(IN);

	return $output_table;
}

=head2
 
 filter_blast_table -- filter blast table according to evalue and identity

 Rules: the coverage of hsp should be more than 0.75 (controlled by coverage param) for query or hit
	for one query, get the best hsp base on evalue; then get the best identity of the best hit. 
        keep the other hits which show identity lower than 5% (controlled by drop_off param) of best identity.

 Usage: $filtered_blast_table = filter_blast_table(input_blast_table, coverage, identity_dropoff, blast_program);

 modified by kentnf
 20140715: combine blast_filter3 and blast_filter5
           change blast_filter to subroutine parse_blast_table
=cut
sub filter_blast_table
{
	my ($input_blast_table, $coverage, $identity_dropoff, $blast_program) = @_;

	# param will be 1 for megablastn/blastn, param will be 5 for tblastx
	my $param;
	if ( $blast_program =~ m/megablast/ || $blast_program =~ m/blastn/ ) { $param = 1; }
	elsif ($blast_program =~ m/tblastx/ || $blast_program =~ m/blastx/) { $param = 3; }
	else  { die "[ERR]blast program: $blast_program\n"; }

	chomp($input_blast_table);
	my @a = split("\n", $input_blast_table);

	# define var for last query hit
	my $last_query = '';
	my $last_hit = '';
	my $last_query_len;
	my $last_hit_len;
	my $current_query;
	my $current_hit;
	#my $high_identity;
	#my $current_identity;

	# defined output 
	my $output_blast_table = '';

	# define array to store the multiple hsp for one pairs of query and hit
	# [$query_start, $query_end, $hit_start, $hit_end, $identity, $line], 
	# [$query_start, $query_end, $hit_start, $hit_end, $identity, $line],
	#  ....
	my @hsp_region = ();
	
	# defined hash to store highest_iden for each contigs;
	# key: query ID
	# value: identity;
	my %high_iden;

	foreach my $line ( @a )
	{
		if ($line =~ m/^#/) { $output_blast_table.=$line."\n"; next; }
		my @cols = split(/\t/, $line);
		# $query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end, $identity2, $aligned_query, $aligned_hit, $aligned_string#
		next if scalar @cols < 13;
		my $query_id	= $cols[0];
		my $query_len	= $cols[1];
		my $hit_id		= $cols[2];
		my $hit_len		= $cols[3];
		my $hsp_len		= $cols[4] * $param;
		my $identity	= $cols[5];
		my $query_start = $cols[9];
		my $query_end	= $cols[10];
		my $hit_start	= $cols[11];
		my $hit_end		= $cols[12];

		$current_query = $query_id;
		$current_hit = $hit_id;

		if ($current_query eq $last_query) 
		{
			if ($current_hit eq $last_hit) {
				# current and hit are same as previous, pub the his position to array
				push(@hsp_region, [$query_start, $query_end, $hit_start, $hit_end, $identity, $line]);
			} 
			else 
			{
				# same contig hit to different vrl
				my ($m_cov, $m_iden, $m_align) = get_iden_cov(\@hsp_region, $last_query_len, $last_hit_len);
				if ($m_cov >= $coverage) {
					$high_iden{$last_query} = $m_iden unless defined $high_iden{$last_query};
					if ($m_iden >= $high_iden{$last_query} - $identity_dropoff) {
						$output_blast_table.= $m_align;
					}
				}

				# start new for next
				@hsp_region = ();
				push(@hsp_region, [$query_start, $query_end, $hit_start, $hit_end, $identity, $line]);
			}
		} 
		else 
		{
			# process previous query contigs
			if (@hsp_region > 0) { # skip 1 line of blast
				my ($m_cov, $m_iden, $m_align) = get_iden_cov(\@hsp_region, $last_query_len, $last_hit_len);
				if ($m_cov >= $coverage) {
					# define the highest identity if the coverage meet requirement
					$high_iden{$last_query} = $m_iden unless defined $high_iden{$last_query};
					if ($m_iden >= $high_iden{$last_query} - $identity_dropoff) {
   		             	$output_blast_table.= $m_align;
					}
				}
			}

			# start new for next
			@hsp_region = (); 
			push(@hsp_region, [$query_start, $query_end, $hit_start, $hit_end, $identity, $line]);

		}

		$last_query = $query_id;
		$last_hit = $hit_id;
		$last_query_len = $query_len;
		$last_hit_len = $hit_len;
	}

	# parse the last record
	if (scalar @hsp_region > 0) {
		my ($m_cov, $m_iden, $m_align) = get_iden_cov(\@hsp_region, $last_query_len, $last_hit_len);
		if ($m_cov >= $coverage) {
			# define the highest identity if the coverage meet requirement
			$high_iden{$last_query} = $m_iden unless defined $high_iden{$last_query};
			if ($m_iden >= $high_iden{$last_query} - $identity_dropoff) {
				$output_blast_table.= $m_align;
			}
		}
	}

	return $output_blast_table;
}

# sidekick of filter blast table
sub get_iden_cov
{
	my ($hsp_region, $query_len, $hit_len) = @_;

	my $identity = 0;
	my $coverage = 0;
	my $alignment = '';

	# compute query coverage
	# for store query overlap
	my @qo;	
	foreach my $hr (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$hsp_region}) {
		# $hr (hsp region) is array of:
		# [$query_start, $query_end, $hit_start, $hit_end, $identiy, $alignment]
		if (scalar(@qo) == 0) {   
			push(@qo, [$hr->[0], $hr->[1]]);
		}
		elsif ($qo[-1][1] >= $hr->[0]-1) {  # combine overlapped block  
			$qo[-1][1] < $hr->[1] and $qo[-1][1] = $hr->[1];
		}
		else {   
			push(@qo, [$hr->[0], $hr->[1]]);
		}

		$identity = $hr->[4] if $identity == 0; # use the best hit for identity;
		$alignment .= $hr->[5]."\n";
	}

	my $query_cov_len = 0;
	foreach my $qr (@qo) {
		 $query_cov_len += abs($qr->[1] - $qr->[0] + 1);
	}

	unless (defined $query_len) {
		print "DEBUG\n";
		print $alignment."\n";
		foreach my $hr (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$hsp_region}) {
			print $hr->[5]."\n";
		}
		print $hit_len."\n";
		print "DEBUG\n";
		exit;
	} 
	
	$coverage = $query_cov_len/$query_len if $query_cov_len/$query_len > $coverage;

	# compute hit coverage
	# for store hit coverage
	my @ho;
	foreach my $hr (sort { $a->[2]<=>$b->[2] || $a->[3]<=>$b->[3];} @{$hsp_region}) {
		# $hr (hsp region) is array of:
		# [$query_start, $query_end, $hit_start, $hit_end, $identiy, $alignment]
		if (scalar(@ho) == 0) {
			push(@ho, [$hr->[2], $hr->[3]]);
		}
		elsif ($ho[-1][1] >= $hr->[2]-1) {  # combine overlapped block
			$ho[-1][1] < $hr->[3] and $ho[-1][1] = $hr->[3];
		}
		else {
			push(@ho, [$hr->[2], $hr->[3]]);
		}
	}

	my $hit_cov_len = 0;
	foreach my $hr (@ho) {
		$hit_cov_len += abs($hr->[1] - $hr->[0] + 1);
	}
	$coverage = $hit_cov_len/$hit_len if $hit_cov_len/$hit_len > $coverage;

	return ($coverage, $identity, $alignment);
}

=head2

 find_known_contig -- label the contig with 'known' according to blast result

 modified by kentnf
 20140715: change output file to hash, make it subroutine as modular

 usage: %known = find_known_contig($blast_table, identity, min_ratio)

=cut
sub find_known_contig
{
	my ($input_blast_table, $identity, $min_ratio) = @_;
	
	# put blast info to hash
	# %blk ->
	# key: query ID
	# value: [Query start, Query End], [Query start, Query End] ....
	# %query_length -> keu: query ID, value: query length
	# %known -> key: virus ref ID, value: 1
	my %blk;
	my %query_len;
	my %hit_len;
	my %known;	

	chomp($input_blast_table);
	my @a = split(/\n/, $input_blast_table);
	foreach my $line (@a)
	{
		chomp($line);
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		
		if ($ta[5]>=$identity)
		{
			push(@{$blk{$ta[0]}}, [@ta[9,10]]); 				# create blk hash for query ID and [query_start,query_end], query is contig
			defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];
			defined $hit_len{$ta[2]} or $hit_len{$ta[2]} = $ta[3];
		}
	}

	foreach my $tk (sort keys %blk) 
	{
		my @o; # put non-overlap blocks to array @o for one query
		for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}}) 
		{
			if (scalar(@o) == 0) 
			{
				push(@o, [@$ar]);
			} 
			elsif ($o[-1][1] >= $ar->[0]-1) # combine overlapped block
			{ 
				$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1];
			}
			else
			{
				push(@o, [@$ar]);
			}
		}

		my $total_cov = 0; # compute total coverage according to non-overlap blocks for one query 
		for my $ar ( @o ) 
		{
			# print information for error checking
			# #query_name query_len block_start block_end block_length
			# print join("\t", $tk, $query_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1)."\n"; 
			$total_cov += ($ar->[1]-$ar->[0]+1);
		}

		my $ratio_q=int($total_cov/$query_len{$tk}*10000+0.5)/100;
		#my $ratio_s=int($total_cov/$hit_len{$tk}*10000+0.5)/100;

		if($ratio_q >= $min_ratio) {
			defined $known{$tk} or $known{$tk} = 1;
		}
	}

	# filter blast_table to generante known_blast_table
	my $known_blast_table = '';
	foreach my $line (@a)
	{
		chomp($line);
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);

		if ( defined $known{$ta[0]} ) 
		{
			$known_blast_table.=$line."\n";
		}
	}

	return (\%known, $known_blast_table);
}

=head2
 get_hit_coverage: compute coverge of hit (vrl reference) according to blast table
 usage: hit_cov(input, output1, output2, identify, query_cov, ???param???);
 $sample.known.table", "$sample.known.cov", "$sample.known.block" 60, 0.5, 1

 # input is blast table informat : 
 # query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t 
 # strand \t query_start \t query_end \t hit_start \t hit_end ... \n
 # index 2 3 4 5 11 12 are required
 
 # for each pair of [query, hit] should meet the requirement of identity and query_cov
 # ratio = all query coverage / hit_length
 # hit: query1, query2, ... queryN

=cut
sub get_hit_coverage
{
	my ($input_blast_table, $cutoff_identity, $query_cov, $param) = @_;

	# put hit info to hash hash
	# hit is virus reference 
	# query is assembled contigs
	#
	# %bkl
	# key: hit_name, 
	# value: arrays of block and query
	#        for each element in array is another array include three element [hit_start, hit_end, hash_of_query]
	#       [hit_start, hit_end, hash_of_query], 
	#       [hit_start, hit_end, hash_of_query], 
	#       ....
	#
	#        % hash_of_query
	#        key: queryID
	#        value: 1
	#
	# %hit_len
	# key: hit_name
	# value: hit_length
	my %blk;
	my %hit_len;
	
	chomp($input_blast_table);
	my @a = split("\n", $input_blast_table);
	foreach my $line (@a) 
	{
		chomp($line);
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		die "Error in line $_" if scalar @ta < 13;
		my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue,
		$score, $strand, $query_start, $query_end, $hit_start, $hit_end) = (@ta[0..12]);

		# remove the identity and query check, because the blast table has been filtered
		#   using query identity, and (query | hit) coverage
		#my $query_covered= $param * $hsp_length / $query_length;
		#if ( $identity >= $cutoff_identity && $query_covered >= $query_cov)
		#{
			# put hit and (hit_start,hit_end) to hash

			# change for blastX result
			if ($hit_start > $hit_end) {
				my $hit_start_temp = $hit_start;
				$hit_start = $hit_end;
				$hit_end = $hit_start_temp;
			}

			push(@{$blk{$hit_name}}, [$hit_start, $hit_end]);
			$blk{$ta[2]}[-1][2]{$ta[0]}=1;
			defined $hit_len{$hit_name} or $hit_len{$hit_name} = $hit_length;
		#}
	}

	my %known_hit_cov;
	my $known_hit_cov = '';
	my $known_block	= "#hit_name\tblock_len\tblock_start\tblock_end\tcovered_len\tquery_names\n";

	foreach my $tk (sort keys %blk) # sort by hit name 
	{
		# structure of array '@o' and $ar
		# [array]               , [array]                 ... [array]
		# [start, end, %contigs], [start, end, %contigs], ... [start, end, %contigs]

		my @o; # array for hit of non-overlap blocks

		# sort by hit start, hit end, and query name
		foreach my $ar (sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $a->[2] cmp $b->[2];} @{$blk{$tk}})
		{
			# if ($tk eq 'D10663')
			# {
			#	print "$$ar[0]\t$$ar[1]\n";
			#	my %sss = %{$$ar[2]};
			#	foreach my $k (sort keys %sss) {
			#		print $k."\n";
			#	}
			#	die;
			# }

			if (scalar(@o) == 0) 
			{
				push(@o, [@$ar]);
			}
			elsif ($o[-1][1] >= $ar->[0]-1) 
			{ # if the hit start < previous hit end, should combine two hits
				for my $qname (keys %{$ar->[2]}) 
				{
					$o[-1][2]{$qname}=1;
				}
				$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1];
			}
			else
			{
				push(@o, [@$ar]);
			}
		}

		my $total_cov = 0; # hit totoal coverage 
		my %aa;
		for my $ar (@o) 
		{
			my @query_names = sort keys %{$ar->[2]};
	
			# output to file
			# hit_name \t hit_len \t hit_block_start \t hit_block_end \t hit_block_length \t query names
			$known_block.=join("\t", $tk, $hit_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1, join(",", @query_names))."\n";
			$total_cov += ($ar->[1]-$ar->[0]+1);            # get total length of all non-overlap block
			@aa{@query_names} = (1) x scalar(@query_names); # put all the contigs into hash table
		}

		# output file: known.cov
		# format:
		# hit_id \t hit_length \t hit_converage_len \t hit_coverage_% \t contigs_name \t contigs_num

		$known_hit_cov{$tk}{'length'} = $hit_len{$tk};
		$known_hit_cov{$tk}{'cov'} = $total_cov;
		$known_hit_cov{$tk}{'contig'} = join(",", sort keys %aa);
		$known_hit_cov{$tk}{'contig_num'} = scalar(keys %aa);

		my $bbb =join("\t", $tk, $hit_len{$tk}, $total_cov, $total_cov/$hit_len{$tk}*1.0, join(",", sort keys %aa), scalar(keys %aa))."\n"; 
		$known_hit_cov.=$bbb;
	}

	return ($known_hit_cov, $known_block);
}

=head2

 remove_redundancy_hit -- this script will filter the hit record (virus reference) accroding contig supporting

 Rule: 
 The hit record will be treated as redundancy (not print) when meet below condtions at same time
 1. the nubmer of diff contigs should be < 25%
 2. length of diff contigs < 100bp
 3. cover of different contigs < 50%

 # format of known coverage
 AB373203        9535    161     0.0168851599370739      CONTIG257,CONTIG55      2

 # format of known blast table
 parse_blast_to_table output

=cut

sub remove_redundancy_hit
{
	my ($known_coverage, $known_blast_table, $diff_ratio, $diff_contig_cover, $diff_contig_length) = @_;

	# load known coverage information to array
	# struc of array:
	# [contig num, covered bp, line],[contig num, covered bp, line], ... ,[contig num, covered bp, line]
	my (@all_data, $name, $sequence);
	$name = ''; $sequence = '';

	chomp($known_coverage);
	my @a = split(/\n/, $known_coverage);
	foreach my $line (@a)
	{
		chomp($line);
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		my @contigs= split(/,/, $ta[4]);	# column 5 has all contig name for each hit

		# shan use bowtie align read to virus ref to get covered bp
		# for simply, I use the contig covered length to instead of covered bp
		# the covered is used for sort the coverage result in next
		# my $covered_bp = $ta[7];		# 
		my $covered_bp = $ta[2];

		if ( $ta[2] eq "") { $covered_bp = 0; }
		push(@all_data, [scalar(@contigs), $covered_bp, $line]);
	}
	
	# sort all dataset by 1) number of contigs for each hit
	# 		      2) covered bp
	# 		      in descending order
	@all_data = sort { -1*($a->[0] <=> $b->[0]) || -1*($a->[1] <=> $b->[1])} @all_data;

	# load contig hit covered start, end to hash
	# hash : contig_position
	# key : contig \t hit
	# value : hit_start \t hit_end
	my %contig_position;
	
	chomp($known_blast_table);
	my @b = split(/\n/, $known_blast_table);
	foreach my $line (@b)
	{
		chomp($line);
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		$contig_position{$ta[0]."\t".$ta[2]} = "$ta[11]\t$ta[12]";
	}

	# story non-redundancy and redundancy lines to @inset, and @restset 
	# @inset : lines of non-redundancy contigs
	# @restset : lines of redundancy contigs
	my @inset;
	my @restset = ();
	my $contig_count = 1;

	foreach my $tr (@all_data)
	{
		if (scalar(@inset)  == 0) 
		{
			push(@inset, $tr->[2]);
		}
		else
		{
			my @aa = split(/\t/, $tr->[2]);
			# format of @aa, $tr->[2]
			# ref_ID \t length \t HSP_len \t Coverage \t Contigs1,...,ContigsN \n;
			my ($hit_id, $hit_length, $hsp_len, $coverage, $query_ctg) = @aa;
			
			# get return string for non-redundancy and redundancy status
			# return_string: n -- non-redundancy
			# return_string: r -- redundancy
			# put non-redundancy and redundancy data to array inset, restset

			my $return_string = ifRedundant(\@inset, \$query_ctg, $hit_id, $hit_length, $diff_ratio, $diff_contig_cover, $diff_contig_length, \%contig_position);

			if ($return_string eq "n")  { push(@inset, $tr->[2]); }
			else { push(@restset, $tr->[2]); }
		}
	}

	my $output_identified = '';
	foreach my $tr (@inset) {  $output_identified.= $tr."\n"; }

	return ($output_identified);
}

# ifRedundant: check if the datase4t is redundant
sub ifRedundant
{
	my ($inset, $query_ctg, $hit_id, $hit_length, $diff_ratio, $diff_contig_cover, $diff_contig_length, $contig_position) = @_;
	
	my @query_contigs = split(/,/, $$query_ctg);
	my $total_contigs = scalar(@query_contigs);
	my $ratio;

	# scan previous stored record (@inset) to check redundancy
	foreach my $tr ( @$inset )
	{
		# put previous stored contig to hash for non-redundancy data
		my @aa = split(/\t/, $tr);
		my @contigs= split(/,/, $aa[4]);
		my %contigs;
		foreach my $each_contig ( @contigs ) { $contigs{$each_contig} = 1; }

		# check query contings in current record
		my @hit_range;
		my $diff_contigs = 0;

		foreach my $each_contig ( @query_contigs )
		{
			unless ( defined $contigs{$each_contig} ) 
			{
				$diff_contigs++;
				if ( defined $$contig_position{$each_contig."\t".$hit_id} )
				{
					push(@hit_range, $$contig_position{$each_contig."\t".$hit_id});
				}
				else
				{
					die "[ERROR]Undef Blast for pair of $each_contig $hit_id\n";
				}
			}
		}

		# ratio = different contig in current record / total contig in previous record
		$ratio = $diff_contigs * 1.0 / $total_contigs;

		# compute length of different contig in current record covered for hit range
		my %pos;
		foreach my $range (@hit_range)	
		{
			my ($start, $end) = split(/\t/, $range);
			for(my $i=$start; $i<=$end; $i++)
			{
				$pos{$i}++;
			}
		}
		my $length_of_diff_contigs_cover = scalar(keys(%pos));

		my $cover_of_diff_contigs = $length_of_diff_contigs_cover / $hit_length;
		
		if( $ratio <= $diff_ratio && $cover_of_diff_contigs <= $diff_contig_cover && $length_of_diff_contigs_cover <= $diff_contig_length )
		{ return "r"; }

	}

	return "n";
}

=head2

 pileup_depth -- get read depth of contig from pileup file

=cut
sub pileup_depth
{
	my $input_pileup = shift;

	# put mean depth, and total coverage to hash
	# key: ref_id, mean, total. cover
	# value: mean, total, cover
	my %depth;

	my $prev_ref = '';	
	my $mean_depth;
	my $total_depth = 0;
	my $total_position = 1;

	open(IN, $input_pileup) || die $!;
	while (<IN>) 
	{
		chomp;
		next if $_ =~ m/^#/;
		# SL2.40ch00      23059   T       1       ^S.     I
		my @a = split(/\t/, $_);
		die "[ERR]pileup line\n" if scalar @a < 4;
		next if scalar @a < 6;
		my ($ref, $pos, $base, $depth, $align, $qual) = @a;

		if ($ref ne $prev_ref) 
		{
                	if ($prev_ref)
			{
				$mean_depth = 1.0*$total_depth/$total_position;
				$depth{$prev_ref}{'mean'} = $mean_depth;
				$depth{$prev_ref}{'total'} = $total_depth;
				$depth{$prev_ref}{'cover'} = $total_position;
			}	
			$total_depth = $depth;		# init total depth
			$total_position = 1; 		# init total position
		}
		else
		{
			$total_depth+=$depth;
			$total_position++;
		}

		$prev_ref = $ref;	
	}
	
	$mean_depth = 1.0*$total_depth/$total_position;
	$depth{$prev_ref}{'mean'} = $mean_depth;
	$depth{$prev_ref}{'total'} = $total_depth;	
	$depth{$prev_ref}{'cover'} = $total_position;

	close(IN);

	return %depth;
}

=head2

 blast_table_to_sam -- convert blast table information to sam format

=cut
sub blast_table_to_sam
{
	my $blast_table = shift;

	# put all hit name and length to hash
	my %name_len;

	chomp($blast_table);
	my @a = split(/\n/, $blast_table);
	foreach my $line (@a)
	{
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		my ($hit_name, $hit_len) = ($ta[2], $ta[3]);
		defined $name_len{$hit_name} or $name_len{$hit_name} = $hit_len;
	}

	my $output_sam = '';

	# generate header
	foreach my $key (keys(%name_len))
	{
		$output_sam.="\@SQ\tSN:".$key."\tLN:".$name_len{$key}."\n";
	}

	my ($qlen, $slen, $q, $s, $qbeg, $qend, @sam, @cigar, @cmaux);
	@sam = (); 
	@sam[4,6,7,8,10] = (255, '*', 0, 0, '*'); # fixed number in sam output

	foreach my $line (@a)
	{
		next if $line =~ m/^#/;
		my @cols= split(/\t/,$line);

		$sam[0] = $cols[0]; 							# query ID = read ID
		if($cols[8] eq "-1"){ $sam[1] = 0x10; } else{ $sam[1] = 0x00; }		# strand
		$sam[2] = $cols[2]; 							# hit ID = ref ID
		$qlen = $cols[1]; 
		$qend = $cols[10]; 
		$sam[3] = $cols[11]; 							# hit start = pos

		@cigar = ();
		@cmaux = (0, 0);
		my @q_seqs= split(/,/,$cols[14]); 	# query seq
		my @s_seqs= split(/,/,$cols[15]); 	# hit seq
		my $l = scalar(@q_seqs);		# 

		$sam[9] = '';				# init query seq

		for (my $i = 0; $i < $l; ++$i)
		{
			$q_seqs[$i] =~ /Query\:\s(\d+)\s*(\S+)\s(\d+)/;
			$q = $2;
			my $x = $q;
			$x =~ s/-//g; 
			$sam[9] .= $x;				# connect query seq from blast result
		
			$s_seqs[$i] =~ /Sbjct\:\s(\d+)\s*(\S+)\s(\d+)/;
			$s = $2;
			aln2cm(\@cigar, \$q, \$s, \@cmaux);	# get cigar from query and hit seq
		}	

		# additional info
		my ($as, $ev) = (int($cols[7] + .499), $cols[6]);
		$ev = "1$ev" if ($ev =~ /^e/);
		@sam[11,12] = ("AS:i:$as", "EV:Z:$ev");
		my $sam_info = blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
		$output_sam .= $sam_info;
	}

	return $output_sam;
}

sub blast_print_sam 
{
	my ($sam, $cigar, $cmaux, $qrest) = @_;
  
	push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
  	push(@$cigar, $qrest . 'H') if ($qrest);
  	if ($sam->[1] & 0x10) 
	{
		@$cigar = reverse(@$cigar);
		$sam->[9] = reverse($sam->[9]);
		$sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
	}
	$sam->[9] = '*' if (!$sam->[9]);
	$sam->[5] = join('', @$cigar);
	my $sam_info = join("\t", @$sam)."\n";
	return $sam_info;
}

sub aln2cm 
{
	my ($cigar, $q, $s, $cmaux) = @_;
  	my $l = length($$q);
  
	for (my $i = 0; $i < $l; ++$i) 
	{
		my $op;
		# set $op
		if (substr($$q, $i, 1) eq '-') { $op = 2; }	# Query '-' should be "D"
		elsif (substr($$s, $i, 1) eq '-') { $op = 1; }	# Hit '-' should be "I"
		else { $op = 0; }				# Otherwise "M"

		# for CIGAR
		if ($cmaux->[0] == $op) {
	  		++$cmaux->[1];	# MID count
		} else {
	  		push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
	  		$cmaux->[0] = $op; $cmaux->[1] = 1;
		}
  	}
}

=head2
 plot_select -- plot the select contigs which sRNA enriched in 21-22nt
=cut
sub plot_select
{
	my ($select_ctg, $label_ctg, $contig_best_blast, $virus_info, $map_sRNA_len_stat, $output_dir, $output_prefix) = @_;

	foreach my $cc (sort keys %$contig_best_blast) {
		print $cc."\t".$$contig_best_blast{$cc}."\n";
	}

	my $output_file = "$output_dir/$output_prefix.html";

	my $out_table;
	$out_table = qq'<!-- Top 6 lines shoud be removed when load html into virome database -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
<script src="http://libs.baidu.com/jquery/2.1.3/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
<script> \$(function() { \$(\'a\').tooltip(); }) </script>
<br><br><br><br><br>

<table class="table" style="font-size:12px; width: 780px;" align=center>
 <tr><td>* potential novel viruses do not have sequence similarity with known viruses</td></tr>
</table>

<table class="table table-bordered table-condensed" style="font-size:12px; width: 780px;" align=center>
  <tr bgcolor=#e2e8ec>
    <th>Contig ID</th>
    <th>Length</th>
';
	for (15 .. 40) {
		$out_table.= qq'<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="sRNA length">$_</a></th>\n';
	}
	$out_table.= "<th>BLAST</th>";
	$out_table.= "</tr>\n";

	my ($tr_style, $td_style);

	foreach my $cid (sort {$$select_ctg{$b}<=>$$select_ctg{$a}} keys %$select_ctg) {
		my $ctg_len = $$select_ctg{$cid};

		if (defined $$label_ctg{$cid}) {
			$tr_style = 'style="font-weight:bold"';
			$td_style = '';
		} else {
			$tr_style = '';
			$td_style = '';
		}

		$out_table.= "<tr $tr_style><td>$cid";
		$out_table.= "*" if defined $$label_ctg{$cid};
		$out_table.= "</td><td>$ctg_len</td>";
		for(15 .. 40) {
			my $n = 0;
			$n = $$map_sRNA_len_stat{$cid}{$_} if defined $$map_sRNA_len_stat{$cid}{$_};
			if ($_ == 21 || $_ == 22) {
				$out_table.= "<td bgcolor=\"#e6ccb3\">$n</td>";
			} else {
				$out_table.= "<td>$n</td>";
			}
		}

		# add virus information
		if (defined $$contig_best_blast{$cid}) {
			my $vid = $$contig_best_blast{$cid};
			my $genus = $$virus_info{$vid}{'genus'};
			my $desc = $$virus_info{$vid}{'genus'};
			my $vir_desc = "[$vid] $desc";
			$out_table.= qq'<td><a title="" data-placement="top" data-toggle="tooltip" data-original-title="$vir_desc">$genus</a></td>\n';
		} else {
			$out_table.= "<td>NA</td>\n";
		}

		$out_table.= "</tr>\n";
	}
	$out_table.= qq'</table>';

    # output table to file
    save_file($out_table, $output_file);
}
=head2
 plot_result -- plot the result
=cut
sub plot_result
{
	my ($known_identified, $known_contig_blast_table, $output, $type) = @_;

	# parse blast result to get contig_length and match length
	# push contig and match info to hash
	# key 1: contig # reference virus ID
	# 	key 2: length
	# 		   match
	# 		   iden
	my %contig_match;
	chomp($known_contig_blast_table);
	my @bs = split(/\n/, $known_contig_blast_table);
	foreach my $bs (@bs) {
		my @line = split(/\t/, $bs);
		if ($line[13] =~ m/(\d+)\/(\d+)\(\d+%\)/) {
			$contig_match{$line[0]."#".$line[2]}{'match'}  = $1;
			$contig_match{$line[0]."#".$line[2]}{'length'} = $2;
			$contig_match{$line[0]."#".$line[2]}{'iden'}   = $1/$2;
		}
		else {
			die "[ERR]CONTIG $bs\n";
		}
	}

	my $OUTPUT_DIR  = "$output/$type"."_references";
	mkdir $OUTPUT_DIR unless -e $OUTPUT_DIR;
	my $output_file = "$output/$type.html"; 	

	my %all_hits;
	my %out;
	my %index;

	my $out_table = '';
	my $line_number=0;

	# save contig and virus to hash
	# key: contig ID
	# value: virus accession ID
	my %contig_virus;
	chomp($known_identified);
	my @a = split(/\n/, $known_identified);
	foreach my $line (@a)
	{
		next if $line =~ m/^#/;
		$line_number++;

		my @ta = split(/\t/, $line);
		# 0 [col1] - virus seq ID
		# 1 [col2] - virus seq length
		# 2 [col3] - virus seq cover length
		# 3 [col4] - coverage (cover_length/seq_length: col3/col2)
		# 4 [col5] - contigs
		# 5 [col6] - contig number
		# 6 [col7] - raw depth (read depth, from pileup file)
		# 7 [col8] - normalized depth
		# 8 [col9] - genus
		# 9 [col10]- desc

		die "[ERR]Num of cols: $line\n" unless scalar @ta == 10;
		my ($ref, $len, $cov_base, $coverage, $contig, $contig_num, $raw_depth, $norm_depth, $genus, $desc) = @ta;

		my ($match, $len_c, $max_iden, $min_iden) = (0, 0, 0, 100);
		my @c = split(/,/, $ta[4]);
		foreach my $cid (@c) {
			my $key = $cid."#".$ref;
			$match += $contig_match{$key}{'match'};
			$len_c += $contig_match{$key}{'length'};
			my $iden_c = $contig_match{$key}{'iden'};
			$max_iden = $iden_c if $iden_c > $max_iden;
			$min_iden = $iden_c if $iden_c < $min_iden;
		}
		my $iden = sprintf("%.2f", $match/$len_c*100);
		$max_iden = sprintf("%.2f", $max_iden*100);
		$min_iden = sprintf("%.2f", $min_iden*100);
		$iden = 100 if $iden == 100;
		$max_iden = 100 if $max_iden == 100;
		$min_iden = 100 if $min_iden == 100;

		$all_hits{$ref} = 1; 
		$index{$ref} = $line_number;

		my $link = $type."_references/".$ref.".html";
		$coverage = sprintf("%.3f", $coverage) * 100;
		$raw_depth = sprintf("%.1f", $raw_depth);
		$norm_depth = sprintf("%.1f", $norm_depth);

		$out_table .= qq' 
  <tr>
	<td><a href="$link">$ref</a></td>
	<td>$len</td>
	<td>$cov_base ($coverage)</td>
	<td>$contig_num</td>	
	<td>$raw_depth</td>
	<td>$norm_depth</td>
	<td>$iden</td>
	<td>$max_iden</td>
	<td>$min_iden</td>
	<td>$genus</td>	
	<td width="50%">$desc</td>	
  </tr>';	
	}

	$out_table = qq'<!-- Top 6 lines shoud be removed when load html into virome database -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
<script src="http://libs.baidu.com/jquery/2.1.3/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
<script> \$(function() { \$(\'a\').tooltip(); }) </script>
<br><br><br><br><br>

<table class="table table-bordered table-condensed" style="font-size:12px; width: 780px;" align=center>
  <tr bgcolor=#e2e8ec>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="GenBank accession number of reference sequence to which 
significant similarity was found with contigs assembled from the sample">Reference</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="Length of the reference sequence identified">Length</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="percentage of the GenBank reference sequence that is
covered by alignment of contigs assembled from the sample">Coverage (%)</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="Number of contigs, assembled from the sample, that can be
aligned to the reference sequence">#contig</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="The average total number of times each nucleotide of the
reference sequence is covered by sequences identified in the sample">Depth</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="The average number of times each nucleotide of the GenBank
reference sequence is covered by aligned sequences identified in the sample per million of total sequenced reads (small RNAs)">Depth (Norm)</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="The average percentage of nucleotide identity to the GenBank
reference sequence of all contigs aligned to that sequence">%Identity</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="Maximum percentage of sequence identity of the
contigs with significant similarity to the GenBank reference sequence">%Iden Max</a></th>
	<th><a title="" data-placement="top" data-toggle="tooltip" data-original-title="Minimum percentage of sequence identity of the
contigs with significant similarity to the GenBank reference sequence">%Iden Min</a></th>
	<th>Genus</th>
	<th>Description</th>
  </tr>
$out_table
</table>';

	# output table to file
	save_file($out_table, $output_file);

	# put blast information to hash
	# key: hit ID
	# value: line of contig blast to array
	chomp($known_contig_blast_table);
	my @b = split(/\n/, $known_contig_blast_table);
	foreach my $line (@b)
	{
		next if $line =~ m/^#/;
		my @ta = split(/\t/, $line);
		die "[ERR]Undef hit: $ta[2]\n" unless defined $all_hits{$ta[2]};
		push(@{$out{$ta[2]}}, $line);
	}
	
	foreach my $hit (sort { $index{$a} <=> $index{$b} } keys %out)
	{
		my @lines = @{$out{$hit}};
		my $out_hit = $OUTPUT_DIR."/$hit.html";

		open(OUT1, ">$out_hit") or print "you forgot to make $OUTPUT_DIR";
		print OUT1 qq'<!-- Top 5 lines shoud be removed when load html into virome database -->
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
<br><br><br>
		
<div align=center>
';
		draw_img(\$hit, \@lines, \*OUT1, $OUTPUT_DIR);
		print OUT1 qq'</div><br>';
	
		my $table_info = output_aligments(\@lines);
		print OUT1 $table_info;
		close(OUT1);

		# code for checking
		# foreach my $each_line (@lines){ print $each_line."\n"; }
	}
}

# subroutine
sub draw_img
{
	my ($hit_id, $align_info, $out, $OUTPUT_DIR) = @_;
	my @cols = split(/\t/, @$align_info[0]);
	my $full_length = $cols[3];

	my $panel = Bio::Graphics::Panel->new(
		-length => $full_length,
		-width  => 610,
		-pad_left => 25,
		-pad_right => 25,
	);
       
	my $full_length_seq = Bio::Graphics::Feature->new(
		-start => 1,
		-end => $full_length,
		-display_name => "$$hit_id",
	);
       
	$panel->add_track( 
		$full_length_seq,
		-glyph => 'arrow',
		-tick => 2,
		-fgcolor => 'black',
		-double => 1,
	);

	$panel->add_track(
		$full_length_seq,
		-glyph => 'generic',
		-bgcolor => 'blue',
		-label => 1,
		-link => "http://www.ncbi.nlm.nih.gov/nuccore/$$hit_id",
		-target => '_blank'
	);

	my $track = $panel->add_track(
		-glyph => 'graded_segments',
		-label => 1,
		-bgcolor => 'red',
		-sort_order  => 'high_score',
		-min_score   => 80,
		-max_score   => 100, # identity for score
		-key         => 'domains',
		-link        => '#$id',
		-description => sub {
			my $feature = shift;
			return unless $feature->has_tag('description');
			my ($description) = $feature->each_tag_value('description');
			"$description";
		},
	);

	# update for multiple HSP
	# 1. add fist record to $last_query_id and hsp
	my @align_info = @$align_info;
	my $first_hsp = shift @align_info;
	my @fmm = split(/\t/, $first_hsp);
	my $last_query_id = $fmm[0];
	my @hsp = ();
	push(@hsp, [$fmm[11], $fmm[12], $fmm[5]]);

	# 2. scan other hsp from 2nd
	foreach my $each_hsp (@align_info) 
	{
		my @mm = split(/\t/, $each_hsp);
		my $query_id = $mm[0];
		my $score = $mm[5];
		my $start = $mm[11];
		my $end   = $mm[12];

		if ($query_id ne $last_query_id)
		{
			# add hit feature
			my $feature = Bio::SeqFeature::Generic->new(
				-seq_id => $last_query_id,
				#-name	=> $last_query_id,
				#-score	=> $mm[5], # identify for score
				#-start	=> $mm[11],
				#-end	=> $mm[12],
				-primary_id	=> $last_query_id,
			);

			# add hsp
			foreach my $sub (@hsp) {
				my $subfeature = Bio::SeqFeature::Generic->new(
					-start 	=> $sub->[0],
					-end	=> $sub->[1], 
					-score 	=> $sub->[2]
				);
				$feature->add_sub_SeqFeature($subfeature,'EXPAND');
			}
			$track->add_feature($feature);
	
			# clean and add hsp
			@hsp = ();
			push(@hsp, [$start, $end, $score]);
			$last_query_id = $query_id;
		} 
		else 
		{
			push(@hsp, [$start, $end, $score]);
		}
	}

	# 3. add last hsp
	if (@hsp > 0 && $last_query_id) {
		# add hit feature
		my $feature = Bio::SeqFeature::Generic->new(
			-seq_id => $last_query_id,
			#-name   => $last_query_id,
			#-score => $mm[5], # identify for score
			#-start => $mm[11],
			#-end   => $mm[12],
			-primary_id => $last_query_id,
		);
            
		# add hsp
		foreach my $sub (@hsp) {
			my $subfeature = Bio::SeqFeature::Generic->new(
				-start  => $sub->[0], 
				-end    => $sub->[1],
				-score  => $sub->[2]
			);
			$feature->add_sub_SeqFeature($subfeature,'EXPAND');
		}
		$track->add_feature($feature);
	}
	
	my ($url,$map,$mapname) = $panel->image_and_map(
		-root => $OUTPUT_DIR,
		-url=> "$$hit_id",
	);
						
	print $out qq(<img src="$url" border=0 usemap="#$mapname">),"\n";
	print $out $map;
}

sub output_aligments
{
	my $align_info = shift;
	# produce domains alignment for one protein
	my $break_length = 95;
	
	# table title
	my $out_table = qq'
<table class="table table-bordered table-condensed" style="font-size:12px; width: 700px;" align=center>
  <tr bgcolor=#e2e8ec>
    <th>Order</th>
    <th>Query ID</th>
    <th>Query Start</th>
    <th>Query End</th>
    <th>Subjct Start</th>
    <th>Subjct End</th>
    <th>Identity</th>
    <th>E value</th>
    <th>Strand</th>
  </tr>
';

	my $k = 0;

	foreach my $each_contig (@$align_info)
	{
		$k++;

		my @cols = split(/\t/, $each_contig); 
		#print $each_contig."\n";		
		my @aligned_query = split(/,/, $cols[14]);
		my @aligned_hit = split(/,/, $cols[15]);
		my @aligned_string = split(/,/, $cols[16]);
		my $i=0;
		my $align;
		foreach my $each_seg (@aligned_query)
		{
			#print $each_seg."\n";			
			$align.="<span style=\"background-color:#99FFFF\";>$each_seg</span><br>";
			$align.="<span style=\"background-color:#99FFFF\";>$aligned_string[$i]</span><br>";
			$align.="<span style=\"background-color:#99FFFF\";>$aligned_hit[$i]</span><br> <br>";
			$i++;
		}

		$out_table .= qq' 
   <tr>
	<td>$k</td>
        <td><a name=$cols[0]>$cols[0]</a></td>
        <td>$cols[9]</td>
        <td>$cols[10]</td>
        <td>$cols[11]</td>
        <td>$cols[12]</td>
        <td>$cols[13]</td>
	<td>$cols[6]</td>
	<td>$cols[8]</td>
   </tr>
   <tr>
        <td colspan=9>
	<div style="padding:5px;">
        <b>Alignment:   </b><br>
	<span style="font-size:12px; color:#000; font-face: courier">
<pre>
$align
</pre>
	</span>
        </td>
	</div>
   </tr>';
	}

	$out_table.= qq'</table>';

	return $out_table;
}

1;
