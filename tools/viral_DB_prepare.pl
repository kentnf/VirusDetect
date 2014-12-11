#!/usr/bin/perl

=head
 viral_DB_prepare.pl -- create viral database from genbank viral sequnece
 ICTV virus list (http://talk.ictvonline.org/files/tags/default.aspx)
 virus and host (http://www.mcb.uct.ac.za/tutorial/ICTV%20Species%20Lists%20by%20host.htm)
 then get relationship between host and virus: host.info, Algae, Archaea, Bacteria, Fungi, Invertebrata, Protozoa, Plant, Protozoa, Vertebrata
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Net::FTP;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }
if	($options{'t'} eq 'download')	{ vrl_download(\%options, \@ARGV); }    # parse multi dataset
elsif	($options{'t'} eq 'category')	{ vrl_category(\%options, \@ARGV); }    # parse multi dataset
else	{ usage($version); }

sub usage
{
	my $version = shift;
	my $usage = qq'
PIPELINE (version $version):
1. download viral sequnece from genbank ftp (ftp://ftp.ncbi.nih.gov/genbank/)
   \$ $0 -t download
2. run viral_DB_prepare.pl script to category virus 
   \$ $0 -t category gbvr*.gz

';
	print $usage;
	exit;
}

=head2
 vrl_download : download vrl database, just generate download command using wget
=cut
sub vrl_download
{
	my ($options, $files) = @_;
	my $download_cmd = '';
	my $download_cmd_all = '';
	my $genbank_ftp = "ftp.ncbi.nih.gov";
	my $ftp = Net::FTP->new($genbank_ftp, Debug=>0) || die "[ERR]Cannot connect to $genbank_ftp $@\n";
	$ftp->login("anonymous",'-anonymous@') || die "[ERR]Cannot login ", $ftp->message, "\n";
	$ftp->cwd("/genbank") || die "[ERR]Cannot change working directory ", $ftp->message;
	my @vrl_files;
	my @files = $ftp->ls();
	foreach my $f (@files) {
		if ($f =~ m/gbvrl\d+\.seq\.gz/) {
			#$ftp->get($f) unless -s $f;
			#push(@vrl_files, $f);
			if ( -s $f ) {
				print "# == $f has been download from genbank ==\n";
			} else {
				$download_cmd.="wget ftp://ftp.ncbi.nih.gov/genbank/$f\n";
			}
			$download_cmd_all.="wget ftp://ftp.ncbi.nih.gov/genbank/$f\n";
		}
	}
	$ftp->quit;

	print "# cmd for download vrl database\n$download_cmd\n\n";
	print "# cmd for download all vrl database files\n$download_cmd_all\n\n";
	exit;
}

=head2
 vrl_category : 
=cut
sub vrl_category
{
	my ($options, $files) = @_;	

	my $usage = qq'
USAGE: $0 -t category input_file1 ... input_fileN

';

	print $usage and exit unless defined $$files[0];
	my @vrl_files = @$files;

	# load taxonomy information to hash for 
	# # 1. tracking host to division
	# # 2. tracking genus of virus
	# # description of hash
	# # hash division: key: taxon id, value division name
	# # hash node_table: key: taxon id, parent, rank, division; value: parent taxon id, rank, and division
	# # hash name_table: key: name; value: taxon id;
	# # hash virus_genus_taxon_id: key: taxon_id; value: genus scientific name

	my $node_file = 'nodes.dmp';
	my $name_file = 'names.dmp';
	my %division  = get_division();
	my ($node_table, $virus_genus_taxon_id) = load_taxon_node($node_file);
	my ($name_table, $virus_genus_taxon) = load_taxon_name($name_file, $virus_genus_taxon_id);
	my %abnormal_host_name = load_abnormal_host_name('abnormal_host_name.txt');
	print "== ", scalar(keys(%$virus_genus_taxon_id)), " virus genus taxon has been load into hash ==\n";
	print "== ", scalar(keys(%$virus_genus_taxon)), " virus genus taxon has been load into hash ==\n";
	print "== taxon info has been load into hash ==\n";

	# debug: generate virus genus id and name
	# foreach my $tid (sort keys %$virus_genus_taxon) { print $tid."\t".$$virus_genus_taxon{$tid}."\n"; } exit;

	# load host info to hash (ICTV-Master-Species-List-2013_v2: http://www.mcb.uct.ac.za/tutorial/ICTV%20Species%20Lists%20by%20host.htm )
	my $host_info_file = 'host_info.txt';
	my %virus_hostdiv; # key: virus genus name, value: host divition num
	my $fh2 = IO::File->new($host_info_file) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		next if $a[1] eq 'NA';
		# think about virus may infected more than 1 divison
		my $gname = lc($a[0]);
		$virus_hostdiv{$gname} = $a[1];
	}
	$fh2->close;
	print "== ", scalar(keys(%virus_hostdiv)), " virus genus host info has been load into hash ==\n";

	# parse viral sequence to database file
	my $vrl_info_classified   = "vrl_genbank.txt";
	my $vrl_info_unclassified = "vrl_genbank_unclassified.txt";
	my $vrl_seq_unclassified  = "vrl_genbank_unclassified.fasta"; 
	my $vrl_origin = "vrl_origin";

	# generate viral sequence with fasta format
	my %vrl_seq_info;
	
	#print "[ERR]viral original sequence file exist $vrl_info_classified\n" and exit if -s $vrl_info_classified;
	#exit;

	my $out1 = IO::File->new(">".$vrl_info_classified) || die $!;
	my $out2 = IO::File->new(">".$vrl_info_unclassified) || die $!;
	my $out3 = IO::File->new(">".$vrl_seq_unclassified) || die $!;

	#@vrl_files = ('sequence.gb.gz');
	foreach my $f (@vrl_files)
	{
		print "parsing file $f start ......\n";
		my $f_head = `less $f | head`;
		my $seq_num = 70000;
		if ($f_head =~ m/(\d+) loci/) { $seq_num = $1; }
		my $parse_seq_num = 0; 
		my $pre_parse_percent = 0;

		my $seqin = Bio::SeqIO->new(-format => 'GenBank', -file => "gunzip -c $f |" );
		while(my $inseq = $seqin->next_seq)
		{
			# print $inseq->id."\n".$inseq->seq_version."\n".$inseq->desc."\n";
			my $sid = $inseq->id;
			$vrl_seq_info{$sid}{'ver'} = $inseq->seq_version;
			$vrl_seq_info{$sid}{'seq'} = $inseq->seq;
			$vrl_seq_info{$sid}{'des'} = $inseq->desc;
						
			# get org taxon id or host name 
			my ($org_taxon, $org_name, $genus_taxon, $genus_name, $genus_div, $host_taxon, $host_name, $host_div);
			my %source; # key: org_taxon; value: host_name
			foreach my $feat_object ($inseq->get_SeqFeatures) 
			{ 
				my $primary_tag = $feat_object->primary_tag;
				next unless $primary_tag eq "source";

				$org_taxon = 'NA'; $host_name = 'NA';
				foreach my $tag ($feat_object->get_all_tags) { 
 					foreach my $value ($feat_object->get_tag_values($tag)) {
						if ($tag eq 'db_xref' && $value =~ m/taxon:(\d+)/) {
							if ( $org_taxon eq 'NA' ) { $org_taxon = $1; } else { print "[ERR]Repeat organism taxon id $sid\n"; exit; }
						}
 
						if ($tag eq 'host' ) { 
							my $normal_host_name = $value;
							$normal_host_name = $abnormal_host_name{$value} if defined $abnormal_host_name{$value};
							if ( $host_name eq 'NA' ) { $host_name = $normal_host_name; } else { $host_name.= "\t".$normal_host_name; } # one org <=> more host
						}      
      					}          
   				}
			
				unless (defined $$node_table{$org_taxon}{'division'} && $$node_table{$org_taxon}{'division'} eq '9') {
					my $changed_org_taxon = correct_org_taxon_division($sid);
					$org_taxon = $changed_org_taxon unless $changed_org_taxon eq 'NA';
				} 

				next unless $$node_table{$org_taxon}{'division'} eq '9'; # example EF710638
				$source{$org_taxon} = $host_name;
			}

			#$vrl_seq_info{$sid}{'taxon'} = $org_taxon;
			#$vrl_seq_info{$sid}{'host_name'} = $host_name;

			if (scalar(keys(%source)) < 1 ) {  
				print "STAT$f\t$sid\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
			}

			my %host_division = (); # key: division id; value: source
			foreach my $otid (sort keys %source) 
			{

				#=== find genus taxon id according organism taxon id, then genus name, then division ===
				my $genus_taxon_tracking = $otid;
                        	while(1) {
					if (defined $$node_table{$genus_taxon_tracking}{'parent'}) {
						$genus_taxon_tracking = $$node_table{$genus_taxon_tracking}{'parent'};
						if (defined $$node_table{$genus_taxon_tracking}{'rank'} && $$node_table{$genus_taxon_tracking}{'rank'} eq 'genus') { last; }
							last if $genus_taxon_tracking == 1;
					} else {
						last;
					}
				}

				if ($$node_table{$genus_taxon_tracking}{'rank'} eq 'genus' ) {
					$genus_taxon = $genus_taxon_tracking;
				} else {
					$genus_taxon = "NA";
				}

				$genus_name = 'NA';
				$genus_name = $$virus_genus_taxon{$genus_taxon} if defined $$virus_genus_taxon{$genus_taxon}; 

				$genus_div = 'NA';	# the host div tracking by org_genus -- org_taxon_id
				if ( defined $virus_hostdiv{$genus_name}) {
					my @gdiv = split(/,/, $virus_hostdiv{$genus_name});
					foreach $genus_div (@gdiv) { $host_division{"$genus_div"} = 'ICTV'; } # on organism may categorized into two division;
				}

				#=== find host name and host id, tracking host division ===
				my @hname = split(/\t/, $source{$otid}); 
				foreach my $hname (@hname)
				{
					$host_taxon = "NA"; $host_div = 'NA';
					# find host taxon according to host name
					my $host_taxon_fixed = "NA";
					my $lchname = lc($hname);
					if (defined $$name_table{$lchname}) {
						$host_taxon = $$name_table{$lchname};
					} else {
						# fixed the host name for tracking host division 
						my @hn = split(" ", $lchname);
						while(scalar(@hn) > 1) {
							pop @hn;
							my $sname = join(" ", @hn);
							chomp($sname);
							$sname=~ s/ $//;
							$sname=~ s/;$//;
							$sname=~ s/,$//;
							if (defined $$name_table{$sname}) {
								$host_taxon = $$name_table{$sname};
								$host_taxon_fixed = 'fixed';
								last;
							}
						}
					}

					# get host division, and put it to hash
					if ( defined $$node_table{$host_taxon}{'division'} ) {
						$host_div = "G".$$node_table{$host_taxon}{'division'};
						$host_div = "G10" if ($host_div eq "G2" || $host_div eq "G5" || $host_div eq "G6");
						if (  defined $host_division{$host_div} ) {
							if ( $host_division{$host_div} ne 'genbank' ) {
								$host_division{$host_div}.=",genbank";
							}
						} else {
							$host_division{$host_div} = "genbank";
						}
					}

					print "STAT$f\t$sid\t$otid\t$genus_taxon\t$genus_name\t$genus_div\t$hname\t$host_taxon\t$host_taxon_fixed\t$host_div\n";
				}
			}

			##### tracking host to high level #####
			#my $round = 0;
			#my $this_tax = $host_taxon;
			#while(1) {
			#       if (defined $node_table{$this_tax}) {
			#               $this_tax = $node_table{$this_tax};
			#               $round++;
			#               last if $round > 200;
			#               last if $this_tax == 33208;
			#       } else {
			#               last;
			#       }       
			#}
			#$vrl_seq_info{$inseq->id}{'host_kingdom'} = $this_tax;
			# print "$host_taxon\t$value\t$this_tax\n";

			if ( scalar(keys(%host_division)) > 0 ) {	
				$vrl_seq_info{$sid}{'host_div'} = \%host_division;

				foreach my $host_div_id (sort keys %host_division) {

					unless ( defined $division{$host_div_id} ) {
						print $host_div_id,"\n";
						print $sid,"\n";
						foreach my $aa (sort keys %host_division) {
							print $aa."\t".$host_division{$aa}."\n";
						}
						exit;
					}

					my $host_div_name = $division{$host_div_id};
					my $host_source = $host_division{$host_div_id};
					print $out1 $sid,"\t",$inseq->length,"\t",$inseq->desc,"\t",$inseq->version,"\t",$host_div_name,"\t",$host_source,"\n";
				}

			} else {
				$host_name =~ s/\t/,/ig;
				print $out2 $sid,"\t",$inseq->length,"\t",$inseq->desc,"\t",$inseq->version,"\tNA\t$host_name\n";
				print $out3 ">",$sid,"\n",$inseq->seq,"\n";
			}


			

			#$parse_seq_num++;
			#my $parse_percent = int(($parse_seq_num/$seq_num) * 100);
			#if ($parse_percent % 5 == 0 && $parse_percent ne $pre_parse_percent) {
			#	print "......$parse_percent......\n";
			#}
			#$pre_parse_percent = $parse_percent;
		}

		print print "parsing file $f end ......\n\n";
	}

	$out1->close;
	$out2->close;

	# output sequence for each division and uniq them by uclust
	foreach my $div_id (sort keys %division)
	{
		my $div_name = $division{$div_id};
		$div_name =~ s/ /_/ig;
		my $vrl_seq_file = "vrl_".$div_name."_all.fasta";
		my $vrl_seq = '';
		foreach my $id (sort keys %vrl_seq_info)
		{
			my $host_division = $vrl_seq_info{$id}{'host_div'} if defined $vrl_seq_info{$id}{'host_div'};
			foreach my $host_div_id (sort keys %$host_division)
			{
				if ($host_div_id eq $div_id)
				{
					$vrl_seq.=">$id\n$vrl_seq_info{$id}{'seq'}\n";
				}
			}
		}

		if ($vrl_seq) {
			my $fhs = IO::File->new(">".$vrl_seq_file) || die $!;
			print $fhs $vrl_seq;
			$fhs->close;
		}
		# cluster by uclust
	}

}

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 load_taxon_node
=cut
sub load_taxon_node
{
	my $node_file = shift;
	my %node_table;		  # key: taxon_id, parent, division, rank; value: parent, division, rank
	my %virus_genus_taxon_id; # key: taxon_id; value: 1
	my ($tax_id, $parent_tax_id, $rank, $division);
	my $fh = IO::File->new($node_file) || die $!;
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t\|\t/, $_);
		$tax_id = $a[0];
		$parent_tax_id = $a[1];
		$rank = $a[2];
		$division = $a[4];
		#print "s $tax_id es $parent_tax_id es $rank se $division e\n"; exit;
		$node_table{$tax_id}{'parent'} = $parent_tax_id;
		$node_table{$tax_id}{'division'} = $division;
		$node_table{$tax_id}{'rank'} = $rank;

		# get virus genus taxon id and 
		if ($division eq 9 && $rank eq 'genus') {
			$virus_genus_taxon_id{$tax_id} = 1;
		}
	}
	$fh->close;

	return (\%node_table, \%virus_genus_taxon_id);
}

=head2
 load_taxon_name
=cut
sub load_taxon_name
{
	my ($name_file, $virus_genus_taxon_id) = @_;
	my %name_table;		# key: name; value: taxon_id
	my %virus_genus_taxon;  # key: taxon_id; value: name;

	my ($tax_id, $name, $name_class);
	my $fh = IO::File->new($name_file) || die $!;
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t\|\t/, $_);
		$tax_id = $a[0];
		$name = lc($a[1]);
		$name_class = $a[3];
		$name_class =~ s/\t.*//;
		# print "s $tax_id es $name es $name_class e\n"; exit;
		$name_table{$name} = $tax_id;

		# notice, some virus genus using synonym name in ICTV table		
		if ( defined $$virus_genus_taxon_id{$tax_id} && $name_class eq 'scientific name' ) {
			if ( defined $virus_genus_taxon{$tax_id} ) {
				print "[ERR]Repeat virus genus: $tax_id $name\n";
			} else {
				$virus_genus_taxon{$tax_id} = $name;
			}
		}
	}
	$fh->close;
	return (\%name_table, \%virus_genus_taxon);
}

=head2
 get_division: division for both ICTV and genbank
=cut
sub get_division
{
	my %division = (
		'G0' =>  'Bacteria',
		'G1' =>  'Invertebrates',
		'G2' =>  'Vertebrates',
		'G3' =>  'Phages',
		'G4' =>  'Plants',
		'G5' =>  'Vertebrates',
		'G6' =>  'Vertebrates',
		'G7' =>  'Synthetic',
		'G8' =>  'Unassigned',
		'G9' =>  'Viruses',
		'G10' => 'Vertebrates',
		'G11' => 'EnvironmentalSamples',
		'X1' => 'Algae',
		'X2' => 'Archaea',
		'X3' => 'Fungi',
		'X4' => 'Protozoa'
	);
	return %division;
}

=head2
 correct_org_taxon_division:
 will fix it later
=cut
sub correct_org_taxon_division
{
	my $seq_id = shift;
	
	my %sid_org_taxon = (
		'KF178708' => 750074,
		'KF178710' => 750074,
		'KF178712' => 750074,
		'KC167159' => 1449545,
		'KC167160' => 1449545,
		'KC167161' => 1449545,
		'FJ654336' => 569708,
		'EU100683' => 461790,
		'EU100684' => 461791,
		'EU478797' => 1503438,
		'GU052197' => 352522,
		'GU052198' => 352522,
		'GU052199' => 352522,
		'GU052200' => 352522,
		'GU052201' => 352522,
		'GU052202' => 352522
	);

	my $taxon = 'NA';
	$taxon = $sid_org_taxon{$seq_id} if defined $sid_org_taxon{$seq_id};
	return $taxon;
}

=head2
 load_abnormal_host_name
=cut
sub load_abnormal_host_name
{
	my $input = shift;
	my %name;
	my $fh = IO::File->new($input) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		chomp($a[1]);
		$name{$a[0]} = $a[1] if (defined $a[1] && $a[1]=~ /\S+/);
	}
	return %name;
}
