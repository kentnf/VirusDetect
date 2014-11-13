#!/usr/bin/perl

=head
 viral_DB_prepare.pl -- create viral database from genbank viral sequnece
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Net::FTP;
use IO::File;

my $usage = qq'
USAGE: 

1. download viral sequnece from genbank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

';

# download viral sequence file from genbank
my $genbank_ftp = "ftp.ncbi.nih.gov";
my $ftp = Net::FTP->new($genbank_ftp, Debug=>0) || die "[ERR]Cannot connect to $genbank_ftp $@\n";
$ftp->login("anonymous",'-anonymous@') || die "[ERR]Cannot login ", $ftp->message, "\n";
$ftp->cwd("/genbank") || die "[ERR]Cannot change working directory ", $ftp->message;
my @vrl_files;
my @files = $ftp->ls();
foreach my $f (@files) {
	if ($f =~ m/gbvrl11\.seq\.gz/) {
		$ftp->get($f) unless -s $f;
		push(@vrl_files, $f);
		print "== $f has been download from genbank ==\n";
	}
}
$ftp->quit;

# load taxonomy information to hash for 
# 1. tracking host to division
# 2. tracking genus of virus
# description of hash
# hash division: key: taxon id, value division name
# hash node_table: key: taxon id, parent, rank, division; value: parent taxon id, rank, and division
# hash name_table: key: name; value: taxon id;
# hash virus_genus_taxon_id: key: taxon_id; value: genus scientific name
my $node_file = 'nodes.dmp';
my $name_file = 'names.dmp';
my $division_file = 'division.dmp';
my %division;
my $fh1 = IO::File->new($division_file) || die $!;
while(<$fh1>)
{
	chomp;
	next if $_ =~ m/^#/;
	my @a = split(/\t\|\t/, $_);
	$division{$a[0]} = $a[2];
	#print "s $a[0] es $a[2] e\n";
}
$fh1->close;
my ($node_table, $virus_genus_taxon_id) = load_taxon_node($node_file);
my ($name_table, $virus_genus_taxon) = load_taxon_name($name_file, $virus_genus_taxon_id);
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
	# think about virus may infected more than 1 divison
	$virus_hostdiv{$a[1]} = $a[0];
}
$fh2->close;
print "== ", scalar(keys(%virus_hostdiv)), " virus genus host info has been load into hash ==\n";

# parse viral sequence to database file
my $vrl_info = "vrl_genbank.info";
my $vrl_origin = "vrl_origin";

# generate viral sequence with fasta format
my %vrl_seq_info;

if (-s $vrl_origin) 
{
	print "[ERR]viral original sequence file exist $vrl_origin\n";
} 
else
{
	my $out1 = IO::File->new(">".$vrl_info) || die $!;

	foreach my $f (@vrl_files)
	{
		my $f_head = `less $f | head`;
		my $seq_num = 70000;
		if ($f_head =~ m/(\d+) loci/) { $seq_num = $1; }
		my $parse_seq_num = 0; 
		my $pre_parse_percent;

		my $seqin = Bio::SeqIO->new(-format => 'GenBank', -file => "gunzip -c $f |" );
		while(my $inseq = $seqin->next_seq)
		{
			#print $inseq->id."\n";
			#print $inseq->seq_version."\n";
			#print $inseq->desc."\n";
			$vrl_seq_info{$inseq->id}{'ver'} = $inseq->seq_version;
			$vrl_seq_info{$inseq->id}{'seq'} = $inseq->seq;
			$vrl_seq_info{$inseq->id}{'des'} = $inseq->desc;
			
			my %host_division = (); # key: division id; value: source	
			foreach my $feat_object ($inseq->get_SeqFeatures) { 
				my $primary_tag = $feat_object->primary_tag;
				next unless $primary_tag eq "source";
				foreach my $tag ($feat_object->get_all_tags) { 
 					foreach my $value ($feat_object->get_tag_values($tag)) {
						if ($tag eq 'db_xref' && $value =~ m/taxon:(\d+)/) {
							my $org_taxon = $1;
							$vrl_seq_info{$inseq->id}{'taxon'} = $org_taxon;
							# convert the taxon to viral genus, then tracking host through host_info.txt
							# print "$org_taxon\t$$node_table{$org_taxon}{'division'}\t$division{$$node_table{$org_taxon}{'division'}}\n";

							my $round = 0;
							my $genus_taxon = $org_taxon;
							while(1) {
								if (defined $$node_table{$genus_taxon}{'parent'}) {
									$genus_taxon = $$node_table{$genus_taxon}{'parent'};
									if (defined $$node_table{$genus_taxon}{'rank'} && $$node_table{$genus_taxon}{'rank'} eq 'genus') { last; }
									last if $genus_taxon == 1;
								} else {
									last;
								}
							}

							# print $genus_taxon,"\t", $$virus_genus_taxon{$genus_taxon},"\n"; exit;
							my $genus_name = $$virus_genus_taxon{$genus_taxon} if defined $$virus_genus_taxon{$genus_taxon};

							if (defined $virus_hostdiv{$genus_name}) {
								$host_division{$virus_hostdiv{$genus_name}} = 'ICTV';
							}
						} 
						if ($tag eq 'host' && defined $$name_table{$value} ) { 
							my $host_taxon = $$name_table{$value};
							$vrl_seq_info{$inseq->id}{'host_name'} = $value;
							$vrl_seq_info{$inseq->id}{'host_taxon'} = $host_taxon;

							##### get divison of host taxon #####
							my $host_div = "NA";
							$host_div = $$node_table{$host_taxon}{'division'} if defined $$node_table{$host_taxon}{'division'};

							##### use the host div from feature attribute as host division #####
							if ( defined $host_division{$host_div} ) {
								$host_division{$host_div}.= ",genbank";
							} else {
								$host_division{$host_div} = "genbank";
							}
							#print "$host_taxon\t$value\t$host_div\t$division{$host_div}\n";

							##### tracking host to high level #####
							#my $round = 0;
							#my $this_tax = $host_taxon;
							#while(1) {
							#	if (defined $node_table{$this_tax}) {
							#		$this_tax = $node_table{$this_tax};
							#		$round++;
							#		last if $round > 200;
							#		last if $this_tax == 33208;
							#	} else {
							#		last;
							#	}	
							#}
							#$vrl_seq_info{$inseq->id}{'host_kingdom'} = $this_tax;
							# print "$host_taxon\t$value\t$this_tax\n";
						}      
      					}          
   				}       
			}

			$vrl_seq_info{$inseq->id}{'host_div'} = \%host_division;

			foreach my $host_div_id (sort keys %host_division)
			{
				my $host_name = $division{$host_div_id};
				my $host_source = $host_division{$host_div_id};
				print $out1 $inseq->id,"\t",$inseq->length,"\t",$inseq->desc,"\t",$inseq->version,"\t",$host_name,"\t",$host_source,"\n";
			}

			$parse_seq_num++;
			my $parse_percent = int(($parse_seq_num/$seq_num) * 100);
			if ($parse_percent % 5 == 0 && $parse_percent ne $pre_parse_percent) {
				print "......$parse_percent......\n";
			}
			$pre_parse_percent = $parse_percent;
		}
	}

	$out1->close;
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
		$name = $a[1];
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
