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
use FindBin;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }
if	($options{'t'} eq 'download')	{ vrl_download(\%options, \@ARGV); }    # download dataset
elsif	($options{'t'} eq 'category')	{ vrl_category(\%options, \@ARGV); }    # automatically classification
elsif	($options{'t'} eq 'manually')	{ vrl_manually(\%options, \@ARGV); }	# manually classification
elsif	($options{'t'} eq 'patch')	{ vrl_patch(\%options, \@ARGV); }	# combine the result
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
3. generate manually classification file
   \$ $0 -t manually  
4. combine the result (patch the manually to update file)
   \$ $0 -t patch

';
	print $usage;
	exit;
}

=head2
 vrl_patch : patch the manually corrected result to update file 
=cut
sub vrl_patch
{
	# do not need any parameters, just path the result to update files

	# load default version
        my $version_file = $FindBin::RealBin."/default_version";
        my $version = `head -n 1 $version_file`; chomp($version);
        die "[ERR]version $version\n" if $version < 0;
	my $ver_num = $version;
        $version = "v$version";		
	die "[ERR]version mumber:$ver_num\n" if $ver_num < 0;	

	# check input files
	my $host_info_file = $FindBin::RealBin."/ICTV_2013_genus_uniq_host_info_raw.txt";
	my $update_hname = $FindBin::RealBin."/update_hname_table_$version.txt";
	my $update_genus = $FindBin::RealBin."/update_genus_table_$version.txt";
	my $update_desc  = $FindBin::RealBin."/update_desc_table_$version.txt";
		
	my $manual_hname = "manual_hname_table.txt";
        my $manual_genus = "manual_genus_table.txt";
        my $manual_desc  = "manual_desc_table.txt";

	foreach my $f (($host_info_file, $update_hname, $update_genus, $update_desc, $manual_hname, $manual_genus, $manual_desc )) {
		print "[ERR]file not exist $f\n" and exit unless -e $f;
	}	

	my %division  = get_division();
        my %virus_hostdiv = load_host_div($host_info_file, $update_genus);
        my %update_hname = load_update_hname($update_hname);
        my %update_desc  = load_update_desc($update_desc, \%division);

	# set output file
	$ver_num++;
	my $new_version = "v$ver_num";
	my $patch_hname = "update_hname_table_$new_version.txt";
	my $patch_genus = "update_genus_table_$new_version.txt";
	my $patch_desc  = "update_desc_table_$new_version.txt";

	my ($genus_patch_num, $hname_patch_num, $desc_patch_num) = (0, 0, 0);

	# patch genus classification 
	my $genus_cat = `cat $update_genus`;
	$genus_cat.= "# patch $new_version\n";
	my $in1 = IO::File->new($manual_genus) || die $!;
	while(<$in1>)
	{
		chomp;
		next if $_ =~ m/^#/;
		# format
		# genus_name      div,num,freq	  div,num,freq
		# cucumovirus     G4,16,0.94      G1,1,0.06
		my @a = split(/\t/, $_);
		my $genus = shift @a;
		next if defined $virus_hostdiv{$genus};
		my @div_id = ();
		my @div_name = ();
		foreach my $a (@a) {
			my @b = split(/,/, $a);
			die "[ERR]patch genus: $_\n" unless (scalar @b == 3);
			die "[ERR]patch genus div: $_\n" unless (defined $division{$$b[0]});
			push(@div_id, $b[0]);
			push(@div_name, $division{$$b[0]});
		}
		$genus_cat.= $genus. "\t" . join(",", @div_id) . "\t" . join(",", @div_name) . "\n";
		$genus_patch_num++;
	}
	$in1->close;

	my $out1 = IO::File->new(">".$patch_genus) || die $!;
	print $out1 $genus_cat;
	$out1->close;

	# patch host name
	my $host_name = `cat $update_hname`;
	$host_name.= "# patch $new_version\n";
	my $in2 = IO::File->new($manual_hname) || die $!;
	while(<$in2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		next unless defined $a[1];
		next if defined $update_hname{$a[0]};
		$host_name.= $_."\n";
		$hname_patch_num++;
	}
	$in2->close;

	my $out2 = IO::File->new(">".$patch_hname) || die $!;
	print $out2 $host_name;
	$out2->close;

	# patch manually classification
	my $manual_cat = `cat $update_desc`;
	$manual_cat.= "# patch $new_version\n";
	my $in3 = IO::File->new($manual_desc) || die $!;
	while(<$in3>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		next if $a[2] eq 'NA-CAT';
		$manual_cat.= $_."\n";	
		$desc_patch_num++;
	}
	$in3->close;

	my $out3 = IO::File->new(">".$patch_desc) || die $!;
	print $out3 $manual_cat;
	$out1->close;
	
	# report the patch
	print qq'
No. of patched genus: $genus_patch_num
No. of patched host name: $hname_patch_num
No. of manual classifiction: $desc_patch_num

!!! Notice !!!
Please copy patched file into $FindBin::RealBin
Please change version number of default_version
';
	
}

=head2
 vrl_manually : 
=cut
sub vrl_manually
{
	print qq'
Please manually check and update these files: 
	1. manual_hname_table.txt
	   one abnormal name per line, please search the abnormal name using google, genbank taxon database, and iciba 
	   to infer the correct name, this correct name must be include in genbank taxon database (in name.dmp.gz), 
	   then append the correct name to the abnormal name in one line: 
  	   abnormal name [tab] correct name [return]

	2. manual_genus_table.txt
	   one genus & classification perl line. If the genus has one division, that is pretty good result. If genus 
	   has more than one division, make it to keep one or more according frequency of each division, or searching
	   background of this genus to make correct decision.

	3. manual_desc_table.txt
	   The format of manual_desc:
	   1 - ID: GFLRNA3 
	   2 - Desc: Grapevine fanleaf virus satellite RNA (RNA3), complete cds.
	   3 - Same Desc: Grapevine fanleaf virus 
	   4 - Plants  
	   5 - No. of same: 1       
	   6 - Freq of same: 100.00
	   7 - ID of same: GFLRNA1

	   The script seach the same description of unclassified virus against classified virus, then brorrow the info 
	   of classification to the unclassified. Just check the 2nd, 3rd, and 4th column. If it does not make sense, 
	   correct them manually. Please also concern the 5th and 6th column, the error usually happens to the record 
	   without 100% freq. 

	After correction, please save the correction in same file with same format.

';
	exit
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

	# load default version
	my $version_file = $FindBin::RealBin."/default_version";
	my $version = `head -n 1 $version_file`; chomp($version);
	die "[ERR]version $version\n" if $version < 0;
	$version = "v$version";

	# load taxonomy information to hash for 
	# # 1. tracking host to division
	# # 2. tracking genus of virus
	# # description of hash
	# # hash division: key: taxon id, value division name
	# # hash node_table: key: taxon id, parent, rank, division; value: parent taxon id, rank, and division
	# # hash name_table: key: name; value: taxon id;
	# # hash virus_genus_taxon_id: key: taxon_id; value: genus scientific name

	my $node_file = $FindBin::RealBin."/nodes.dmp.gz";
	my $name_file = $FindBin::RealBin."/names.dmp.gz";
	my $host_info_file = $FindBin::RealBin."/ICTV_2013_genus_uniq_host_info_raw.txt";
	my $update_hname = $FindBin::RealBin."/update_hname_table_$version.txt";
	my $update_genus = $FindBin::RealBin."/update_genus_table_$version.txt";
	my $update_desc  = $FindBin::RealBin."/update_desc_table_$version.txt";
	foreach my $f (($node_file, $name_file, $host_info_file, $update_hname, $update_genus, $update_desc)) {
		print "[ERR]file not exist: $f\n" and exit unless -s $f;
	}

	my %division  = get_division();
	my ($node_table, $virus_genus_taxon_id) = load_taxon_node($node_file);
	my ($name_table, $virus_genus_taxon) = load_taxon_name($name_file, $virus_genus_taxon_id);
	my %virus_hostdiv = load_host_div($host_info_file, $update_genus);
	my %update_hname = load_update_hname($update_hname);
	my %update_desc  = load_update_desc($update_desc, \%division);

	print "== For $version ==\n";
	print "== ", scalar(keys(%$virus_genus_taxon_id)), " virus genus taxon has been load into hash ==\n";
	print "== ", scalar(keys(%$virus_genus_taxon)), " virus genus taxon has been load into hash ==\n";
	print "== ", scalar(keys(%virus_hostdiv)), " virus genus host info has been load into hash ==\n";
	print "== ", scalar(keys(%update_hname)), " corrected host name has been load into hash ==\n";
	print "== ", scalar(keys(%update_desc)), " manually classification has been load into hash ==\n";
	print "== taxon info has been load into hash ==\n";

	# debug: generate virus genus id and name
	# foreach my $tid (sort keys %$virus_genus_taxon) { print $tid."\t".$$virus_genus_taxon{$tid}."\n"; } exit;

	# parse viral sequence to database file
	my $vrl_info_classified   = "vrl_genbank.txt";
	my $vrl_info_unclassified = "vrl_genbank_unclassified.txt";
	my $vrl_seq_unclassified  = "vrl_genbank_unclassified.fasta"; 
	my $manual_hname = "manual_hname_table.txt";
	my $manual_genus = "manual_genus_table.txt";
	my $manual_desc  = "manual_desc_table.txt";
	my %manual_hname_check;	# key: hname
	my %manual_genus_check; # key: genus, div; value count

	# generate viral sequence with fasta format
	my %vrl_seq_info;

	my $out1 = IO::File->new(">".$vrl_info_classified) || die $!;
	my $out2 = IO::File->new(">".$vrl_info_unclassified) || die $!;
	my $out3 = IO::File->new(">".$vrl_seq_unclassified) || die $!;

	foreach my $f (@vrl_files)
	{
		print "parsing file $f start ......\n";

		# get number of sequence for each file
		my $f_head = `less $f | head`;
		my $seq_num = 70000;
		if ($f_head =~ m/(\d+) loci/) { $seq_num = $1; }
		my $parse_seq_num = 0; 
		my $pre_parse_percent = 0;

		# parse each sequence
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
							$normal_host_name = $update_hname{$value} if defined $update_hname{$value};	# correct abnormal host name
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
				#print "STAT$f\t$sid\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
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

					# record abnormal host name (for manually check table 1)
					if ($host_taxon eq 'NA') {
						$manual_hname_check{$hname} = 1;
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

						# record genus and classified division (for manually check table 2)
						if ($genus_name ne 'NA') {
							if (defined $manual_genus_check{$genus_name}{$host_div}) {
								$manual_genus_check{$genus_name}{$host_div}++;
							} else {
								$manual_genus_check{$genus_name}{$host_div} = 1;
							}
						}	
					}
					#print "STAT$f\t$sid\t$otid\t$genus_taxon\t$genus_name\t$genus_div\t$hname\t$host_taxon\t$host_taxon_fixed\t$host_div\n";
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
				# check manually classification result
				if ( defined $update_desc{$sid} ) {
					print $out1 $sid,"\t",$inseq->length,"\t",$inseq->desc,"\t",$inseq->version,"\t",$update_desc{$sid},"\tmanual\n";
				} else { 
					$host_name =~ s/\t/,/ig;
					print $out2 $sid,"\t",$inseq->length,"\t",$inseq->desc,"\t",$inseq->version,"\tNA\t$host_name\n";
					print $out3 ">",$sid,"\n",$inseq->seq,"\n";
				}
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

	# output put manually changed file
	my $outx = IO::File->new(">".$manual_hname) || die $!;
	foreach my $hname (sort keys %manual_hname_check) {
		print $outx $hname."\n";
	}
	$outx->close;
		
	my $outy = IO::File->new(">".$manual_genus) || die $!;
	foreach my $genus (sort keys %manual_genus_check) 
	{
		print $outy $genus;
		my $total_vrl = 0; my $best_num = 0; my $best_div;
		foreach my $div (sort keys %{$manual_genus_check{$genus}}) {
			my $num = $manual_genus_check{$genus}{$div};
			$total_vrl = $total_vrl + $num;
			if ($num > $best_num) {
				$best_num = $num;
				$best_div = $div;
			}
		}

		my $best_freq = sprintf("%.2f", ($best_num / $total_vrl));
		print $outy "\t$best_div,$best_num,$best_freq";

		foreach my $div (sort keys %{$manual_genus_check{$genus}}) {
			next if $div eq $best_div;
			my $num = $manual_genus_check{$genus}{$div};
			my $freq = sprintf("%.2f", ($num / $total_vrl));
			print $outy "\t$div,$num,$freq";
		}
		print $outy "\n";
	}
	$outy->close;

	classify_by_classified($vrl_info_unclassified, $vrl_info_classified,  $manual_desc);
}

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 load_host_div : # load host info to hash (ICTV-Master-Species-List-2013_v2: http://www.mcb.uct.ac.za/tutorial/ICTV%20Species%20Lists%20by%20host.htm )
=cut
sub load_host_div
{
	my ($host_info_file, $update_genus) = @_;
        my %virus_hostdiv; # key: virus genus name, value: host divition num
        my $fh1 = IO::File->new($host_info_file) || die $!;
        while(<$fh1>)
        {
                chomp;
		# format
		# 1. genus name		T4likevirus		# make the genus name lower case for each match
		# 2. division ID	X2,G0			# virus may infected more than 1 divison
		# 3. division name	Archaea,Bacteria
                next if $_ =~ m/^#/;
                my @a = split(/\t/, $_);
                next if $a[1] eq 'NA';
                my $gname = lc($a[0]);
                $virus_hostdiv{$gname} = $a[1];
        }
        $fh1->close;

	my $fh2 = IO::File->new($update_genus) || die $!;
	while(<$fh2>)
	{
		chomp;
		# the format is like previous file
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		next if $a[1] eq 'NA';
		my $gname = lc($a[0]);
		if (defined $virus_hostdiv{$gname}) {
			if ($virus_hostdiv{$gname} ne $a[1]) {
				print "[Conflict|Repeat]$gname\t$virus_hostdiv{$gname}\t$a[1]\n";
			}
		} else {
			$virus_hostdiv{$gname} = $a[1];
		}
	}
	$fh2->close;

	return %virus_hostdiv;
}

=head2
 load_taxon_node
=cut
sub load_taxon_node
{
	my $node_file = shift;
	my %node_table;		  # key: taxon_id, parent, division, rank; value: parent, division, rank
	my %virus_genus_taxon_id; # key: taxon_id; value: 1
	my ($tax_id, $parent_tax_id, $rank, $division);

	my $fh;
	if ($node_file =~ m/\.gz$/) {
		$fh = IO::File->new("gunzip -c $node_file |") || die $!;
	} else {
		$fh = IO::File->new($node_file) || die $!;
	}

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

	my $fh;
	if ($name_file =~ m/\.gz$/) {
		$fh = IO::File->new("gunzip -c $name_file |") || die $!;
	} else {
		$fh = IO::File->new($name_file) || die $!;
	}

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
 load_update_hname : correct abnormal host name 
=cut
sub load_update_hname
{
	my $input = shift;
	my %name;
	my $fh = IO::File->new($input) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		if (defined $a[1] && $a[1]=~ /\S+/) {
			chomp($a[1]);
			$name{$a[0]} = $a[1];
		}
	}
	return %name;
}

=head2
 load_update_desc: load manually classified virus to hash
=cut
sub load_update_desc
{
	my ($input_file, $division) = @_;

	# reverse division
	my %division = reverse %$division;

	my %update_desc;	
	my $fh = IO::File->new($input_file) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		# format
		# 0-ID			SRU26458
		# 1-Description		Snakehead retrovirus (SnRV), complete genome.
		# 2-sub Description	Snakehead
		# 3-Division		Invertebrates,Vertebrates
		# 4-total		2
		# 5-%div		100
		# 6-top			AF147498,AF147498
		my @a = split(/\t/, $_);
		my @b = split(/,/, $a[3]);
		foreach my $b (@b) {
			die "[ERR]undef div $b:$_\n" unless defined $division{$b};
		}
		$update_desc{$a[0]} = $a[3];
	}
	$fh->close;
	return %update_desc;
}

=head2
 classify_by_classified : classify by classified 
=cut
sub classify_by_classified
{
        # set input and output file;
        my ($unclassified, $classified, $manual_desc) = @_;
        print "[ERR]input file not exist\n" and exit unless (-s $unclassified && -s $classified);

        ####################################################
        # ==== generate classification by description ==== #
        ####################################################
        # load classified to hash
        my %vrl_cat;
        my $fh1 = IO::File->new($classified) || die $!;
        while(<$fh1>)
        {
                # CY151054        2270    Influenza B virus (B/Brisbane/11/2002) polymerase PA (PA) gene, complete cds.   1       Vertebrates     ICTV,genbank
                chomp;
                my @a = split(/\t/, $_);
                if (defined $vrl_cat{$a[0]}) {
                        $vrl_cat{$a[0]}.="\t".$a[4];
                } else {
                        $vrl_cat{$a[0]} = $a[4];
                }
        }
        $fh1->close;

        # temp hash for saving time
        my %temp_cat;

        # check unclassfied
        my $out3 = IO::File->new(">".$manual_desc) || die $!;
        my $fh2 = IO::File->new($unclassified) || die $!;
        while(<$fh2>)
        {
                # MTSCA   322     Lucerne transient streak virus satellite RNA, complete sequence.        1       NA      NA
                chomp;
                my @a = split(/\t/, $_);
                my ($id, $desc) = ($a[0], $a[2]);
                my @b = split(/ /, $desc);

                my $length = scalar(@b);
                my ($subdesc, $best_cat, $best_freq, $best_ratio, $cat_ids);
                for(my $i=0; $i<$length-1; $i++) {
                        pop @b;
                        $subdesc = join(" ", @b);
                        last if defined $temp_cat{$subdesc}; # check previous result

                        # get the best match of desc, the find the best classification
                        my $match_content = `grep \"$subdesc\" $classified`;
                        my %freq_cat = ();
                        if (defined $match_content && $match_content =~ m/\S+/)
                        {
                                chomp($match_content);
                                my @c = split(/\n/, $match_content);
                                my @ids;
                                foreach my $c (@c) {
                                        my @d = split(/\t/, $c);
                                        my $cat = $vrl_cat{$d[0]};
                                        push(@ids, $d[0]);
                                        $best_cat = $cat unless defined $best_cat;
                                        if (defined $freq_cat{$cat}) {
                                                $freq_cat{$cat}++;
                                        } else {
                                                $freq_cat{$cat} = 1;
                                        }
                                }

                                my $n = 4;
                                $n = scalar(@c) - 1 if @c < 5;
                                $cat_ids = join(",", @ids[0 .. $n]);
                                $best_freq = $freq_cat{$best_cat};
                                if (scalar(keys(%freq_cat)) == 1) {
                                        $best_ratio = 1;
                                        last;
                                }

                                foreach my $cname (sort keys %freq_cat) {
                                        if ($freq_cat{$cname} > $best_freq) {
                                                $best_freq = $freq_cat{$cname};
                                                $best_cat = $cname;
                                        }
                                }
                                $best_ratio = $best_freq / scalar(@c);
                                last;
                        }
                }

                # output result
                if (defined $temp_cat{$subdesc}) {
                        print $out3 "$id\t$desc\t$subdesc\t$temp_cat{$subdesc}\n";
                } elsif (defined $best_cat) {
                        $best_ratio = sprintf("%.2f", $best_ratio * 100);
                        print $out3 "$id\t$desc\t$subdesc\t$best_cat\t$best_freq\t$best_ratio\t$cat_ids\n";
                        $temp_cat{$subdesc} = "$best_cat\t$best_freq\t$best_ratio\t$cat_ids" unless defined $temp_cat{$subdesc};
                } else {
                        print $out3 "$id\t$desc\tNA-CAT\n";
                        #exit;
                }
        }
        $fh2->close;
        $out3->close;
}

