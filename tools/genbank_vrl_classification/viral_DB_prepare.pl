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
use Bio::DB::GenBank;
use Net::FTP;
use IO::File;
use Getopt::Std;
use FindBin;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }
if		($options{'t'} eq 'download')	{ vrl_download(\%options, \@ARGV); }    # download dataset
elsif	($options{'t'} eq 'category')	{ vrl_category(\%options, \@ARGV); }    # automatically classification
elsif	($options{'t'} eq 'manually')	{ vrl_manually(\%options, \@ARGV); }	# manually classification
elsif	($options{'t'} eq 'patch')		{ vrl_patch(\%options, \@ARGV); }		# combine the result
elsif	($options{'t'} eq 'unique')		{ vrl_unique(\%options, \@ARGV); }		# unique the classified virus
elsif	($options{'t'} eq 'refine')		{ vrl_refine(); }						# revine manually changed classification
elsif	($options{'t'} eq 'extProt')	{ vrl_extract_protein(\%options, \@ARGV); }	# extract protein sequence
elsif   ($options{'t'} eq 'genProt')    { vrl_genProt(\%options, \@ARGV); }		# generate protein files
elsif	($options{'t'} eq 'rmDup')		{ vrl_rmDup(\%options, \@ARGV);	}		# remove duplication
else	{ usage($version); }

sub usage
{
	my $version = shift;
	my $usage = qq'
Virus Classification Pipeline (VCP; version $version)

USAGE: $0 -t [tool] [options]

	[tool]		[description]
	download	Generate the download script
	category	Classify virus sequences into different kingdoms
	path		Append manually corrected files to the classification
	extProt		Extract virus protein sequences
	unique		Remove redundancy in virus sequences
	genProt		Retrieve proteins for each division

';
	print $usage;
	exit;
}

=head2
 vrl_genProt: generate protein for database file
=cut
sub vrl_genProt
{
        my ($options, $files) = @_;

        my $usage = qq'
USAGE: $0 -t genProt input_seq vrl_genbank_prot vrl_genbank_tab
       
		* the output file is protein sequence of input_seq
		* input_seq_prot

';

        print $usage and exit if @$files < 2;
        my $input_seq = $$files[0];
		my $vrl_prot  = $$files[1];
		my $vrl_tab   = $$files[2];

		# load protein seq to hash
		# key: prot id
		# value: prot seq
		my %prot_seq;
		my $inp = Bio::SeqIO->new(-format=>'fasta', -file=>$vrl_prot);
		while(my $inseq = $inp->next_seq) {
			my $pid = $inseq->id;
			my $seq = $inseq->seq;
			$prot_seq{$pid} = $seq;
		}

		# load protein nt id mapping to hash
		# key: nt id
		# value: array of prot id
		my %nt_prot;
		my $fhi = IO::File->new($vrl_tab) || die $!;
		while(<$fhi>) {
			chomp;
			next if $_ =~ m/^#/;
			my @a = split(/\t/, $_);
			push(@{$nt_prot{$a[0]}}, $a[1]);
		}
		$fhi->close;

		# generate protein seq for input nt seq
		my $fho = IO::File->new(">".$input_seq."_prot") || die $!;
		my $inn = Bio::SeqIO->new(-format=>'fasta', -file=>$input_seq);
		while(my $inseq = $inn->next_seq) {
			my $id = $inseq->id;
			next unless defined $nt_prot{$id};
			foreach my $pid (@{$nt_prot{$id}}) {
				die "[ERR]undef seq for $pid\n" unless defined $prot_seq{$pid};
				my $pseq = $prot_seq{$pid};
				print $fho ">$pid\n$pseq\n";
			}
		}
		$fho->close;
}

=head
 extract_protein: extract protein sequence 
=cut
sub vrl_extract_protein
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t extProt [options] gbvrl1.seq.gz gbvrl2.seq.gz ... gbvrlN.seq.gz 
	-o output prefix (default: vrl_genbank)

	* two file will be generate
	* 1. vrl_genbank_prot: protein sequences used for blastx
	* 2. vrl_genbank_tab: id_mapping file for virus detection

';
	print $usage and exit unless defined $$files[0];
	foreach my $f (@$files) { 
		die "[ERR]file not exist: $f\n" unless -s $f;
	}

	my $output_prefix = 'vrl_genbank';
	my $output_protein = $output_prefix."_prot";
	my $output_table   = $output_prefix."_tab";

	# extract_protein_seq
	my $out1 = IO::File->new(">".$output_protein) || die $!;
	my $out2 = IO::File->new(">".$output_table) || die $!;

	foreach my $f (@$files)
	{
		my $seqin = Bio::SeqIO->new(-format => 'GenBank', -file => "gunzip -c $f |" );
		while(my $inseq = $seqin->next_seq)
		{
			my $sid = $inseq->id;
			my $acc = $inseq->accession;
			my $pid;
			foreach my $feat_object ($inseq->get_SeqFeatures) 
			{
				my $primary_tag = $feat_object->primary_tag;
				next unless $primary_tag eq "CDS";
				foreach my $tag ($feat_object->get_all_tags) {
					foreach my $value ($feat_object->get_tag_values($tag)) {
						$pid = $value if ($tag eq 'protein_id');
						if ($tag eq 'translation' && length($value) > 0) {
							print $out1 ">".$pid."\n".$value."\n";
							print $out2 "$acc\t$pid\n";
						}
					}
				}
			}
		}
	}
	$out1->close;
	$out2->close;
}

sub vrl_refine
{
	my @fasta = qw/vrl_Algae_all.fasta vrl_Archaea_all.fasta vrl_Bacteria_all.fasta vrl_Fungi_all.fasta vrl_Invertebrates_all.fasta vrl_Plants_all.fasta vrl_Protozoa_all.fasta vrl_Unassigned_all.fasta vrl_Vertebrates_all.fasta vrl_Viruses_all.fasta/;
	my $input_file = 'vrl_genbank.txt';

	# load seq to hash	
	my %seq_info;
	my %div_fh;
	my %div_file;
	foreach my $f (@fasta) {
		die "[ERR]file not exist $f\n" unless -s $f;
		my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$f);
		while(my $inseq = $in->next_seq) {
			next if defined $seq_info{$inseq->id};
			$seq_info{$inseq->id} = $inseq->seq;
		}

		my $new_file = $f.".new";
		my $div = $f;
		$div =~ s/vrl_//; $div =~ s/_all\.fasta//;
		my $fh = IO::File->new(">".$new_file) || die $!;
		$div_fh{$div} = $fh;
		$div_file{$div} = $new_file;
	}

	# load classification to hash
	my $in = IO::File->new($input_file) || die $!;
	while(<$in>)
	{
		chomp;
		my @a = split(/\t/, $_);
		my @b = split(/,/, $a[5]);
		foreach my $b (@b) {
			die "[ERR]undef div $b\n" unless defined $div_fh{$b};
			print {$div_fh{$b}} ">".$a[0]."\n".$seq_info{$a[0]}."\n";
		}
	}
	$in->close;

	foreach my $div (sort keys %div_fh) {
		$div_fh{$div}->close;
	}
}

=head2
 vrl_unique : unique the classified virus with different sequence similarity
=cut
sub vrl_unique
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t unique [option] input_file

	-m	max sequence length (default: 40,000)
	-n	min sequence length (default: 100)
	-s	sequence similarity (default: 100)
		100 for 100%, 97 for 97%, 95 for 95%
	-p	number of CPU (default: 1)


* remove sequence longer than 40kb/shorter than 100bp
* generate cluster using cd-hit-est with sequence similarity
* extract representative seq
'; 

	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]file not exist $input_file\n" unless -s $input_file;

	my $prefix = $input_file;
	$prefix =~ s/\.fasta$//ig; $prefix =~ s/\.fa$//ig;

	# parameters
	my $identity = 100;
    $identity = $$options{'s'} if (defined $$options{'s'} && $$options{'s'} > 90);
	my $cls_fasta  = $prefix."_U$identity.fasta";
	$identity = $identity/100;

	my $cpu = 1;
	$cpu = $$options{'p'} if (defined $$options{'p'} && $$options{'p'} > 1);

	my $max_len = 40000;
	my $min_len = 100;
	$max_len = $$options{'m'} if (defined $$options{'m'} && $$options{'m'} >= 1000);
	$min_len = $$options{'n'} if (defined $$options{'n'} && $$options{'n'} >= 30 && $$options{'n'} < $max_len);

	# filter reads by length, put sequence to hash
    my %seq_hash;	# key: id, value: seq
	my $filter_file = $prefix.".filter.fasta";
	my $out1 = IO::File->new(">".$filter_file) || die $!;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_file);
	while(my $inseq = $in->next_seq){   
		my $id  = $inseq->id;
		my $seq = $inseq->seq;
		my $len = $inseq->length;
		if ($len >= $min_len && $len <= $max_len) {
			print $out1 ">$id\n$seq\n";
			$seq_hash{$id} = $seq;
		}
	}
	$out1->close;

	# cd-hit
	my $cdhit_bin = $FindBin::RealBin."/../../bin/cd-hit-est";
	my $cls_file = $prefix.".cls$identity";
	run_cmd("$cdhit_bin -i $filter_file -o $cls_file -c $identity -M 0 -T $cpu");
	
	########
	my $cls_result = $prefix.".cls$identity.clstr";
	my $fho = IO::File->new(">".$cls_fasta) || die $!;
	my $fhi = IO::File->new($cls_result) || die $!;
	while(<$fhi>) {
		chomp;
		if ($_ =~ m/>(\S+)\.\.\. \*/) {
			my $id = $1;
			my $seq = $seq_hash{$id};
			print $fho ">".$id."\n".$seq."\n";
		}
	}
	$fhi->close;
	$fho->close;
}

=head2
 vrl_unique : unique the classified virus 
=cut
sub vrl_rmDup
{
	my ($options, $files) = @_;
	
	my $usage = qq'
USAGE: $0 -t rmDup input_file

* remove sequence longer than 40kb
* remove duplication virus (shorter one)	

';

	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]file not exist $input_file\n" unless -s $input_file;

	# parameters
	my $max_len_cutoff = 40000;
	my $min_len_cutoff = 100;

	# remove seq longer than 40k, the left is ordered by length 
	my $temp_seq = $input_file.".temp.fa";
	my $out = IO::File->new(">".$temp_seq) || die $!;
	my %seq_info;
	my %uniq_seq;
	my %removed;
	my ($id, $seq, $length, $revseq);
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_file);
	while(my $inseq = $in->next_seq)
	{
		$id = $inseq->id;
		$seq = $inseq->seq;
		$length = $inseq->length;
		$revseq = reverse($seq);
		$revseq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

		$seq_info{$id}{'seq'} = $seq;
		$seq_info{$id}{'len'} = $length;

		if ($length > $max_len_cutoff || $length < $min_len_cutoff) {
			$removed{$id} = $length;
			next;
		}

		if (defined $uniq_seq{$seq}) {
			$removed{$id} = "Dup";
			next
		}

		$uniq_seq{$seq} = 1;
		$uniq_seq{$revseq} = 1;

		print $out ">$id\n$seq\n";	
	}
	$out->close;

	# self blast
	my $cmd_format = "formatdb -i $temp_seq -p F";
	my $cmd_blast = "blastall -p blastn -i $temp_seq -d $temp_seq -a 24 -v 2 -b 2 -m 8 -o $temp_seq.blast";
	my $blast_tool = "blastTool.pl";
	# die "[ERR]no blastTool.pl\n" unless -s $blast_tool;
	my $cmd_unique = "$blast_tool -t unique $temp_seq $temp_seq.blast";
	run_cmd($cmd_format);
	run_cmd($cmd_blast);
	run_cmd($cmd_unique);
	
	# check the output files
	foreach my $f (($temp_seq.".removed", $temp_seq.".unique", $temp_seq.".rmDup_table.txt")) {
		die "[ERR]file not exist $f\n" unless -e $f;
		#run_cmd("mv $temp_seq.removed $input_file.removed");
		#run_cmd("mv $temp_seq.unique  $input_file.unique");
		#run_cmd("mv $temp_seq.rmDup_table.txt $input_file.rmDup_table.txt");
	} 

	# update with pre removed sequences 
	my $out1 = IO::File->new(">>"."$input_file.removed") || die $!;
	my $out2 = IO::File->new(">>"."$input_file.rmDup_table.txt") || die $!;
	
	foreach my $id (sort keys %removed)
	{
		print $out1 ">".$id."\n".$seq_info{$id}{'seq'}."\n";
		print $out2 $id."\t".$removed{$id}."\n";
	}
	
	$out1->close;
	$out2->close;
}

=head2
 run_cmd: run other program using tophat
=cut
sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}


=head2
 vrl_patch : patch the manually corrected result to update file 
=cut
sub vrl_patch
{
	# do not need any parameters, just path the result to update files
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t path [options]
	-v version( default +1 of current )
	   if you input 211, the version should be updated to 211

';

	# load default version
	my $version_file = $FindBin::RealBin."/current_genbank_version";
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
	$ver_num = $$options{'v'} if (defined $$options{'v'} && $$options{'v'} > $ver_num);

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
			die "[ERR]patch genus div: $_\t$b[0]\n" unless (defined $division{$b[0]});
			push(@div_id, $b[0]);
			push(@div_name, $division{$b[0]});
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
		next unless $a[3];
		next if $a[3] eq 'NA-CAT';
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
Please change version number of current_genbank_version
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

	>> Steps
        A. open file manual_hname_table.txt 
        B. search the abnormal name using google and genbank
        C. add the correct to behind the abnormal name

	2. manual_genus_table.txt
	   one genus & classification perl line. If the genus has one division, that is pretty good result. If genus 
	   has more than one division, make it to keep one or more according frequency of each division, or searching
	   background of this genus to make correct decision.
	
	>> usually ignore this steps 

	3. manual_desc_table.txt
	   The format of manual_desc:
	   1 - ID: GFLRNA3 
	   2 - Desc: Grapevine fanleaf virus satellite RNA (RNA3), complete cds.
	   3 - Same Desc: Grapevine fanleaf virus 
	   4 - Plants  
	   5 - No. of same: 1       
	   6 - Freq of same: 100.00
	   7 - ID of same: GFLRNA1
	   8 - Div name by blast
	   9 - best hit of blast
	   10- match length of blast
	   11- percentage identity
	   12- match score

	   The script seach the same description of unclassified virus against classified virus, then brorrow the info 
	   of classification to the unclassified. Just check the 2nd, 3rd, and 4th column. If it does not make sense, 
	   correct them manually. Please also concern the 5th and 6th column, the error usually happens to the record 
	   without 100% freq. 
	   
	   Beside correction with word search method, we also blast the unclassified virus against the classified virus
	   (col 8-12). Most of word search result and blast search result are same, less the 100 of them are diff. So 
	   manually corret the different part will saving a lot of time
	
	>> Steps
	A. sort to find different in word search and blast search (col 4 and 8).
	B. check the blast match length, identify. Lower than 100 match base, 90% identity show low blast clue for classification
	C. check the word search column 3, 5, 6, and 7. 
	D. assign a correct div to column 4 for the virus 

	After correction, please save the correction in same file with same format.
   	then patch the changed files to previous file using:
      	  \$ perl $0 -t patch

        check the updated file in current folder
        1. update_genus_table_v205.txt
                old :   
                        Alphapermutotetravirus  NA      NA
                update: # patch v205
                        Alphapermutotetravirus  G1      Invertebrates
                Then remove the line on old file
        2. update_hname_table_v205.txt
                The manually changed name should be in the end of this file

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
    #my $genbank_ftp = "ftp.ncbi.nih.gov";
    #my $ftp = Net::FTP->new($genbank_ftp, Debug=>0) || die "[ERR]Cannot connect to $genbank_ftp $@\n";
    #$ftp->login("anonymous",'-anonymous@') || die "[ERR]Cannot login ", $ftp->message, "\n";
    #$ftp->cwd("/genbank") || die "[ERR]Cannot change working directory ", $ftp->message;

    my $doc=get("https://ftp.ncbi.nih.gov/genbank/");
    my @lines = split(/\n/, $doc);

    my @vrl_files;
    #my @files = $ftp->ls();
    my @files;
    foreach my $f (@lines) {
        if ($f =~ m/gbvrl(\d+)\.seq\.gz/) {
            my $file = "gbvrl$1.seq\.gz";
            #$ftp->get($f) unless -s $f;
            #push(@vrl_files, $f);
            #print $f."\t$file\n";
            #if ( -s $f ) {
            #   print "# == $f has been download from genbank ==\n";
            #} else {
            #   $download_cmd.="wget http://ftp.ncbi.nih.gov/genbank/$f\n";
            #}
            $download_cmd_all.="wget https://ftp.ncbi.nih.gov/genbank/$file\n";
        }
    }
    #$ftp->quit;

    #print "# cmd for download vrl database\n$download_cmd\n\n";
    print "# cmd for download all vrl database files\n$download_cmd_all\n\n";
    print "cat gbvrl*.seq.gz > GB_virus.gz\nrm gbvrl*.seq.gz\n";
    print "# cmd for dwonload taxonomy database\nwget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n";
    my @del_f = qw/citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp readme.txt taxdump.tar.gz/;
    my @zip_f = qw/names.dmp nodes.dmp/;
    print "tar -xzvf taxdump.tar.gz\n";
    print "rm ".join(" ", @del_f)."\n";
    print "gzip ".join(" ", @zip_f)."\n";
    exit;
}

=head2
 vrl_category : 
=cut
sub vrl_category
{
	my ($options, $files) = @_;	

	my $usage = qq'
USAGE: $0 -t category [options] input_file1 ... input_fileN

	-b previous classification info
	-c 1 make classify by classification enable

';
	print $usage and exit unless defined $$files[0];
	my @vrl_files = @$files;

	my $cbc = 0;
	$cbc = 1 if (defined $$options{'c'} && $$options{'c'} > 0);

	my $pre_cat;
	if (defined $$options{'b'}) {
		$pre_cat = $$options{'b'};
		die "[ERR]file not exit $pre_cat\n" unless -s $pre_cat;
	}

	# load default version
	my $version_file = $FindBin::RealBin."/current_genbank_version";
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

	my $node_file = "nodes.dmp.gz";
	my $name_file = "names.dmp.gz";
	die "[ERR]no taxonomy file\n" unless (-s $node_file && -s $name_file);
	my $host_info_file = $FindBin::RealBin."/ICTV_2013_genus_uniq_host_info_raw.txt";
	my $update_hname = $FindBin::RealBin."/update_hname_table_$version.txt";
	my $update_genus = $FindBin::RealBin."/update_genus_table_$version.txt";
	my $update_desc  = $FindBin::RealBin."/update_desc_table_$version.txt";
	foreach my $f (($node_file, $name_file, $host_info_file, $update_hname, $update_genus, $update_desc)) {
		print "[ERR]file not exist: $f\n" and exit unless -s $f;
	}

	my %division  = get_division();
	my %rev_division = reverse %division;
	my %pre_division = load_pre_classification($pre_cat) if defined $pre_cat;
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
	my $manual_hname = "manual_hname_table.txt";
	my $manual_genus = "manual_genus_table.txt";
	my $manual_desc  = "manual_desc_table.txt";
	my %manual_hname_check;	# key: hname
	my %manual_genus_check; # key: genus, div; value count

	# generate viral sequence with fasta format
	my %vrl_seq_info;

	my $out1 = IO::File->new(">".$vrl_info_classified) || die $!;

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
			#my $sid = $inseq->id;
			my $sid = $inseq->accession;
			$vrl_seq_info{$sid}{'ver'} = $inseq->seq_version;
			$vrl_seq_info{$sid}{'seq'} = $inseq->seq;
			$vrl_seq_info{$sid}{'des'} = $inseq->desc;

			my %host_division = ();	# key: divID, value, classification method

			# check if the virus was classified before
			if (defined $pre_division{$sid}) {
				print $out1 $pre_division{$sid}."\n";

				my @m = split(/\t/, $pre_division{$sid});
				my @dname = split(/,/, $m[5]);
                               	foreach my $dname (@dname) {
                                	my $div_id = $rev_division{$dname};
                             		$host_division{$div_id} = 'manual';
                               	}
				$vrl_seq_info{$sid}{'host_div'} = \%host_division; # for sequence classification
				next;
			}

			# get org taxon id or host name 
			# org taxon id -- the taxon id for the virus, which could be find on the db_xref tag of feature.
			# host name -- could be find in host tag 
			my ($org_taxon, $org_name, $genus_taxon, $genus_name, $genus_div, $host_taxon, $host_name, $host_div);
			my %source; # key: org_taxon; value: host_name
			foreach my $feat_object ($inseq->get_SeqFeatures) 
			{ 
				my $primary_tag = $feat_object->primary_tag;
				next unless $primary_tag eq "source";

				# get org taxon name and host name, they are two clues for classification
				$org_taxon = 'NA'; $host_name = 'NA';
				foreach my $tag ($feat_object->get_all_tags) { 
 					foreach my $value ($feat_object->get_tag_values($tag)) {
						if ($tag eq 'db_xref' && $value =~ m/taxon:(\d+)/) {
							if ( $org_taxon eq 'NA' ) { $org_taxon = $1; } else { die "[ERR]Repeat organism taxon id $sid\n"; }
						}
 
						if ($tag eq 'host' ) { 
							my $normal_host_name = $value;
							$normal_host_name = $update_hname{$value} if defined $update_hname{$value};	# correct abnormal host name
							if ( $host_name eq 'NA' ) { $host_name = $normal_host_name; } else { $host_name.= "\t".$normal_host_name; } # one org <=> more host
						}      
      					}          
   				}

				# find the organism that not linked to virus through the division number is not equal to 9
				# if get warn on screen, check the corrected virus taxon id in genbank, then change it manally in sub correct_org_taxon_division
				unless ( defined $$node_table{$org_taxon}{'division'} ) {
					my $changed_org_taxon = correct_org_taxon_division($sid);
					$org_taxon = $changed_org_taxon unless $changed_org_taxon eq 'NA';
					warn "[WARN]organism is not virus(missing): $sid\n" if $changed_org_taxon eq 'NA';
				}

				unless ( $$node_table{$org_taxon}{'division'} eq '9' ) {
					my $changed_org_taxon = correct_org_taxon_division($sid);
					$org_taxon = $changed_org_taxon unless $changed_org_taxon eq 'NA';
					warn "[WARN]organism is not virus(mistake): $sid\n" if $changed_org_taxon eq 'NA';
				}

				# try next feature of the organism is not virus
				next unless $$node_table{$org_taxon}{'division'} eq '9'; # example EF710638
				$source{$org_taxon} = $host_name;
			}

			# find the organism that not linked to virus through the division number is not equal to 9
			# if get warn on screen, check the corrected virus taxon id in genbank, then change it manally in sub correct_org_taxon_division

			#if (scalar(keys(%source)) < 1 ) {  
			#	#print "STAT$f\t$sid\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
			#}

			# ==== classification ====
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
					# === find host taxon ID according to host name ===
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

					# record abnormal host name (for manually chek)
					# this will get all abnormal name, some of them could be classified by genus
					# if ($host_taxon eq 'NA') { $manual_hname_check{$hname} = 1; }

					# get host division, and put it to hash
					if ( defined $$node_table{$host_taxon}{'division'} ) {
						next if $$node_table{$host_taxon}{'division'} == 8; # skip unsigned
						next if $$node_table{$host_taxon}{'division'} == 9; # skip virus
						$host_div = "G".$$node_table{$host_taxon}{'division'};
						$host_div = "G10" if ($host_div eq "G2" || $host_div eq "G5" || $host_div eq "G6");
				
						if (  defined $host_division{$host_div} ) {
							# notice -- disable classification using host name if the virus has already been 
							#           classified by virus taxon id
							#
							# if ( $host_division{$host_div} ne 'genbank' ) {
							# 	$host_division{$host_div}.=",genbank";
							# }
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
				# output classification to file
				my $host_div_name = '';
				my %host_source = ();
				foreach my $host_div_id (sort keys %host_division) {

					# check the undefined host_div_id if it do not have host div name
					unless ( defined $division{$host_div_id} ) {
						print $host_div_id,"\n";
						print $sid,"\n";
						foreach my $aa (sort keys %host_division) {
							print $aa."\t".$host_division{$aa}."\n";
						}
						exit;
					}

					$host_div_name.=",".$division{$host_div_id};
					$host_source{$host_division{$host_div_id}} = 1;
				}
				$host_div_name =~ s/^,//;
				my $output_source = join(",", sort keys %host_source);
				print $out1 $sid,"\t",$inseq->length,"\t",$genus_name,"\t",$inseq->desc,"\t",$inseq->version,"\t",$host_div_name,"\t",$output_source,"\n";

			} else {
				# check manually classification result
				if ( defined $update_desc{$sid} ) {
					# put classification to hash for sequence classification
					my @dname = split(/,/, $update_desc{$sid});
					foreach my $dname (@dname) {
						my $div_id = $rev_division{$dname};
						$host_division{$div_id} = 'manual';
					}

					# output classification to file
					print $out1 $sid,"\t",$inseq->length,"\t",$genus_name,"\t",$inseq->desc,"\t",$inseq->version,"\t",$update_desc{$sid},"\tmanual\n";

				} else {
			
					# record abnormal host name (for manually chek)	 
					my @hname = split(/\t/, $host_name);
					foreach my $hname (@hname) { $manual_hname_check{$hname} = 1; }

					# put classification to hash for sequence classification
					$host_division{'G8'} = 'manual';

					# output classification to file
					$host_name =~ s/\t/,/ig;
					print $out1 $sid,"\t",$inseq->length,"\t",$genus_name,"\t",$inseq->desc,"\t",$inseq->version,"\tUnassigned\tNA\n";
				}
			}
			$vrl_seq_info{$sid}{'host_div'} = \%host_division; # for sequence classification

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

	# output sequence for each division
	my %out_div_seq;	
	foreach my $id (sort keys %vrl_seq_info) {
		my $host_division = $vrl_seq_info{$id}{'host_div'} if defined $vrl_seq_info{$id}{'host_div'};
		foreach my $host_div_id (sort keys %$host_division) {
			die "[ERR]undef div id $host_div_id\n" unless defined $division{$host_div_id};
			my $div_name = $division{$host_div_id};
			if (defined $out_div_seq{$div_name}) {
				$out_div_seq{$div_name}.= ">$id\n$vrl_seq_info{$id}{'seq'}\n";
			} else {
				$out_div_seq{$div_name} = ">$id\n$vrl_seq_info{$id}{'seq'}\n";
			}
		}
	}

	foreach my $div_name (sort keys %out_div_seq) {
		$div_name =~ s/ /_/ig;
		my $vrl_seq_file = "vrl_".$div_name."_all.fasta";
		my $vrl_seq = $out_div_seq{$div_name};
		if ($vrl_seq) {
			my $fhs = IO::File->new(">".$vrl_seq_file) || die $!;
			print $fhs $vrl_seq;
			$fhs->close;
		}
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

	classify_by_classified($vrl_info_classified,  $manual_desc) if $cbc;
}

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 load_pre_classification
=cut
sub load_pre_classification
{
	my $pre_cat = shift;

	my %pre_division;
	
	my $fh = IO::File->new($pre_cat) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$pre_division{$a[0]} = $_;
	}
	$fh->close;
	
	return %pre_division;
}

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

	# create blank add_taxid.txt file if there is no info
	system("touch add_taxid.txt") unless -e 'add_taxid.txt';

	# check the exist record for accession ID and taxon id, save to hash
	my %sid_taxon;
	open(FH, 'add_taxid.txt') || die $!;
	while(<FH>) {
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$sid_taxon{$a[0]} = $a[1];
	}
	close(FH);

	# return the taxon id if the assession were searched previously
	if (defined $sid_taxon{$seq_id}) {
		return $sid_taxon{$seq_id};
	}
	
	# if there is no record for accesion, find the taxon in GenBank automatically
	my $taxon = 'NA';
	my $gb = Bio::DB::GenBank->new();
	my $seq = $gb->get_Seq_by_acc($seq_id); # Unique ID, *not always the LOCUS ID*
	foreach my $feat_object ($seq->get_SeqFeatures)
	{
		my $primary_tag = $feat_object->primary_tag;
		next unless $primary_tag eq "source";

		# get org taxon id
		foreach my $tag ($feat_object->get_all_tags) {
			foreach my $value ($feat_object->get_tag_values($tag)) {
				if ($tag eq 'db_xref' && $value =~ m/taxon:(\d+)/) {
					if ( $taxon eq 'NA' ) { 
						$taxon = $1; 

						# update seq_id and taxon id
						$sid_taxon{$seq_id} = $taxon;

						open(UP, ">add_taxid.txt") || die $!;
						foreach my $sid (sort keys %sid_taxon) {
							print UP $sid."\t".$sid_taxon{$sid}."\n";
						}
						close(UP);

						return $taxon; 
					} 
					# else { die "[ERR]Repeat organism taxon id $seq_id\n"; }
				}
			}
		}
	}

	# return taxon with NA
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
        my ($classified, $manual_desc) = @_;
        print "[ERR]input file not exist\n" and exit unless -s $classified;

	####################################################
	# ==== generate classification by blast	========== #
	####################################################
	
	# create hash for store classification by blast
	# key: vrl_i, div, identify, match, mismatch, gap
	# value: divname, 
	my %vid_cat;

	my $s_file = 'vrl_Unassigned_all.fasta';
	if (-s $s_file) 
	{
		# format s_file database;
		my $seq_num = 0;
		my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$s_file); 
		while(my $inseq = $in->next_seq) { $seq_num++;}
		return 1 if $seq_num == 0;
	
		my $format = 0;
		foreach my $f ( ("$s_file.nhr", "$s_file.nhr", "$s_file.nhr") ) {
			$format = 1 unless -s $f;
		}
		run_cmd("formatdb -i $s_file -p F") if $format == 1;
	
		# hash for query
		# key: file, value: division
		my %q_files;
		$q_files{'vrl_Algae_all.fasta'} = 'Algae';
		$q_files{'vrl_Archaea_all.fasta'} = 'Archaea';
		$q_files{'vrl_Bacteria_all.fasta'} = 'Bacteria';
		$q_files{'vrl_Fungi_all.fasta'} = 'Fungi';
		$q_files{'vrl_Invertebrates_all.fasta'} = 'Invertebrates';
		$q_files{'vrl_Plants_all.fasta'} = 'Plants';
		$q_files{'vrl_Protozoa_all.fasta'} = 'Protozoa';
		$q_files{'vrl_Vertebrates_all.fasta'} = 'Vertebrates';

		# bast to get best div
		foreach my $q_f (sort keys %q_files) {
			next unless -s $q_f;
			my $div = $q_files{$q_f};
			my $blast_file = "vrl_".$div.".blast";
			print "blastall -i $q_f -d $s_file -p blastn -m 8 -e 1e-5 -F F -a 24 -v 5 -b 5 -o $blast_file\n";
			run_cmd("blastall -i $q_f -d $s_file -p blastn -m 8 -e 1e-5 -F F -a 24 -v 5 -b 5 -o $blast_file");

			my $bfh = IO::File->new($blast_file) || die $!;
			while(<$bfh>) {
				chomp;
				next if $_ =~ m/^#/;
				my @a = split(/\t/, $_);
				next if scalar(@a) < 12;
				# AB000709        SLLCP   100.00  33      0       0       6219    6251    919     951     1e-09   65.9
				my ($qid, $sid, $iden, $match, $mismatch, $gap, $evalue, $score) = @a[0,1,2,3,4,5,10,11];
				$match = $match - $mismatch - $gap;
				if (defined $vid_cat{$sid}) {
					if ( $match > $vid_cat{$sid}{'match'} && $iden > $vid_cat{$sid}{'iden'} && $score > $vid_cat{$sid}{'score'} ) {
						$vid_cat{$sid}{'qid'} = $qid;
						$vid_cat{$sid}{'div'} = $div;
						$vid_cat{$sid}{'match'} = $match;
						$vid_cat{$sid}{'iden'} = $iden;
						$vid_cat{$sid}{'score'} = $score;
					}
				} else {
					$vid_cat{$sid}{'qid'} = $qid;
					$vid_cat{$sid}{'div'} = $div;
					$vid_cat{$sid}{'match'} = $match;
					$vid_cat{$sid}{'iden'} = $iden;
					$vid_cat{$sid}{'score'} = $score;
				}
			}
		} 
	}

        ####################################################
        # ==== generate classification by description ==== #
        ####################################################
        # load classified to hash, unclassified to array
	my @unclassified = ();
        my %vrl_cat;
        my $fh1 = IO::File->new($classified) || die $!;
        while(<$fh1>)
        {
                # CY151054        2270	influenza b	Influenza B virus (B/Brisbane/11/2002) polymerase PA (PA) gene, complete cds.   1       Vertebrates     ICTV,genbank
                chomp;
                my @a = split(/\t/, $_);

		if ($a[5] eq 'Unassigned') {
			push(@unclassified, $_);
		} else {
                        $vrl_cat{$a[0]} = $a[5];
		}
        }
        $fh1->close;

        # temp hash for saving time
        my %temp_cat;

        # check unclassfied
        my $out3 = IO::File->new(">".$manual_desc) || die $!;
	foreach my $line ( @unclassified )
	{
                # MTSCA   322     Lucerne transient streak virus satellite RNA, complete sequence.        1       NA      NA
                chomp($line);
                my @a = split(/\t/, $line);
                my ($id, $desc) = ($a[0], $a[3]);
                my @b = split(/ /, $desc);

                my $length = scalar(@b);
                my ($subdesc, $best_cat, $best_freq, $best_ratio, $cat_ids);
                for(my $i=0; $i<$length-1; $i++) {
                        pop @b;
                        $subdesc = join(" ", @b);
                        last if defined $temp_cat{$subdesc}; # get div using previous result

                        # get the best match of desc, the find the best classification
                        my $match_content = `grep \"$subdesc\" $classified | grep -v Unassigned`;
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

		#warn "[WARN]undef blast for $id\n" unless $vid_cat{$id}{'qid'};
		my $blast_info = 'NA-BLAST';
		if ( defined $vid_cat{$id}{'qid'} ) {
			$blast_info = $vid_cat{$id}{'div'}."\t".$vid_cat{$id}{'qid'}."\t".$vid_cat{$id}{'match'}."\t".$vid_cat{$id}{'iden'}."\t".$vid_cat{$id}{'score'};
		}

                # output result
                if (defined $temp_cat{$subdesc}) {
                        print $out3 "$id\t$desc\t$subdesc\t$temp_cat{$subdesc}\t$blast_info\n";
                } elsif (defined $best_cat) {
                        $best_ratio = sprintf("%.2f", $best_ratio * 100);
                        print $out3 "$id\t$desc\t$subdesc\t$best_cat\t$best_freq\t$best_ratio\t$cat_ids\t$blast_info\n";
                        $temp_cat{$subdesc} = "$best_cat\t$best_freq\t$best_ratio\t$cat_ids" unless defined $temp_cat{$subdesc};
                } else {
                        print $out3 "$id\t$desc\tNA-CAT\t\t\t\t\t$blast_info\n";
                        #exit;
                }
        }
        $out3->close;
}


