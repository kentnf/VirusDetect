#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use Bio::SeqIO;
use FindBin;

my $usage = qq'
USAGE: $0 datbase_file (default: vrl_plant.u95.fasta.gz)

';

my $db_file = "vrl_plant.u95.fasta.gz";
$db_file = $ARGV[0] if defined $ARGV[0];
$db_file = $FindBin::RealBin."/databases/$db_file";
die $usage unless -s $db_file;

print $usage;
print "[MSG]selecting $db_file as default database\n";
print "[MSG]the script will create plant virus database for you automatically ...\n";
sleep(5);

my $db = $FindBin::RealBin."/databases/vrl_plant";
my $samtools_bin = $FindBin::RealBin."/bin/samtools";
my $bwa_bin      = $FindBin::RealBin."/bin/bwa";
my $formatdb_bin = $FindBin::RealBin."/bin/formatdb";

print "[MSG]building plant virus database ...\n";
run_cmd("gzip -cd $db_file > $db");
run_cmd("$samtools_bin faidx $db");
run_cmd("$bwa_bin index $db");
run_cmd("$formatdb_bin -i $db -p F");

print "[MSG]building plant virus protein database ...\n";
my $db_prot_s  = $FindBin::RealBin."/databases/vrl_plant_prot.u100.fasta.gz";
my $db_table_s = $FindBin::RealBin."/databases/vrl_plant_prot_table.u100.gz"; 
die "[ERR]no file $db_prot_s\n"  unless -s $db_prot_s;
die "[ERR]no file $db_table_s\n" unless -s $db_table_s;

my %id; # save plant virus id to hash
my $in1 = Bio::SeqIO->new(-format=>'fasta', -file=>$db);
while(my $inseq = $in1->next_seq) {
        my $sid = $inseq->id;
        $id{$sid} = 1;
}

my %prot; # save plant protein seq to hash
my $in2;
if ($in2 =~ m/\.gz$/) {
	$in2 = Bio::SeqIO->new(-format=>'fasta', -file=>"gzip -cd $db_prot_s |");
} else {
	$in2 = Bio::SeqIO->new(-format=>'fasta', -file=>"$db_prot_s");
}
while(my $inseq = $in2->next_seq) {
        my $sid = $inseq->id;
        my $seq = $inseq->seq;
        $prot{$sid} = $seq;
}

my $db_prot_d  = $FindBin::RealBin."/databases/vrl_plant_prot";
my $db_table_d = $FindBin::RealBin."/databases/vrl_plant_prot_table";

open(OUT1, ">".$db_prot_d) || die $!;
open(OUT2, ">".$db_table_d) || die $!;

my $fh2;
if ($db_table_s =~ m/\.gz$/) {
	open($fh2, "-|", "gzip -cd $db_table_s") || die $!;
} else {
	open($fh2, $db_table_s) || die $!;
}
while(<$fh2>)
{
        chomp;
        my @a = split(/\t/, $_);
        if (defined $id{$a[0]}) {
                die "[ERR]undef prot seq $a[1]\n" unless defined $prot{$a[1]};
                my $prot_seq = $prot{$a[1]};
                print OUT1 ">$a[1]\n$prot_seq\n";
                print OUT2 $_."\n";
        }
}
close($fh2);

close(OUT1);
close(OUT2);

run_cmd("$formatdb_bin -i $db_prot_d -p T");

print "[MSG]done\n";

sub run_cmd {
	my $cmd = shift;
	print "[CMD]".$cmd."\n";
	system($cmd) && die "[ERR][CMD]$cmd\n";
}

