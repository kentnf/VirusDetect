#!/usr/bin/perl -w 
use strict; 

our $max_Xratio = 80/100; 
our $min_BaseNum = 40; 
#将fasta文件中的所有序列满足，[XxNn]碱基比例小于一个ratio，或者非[XxNn]碱基总数大于一个绝对值的序列输出到output
 
if (@ARGV < 2)
{
  print "usage: trm_XNseq1.pl input1 input2 max_Xratio min_BaseNum > output\n";
  exit(0);
}

our $input1 = $ARGV[0];
our $input2 = $ARGV[1];
$max_Xratio = $ARGV[2];
$min_BaseNum = $ARGV[3];

my %seq_required;
my ($pre_name, $seq) = ('', ''); 
open(IN1, "$input1");
while (<IN1>) {
	if (/^>(\S+)/) {
		my $current_name = $1; 
		if ($pre_name ne '') {
			my $is_out = 1; #表示该序列应该输出，这个变量调试用
			$seq =~ s/\s//g; 
			my $xnnum = ($seq =~ tr/XxNn/XxNn/); #得到非4种常见碱基[XxNn]的碱基总数
			my $seqlen = length($seq); 
			if ($xnnum/$seqlen >= $max_Xratio) {#[XxNn]碱基比例大于一个ratio
				$is_out = 0; #表示该序列不应该输出
			}elsif ($seqlen-$xnnum < $min_BaseNum) {#或者[ATCG]碱基总数不够一定数量
				$is_out = 0; #表示该序列不应该输出
			}else{
				#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
				#print ">$pre_name\n$seq\n"; #都不满足的就输出到标准输出
				defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#建立query和query_length之间的映射
			}
			#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; #调试用
		}
		$pre_name=$current_name; $seq = ''; 
	}else{
		$seq .= $_; 
	}
}
#不要忘记处理剩下的
if ($pre_name ne '') {
	my $is_out = 1; 
	$seq =~ s/\s//g; 
	my $xnnum = ($seq =~ tr/XxNn/XxNn/);
	my $seqlen = length($seq);
	if ($xnnum/$seqlen >= $max_Xratio) {
		$is_out = 0;
	}elsif ($seqlen-$xnnum < $min_BaseNum) {
		$is_out = 0;
	}else{
		#$seq =~ s/(.{50})/$1\n/g; chomp($seq); 
		#print ">$pre_name\n$seq\n"; #都不满足的就输出到标准输出
		defined $seq_required{$pre_name} or $seq_required{$pre_name} = 1;#建立query和query_length之间的映射
	}
	#$is_out == 0 and warn "[Record] [$pre_name] dropped.\n"; 
}
close(IN1);

#从input2中把在%seq_required中的序列提取出来
open(IN2, $input2);

my $flag = "off";
while(<IN2>) {
	if($_ =~ m/^>/) {
		my $head = $_;
		chomp($head);
		$head=~s/>//;

		if(defined $seq_required{$head}) {#如果包括这个name
			print $_;#就输出
			$flag = "on";#同时改变标志，表示后面的序列需要继续向OUT1输出
		}
		else {#如果不包括这个name
			#print OUT $_;#就输出到OUT
			$flag = "off";#同时改变标志，表示后面的序列需要继续向OUT2输出
		}
	}
	else {
		if($flag eq "on") {#表示为"on"
			print $_;#后面的序列需要继续输出
		}
	}
}
close(IN2);
