#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0];
my $idx = $ARGV[1]; # column number of key indicator
my @idx;
my $weight = $ARGV[2]; # weight of indicators
my @weight;
my %hash;
my $head;
my $suffix;
my ($num, $num2);
my $PATH = "./";

open (IN, $PATH.$file.".txt") or die "$!";

for (my $i=0;$i<=0;$i++){ # skipping head
	$_ = <IN>;
	$_ =~ s/^\s+|\s+$//g;
	$head = $_;
	my @array = split("\t",$_);
	@idx = split(",",$idx);
	for (my $i=0;$i<=$#idx;$i++){
		$idx[$i] -= 1;
		$suffix .= "_".$array[$idx[$i]];
	}
}

while (<IN>){
	$_ =~ s/^\s+|\s+$//g;
	my @array = split("\t",$_);
	for (my $i=$idx[0];$i<=$idx[$#idx];$i++){
		if ($array[$i] ne ""){
			$num2 += 1;
		}
	}
	@weight = split(",",$weight);
	if ($num2 == $#idx + 1){
		for (my $j=0;$j<=$#array;$j++){
			if ($array[$j] ne "" & $array[$j] ne "不配合" & $array[$j] ne "不能配合" & $array[$j] ne "不能完成" & $array[$j] ne "未测" & $array[$j] ne "-" & $array[$j] ne "—"){ # please use your own key words for missing data
				$num += $weight[$j];
			}
		}
		$hash{$array[0]}{$num} = join("\t",@array);
		$num = 0;
	}
	$num2 = 0;
}

close IN or die "$!";

#$file =~ s/.txt/./;
open (OUT, "> ".$PATH.$file.$suffix) or die "$!";

print OUT $head."\n";
foreach my $key (sort keys %hash){
	foreach my $key2 (sort {$b <=> $a} keys %{$hash{$key}}){
		print OUT $hash{$key}{$key2}."\n";
		last;
	}
}

close OUT or die "$!";

open (IN, $PATH.$file.$suffix) or die "$!";
open (OUT, "> ".$PATH.$file.$suffix.".rdata") or die "$!";

$_ = <IN>;
$_ =~ s/^\s+|\s+$//g;
my @array = split("\t",$_);
my $len = $#array;
print OUT join("\t",@array)."\n";

while (<IN>){
	$_ =~ s/^\s+|\s+$//g;
	my @array = split("\t",$_);
	if ($#array < $len){
		print OUT join("\t",@array,"\t"x($len - $#array - 1))."\n";
	}else{
		print OUT join("\t",@array)."\n";
	}
}

close IN or die "$!";
close OUT or die "$!";
