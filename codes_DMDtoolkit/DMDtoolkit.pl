################################
# 
#       by Jiapeng Zhou
#      zhoujp111@126.com
#
################################

#!/usr/bin/perl
#use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/; 

my $PATH = "./";
my $file = $ARGV[0];

=pod

### gene seq ###

sub gene_seq{

	open (IN, $PATH."Dp427m CDs.txt") or die "$!";
	open (IN2, $PATH."DMD gene.fa") or die "$!";
	open (OUT, "> ".$PATH."Dp427m CDs.fa") or die "$!";
	
	my $seq = <IN2>;
	my $subseq = ();
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		my $subseq .= substr($seq, 2220391 + 31119219 - $array[1] - 1, $array[1] - $array[0] + 1);
		print OUT $subseq;
	}
	
	close IN or die "$!";
	close IN2 or die "$!";
	close OUT or die "$!";
}

&gene_seq;

### normal protein seq ###

sub pro_seq{

	my %hash = ();
	
	open (IN, $PATH."codon list.txt") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		$hash{$array[0]} = $array[2];
	}
	
	close IN or die "$!";
	
	open (IN, $PATH."Dp427m CDs.fa") or die "$!";
	open (OUT, "> ".$PATH."Dp427m protein.fa") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		for (my $i=0; $i<length($_) - 5; $i+=3){
			my $codon = substr($_, $i, 3);
			if (exists $hash{$codon}){
				print OUT $hash{$codon};
			}
		}
	}
	
	close IN or die "$!";
	close OUT or die "$!";
	
}

&pro_seq;

### protein domains ###

sub pro_dom{

	my %hash = ();
	
	open (IN, $PATH."Dp427m CDs.txt") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		$hash{$array[0]}{"start"} = $array[3];
		$hash{$array[0]}{"end"} = $array[4];
	}
	
	close IN or die "$!";
	
	open (IN, $PATH."Dp427m Domains.txt") or die "$!";
	open (OUT, "> ".$PATH."Dp427m Domains.rdata") or die "$!";
	
	$_ = <IN>;
	$_ =~ s/^\s+|\s+$//g;
	my @array = split("\t",$_);
	print OUT join("\t",@array,"exon_start","exon_end")."\n";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split("\t",$_);
		print OUT join("\t",@array);
		for (my $i=1; $i<=79; $i++){
			if ($hash{$i}{"start"} <= $array[1]*3 & $hash{$i}{"end"} >= $array[1]*3){
				print OUT "\t".$i;
				last;
			}
		}
		for (my $i=1; $i<=79; $i++){
			if ($hash{$i}{"start"} <= $array[2]*3 & $hash{$i}{"end"} >= $array[2]*3){
				print OUT "\t".$i."\n";
				last;
			}
		}
	}
	
	close IN or die "$!";
	close OUT or die "$!";
	
}

&pro_dom;

=cut

### mutation to disordered protein seq ###

sub main_pros{

	my ($key1, $key2);
	my %hash = ();
	my %hash2 = ();
	my %hash3 = ();
	my $len;
	my (@tmp, @tmp2);
	my ($pos, $codon);
	
	open (IN, $PATH."Dp427m CDs.txt") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		$hash{$array[0]} = join("\t",@array[3,4]);
	}
	
	close IN or die "$!";
	
	open (IN, $PATH."codon list.txt") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		$hash3{$array[0]} = $array[2];
	}
	
	close IN or die "$!";
	
	open (IN, $PATH.$file.".txt") or die "$!";
	open (IN2, $PATH."Dp427m CDs.fa") or die "$!";
	open (OUT, ">".$PATH.$file.".Dp427m.pros") or die "$!";
	open (OUT2, ">".$PATH.$file.".Dp427m.rdata") or die "$!";
	
	my $head = <IN>; 
	chomp($head);
	my $CDs = <IN2>;
	print OUT join("\t", $head, "cds", "pro1", "pro2")."\n";
	print OUT2 join("\t", $head, "exon_mut", "pos_mut", "cds_range")."\n";
	
	while (<IN>){
	
		%hash2 = ();
		
		$_ =~ s/^\s+|\s+$//g;
		print OUT "$_\t";
		my @array = split(/\t/,$_);
		my @array2 = split(/, |,|£¬|; |;/,$array[4]);
		print OUT2 join("\t",@array[0..3],join("; ",@array2))."\t";
		
		for (my $i=0; $i<=$#array2; $i++){
				#exonic del/dup
			if ($array2[$i] =~ /exon([0-9]+)\-([0-9]+)(dup|del)/){
				@tmp = split(/\t/,$hash{$1});
				@tmp2 = split(/\t/,$hash{$2});
				$hash2{"exons"}{$1."-".$2.$3} = [$3,$tmp[0],$tmp2[1]];
				print OUT2 join(" ",$1,$2,"exon_".$3)."\t".join(" ",$tmp[0],$tmp2[1])."\t";
				#exonic del/dup
			}elsif ($array2[$i] =~ /exon([0-9]+)(dup|del)/){
				@tmp = split(/\t/,$hash{$1});
				$hash2{"exon"}{$1.$2} = [$2,$tmp[0],$tmp[1]];
				print OUT2 join(" ",$1,$1,"exon_".$2)."\t".join(" ",$tmp[0],$tmp[1])."\t";
				#short fragment del/dup
			}elsif ($array2[$i] =~ /c\.([0-9]+)_([0-9]+)(delins)([ACGT]+)/){
				$hash2{"indel"}{$3} = [$1,$2,$4];
				foreach my $key (sort keys %hash){
					my @array3 = split(/\t/,$hash{$key});
					if (($1 >= $array3[0]) & ($2 <= $array3[1])){
						print OUT2 join(" ",$key,$key,$3)."\t".join(" ",$1,$2)."\t";
					}
				}
				#short fragment del/dup
			}elsif ($array2[$i] =~ /c\.([0-9]+)_([0-9]+)(ins|del|dup)([ACGT]*)/){
				$hash2{"indel"}{$3} = [$1,$2,$4];
				foreach my $key (sort keys %hash){
					my @array3 = split(/\t/,$hash{$key});
					if (($1 >= $array3[0]) & ($2 <= $array3[1]) & (length($4) == 1)){
						print OUT2 join(" ",$key,$key,$3)."\t".join(" ",$1,$1)."\t";
					}elsif (($1 >= $array3[0]) & ($2 <= $array3[1]) & (length($4) != 1)){
						print OUT2 join(" ",$key,$key,$3)."\t".join(" ",$1,$2)."\t";
					}
				}
				#short fragment del/dup
			}elsif ($array2[$i] =~ /c\.([0-9]+)(dup|del)([ACGT]*)/){
				$hash2{"indel"}{$2} = [$1,$1];
				foreach my $key (sort keys %hash){
					my @array3 = split(/\t/,$hash{$key});
					if (($1 >= $array3[0]) & ($1 <= $array3[1]) & ($2 eq "dup")){
						print OUT2 join(" ",$key,$key,$2)."\t".join(" ",$1+1,$1+1)."\t";
					}elsif (($1 >= $array3[0]) & ($1 <= $array3[1]) & ($2 eq "del")){
						print OUT2 join(" ",$key,$key,$2)."\t".join(" ",$1,$1)."\t";
					}
				}
				#splice site
			}elsif ($array2[$i] =~ /c\.([0-9]+)(\+|\-)[0-9]{1}[ACGT]\>[ACGT]/){
				$hash2{"as"}{$1} = $1;
				foreach my $key (sort keys %hash){
					my @array3 = split(/\t/,$hash{$key});
					if (($1 == $array3[0]) | ($1 == $array3[1])){
						print OUT2 join(" ",$key,$key,"as")."\t".join(" ",$array3[0],$array3[1])."\t";
					}
				}
				#point mutation
			}elsif ($array2[$i] =~ /c\.([0-9]+)([ACGT])\>([ACGT])/){
				$hash2{">"}{$1} = $3;	
				foreach my $key (sort keys %hash){
					my @array3 = split(/\t/,$hash{$key});
					if (($1 >= $array3[0]) & ($1 <= $array3[1])){
						print OUT2 join(" ",$key,$key,"point")."\t".join(" ",$1,$1)."\t";
					}
				}
				#promoter site
			}elsif ($array2[$i] =~ /(Dp[0-9]+).*(del)/){
				$hash2{"promoter"}{$1} = [$1,$2];
				print OUT2 join(" ",$1,$1,"promoter")."\t".join(" ",$1,$2)."\t";
			}
		}
		
		my $CDS = $CDs;
		
		foreach $key1 (keys %hash2){
			foreach $key2 (reverse sort keys %{$hash2{$key1}}){
				#point mutation
				if ($key1 eq ">"){
					substr($CDS, $key2 - 1, 1) = $hash2{$key1}{$key2};
				}
			}
		}
		
		foreach $key1 (keys %hash2){
			foreach $key2 (reverse sort keys %{$hash2{$key1}}){
				#exonic del/dup
				if ($key1 eq "exons" & $hash2{$key1}{$key2}[0] eq "del"){
					substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1) = "";
				}elsif ($key1 eq "exons" & $hash2{$key1}{$key2}[0] eq "dup"){
					substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1) = substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1)x2;
				}elsif ($key1 eq "exon" & $hash2{$key1}{$key2}[0] eq "del"){
					substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1) = "";
				}elsif ($key1 eq "exon" & $hash2{$key1}{$key2}[0] eq "dup"){
					substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1) = substr($CDS, $hash2{$key1}{$key2}[1] - 1, $hash2{$key1}{$key2}[2] - $hash2{$key1}{$key2}[1] + 1)x2;
				}
				#short fragment del/dup
				if ($key1 eq "indel" & $key2 eq "del"){
					substr($CDS, $hash2{$key1}{$key2}[0] - 1, $hash2{$key1}{$key2}[1] - $hash2{$key1}{$key2}[0] + 1) = "";
				}elsif ($key1 eq "indel" & $key2 eq "ins"){
					substr($CDS, $hash2{$key1}{$key2}[0] - 1, 1) = substr($CDS, $hash2{$key1}{$key2}[0] - 1, 1).$hash2{$key1}{$key2}[2];
				}elsif ($key1 eq "indel" & $key2 eq "delins"){
					substr($CDS, $hash2{$key1}{$key2}[0] - 1, $hash2{$key1}{$key2}[1] - $hash2{$key1}{$key2}[0] + 1) = $hash2{$key1}{$key2}[2];
				}elsif ($key1 eq "indel" & $key2 eq "dup"){
					substr($CDS, $hash2{$key1}{$key2}[0] - 1, $hash2{$key1}{$key2}[1] - $hash2{$key1}{$key2}[0] + 1) = substr($CDS, $hash2{$key1}{$key2}[0] - 1, $hash2{$key1}{$key2}[1] - $hash2{$key1}{$key2}[0] + 1)x2;
				}
				#splice site
				if ($key1 eq "as"){
					for (my $j=1; $j<= keys %hash; $j++){
						@tmp = split(/\t/,$hash{$j});
						if ($tmp[0] eq $hash2{"as"}{$key2} | $tmp[1] eq $hash2{"as"}{$key2}){
							substr($CDS, $tmp[0] - 1, $tmp[1] - $tmp[0] + 1) = "";
						}
					}
				}
			}
		}
		
		foreach $key1 (keys %hash2){
			foreach $key2 (reverse sort keys %{$hash2{$key1}}){
				#promoter site
				if ($key1 eq "promoter"){
					$CDS = "";
				}
			}
		}
		
		print OUT $CDS."\t";
		
		if ($CDS =~ m/^$/){
			print OUT2 "0 0";
		}elsif ($CDS !~ m/^ATG/){
			$pos = index($CDS, "ATG");
			print OUT2 $pos+1;
		}elsif ($CDS =~ m/^ATG/){
			$pos = 0;
			print OUT2 $pos+1;
		}
		
		for (my $i = $pos; $i < length($CDS) - 2; $i += 3){
			$codon = substr($CDS, $i, 3);
			if (exists $hash3{$codon} & $hash3{$codon} ne "X"){
				print OUT $hash3{$codon};
			}elsif ($hash3{$codon} eq "X"){
				print OUT2 " ".$i;
				last;
			}
		}
		
		print OUT "\t";
		print OUT2 "\n";
		
		for (my $i = $pos; $i < length($CDS) - 2; $i += 3){
			$codon = substr($CDS, $i, 3);
			if (exists $hash3{$codon}){
				print OUT $hash3{$codon};
			}
		}
		
		print OUT "\n";
		
	}
	
	close IN or die "$!";
	close IN2 or die "$!";
	close OUT or die "$!";
	close OUT2 or die "$!";

}

&main_pros;	

### statistics for disordered proteins ###

sub pros_stats{

	open (IN, $PATH.$file.".Dp427m.pros") or die "$!";
	open (OUT, "> ".$PATH.$file.".Dp427m.pros.stats") or die "$!";
	
	$_ = <IN>;
	$_ =~ s/^\s+|\s+$//g;
	my @array = split(/\t/,$_);
	print OUT join("\t", @array[0..4], "pro_length", "codon_stops")."\n";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split("\t",$_);
		my $len = length($array[6]);
		my $Xs = 0;
		++$Xs while ($array[7] =~ m/X/g);
		print OUT join("\t", @array[0..4],$len,$Xs)."\n";
	}
	
	close IN or die "$!";
	close OUT or die "$!";

}

&pros_stats;

### reading-frame rule ###
	
sub fs_rule{

# restricted to exon deletions/duplications

	my %hash = ();

	open (IN, $PATH."Dp427m CDs.txt") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		my $tmp = ($array[3] - 1) % 3;
		my $tmp2 = $array[4] % 3;
		$hash{$array[0]} = join("\t", $tmp, $tmp2);
	}
	
	close IN or die "$!";
	
	open (IN, $PATH.$file.".Dp427m.pros.stats") or die "$!";
	open (OUT, "> ".$PATH.$file.".Dp427m.pros.diag") or die "$!";
	
	$_ = <IN>;
	$_ =~ s/^\s+|\s+$//g;
	print OUT join("\t", $_, "fs_rule")."\n";
	
	while (<IN>){
	
		$_ =~ s/^\s+|\s+$//g;
		my @array = split("\t", $_);
		
		if ($array[4] =~ m/^exon(\d+)(del|dup)/){
			if ($1 == 1){
				print OUT join("\t", @array, "DMD")."\n";
			}elsif ($2 eq "del"){
				my @array1 = split("\t", $hash{$1 - 1});
				my @array2 = split("\t", $hash{$1 + 1});
				if ($array1[1] == $array2[0]){
					print OUT join("\t", @array, "BMD")."\n";
				}else{
					print OUT join("\t", @array, "DMD")."\n";
				}
			}elsif ($2 eq "dup"){
				my @array1 = split("\t", $hash{$1});
				if ($array1[1] == $array1[0]){
					print OUT join("\t", @array, "BMD")."\n";
				}else{
					print OUT join("\t", @array, "DMD")."\n";
				}
			}
		}elsif ($array[4] =~ m/^exon(\d+)\-(\d+)(del|dup)/){
			if ($1 == 1){
				print OUT join("\t", @array, "DMD")."\n";
			}elsif ($3 eq "del"){
				my @array1 = split("\t", $hash{$1 - 1});
				my @array2 = split("\t", $hash{$2 + 1});
				if ($array1[1] == $array2[0]){
					print OUT join("\t", @array, "BMD")."\n";
				}else{
					print OUT join("\t", @array, "DMD")."\n";
				}
			}elsif ($3 eq "dup"){
				my @array1 = split("\t", $hash{$1});
				my @array2 = split("\t", $hash{$2});
				if ($array2[1] == $array1[0]){
					print OUT join("\t", @array, "BMD")."\n";
				}else{
					print OUT join("\t", @array, "DMD")."\n";
				}
			}
		}else{
			print OUT join("\t", @array, "-")."\n";
		}
	}	
	close IN or die "$!";
	close OUT or die "$!";

# expand to small deletions/duplications and splice sites 

	my %hash = ();
	
	open (IN, $PATH.$file.".Dp427m.pros") or die "$!";
	
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		if (length($array[5]) == 11058){
			$hash{$array[4]} = "-";
		}elsif ((length($array[5]) > 0) & (length($array[5]) % 3 == 0)){
			$hash{$array[4]} = "BMD";
		}else{
			$hash{$array[4]} = "DMD";
		}
	}
		
	close IN or die "$!";
	
	open (IN, $PATH.$file.".Dp427m.pros.stats") or die "$!";
	open (OUT, "> ".$PATH.$file.".Dp427m.pros.diag2") or die "$!";
	
	$_ = <IN>;
	$_ =~ s/^\s+|\s+$//g;
	print OUT join("\t", $_, "fs_rule")."\n";
	
	while (<IN>){
	
		$_ =~ s/^\s+|\s+$//g;
		my @array = split("\t", $_);
		
		if (exists $hash{$array[4]}){
			print OUT join("\t", @array, $hash{$array[4]})."\n";
		}
	}
	
	close IN or die "$!";
	close OUT or die "$!";

# apply ESE rule to nonsense mutations, and apply size and location info to in-frame deletions
=pod
	my %hash = ();
	open (IN, $PATH."matrices.txt") or die "$!";
	while (<IN>){
		if ($_ =~ m/^SRSF1/){
			for (my $i=0; $i<=6; $i++){
				$_ = <IN>;
				$_ =~ s/^\s+|\s+$//g;
				my @array = split("\t", $_);
				$hash{"SRSF1"}{$i}{"A"} = $array[1];
				$hash{"SRSF1"}{$i}{"C"} = $array[2];
				$hash{"SRSF1"}{$i}{"G"} = $array[3];
				$hash{"SRSF1"}{$i}{"T"} = $array[4];
			}
		}elsif ($_ =~ m/^SRSF2/){
			for (my $i=0; $i<=7; $i++){
				$_ = <IN>;
				$_ =~ s/^\s+|\s+$//g;
				my @array = split("\t", $_);
				$hash{"SRSF2"}{$i}{"A"} = $array[1];
				$hash{"SRSF2"}{$i}{"C"} = $array[2];
				$hash{"SRSF2"}{$i}{"G"} = $array[3];
				$hash{"SRSF2"}{$i}{"T"} = $array[4];
			}
		}elsif ($_ =~ m/^SRSF5/){
			for (my $i=0; $i<=6; $i++){
				$_ = <IN>;
				$_ =~ s/^\s+|\s+$//g;
				my @array = split("\t", $_);
				$hash{"SRSF5"}{$i}{"A"} = $array[1];
				$hash{"SRSF5"}{$i}{"C"} = $array[2];
				$hash{"SRSF5"}{$i}{"G"} = $array[3];
				$hash{"SRSF5"}{$i}{"T"} = $array[4];
			}
		}elsif ($_ =~ m/^SRSF6/){
			for (my $i=0; $i<=5; $i++){
				$_ = <IN>;
				$_ =~ s/^\s+|\s+$//g;
				my @array = split("\t", $_);
				$hash{"SRSF6"}{$i}{"A"} = $array[1];
				$hash{"SRSF6"}{$i}{"C"} = $array[2];
				$hash{"SRSF6"}{$i}{"G"} = $array[3];
				$hash{"SRSF6"}{$i}{"T"} = $array[4];
			}
		}
	}
	close IN or die "$!";

	my %hash2 = ();	
	open (IN, $PATH."Dp427m CDs.txt") or die "$!";
	<IN>;
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
		$hash2{$array[0]}{"start"} = $array[3];
		$hash2{$array[0]}{"end"} = $array[4];
	}
	close IN or die "$!";
=cut
	open (IN, $PATH.$file.".Dp427m.pros.diag2") or die "$!";
	open (IN2, $PATH."Dp427m CDs.fa") or die "$!";
	open (OUT, "> ".$PATH.$file.".Dp427m.pros.diag3") or die "$!";
	#my $seq = <IN2>;
	#my $subseq;
	#my @threshold = (1.96, 2.38, 2.67, 2.68);
	while (<IN>){
		$_ =~ s/^\s+|\s+$//g;
		my @array = split(/\t/, $_);
=pod
		if (($array[4] =~ /c\.([0-9]+)([ACGT])\>([ACGT])/) & ($array[$#array - 1] eq "2")){
			for (my $i=$1 - 7; $i<=$1 - 1; $i++){
				my @score = (0, 0, 0, 0);
				my @score2 = (0, 0, 0, 0);
				$subseq = substr($seq, $i, 7);
				my @array2 = split("", $subseq);
				for (my $j=0; $j<=$#array2; $j++){
					$score[0] += $hash{"SRSF1"}{$j}{$array2[$j]};
					$score[2] += $hash{"SRSF5"}{$j}{$array2[$j]};
				}
				$score2[0] = $score[0] + $hash{"SRSF1"}{$1 - $i -1}{$3} - $hash{"SRSF1"}{$1 - $i -1}{$2};
				$score2[2] = $score[2] + $hash{"SRSF5"}{$1 - $i -1}{$3} - $hash{"SRSF5"}{$1 - $i -1}{$2};
				#print OUT "SRSF1\t".$threshold[0]."\t".$score[0]."\t".$score2[0]."\n";
				#print OUT "SRSF5\t".$threshold[2]."\t".$score[2]."\t".$score2[2]."\n";
				for (my $k=0; $k<=$#score; $k++){
					if (($score[$k] - $threshold[$k] >= 0.1) & ($threshold[$k] - $score2[$k] >= 0.1)){
						foreach my $key (keys %hash2){
							if (($hash2{$key}{"start"} <= $1) & ($hash2{$key}{"end"} >= $1)){
								if (($hash2{$key}{"end"} - $hash2{$key}{"start"} + 1) % 3 == 0){
									$array[$#array - 1] = $array[$#array - 1].", "."BMD";
								}
							}
						}
					}
				}
			}
			for (my $i=$1 - 8; $i<=$1 - 1; $i++){
				my @score = (0, 0, 0, 0);
				my @score2 = (0, 0, 0, 0);
				$subseq = substr($seq, $i, 8);
				my @array2 = split("", $subseq);
				for (my $j=0; $j<=$#array2; $j++){
					$score[1] += $hash{"SRSF2"}{$j}{$array2[$j]};
				}
				$score2[1] = $score[1] + $hash{"SRSF2"}{$1 - $i -1}{$3} - $hash{"SRSF2"}{$1 - $i -1}{$2};
				#print OUT "SRSF2\t".$threshold[1]."\t".$score[1]."\t".$score2[1]."\n";
				for (my $k=0; $k<=$#score; $k++){
					if (($score[$k] - $threshold[$k] >= 0.1) & ($threshold[$k] - $score2[$k] >= 0.1)){
						foreach my $key (keys %hash2){
							if (($hash2{$key}{"start"} <= $1) & ($hash2{$key}{"end"} >= $1)){
								if (($hash2{$key}{"end"} - $hash2{$key}{"start"} + 1) % 3 == 0){
									$array[$#array - 1] = $array[$#array - 1].", "."BMD";
								}
							}
						}
					}
				}
			}
			for (my $i=$1 - 6; $i<=$1 - 1; $i++){
				my @score = (0, 0, 0, 0);
				my @score2 = (0, 0, 0, 0);
				$subseq = substr($seq, $i, 6);
				my @array2 = split("", $subseq);
				for (my $j=0; $j<=$#array2; $j++){
					$score[3] += $hash{"SRSF6"}{$j}{$array2[$j]};
				}
				$score2[3] = $score[3] + $hash{"SRSF6"}{$1 - $i -1}{$3} - $hash{"SRSF6"}{$1 - $i -1}{$2};
				#print OUT "SRSF6\t".$threshold[3]."\t".$score[3]."\t".$score2[3]."\n";
				for (my $k=0; $k<=$#score; $k++){
					if (($score[$k] - $threshold[$k] >= 0.1) & ($threshold[$k] - $score2[$k] >= 0.1)){
						foreach my $key (keys %hash2){
							if (($hash2{$key}{"start"} <= $1) & ($hash2{$key}{"end"} >= $1)){
								if (($hash2{$key}{"end"} - $hash2{$key}{"start"} + 1) % 3 == 0){
									$array[$#array - 1] = $array[$#array - 1].", "."BMD";
								}
							}
						}
					}
				}
			}
		}
=cut
		if (($array[4] =~ /^exon([0-9]+)\-([0-9]+)del$/) & ($array[$#array] eq "BMD")){
			if (($1 <= 8) & ($2 >= 10)){
				$array[$#array] = "DMD";
			}elsif (($1 >= 10) & ($2 <= 61) & ($2 - $1 >= 36)){
				$array[$#array] = "DMD";
			}elsif (($1 >= 62) & ($2 <= 69)){
				$array[$#array] = "DMD";
			}elsif ($1 >= 74){
				$array[$#array] = "DMD";
			}
		}
		print OUT join("\t", @array)."\n";
	}
	close IN or die "$!";
	close IN2 or die "$!";
	close OUT or die "$!";
	
}

&fs_rule;
