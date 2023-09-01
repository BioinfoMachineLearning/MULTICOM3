#!/usr/bin/perl -w
#define hash

$numArgs = @ARGV;
if($numArgs != 2)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$infolder	= "$ARGV[0]";
$outf		= "$ARGV[1]";

%h1=();
%h2=();

open(OUT, ">$outf/zscore_1.txt") || die("Couldn't open file $outf/zscore_1.txt\n");
open(OUT2, ">$outf/z_1.log") || die("Couldn't open file $outf/z_1.log\n");

opendir ( DIR3, "$infolder" ) || die "Error in opening dir $infolder\n";
while(($ln3 = readdir(DIR3))){

	if (length($ln3)>2 && substr($ln3,0,2) ne "TS") {
		print $ln3;
		open(IN1, "$infolder/$ln3") || die("Couldn't open file $infolder/$ln3\n"); 
		@a=<IN1>;
		close IN1;

		$name=substr($ln3,0,index($ln3,"."))."_zscore_1.txt";
		open(OUT1, ">$outf/$name") || die("Couldn't open file $outf/$name\n");

		@id1=();
		@zscore=();
		$cc9=0;

		$avg=0;
		$c=0;
		foreach $line (@a){
			chomp($line);
			$c++;
			$ss=substr($line,20);
			@i1=split(/ /,$ss);
			$avg+=$i1[0];
		}

		$avg=$avg/$c;
		print OUT2 substr($ln3,0,index($ln3,"."))."\n";
		print OUT2 "Avg: $avg\n";

		$sum=0;
		$stdev=0;
		foreach $line (@a){
			chomp($line);
			$ss=substr($line,20);
			@i1=split(/ /,$ss);
			$sum+=(($i1[0]-$avg)*($i1[0]-$avg));
		}

		$stdev=sqrt($sum/($c-1));
		print OUT2 "Stdev: $stdev\n";

		foreach $line (@a){
			chomp($line);
			#$id=substr($line,0,index($line," "));
			$id=substr($line,0,20);
			$id =~ s/\s//g;
			$ss=substr($line,20);
			@i1=split(/ /,$ss);
			$z=($i1[0]-$avg)/$stdev;
			if ($z<0) {$z=0;}

			print OUT1 "$id\t$z\n";

			$id1[$cc9]=$id;
			$zscore[$cc9]=$z;
			$cc9++;

			if (exists($h1{$id})) {
				$h1{$id}+=$z;
				$h2{$id}++;
			}
			else{
				$h1{$id}=$z;
				$h2{$id}=1;
			}
		}
		close OUT1;

		$name1=substr($name,0,index($name,"."))."_sorted.txt";
		open(OUT1, ">$outf/$name1") || die("Couldn't open file $outf/$name1\n");

		for ($i=0; $i<$cc9-1; $i++) {
			for ($j=$i+1; $j<$cc9; $j++) {
				if ($zscore[$i]<$zscore[$j]) {
					$t=$zscore[$j];
					$zscore[$j]=$zscore[$i];
					$zscore[$i]=$t;

					$t=$id1[$j];
					$id1[$j]=$id1[$i];
					$id1[$i]=$t;
				}
			}
		}

		for ($i=0; $i<$cc9; $i++) {
			print OUT1 "$id1[$i]\t$zscore[$i]\n";
		}
		close OUT1;
	}
}
closedir(DIR3);

@id=();
@zscore=();
@number=();

$cou=0;
while(my ($key, $val) = each(%h1)) { 
	$id[$cou]=$key;
	$zscore[$cou]=$val;
	$number[$cou]=$h2{$key};
	$cou++;
}

for ($i=0; $i<$cou-1; $i++) {
	for ($j=$i+1; $j<$cou; $j++) {
		if ($zscore[$i]<$zscore[$j]) {
			$t=$zscore[$j];
			$zscore[$j]=$zscore[$i];
			$zscore[$i]=$t;

			$t=$id[$j];
			$id[$j]=$id[$i];
			$id[$i]=$t;

			$t=$number[$j];
			$number[$j]=$number[$i];
			$number[$i]=$t;
		}
	}
}

for ($i=0; $i<$cou; $i++) {
	print OUT $id[$i]." ";
	print OUT $zscore[$i]." ";
	print OUT $number[$i]."\n";
}
close OUT;
close OUT2;

