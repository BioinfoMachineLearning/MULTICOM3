$num = @ARGV;
if($num !=5)
{
	die "Parameters are not correct!\n";
}

$fastafile = $ARGV[0];
$distfile = $ARGV[1];
$outfile = $ARGV[2];
$start = $ARGV[3];
$end = $ARGV[4];


open(IN,$fastafile) || die "Failed to open file $fastafile\n";
@content = <IN>;
close IN;
shift @content;
$seq = shift @content;
chomp $seq;

open(IN,$distfile) || die "Failed to open file $distfile\n";
@mapcontent = <IN>;
close IN;

if(@mapcontent != length($seq))
{
	die "The length is not matching!\n";
}
open(OUT,">$outfile")  || die "Failed to open file $outfile\n";

$c=0;
foreach (@mapcontent)
{
	$line=$_;
	chomp $line;
	$c++;
	if($c>=$start and $c<=$end)
	{
		@tmp = split(/\s+/,$line);
		for($i = $start;$i<=$end;$i++)
		{
			if($i==$start)
			{
				print OUT $tmp[$start-1];
			}else{
				print OUT " ".$tmp[$i-1];
			}
		}
		print OUT "\n";
	}
}
close OUT;
