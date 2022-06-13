$num = @ARGV;
if($num !=5)
{
	die "Parameters are not correct!\n";
}

$fastafile = $ARGV[0];
$rrfile = $ARGV[1];
$outfile = $ARGV[2];
$start = $ARGV[3];
$end = $ARGV[4];


open(IN,$fastafile) || die "Failed to open file $fastafile\n";
@content = <IN>;
close IN;
shift @content;
$seq = shift @content;
chomp $seq;

open(IN,$rrfile) || die "Failed to open file $rrfile\n";
@rrcontent = <IN>;
close IN;

$seq_rr = shift @rrcontent;

if($seq_rr != $seq)
{
	die "Seuqence is not matching!\n";
}

$domain_seq = substr($seq, $start-1, $end-$start+1);

$out_tmp = $outfile.".ori";

open(OUT,">$out_tmp")  || die "Failed to open file $out_tmp\n";
open(OUT2,">$outfile")  || die "Failed to open file $outfile\n";

print OUT $domain_seq, "\n";
print OUT2 $domain_seq, "\n";

foreach (@rrcontent)
{
	$line=$_;
	chomp $line;
	
	my @fields = split(/\s+/, $line);

	$x = $fields[0];
	$y = $fields[1];

	if(@fields == 5)
    {
		$prob = $fields[4];
	}
	elsif(@fields == 3)
	{
		$prob = $fields[2];
	}

	if($x >= $start and $x <= $end && $y >= $start && $y <= $end)
	{
		print OUT $x, " ", $y, " ", "0 8 ", $prob, "\n";
		###reindex
		print OUT2 $x-$start+1, " ", $y-$start+1, " ", "0 8 ", $prob, "\n";
	}
}
close OUT;
close OUT2;

