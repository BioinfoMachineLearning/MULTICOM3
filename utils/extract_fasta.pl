if(@ARGV != 5)
{
        print "the number of parameters is not correct!\n";
        exit(1);
}

$infile = $ARGV[0];
$outfile = $ARGV[1];
$targetname = $ARGV[2];
$start  = "$ARGV[3]";
$end = "$ARGV[4]";

if($start > $end)
{
        die "wrong index in extract_domain.pl <start:$start,  end:$end>\n";
}

$start = $start-1;
$end = $end-1;

open INPUTFASTA, $infile or die "ERROR! Could not open $infile";
my @lines_fasta = <INPUTFASTA>;
close INPUTFASTA;

open OUTFASTA, ">$outfile" or die "ERROR! Could not open $outfile";
print OUTFASTA ">$targetname\n";
my @new_fasta;
foreach (@lines_fasta) {
        if (substr($_, 0, 1) eq '>')
	{
		next;
	} 
	$new_seq = substr($_, $start, $end-$start+1);
        print OUTFASTA $new_seq;
}
print OUTFASTA "\n";
close OUTFASTA;

print length($new_seq);
