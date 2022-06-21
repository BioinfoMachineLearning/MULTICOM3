if(@ARGV != 5)
{
        print "the number of parameters is not correct!\n";
        exit(1);
}

$infile = $ARGV[0];
$outfile = $ARGV[1];
$targetname = $ARGV[2];
$rstart  = "$ARGV[3]";
$rend = "$ARGV[4]";

if($rstart > $rend)
{
        die "wrong index in extract_domain.pl <start:$start,  end:$end>\n";
}

$rstart = $rstart-1;
$rend = $rend-1;

open INPUTFASTA, $infile or die "ERROR! Could not open $infile";
my @lines_fasta = <INPUTFASTA>;
close INPUTFASTA;

open OUTFASTA, ">$outfile" or die "ERROR! Could not open $outfile";
print OUTFASTA ">$targetname\n";
my @new_fasta;
foreach (@lines_fasta) {
	$line = $_;
	chomp $line;
        if (substr($line, 0, 1) eq '>')
	{
		next;
	} 
	$new_seq1 = substr($line, 0, $rstart);
	$new_seq2 = substr($line, $rend + 1);
	$new_seq = "$new_seq1"."$new_seq2";
	#$new_seq = substr($_, $start, $end-$start+1);
        print OUTFASTA $new_seq;
}
print OUTFASTA "\n";
close OUTFASTA;

print length($new_seq);
