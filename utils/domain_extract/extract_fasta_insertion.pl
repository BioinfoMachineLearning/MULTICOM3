if(@ARGV != 3)
{
        print "the number of parameters is not correct!\n";
        exit(1);
}

$fasta_file = $ARGV[0];
$domain_info = $ARGV[1];
$output_dir = $ARGV[2];

open FASTA, $fasta_file or die "ERROR! Could not open $fasta_file";
my @lines = <FASTA>;
close FASTA;

open(IN,"$domain_info");
@content = <IN>;
close IN;
$domain_num = @content;

$i = 0;
foreach $info (@content)
{
        @domain_start = ();
        @domain_end   = (); 
        #Domain 1:1-24,92-139,213-232
        chomp $info;
	@tmp = split(':',$info);
	$range = $tmp[1];
        @domain_infos = split(',',$range);
	foreach $domain_info (@domain_infos)
        {
                @domain_range = split('-',$domain_info);
                $start = $domain_range[0];
                $end = $domain_range[1];
                push @domain_start, $start; 
	        push @domain_end, $end; 
                print "domain$i: $start-$end\n";
        }
        $outfile = "$output_dir/domain$i.fasta";
        open OUT, ">$outfile" or die "ERROR! Could not open $outfile";
        print OUT ">domain$i\n";
        foreach (@lines) 
        {
                if (substr($_, 0, 1) eq '>')
                {
                        next;
                }
                $new_seq = "";
                for ($j = 0; $j < @domain_start; $j++)
		{
                        $new_seq = $new_seq.substr($_, $domain_start[$j]-1, $domain_end[$j]-$domain_start[$j]+1);
                }    
                print OUT $new_seq;
        }
        print OUT "\n";
        close OUT;
        $i++;
}
