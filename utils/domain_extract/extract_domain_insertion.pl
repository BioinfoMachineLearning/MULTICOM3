$numArgs = @ARGV;
if($numArgs != 3)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$pdbfile	= $ARGV[0];
$domain_info = $ARGV[1];
$output_file = $ARGV[2];

open(IN,"$domain_info");
@content = <IN>;
close IN;
$domain_num = @content;

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
		open INPUTPDB, $pdbfile or die "ERROR! Could not open $pdbfile";
		my @lines_PDB = <INPUTPDB>;
		close INPUTPDB;

		my @new_PDBlines;
		foreach (@lines_PDB) {
			next if $_ !~ m/^ATOM/;
			my $this_rnum = parse_pdb_row($_,"rnum");
			for ($j = 0; $j < @domain_start; $j++)
			{
				if($this_rnum>=($domain_start[$j]) and $this_rnum<=($domain_end[$j]))
				{
					push @new_PDBlines, $_;
				}
            } 
		}
		# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
		my $resCounter = 0;
		my $atomCounter = 0;
		my $prevrNum = "XX";
		open OUTPDB, ">$output_file" or die "ERROR! Could not open $output_file";
		foreach (@new_PDBlines) {
			next if $_ !~ m/^ATOM/;
			my $this_rnum = parse_pdb_row($_,"rnum");
			if ($prevrNum ne $this_rnum) {
				$prevrNum = $this_rnum;
				$resCounter++;
			}
			$atomCounter++;
			my $rnum_string = sprintf("%4s", $resCounter);
			my $anum_string = sprintf("%5s", $atomCounter);
		#		my $row = substr($_,0,6).$anum_string.substr($_,11,10)." ".$rnum_string." ".substr($_,27);
			my $row = substr($_,0,6).$anum_string.substr($_,11,5)." ".substr($_,17,3)." "." ".$rnum_string." ".substr($_,27);
			print OUTPDB $row;
		}
		print OUTPDB "END\n";
		close OUTPDB;
}

sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	print "Warn: Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}
