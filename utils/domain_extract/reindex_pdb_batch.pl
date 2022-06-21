$numArgs = @ARGV;
if($numArgs != 2)
{   
	print "the number of parameters is not correct!\n";
	exit(1);
}

$pdb_folder	= "$ARGV[0]";
$newpdb_folder	= "$ARGV[1]";

opendir(DIR,"$pdb_folder");
@files = readdir(DIR);
closedir(DIR);

foreach $pdb (@files)
{
	chomp $pdb;
	#if(substr($pdb,length($pdb)-4) ne '.pdb')
	#{
	#	next;
	#}
	if($pdb eq "." || $pdb eq "..")
	{
		next;
	}
	$pdbfile = "$pdb_folder/$pdb";
	$newpdb = "$newpdb_folder/$pdb";
	open INPUTPDB, $pdbfile or die "ERROR! Could not open $pdbfile";
	my @lines_PDB = <INPUTPDB>;
	close INPUTPDB;

	# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
	my $resCounter = 0;
	my $atomCounter = 0;
	my $prevrNum = "XX";
	open OUTPDB, ">$newpdb" or die "ERROR! Could not open $newpdb";
	foreach (@lines_PDB) {
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
