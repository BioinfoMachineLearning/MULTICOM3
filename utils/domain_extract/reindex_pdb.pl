#!/usr/bin/perl -w

##########################################################################
# Author: Jie Hou
# Date: 06/16/2015
##########################################################################

use Carp;
use File::Basename;
$num = @ARGV;
if($num != 2)
{
	print "The number of parameter is not correct!\n";
	exit(1);
}

$file_PDB = $ARGV[0];
$outputPDB = $ARGV[1];

	
# Load PDB
open INPUTPDB, $file_PDB or die "ERROR! Could not open $file_PDB";
my @lines_PDB = <INPUTPDB>;
close INPUTPDB;

# (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
my $resCounter = 0;
my $atomCounter = 0;
my $prevrNum = "XX";
open OUTPDB, ">${outputPDB}" or die "ERROR! Could not open ${outputPDB}";
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
	die "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}