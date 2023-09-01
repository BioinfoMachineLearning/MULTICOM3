#!/usr/bin/perl -w
#################################################################
#Extract sequence from ATOM records of pdb file
#Input: pdb file
#Output: sequence in the pdb file
#Author: Jianlin Cheng
#Date: 12/30/2005
#################################################################

##############standard Amino Acids (3 letter <-> 1 letter)#######
%amino=();
$amino{"ALA"} = 'A';
$amino{"CYS"} = 'C';
$amino{"ASP"} = 'D';
$amino{"GLU"} = 'E';
$amino{"PHE"} = 'F';
$amino{"GLY"} = 'G';
$amino{"HIS"} = 'H';
$amino{"ILE"} = 'I';
$amino{"LYS"} = 'K';
$amino{"LEU"} = 'L';
$amino{"MET"} = 'M';
$amino{"ASN"} = 'N';
$amino{"PRO"} = 'P';
$amino{"GLN"} = 'Q';
$amino{"ARG"} = 'R';
$amino{"SER"} = 'S';
$amino{"THR"} = 'T';
$amino{"VAL"} = 'V';
$amino{"TRP"} = 'W';
$amino{"TYR"} = 'Y';
###################################################################

#read sequence from atom file
sub get_seq_from_atom
{
	#assume the atom file exists
	my $file = $_[0];
	#sleep 5;
	open(ATOM, $file) || die "can't read atom file: $file\n";
	my @atoms = <ATOM>;
	close ATOM; 
	my $prev = -1;
	my $seq = ""; 
	my %chain=();
	$chain{'A'} = 1;
	while (@atoms)
	{
		my $text = shift @atoms;
		if ($text =~ /^ATOM/)
		{
			$ch = substr($text, 21, 1);
			if(!exists($chain{$ch}))
			{
				die "\n\n !!!!!!!!!!! Warning: Pdb file $file has multiple chains, please remove it!\n\n\n";
			}
		}
	}
	return $seq; 
}

if (@ARGV != 1)
{
	die "need two parameters: pdb file, and the chain ID.\n"; 
}

$pdb_file = shift @ARGV;

#sleep 5;
 
$seq = &get_seq_from_atom($pdb_file);
#print "$seq\n";

