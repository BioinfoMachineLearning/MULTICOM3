#!/usr/bin/perl -w
##################################################################################################
#Check the Ca-Ca distance clash in a single-chain pdb model
#Modified from get_atom.pl in prosys/script
#Author: Jianlin Cheng, 5/14/2006 
#################################################################################################

##############standard Amino Acids (3 letter <-> 1 letter)#######################################
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
###################################################################################################

#parse parameters
if (@ARGV != 1)
{
	die "Need four parameters:sequence file(fasta), pdb file\n"
}

$pdb_file =  shift @ARGV;
-f $pdb_file || die "pdb file doesn't exist.\n"; 


open(PDB, "$pdb_file") || die "can't read pdb file.\n";
@content = <PDB>;
close PDB;
##########################################################################

###########################Extract Information from ATOM RECORD###########
#Extract all ATOMS RECORDS
@records = ();
foreach $text(@content)
{
	if ($text =~ /^ATOM\s+/)
	{
		push @records, $text;
	}
}

@x_coords = ();
@y_coords = ();
@z_coords = ();

$i = 0;
while (@records)
	{
		$text = shift @records;
		#print "$text\n";
		#<STDIN>;
		$res = substr($text, 17, 3); 

		$atom_name = substr($text, 12, 4);

		if ($atom_name !~ /CA/)
		{
			next;
		}
		#$org_aa = substr($seq,$i,1);
		$i++;

		#extract the xyz coordinates of CA atom
		$xc = substr($text, 30, 8);
		$xc =~ s/\s//g;
		$yc = substr($text, 38, 8);
		$yc =~ s/\s//g;
		$zc = substr($text, 46, 8); 
		$zc =~ s/\s//g;

		push @x_coords, $xc;
		push @y_coords, $yc;
		push @z_coords, $zc;

		#conver to one letter
		$res = uc($res); 
		
		#check the residue length (3-letter or 1-letter)
		$res =~ s/\s+//g;
		if (length($res) == 3)
		{
			if (exists($amino{$res}) )
			{
				$aa = $amino{$res}; 
			}
			else
			{
				$aa = "X"; 
			}
		}
		else
		{
			$aa = $res;
		}

		# #check if the information match
		# if ($org_aa ne $aa)
		# {
		# 	die "amino acid doesn't match: $org_aa != $aa at pos = $i\n" ;
		# }
}


#check clash
for ($i = 0; $i < @x_coords; $i++)
{
	for ($j = $i + 1; $j < @x_coords; $j++)
	{
		$dist = sqrt(  ($x_coords[$i]-$x_coords[$j])*($x_coords[$i]-$x_coords[$j]) + ($y_coords[$i]-$y_coords[$j])*($y_coords[$i]-$y_coords[$j]) + ($z_coords[$i]-$z_coords[$j])*($z_coords[$i]-$z_coords[$j]) );

		$pos1 = $i + 1;
		$pos2 = $j + 1;
		# $aa1 = substr($seq, $i, 1);
		# $aa2 = substr($seq, $j, 1);

		if ($dist < 0.1)
		{
			print "overlap: ($pos1) <-> ($pos2), dist = $dist\n";
		}
		elsif ($dist < 1.9)
		{
			print "servere clash: ($pos1) <-> ($pos2), dist = $dist\n";
		}
		elsif ($dist < 3.5)
		#elsif ($dist <= 3.5)
		{
			print "clash: ($pos1) <-> ($pos2), dist = $dist\n";
		}
		elsif ($dist > 4.5 && abs($pos1-$pos2) == 1)
		{
			print "chain broken: ($pos1) <-> ($pos2), dist = $dist\n";
		}
	}
}

