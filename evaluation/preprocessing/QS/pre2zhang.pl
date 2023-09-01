#!/usr/bin/perl -w
##################################################################################################
#convert CASP predicted model to model matching with the final true structure file (zhang's target)
#Author: Jianlin Cheng, 9/5/2006, MODIFIED by Zheng Wang, July 21, 2010, after CASP9 
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
$amino{"MSE"} = 'M';
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
if (@ARGV != 3)
{
	die "Need three parameters: casp prediction file, zhang target file, output file\n"
}

$casp_file = shift @ARGV;
$target_file = shift @ARGV;
$out_file = shift @ARGV;

#read target file  TEMPLATE 
open(TARGET, $target_file) || die "can't read target file.\n";
@target = <TARGET>;
close TARGET;


@order = ();
while (@target)
{
	$record = shift @target;
	chomp $record;
	if ($record !~ /^ATOM/)
	{
		next;
	}
	@fields = split(/\s+/, $record);
	#$aa = $fields[3];
	$aa = substr($record, 17, 3);
	$aa =~ s/ //g;
	if (defined($amino{$aa}))
	{
		$aa = $amino{$aa}; 
	}
	else
	{
		die "$aa is not found\n";
	}

#	push @atoms, {
#		new_idx => $fields[1],
#		aa => $aa,
#		org_idx => $fields[4],
#		x => $fields[5],
#		y => $fields[6],
#		z => $fields[7]
#	};
#	push @order, $fields[4];

	my $new_idx_z = substr($record, 6, 5);    #z represent Zheng;
        my $org_idx_z = substr($record, 22, 4);
        my $x_z = substr($record, 30, 8);
        my $y_z = substr($record, 38, 8);
        my $z_z = substr($record, 46, 8);

       push @atoms, {
	       new_idx => $new_idx_z,
               aa => $aa,
               org_idx => $org_idx_z,
               x => $x_z,
               y => $y_z,
               z => $z_z
       };
       push @order, $org_idx_z;

}

#read CASP predicted model    NON-TEMPLATE, CHANGING EVERYTIME
open(CASP, $casp_file) || die "can't read $casp_file\n";
@casp = <CASP>;
close CASP;

@selected = ();
while (@casp)
{
	$line = shift @casp;
	chomp $line;
	#print "original".$line."\n";
	if ($line =~ /^ATOM/)
	{
		@fields = split(/\s+/, $line);
		#$aa = $fields[3];
	        $aa = substr($line, 17, 3);
	        $aa =~ s/ //g;
		if (defined($amino{$aa}))
		{
			$aa = $amino{$aa}; 
		}
		else
		{
			die "$aa is not found\n";
		}
	        my $new_idx_z = substr($line, 6, 5);    #z represent Zheng;
	        my $org_idx_z = substr($line, 22, 4);
        	my $x_z = substr($line, 30, 8);
        	my $y_z = substr($line, 38, 8);
        	my $z_z = substr($line, 46, 8);
		#$idx = $fields[4];
		$idx = $org_idx_z;

		#check if the amino acid is in the true pdb file
		for ($i = 0; $i < @atoms; $i++)
		{
			$ord = $atoms[$i]{"org_idx"};
			#print "idx is $idx\n ord is $ord\n\n";
			#sleep 1;
			if ($idx == $ord)
			{
				#check if amino acid match
				#print '$aa(NON-template):'.$aa."\n";
				#print 'template PDB:'.$atoms[$i]{"aa"}."\n";
				#print '$i:'.$i."\n";	
				#$aa eq $atoms[$i]{"aa"} || die "amino acid doesn't match.\n";
				
				if($aa eq $atoms[$i]{"aa"}){
					#print "printed:".$line."\n";
					#sleep 2;
					push @selected, $line;
				}
				else{
					@selected = ();
					my $pdb_aa = $atoms[$i]{"aa"};
					my $pdb_aa_long = $atoms[$i];
					#print "AA mot mach $pdb_aa $pdb_aa_long \n";
					#die;
					push(@selected, "AMINO ACID DOES NOT MATCH!");
				}
				last;
			}
		}
	}
}

open(OUT, ">$out_file") || die "can't create output file.\n";

print OUT join("\n", @selected), "\n";

close OUT;

