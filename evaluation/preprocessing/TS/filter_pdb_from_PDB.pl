#!/usr/bin/perl -w
############################################
#Read the PDB IDs and target ids from the 
#two list files, pdb-ids, and target_ids.
#then ftp PDB to get the pdb files.Then unzip, and do
#others listed in the readme.txt file.
#
#Author: Zheng Wang, July 8th, 2008.
############################################


use strict;
use Net::FTP;

if (@ARGV != 6){
        die ("Need five parameters: the pdb_ids file, the target_ids file, the location to save the pdb original files, the path of folder of casp original sequence files, the path of alignment folder, the path of filtered pdb files!");
}

my $pdb_file = $ARGV[0];
my $target_file = $ARGV[1];
my $path = $ARGV[2];
my $casp_orig_seq_path = $ARGV[3];
my $align_path = $ARGV[4];
my $pdb_filtered_path = $ARGV[5];

my $clustalw_tool = "../tools/clustalw1.83/clustalw";

open(PDBID, "<$pdb_file");
open(TID, "<$target_file");
my @pdb_ids;
my @target_ids;

while(<PDBID>){
	my $line = $_;
	$line =~ s/\s+$//;
	push(@pdb_ids, $line);
}


while(<TID>){
        my $line = $_;
        $line =~ s/\s+$//;
	push(@target_ids, $line);
}

my $host = "ftp.wwpdb.org"; 
#my $ftp = Net::FTP->new("ftp://ftp.rcsb.org", Debug => 0)
#my $ftp = Net::FTP->new($host, Debug => 1) #this one doesn't work any more
#      or die "Cannot connect to ftp.wwpdb.org/: $@";

#$ftp->login("anonymous",'-anonymous@')
#      or die "Cannot login ", $ftp->message;

#$ftp->cwd("/pub/pdb/data/structures/all/pdb/")
#      or die "Cannot change working directory ", $ftp->message;

my $length = scalar(@pdb_ids);
my $counter = 0;
while($counter < $length){
        my $file = "pdb$pdb_ids[$counter]".'.ent.gz';
        my $file_new = "$target_ids[$counter]"."pdb$pdb_ids[$counter]".'.ent.gz';
	print $target_ids[$counter]."\n";
	my $file_unzipped = "$target_ids[$counter]"."pdb$pdb_ids[$counter]".'.ent';
	print $file."\n";
    #$ftp->get("$file") or die "cannot get the $file from pdb ftp website!\n", $ftp->message;  #need to uncommentted to properly run the program!! something wrong with PDB server, so use the line below.
	#`cd ../pdb_orig/`;
	#chdir("../pdb_orig/");
	if(-e "$file")
	{
		print "$file is found~!\n";
	}else{
		print "Downloading $file from pdb ftp\n";
		print("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/$file\n\n");   # "##" commented lines should be removed if you want to run the program.
		system("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/$file");   # "##" commented lines should be removed if you want to run the program.
      	}
	#`mv $file $path/$file_new`;   #copy the ftp downloaded pdb original files into ./pdb_orig
    `cp $file $path/$file_new`;   #copy the ftp downloaded pdb original files into ./pdb_orig
	`gzip -d $path/$file_new`;  #unzip the gz files (pdb files);
	my $pdb_orig_seq_file = "$target_ids[$counter]".'_pdb_orig_seq';
	my $status = system("perl check_chain_num.pl $path/$file_unzipped"); #Added by jie to check if pdb file has multiple chain, because in this case we don't know wether casp sequence is matched to which chains
	if($status)
	{
		print "\n\n !!!!!!!!!!! Warning: Pdb file $file_unzipped has multiple chains, please remove it!\n\n\n";
	}

	my $casp_orig_seq_file = $casp_orig_seq_path.'/'.$target_ids[$counter].'.fasta';
	if (!-e "$casp_orig_seq_file")
	{
		print "$casp_orig_seq_file is not found !!!\n";
	}
	my $casp_orig_seq;
	my $pdb_orig_seq;
	open(CASPSEQ, "<$casp_orig_seq_file");
	while(<CASPSEQ>){
		my $line = $_;
		if (substr($line, 0, 1) ne ">"){
			$line =~ s/\s+$//;
			$line =~ s/\n//;
			$casp_orig_seq = $line;
		}
	}
	close(CASPSEQ);

	#if($ratio > 0.85 || length($pdb_orig_seq) == 0){
	my $chains = "ABCDEFGHIJKLMNOPQRSTUVXWYZabcdefghijklmnopqrstuvxwyz";
	my $chain_count = 0;
	my $min_ratio = 1;
	my $min_chain_idx = 0;
	my $casp_align_seq = "";
    my $pdb_align_seq = "";
	my $casp_align_seq_min = "";
    my $pdb_align_seq_min = "";

	while(1)
	{
		my $pdb_orig_seq_file = "$target_ids[$counter]".'_pdb_orig_seq';
		my $chain_id = substr($chains, $chain_count, 1);
        `perl pdb2seq.pl $path/$file_unzipped $chain_id > $path/$pdb_orig_seq_file`; #parse the sequence from the downloaded pdb files;
		print("perl pdb2seq.pl $path/$file_unzipped $chain_id > $path/$pdb_orig_seq_file\n");
		open(PDBSEQ, "<$path/$pdb_orig_seq_file");
        while(<PDBSEQ>){
            my $line = $_;
            $line =~ s/\n//;
            $pdb_orig_seq = $line;
		}
        close(PDBSEQ);
		
		if(length($pdb_orig_seq) < 3){
            last;
        }

		my $align_seq_file = $target_ids[$counter].'_align.seq';
        open(ALIGNSEQ, ">$align_path/$align_seq_file");
        print ALIGNSEQ '%CASP'."\n";
        print ALIGNSEQ "$casp_orig_seq\n";
        print ALIGNSEQ '%PDB'."\n";
        print ALIGNSEQ "$pdb_orig_seq\n";
        close(ALIGNSEQ);

        my $align_out_file = $align_path.'/'.$target_ids[$counter].'_align.out';
        `$clustalw_tool  -MATRIX=BLOSUM -TYPE=PROTEIN -INFILE=$align_path/$align_seq_file -OUTFILE=$align_out_file`;
        
		#Parse the alignment results;
        $casp_align_seq = "";
        $pdb_align_seq = "";
        open(ALIGN, "<$align_out_file");
        while(<ALIGN>){
           	my $line = $_;
           	if($line =~ /CASP/){
               	$line =~ s/\n//;
               	$line =~ s/CASP//;
               	$line =~ s/ //g;
               	$casp_align_seq = $casp_align_seq . $line;
           	}
          	if($line =~ /PDB/){
               	$line =~ s/\n//;
                $line =~ s/PDB//;
                $line =~ s/ //g;
                $pdb_align_seq = $pdb_align_seq . $line;
        	}

        }
	    close(ALIGN);
		print $casp_align_seq."\n";
		print $pdb_align_seq."\n";
		# add by Jie
		#If the identity of two sequences are smaller than a threshold, change to another chain. For example, Chain B.

		my $length_seq = length($casp_align_seq);
		my $counter_seq = 0;
		my $no_diff = 0;
		while($counter_seq < $length_seq){
			if(substr($casp_align_seq, $counter_seq, 1) ne substr($pdb_align_seq, $counter_seq, 1) && substr($casp_align_seq, $counter_seq, 1) ne "-" && substr($pdb_align_seq, $counter_seq, 1) ne "-"){
				$no_diff++;
			}
			$counter_seq++;
		}
		my $ratio = $no_diff / $length_seq;
		print "\n\t\t Ratio $ratio\n\n";	

		if ($ratio < $min_ratio){
			$min_ratio = $ratio;
			$min_chain_idx = $chain_count;
			$casp_align_seq_min = $casp_align_seq;
			$pdb_align_seq_min = $pdb_align_seq;
		}
		$chain_count += 1;
	}

	my $chain = substr($chains, $min_chain_idx, 1);
	print "Chain is:".$chain."\n";
	#sleep 3;

	#Filter the original PDB file, only keep the ATOM line, and not eliminate the ATOM lines with aminio acid code has 4 digits, also eliminate the "chain column" from the orig PDB file.
	my $pdb_orig_filtered_file = $path.'/'.$target_ids[$counter].'_pdb_1st_filtered';
	open(PDBORIG, "<$path/$file_unzipped");
	open(PDB1F, ">$pdb_orig_filtered_file");
	while(<PDBORIG>){
		my $line = $_;
		if(substr($line, 0, 4) eq "ATOM"){
			my @items = split(/\s+/, $line);
			if((length($items[3]) == 3 or length($items[3]) == 4 )&& substr($line, 21, 1) eq $chain){
				#print PDB1F $line;
				my $printed_line = substr($line, 0, 20)."  ".substr($line, 22);
				print PDB1F $printed_line;
				#print $printed_line."\n1\n";
				#sleep 1;
			}	
			if(substr($line, 16, 1) eq $chain && substr($line, 21, 1) eq $chain){
				my $printed_line = substr($line, 0, 16)." ".substr($line, 17, 4)." ".substr($line, 22);
				print PDB1F $printed_line;
				#print $printed_line."\n2\n";
				#sleep 1;
			}
		}
	}
	close(PDB1F);
	close(PDBORIG);
	
	#Filter the PDB file, order based on CASP order;
	my $len_align_seq = length($casp_align_seq);
	my $counter_align_seq = 0;
	my $counter_casp_seq = 0;
	my $counter_pdb_seq = 0;
	my $counter_atom = 1;
	my $pdb_filtered_file = $target_ids[$counter].'_filtered.pdb';
	open(PDBFW, ">$pdb_filtered_path/$pdb_filtered_file");
	while($counter_align_seq < $len_align_seq){
		my $cell_casp = substr($casp_align_seq_min, $counter_align_seq, 1);
		my $cell_pdb = substr($pdb_align_seq_min, $counter_align_seq, 1);
		if($cell_casp ne "-"){
			$counter_casp_seq++;
		}
		if($cell_pdb ne "-"){
			$counter_pdb_seq++;
			#print '$counter_pdb_seq: '.$counter_pdb_seq."\n";
		}
		if($cell_casp ne "-" && $cell_pdb ne "-" && $cell_casp eq $cell_pdb){
			#print $cell_casp." ".$cell_pdb."\n";
			#sleep 2;
			open(PDBFR, "<$pdb_orig_filtered_file");
			#my $pdb_filtered_file = $target_ids[$counter].'_filtered.pdb';
			#open(PDBFW, ">$pdb_filtered_path/$pdb_filtered_file");
			my $counter_pdb_num = 0;
			my $pre_num = -1;  # has to be -1 or less, cause pdb starts from 0 or 1;
			while(<PDBFR>){
				my $line = $_;
				if(substr($line, 0, 4) eq "ATOM"){
					my $num = substr($line, 23, 5);
					#print '$num: '.$num."\n";
					#sleep 1;
					if($num != $pre_num){
						$counter_pdb_num++;
						#print '$counter_pdb_num: '.$counter_pdb_num."\n";
						$pre_num = $num;
					}
					if($counter_pdb_num == $counter_pdb_seq){
						#print $counter_pdb_num."\n";
						#print $counter_pdb_seq."\n";
						my $print_casp_order;
						#print "print\n";
						if(length($counter_casp_seq) == 1){
							$print_casp_order = "   $counter_casp_seq";
						}
                                                if(length($counter_casp_seq) == 2){
                                                        $print_casp_order = "  $counter_casp_seq";
                                                }
                                                if(length($counter_casp_seq) == 3){
                                                        $print_casp_order = " $counter_casp_seq";
                                                }
                                                if(length($counter_casp_seq) == 4){
                                                        $print_casp_order = "$counter_casp_seq";
                                                }
                                          
						#print $counter_casp_seq."\n";

						my $print_atom_order;
						if(length($counter_atom) == 1){
							$print_atom_order = "    $counter_atom";
						}
                                                if(length($counter_atom) == 2){
                                                        $print_atom_order = "   $counter_atom";
                                                }
                                                if(length($counter_atom) == 3){
                                                        $print_atom_order = "  $counter_atom";
                                                }
                                                if(length($counter_atom) == 4){
                                                        $print_atom_order = " $counter_atom";
                                                }
                                                if(length($counter_atom) == 5){
                                                        $print_atom_order = "$counter_atom";
                                                }

						print PDBFW substr($line, 0, 6).$print_atom_order.substr($line, 11, 11).$print_casp_order.substr($line, 26);
						#print substr($line, 0, 6).$print_atom_order.substr($line, 11, 11).$print_casp_order.substr($line, 26);
						$counter_atom++;

					}
				}
			}
			close(PDBFR);
                }
	
		$counter_align_seq++;
	}
	close(PDBFR);


	$counter++;
}
#$ftp->quit;

close(TID);
close(PDBID);


