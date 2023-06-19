###################################################################
#Given the models from CASP, the filtered pdb files, evaluate
#CASP models. Top 1 and Top 5.
#
#Difference between this version 2 and previous version is that
#This one won't run TM-score if the output of TM-score has already
#been generated. Save time.
#
#Author: Zheng Wang, July 18th, 2008.
##################################################################

use strict;
if(@ARGV != 7){
	die "Need 6 parameters: The path of filtered PDB files, the path of CASP models, the temp folder, the path of the tm-score program, the path of the output files, the target list file, and the folder for saving target specific results!";
}

my $pdb_filtered_dir = $ARGV[0];
my $casp_models_dir = $ARGV[1];
my $temp_dir = $ARGV[2];
my $tm_score_path = $ARGV[3];
my $output_folder = $ARGV[4];
my $target_list = $ARGV[5];
my $target_results = $ARGV[6]; 

`rm $target_results/*`;

my $log_file = "log.log";
open(LOG, ">$log_file");

my %rank_1;
my %rank_5;

my %rank_1_no;

my %rank_1_TM;
my %rank_1_MaxSub;

my %rank_5_TM;
my %rank_5_MaxSub;


#Get the group names;
opendir(DIR, "$casp_models_dir") or die("Cannot open dir!");
my @targets = readdir DIR;
closedir DIR;
foreach my $target (@targets){
	if($target eq "." || $target eq ".." || index($target,'.tar.gz')>0){
		next;
	}
	opendir(DIR, "$casp_models_dir/$target") or die ("canoot open $casp_models_dir/$target!");
	my @models = readdir DIR;
	close(DIR);
	foreach my $model (@models){
		if($model eq "." || $model eq ".."){
			next;
		}
		print $model;
		my $group = substr($model, 0, index($model, ".pdb")-2);#substr($model, 0, (length($model) - 4));
		if($rank_1{$group} ne "BLANK"){
			print "group is:".$group."\n";
			$rank_1{$group} = "BLANK";
			$rank_5{$group} = "BLANK";

			$rank_1_no{$group} = 0;

			$rank_1_TM{$group} = "BLANK";
			$rank_1_MaxSub{$group} = "BLANK";

			$rank_5_TM{$group} = "BLANK";
			$rank_5_MaxSub{$group} = "BLANK";

		}
	}
	# last;   #only do that on one target. The purpose is to get the name of the group.
}

#open(TLIST, "<./target_ids") or die ("Cannot find target_ids file! put it in the same folder as this program!");
open(TLIST, "<$target_list") or die ("Cannot find target_ids file! put it in the same folder as this program!");
while(<TLIST>){
	my $line = $_;
	$line =~ s/\s+$//;
	#print $line."\n";
	my $target_id = $line;
	$target_id =~ s/\n//;
	$target_id =~ s/ //;
	if(!$target_id){
		next;
	}
        my $target_results_file_g = "$target_results/$target_id.results";
        open(TRESULTSG, ">>$target_results_file_g");
	my %rank_1_TM_target;
        my %to_print_g_hash;
	print "Evaluating Target ".$target_id."\n\n";
	print LOG "Evaluating Target ".$target_id."\n\n";
	#sleep 2;
	`mkdir $temp_dir/$target_id`;
	my $key; #key is the group name;
	my $value;
	open(GDTS, "> $temp_dir/$target_id/gdt_scores_all_groups");	

	print "Evaluating the first model!\n\n";
	#sleep 2;
	while(($key, $value) = each(%rank_1)) {
		#get the rank of rank_1 first, ie the number 1 model;
		#my @key_models;  #store five models of a group;
		print $key . "\n";
		#sleep 1;
		if ($key =~ /BAKER/ || $key =~ /Zhang/) {
			print LOG "$key \n";
		}
		my $rank_1_no_model = 0;
		#print $casp_models_dir.$target_id."\n";
		opendir(DIR, "$casp_models_dir/$target_id") or print "cannot open dir, predictions may be missing in this target!";
		my @models = readdir DIR;
		closedir DIR;

		my $model = $key . '_0.pdb';
		#the number 1 model;
		#print "3\n";
		#print $temp_dir."\n";
		#print "mkdir $temp_dir/$target_id";
		#sleep 5;
		#`mkdir $temp_dir/$target_id`;
		my $pdb_filtered_file = $pdb_filtered_dir . '/' . $target_id . '_filtered.pdb';
		my $model_filtered = $temp_dir . '/' . $target_id . '/' . $model . '_filtered';
		`cp $casp_models_dir/$target_id/$model $model_filtered`;
		my $tm_result = $model_filtered . '_out';
		#print "$tm_score_path $model_filtered $pdb_filtered_file > $tm_result\n";
		unless (-e $tm_result) {
			`$tm_score_path $pdb_filtered_file $model_filtered > $tm_result`;
		}
		#print "$tm_score_path $model_filtered $pdb_filtered_file > $tm_result\n";
		open(POUT, "<$tm_result");
		while (<POUT>) {
			my $line = $_;
			if ($line =~ /There is no common residues in the input structures/) {
				#this is a partial models, deal with it later;
				print "This is a partial models!\n";
				last;
			}
			if (substr($line, 0, 8) eq "TM-score" and index($line, "Structure_1:") > 0) {
				my $tm_score = substr($line, 10, 8);
				#print $gdt_score."\n";
				if ($tm_score) {
					#print $rank_1{$key}."\n";
					$rank_1_TM{$key} = $rank_1_TM{$key} + $tm_score;
					$rank_1_TM_target{$key} = $tm_score;
					print "TM: " . $rank_1_TM{$key} . "\n";
					last;
				}
			}
		}
		close(POUT);
		$rank_1_no{$key} = $rank_1_no{$key} + $rank_1_no_model;
		print "Target Count: " . $rank_1_no{$key} . "\n";
		#sleep 1;
	}
	print "Evaluating the best models of five!\n\n";
	#sleep 2;
	while (($key, $value) = each(%rank_5)) {
		#get the sum of the best model in 5 models;
		my $i = 0;
		print $key . "\n";
		my $to_print_g;
		my $target_results_file = "$target_results/$key.results";
		open(TRESULTS, ">>$target_results_file");
		my $to_print_t = "$target_id ";
		my $highest_score = 0;
		my $highest_TM = 0;
		my $highest_SUB = 0;
		while ($i < 5) {
			my $model = $key . '_' . $i.'.pdb';

			#print $name."\n";
			#sleep 2;

			my $pdb_filtered_file = $pdb_filtered_dir . '/' . $target_id . '_filtered.pdb';
			my $model_filtered = $temp_dir . '/' . $target_id . '/' . $model . '_filtered';
			`cp $casp_models_dir/$target_id/$model $model_filtered`;

			my $tm_result = $model_filtered . '_out';
			#print "$tm_score_path $model_filtered $pdb_filtered_file > $tm_result\n";
			unless (-e $tm_result) {
				`$tm_score_path $pdb_filtered_file $model_filtered > $tm_result`;
			}
			open(POUT, "<$tm_result");
			while (<POUT>) {
				my $line = $_;
				if ($line =~ /There is no common residues in the input structures/) {
					#this is a partial models, deal with it later;
					print "This is a partial models!\n";
					last;
				}
				if (substr($line, 0, 8) eq "TM-score" and index($line, "Structure_1:") > 0) {
					my $tm_score = substr($line, 10, 8);
					$to_print_t .= $tm_score . " ";
					$to_print_g .= $tm_score . " ";
					if ($tm_score) {
						if ($tm_score > $highest_TM) {
							$highest_TM = $tm_score;
						}
						#print GDTS $name." ".$gdt_score."\n";
					}
					print "SingleTM:" . $tm_score . "\n";
					#print "Highest:".$highest_score."\n";
					last;
				}
			}
			close(POUT);
			$i++;
		}
		#print "Previious:".$rank_5{$key}."\n";
		$to_print_g_hash{$key} = $to_print_g;
		print TRESULTS $to_print_t . "\n";
		$rank_5{$key} = $rank_5{$key} + $highest_score;
		$rank_5_TM{$key} = $rank_5_TM{$key} + $highest_TM;
		$rank_5_MaxSub{$key} = $rank_5_MaxSub{$key} + $highest_SUB;
		print $rank_5{$key} . "\n";
		print $rank_5_TM{$key} . "\n";
		print $rank_5_MaxSub{$key} . "\n";
	}
	close(GDTS);
	#last;  #used in debugging, only run on one target;
	#######following is for the target specific evaluataion resulsts##################
	foreach my $value (sort {$rank_1_TM_target{$b} <=> $rank_1_TM_target{$a}} keys %rank_1_TM_target) {
		#my $group = printf '%-10s', $value;
		#print $value."\n";
		#sleep 1;
		printf TRESULTSG '%-30s', $value;
		print TRESULTSG "$to_print_g_hash{$value}\n";
	}
}
close(TLIST);

open(RANK1, ">$output_folder/rank_1");
foreach my $value (sort {$rank_1{$b} <=> $rank_1{$a}} keys %rank_1){
        print RANK1 "$value $rank_1{$value} $rank_1_no{$value}\n";
}
close(RANK1);

open(RANK5, ">$output_folder/rank_5");
foreach my $value (sort {$rank_5{$b} <=> $rank_5{$a}} keys %rank_5){
        print RANK5 "$value $rank_5{$value} $rank_1_no{$value}\n";
}
close(RANK5);

open(RANK1_TM, ">$output_folder/rank_1_TM");
foreach my $value (sort {$rank_1_TM{$b} <=> $rank_1_TM{$a}} keys %rank_1_TM){
        print RANK1_TM "$value $rank_1_TM{$value} $rank_1_no{$value}\n";
}
close(RANK1_TM);

open(RANK5_TM, ">$output_folder/rank_5_TM");
foreach my $value (sort {$rank_5_TM{$b} <=> $rank_5_TM{$a}} keys %rank_5_TM){
        print RANK5_TM "$value $rank_5_TM{$value} $rank_1_no{$value}\n";
}
close(RANK5_TM);

open(RANK5_MaxSub, ">$output_folder/rank_5_MaxSub");
foreach my $value (sort {$rank_5_MaxSub{$b} <=> $rank_5_MaxSub{$a}} keys %rank_5_MaxSub){
        print RANK5_MaxSub "$value $rank_5_MaxSub{$value} $rank_1_no{$value}\n";
	#print "Content of MaxSub rank5 is: $value $rank_5_MaxSub{$value} $rank_1_no{$value}\n";
}
close(RANK5_MaxSub);

open(RANK1_MaxSub, ">$output_folder/rank_1_MaxSub");
foreach my $value (sort {$rank_1_MaxSub{$b} <=> $rank_1_MaxSub{$a}} keys %rank_1_MaxSub){
        print RANK1_MaxSub "$value $rank_1_MaxSub{$value} $rank_1_no{$value}\n";
}
close(RANK1_MaxSub);

close(LOG);
