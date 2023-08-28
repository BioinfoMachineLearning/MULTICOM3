#/bin/sh
rm ../pdb_orig/*
perl filter_pdb.pl pdb_ids target_ids ../pdb_orig ../sequence ../alignment ../pdb_filtered


