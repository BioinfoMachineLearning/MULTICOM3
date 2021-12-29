## About
This program is designed to search for templates based on structural similarity using TMalign and known dimeric interaction between proteins.
The inputs of this program is monomer atom files, dimer_type, database and output directory. The output of this program is a list of templates that has interaction and their tmscores and also structure based alignment of the templates in "pir" format.  


The following steps are performed to achieve the task:
- The input pdbs are copied to the output directory
- Performs TMalign against all the atom files in the atom database that has dimeric interactions (based on precalculated homodimer or heterodimer list)
- Then monomer top_templates (MAX_TEMPLATE_COUNT) based on TMscore are considered for further filtration
- For filtration templates less than a certain threshold value (THRESHOLD ) are filtered
- Then outputs a list of top pairwise templates "result_contact_<threshold_value>.txt"

## First time setup
Open the template_settings.txt and edit the required information.

1. TMalign 
   * Download TMalign (source: https://zhanggroup.org/TM-align/)
   * Then edit TM_ALIGN_DIR with the directory of the TMalign
     * TM_ALIGN_DIR: directory of the TMalign
     * e.g TM_ALIGN_DIR:/home/rsr3gt/programs/template_search/resources/TMalign
     
2. Edit the number of the shorlisted template:
   * e.g INTERMEDIATE_TEMPLATE_COUNT: 10000
   

3. Edit the minimum quality/TMscore of the templates 
   * e.g  THRESHOLD: 0.3
   
### Running Structure Based Tempalte Search
````
General usage: python template_search.py <ligand_atom_path> <receptor_atom_path> <homodimer_flag> <database directory> <Output_directory>
````    
#### 1. For Homodimer
````
  
e.g python template_search.py /home/rsr3gt/programs/template_search/resources/test_pdb/hetero/15C8H.atom /home/rsr3gt/programs/template_search/resources/test_pdb/hetero/15C8L.atom 1 /home/multicom4s_data/HOMO_STD/pdb/ /home/rsr3gt/programs/template_search/resources/output/test_1/
````

#### 2. For Heterodimer
````
e.g python template_search.py /home/rsr3gt/programs/template_search/resources/test_pdb/hetero/15C8H.atom /home/rsr3gt/programs/template_search/resources/test_pdb/hetero/15C8L.atom 0 /home/multicom4s_data/HETERO_STD/pdb/ /home/rsr3gt/programs/template_search/resources/output/test_1/
````

#NOTE
````
1. HETERO_FILE and HOMO_FILE were generated based on interactions that were found during the generation of the dataset hetero30 and homo30 respectively

2. For a larg dataset the "SLEEP_TIME" inside the template serach_settings should be increased. 
````