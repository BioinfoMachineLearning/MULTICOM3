# MULTICOM3 
MULTICOM3 is an add-on package to improve AlphaFold2- and AlphaFold-Multimer-based prediction of protein tertiary and quaternary structures by diverse multiple sequence alignment sampling, template identification, structural prediction evaluation and structural prediction refinement. It can improve AlphaFold2-based tertiary structure prediction by 8-10% and AlphaFold-Multimer-based quaternary structure prediction by 5-8%. In CASP15, MULTICOM3 used AlphaFold v2.2 as the engine to generate structural predictions. In this release, it is adjusted to run on top of AlphaFold v2.3.2 (https://github.com/deepmind/alphafold/releases/tag/v2.3.2) to leverage the latest improvement on AlphaFold2. You can install MULTICOM3 on top of your AlphaFold2 and AlphaFold-Multimer to improve both the tertiary structure prediction of monomers and the quaternary structure prediction of multimers. 

## Overall workflow for the MULTICOM Protein tertiary structure prediction system
![CASP15_TS pipeline](imgs/CASP15_TS_pipeline.png)

## Overall workflow for the MULTICOM Protein complex structure prediction system 
![CASP15_QS pipeline](imgs/CASP15_QS_pipeline.png)

# Download MULTICOM3 package

```
git clone --recursive https://github.com/BioinfoMachineLearning/MULTICOM3 
```

# Installation (Docker version, modified from [alphafold v2.3.2](https://github.com/google-deepmind/alphafold/blob/v2.3.2))

1.  Install [Docker](https://www.docker.com/).
    *   Install
        [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
        for GPU support.
    *   Setup running
        [Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

2.  Check that AlphaFold will be able to use a GPU by running:

    ```bash
    docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
    ```

    The output of this command should show a list of your GPUs. If it doesn't,
    check if you followed all steps correctly when setting up the
    [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
    or take a look at the following
    [NVIDIA Docker issue](https://github.com/NVIDIA/nvidia-docker/issues/1447#issuecomment-801479573).

    If you wish to run AlphaFold using Singularity (a common containerization
    platform on HPC systems) we recommend using some of the third party Singularity
    setups as linked in https://github.com/deepmind/alphafold/issues/10 or
    https://github.com/deepmind/alphafold/issues/24.

3.  Build the Docker image:

    ```bash
    docker build -f docker/Dockerfile -t multicom3 .
    ```

    If you encounter the following error:

    ```
    W: GPG error: https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
    E: The repository 'https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease' is not signed.
    ```

    use the workaround described in
    https://github.com/deepmind/alphafold/issues/463#issuecomment-1124881779.

4.  Install the `run_docker.py` dependencies. Note: You may optionally wish to
    create a
    [Python Virtual Environment](https://docs.python.org/3/tutorial/venv.html)
    to prevent conflicts with your system's Python environment.

    ```bash
    pip3 install -r docker/requirements.txt
    ```

5.  Make sure that the output directory exists (the default is `/tmp/multicom3`)
    and that you have sufficient permissions to write into it.

6. Download Genetic databases in AlphaFold2/AlphaFold-Multimer 

   ```
   bash $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/scripts/download_all_data.sh <YOUR_ALPHAFOLD_DB_DIR>
   ```

   **Note: The download directory `<YOUR_ALPHAFOLD_DB_DIR>` should *not* be a
subdirectory in the MULTICOM3 repository directory.** If it is, the Docker build
will be slow as the large databases will be copied during the image creation.


7. Download Genetic databases and tools in MULTICOM3

   ```
   python docker/download_database_and_tools.py --dbdir <YOUR_MULTICOM3_DB_DIR>
   ```

   **Note: The download directory `<YOUR_MULTICOM3_DB_DIR>` should *not* be a
subdirectory in the MULTICOM3 repository directory.** If it is, the Docker build
will be slow as the large databases will be copied during the image creation.


# Installation (non Docker version)

## Install AlphaFold/AlphaFold-Multimer and other required third-party packages (modified from [alphafold_non_docker](https://github.com/kalininalab/alphafold_non_docker))

### **Install miniconda**

``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
```

### **Create a new conda environment and update**

``` bash
conda create --name multicom3 python==3.8
conda update -n base conda
```

### **Activate conda environment**

``` bash
conda activate multicom3
```

### **Install dependencies**

- Change `cudatoolkit==11.2.2` version if it is not supported in your system

``` bash
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2.2 pdbfixer
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
```

- Change `jaxlib==0.3.25+cuda11.cudnn805` version if this is not supported in your system

``` bash
pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.3.25 ml-collections==0.1.0 numpy==1.21.6 pandas==1.3.4 protobuf==3.20.1 scipy==1.7.0 tensorflow-cpu==2.9.0

pip install --upgrade --no-cache-dir jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### **Download chemical properties to the common folder**

``` bash
wget -q -P $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

### **Apply OpenMM patch**

``` bash
# $alphafold_path variable is set to the alphafold git repo directory (absolute path)

cd ~/anaconda3/envs/multicom3/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch

# or

cd ~/miniconda3/envs/multicom3/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch
```

### Install other required third-party packages

```
conda install tqdm
conda install -c conda-forge -c bioconda mmseqs2=14.7e284 -y
```

### Download Genetic databases in AlphaFold2/AlphaFold-Multimer

```
bash $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/scripts/download_all_data.sh <YOUR_ALPHAFOLD_DB_DIR>
```

### Install the MULTICOM3 addon system and its databases

```
# Note: here the parameters should be the absolute paths
python setup.py --envdir ~/miniconda3/envs/multicom3 --af_dir $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2 --afdb_dir $YOUR_ALPHAFOLD_DB_DIR

# e.g, 
# python setup.py \
# --envdir ~/miniconda3/envs/multicom3 \
# --af_dir /home/multicom3/tools/alphafold-v2.3.2 \
# --afdb_dir /home/multicom3/tools/alphafold_databases/
```
The setup.py python script will 
* Download the additional databases
* Download the required tools in the system
* Copy the alphafold_addon scripts
* Create the configuration file (bin/db_option) for running the system

# Genetic databases used by MULTICOM3

Assume the following databases have been installed as a part of the AlphaFold2/AlphaFold-Multimer installation
*   [BFD](https://bfd.mmseqs.com/),
*   [MGnify](https://www.ebi.ac.uk/metagenomics/),
*   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
*   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
*   [PDB seqres](https://www.rcsb.org/)
*   [UniRef30](https://uniclust.mmseqs.com/),
*   [UniProt](https://www.uniprot.org/uniprot/),
*   [UniRef90](https://www.uniprot.org/help/uniref).

Additional databases will be installed for the MULTICOM system by setup.py:
*   [AlphaFoldDB](https://alphafold.ebi.ac.uk/): ~53G
*   [ColabFold database](https://colabfold.mmseqs.com/): ~1.7T
*   [Integrated Microbial Genomes (IMG)](https://img.jgi.doe.gov/): ~1.5T
*   [Metaclust](https://metaclust.mmseqs.org/current_release/): ~114G
*   [STRING](https://string-db.org/cgi/download?sessionId=bgV6D67b9gi2): ~129G
*   [pdb_complex](https://www.biorxiv.org/content/10.1101/2023.05.16.541055v1): ~38G
*   [pdb_sort90](https://www.biorxiv.org/content/10.1101/2023.05.01.538929v1): ~48G
*   [Uniclust30](https://uniclust.mmseqs.com/): ~87G

# Important parameter values in the db_option file

- For Docker version, please change the contents in [docker/db_option](docker/db_option)
- For non Docker version, please change the contetns in [bin/db_option](bin/db_option) generated by setup.py 

```
# AlphaFold2 parameters
monomer_num_ensemble = 1
monomer_num_recycle = 3
num_monomer_predictions_per_model = 1

# AlphaFold-Multimer parameters
multimer_num_ensemble = 1
multimer_num_recycle = 3
num_multimer_predictions_per_model = 5
```
Please refer to [AlphaFold2](https://github.com/deepmind/alphafold) to understand the meaning of the parameters. The parameter values stored in bin/db_option file are applied to all the AlphaFold2/AlphaFold-Multimer variants in the MULTICOM3 system to generate predictions. The default bin/db_option file is created automatically by setup.py during the installation. The default parameter values above can be changed if needed. 

# Before running the system for non Docker version

## Activate your python environment and add the MULTICOM3 system path to PYTHONPATH

```bash
conda activate multicom3

# Note: here the parameters should be the absolute paths
export PYTHONPATH=$MULTICOM3_INSTALL_DIR

# e.g, 
# conda activate MULTICOM3
# export PYTHONPATH=/home/multicom3/MULTICOM3

```
Now MULTICOM3 is ready for you to make predictions.

# Running the monomer/tertiary structure prediction pipeline

Say we have a monomer with the sequence `<SEQUENCE>`. The input sequence file should be in the FASTA format as follows:

```fasta
>sequence_name
<SEQUENCE>
```

Note: It is recommended that the name of the sequence file in FASTA format should be the same as the sequence name.

Then run the following command:

```bash
# docker version
python3 docker/run_docker.py \
    --mode=monomer \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --af_db_dir=$YOUR_AF_DB_DIR \ --multicom3_db_dir=$YOUR_MULTICOM3_DB_DIR \
    --output_dir=$OUTDIR

# non docker version
python bin/monomer.py \
    --option_file=bin/db_option \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --output_dir=$OUTDIR
```

option_file is a file in the MULTICOM package to store some key parameter values for AlphaFold2 and AlphaFold-Multimer. fasta_path is the full path of the file storing the input protein sequence(s) in the FASTA format. output_dir specifies where the prediction results are stored. Please be aware that we have included a parameter (--run_img) that allows you to turn off the usage of the IMG database for faster prediction (--run_img=False). In the case of --run_img=True, the program will pause at the monomer prediction generation stage to wait for the IMG alignment to be created. Generating alignments from IMG may take a much longer time, potentially several days, because the database is very large. So run_img is set to false by default. It is advised that run_img is set to true only if other alignments cannot yield good results.

## Output

```
$OUTPUT_DIR/                                   # Your output directory
    N1_monomer_alignments_generation/          # Working directory for generating monomer MSAs
    N1_monomer_alignments_generation_img/      # Working directory for generating IMG MSA
        # Note: the img.running file may use many disk space
    N2_monomer_template_search/                # Working directory for searching monomer templates
    N3_monomer_structure_generation/           # Working directory for generating monomer structural predictions
    N4_monomer_structure_evaluation/           # Working directory for evaluating the monomer structural predictions
        - alphafold_ranking.csv    # AlphaFold2 pLDDT ranking
        - pairwise_ranking.tm      # Pairwise (APOLLO) ranking
        - pairwise_af_avg.ranking  # Average ranking of the two
    N5_monomer_structure_refinement_avg/       # Working directory for monomer structure refinement
    N5_monomer_structure_refinement_avg_final/ # Output directory for the refined monomer predictions
        - final_ranking.csv        # AlphaFold2 pLDDT ranking of the original and refined predictions
```

* The predictions and ranking files are saved in the *N4_monomer_structure_evaluation* folder. You can check the AlphaFold2 pLDDT score ranking file (alphafold_ranking.csv) to look for the structure with the highest pLDDT score. The *pairwise_ranking.tm* and *pairwise_af_avg.ranking* are the other two ranking files. 

* The refined monomer predictions are saved in *N5_monomer_structure_refinement_avg_final*.

# Running the multimer/quaternary structure prediction pipeline

## Folding a homo-multimer

Say we have a homomer with 4 copies of the same sequence
`<SEQUENCE>`. The input file should be in the format as follows:

```fasta
>sequence_1
<SEQUENCE>
>sequence_2
<SEQUENCE>
>sequence_3
<SEQUENCE>
>sequence_4
<SEQUENCE>
```

Then run the following command:

```bash
# docker version
python3 docker/run_docker.py \
    --mode=homomer \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --af_db_dir=$YOUR_AF_DB_DIR \ --multicom3_db_dir=$YOUR_MULTICOM3_DB_DIR \
    --output_dir=$OUTDIR

# non docker version
python bin/homomer.py \
    --option_file=bin/db_option \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --output_dir=$OUTDIR
```

## Folding a hetero-multimer

Say we have an A2B3 heteromer, i.e. with 2 copies of
`<SEQUENCE A>` and 3 copies of `<SEQUENCE B>`. The input file should be in the format as follows (the same sequences should be grouped together):

```fasta
>sequence_1
<SEQUENCE A>
>sequence_2
<SEQUENCE A>
>sequence_3
<SEQUENCE B>
>sequence_4
<SEQUENCE B>
>sequence_5
<SEQUENCE B>
```

Then run the following command:

```bash
# docker version
python3 docker/run_docker.py \
    --mode=heteromer \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --af_db_dir=$YOUR_AF_DB_DIR \ --multicom3_db_dir=$YOUR_MULTICOM3_DB_DIR \
    --output_dir=$OUTDIR

# non docker version
python bin/heteromer.py \
    --option_file=bin/db_option \
    --fasta_path=$YOUR_FASTA \
    --run_img=False \
    --output_dir=$OUTDIR
```

## Output

```
$OUTPUT_DIR/                                   # Your output directory
    N1_monomer_alignments_generation/          # Working directory for generating monomer MSAs
        - Subunit A
        - Subunit B
        - ...
    N1_monomer_alignments_generation_img/      # Working directory for generating IMG MSA
        - Subunit A
        - Subunit B
        - ...
    N2_monomer_template_search/                # Working directory for searching monomer templates
        - Subunit A
        - Subunit B
        - ...
    N3_monomer_structure_generation/           # Working directory for generating monomer structural predictions
        - Subunit A
        - Subunit B
        - ...
    N4_monomer_alignments_concatenation/       # Working directory for concatenating the monomer MSAs
    N5_monomer_templates_search/               # Working directory for concatenating the monomer templates
    N6_multimer_structure_generation/          # Working directory for generating multimer structural predictions
    N7_monomer_structure_evaluation            # Working directory for evaluating monomer structural predictions
        - Subunit A
            # Rankings for all the predictions
            - alphafold_ranking.csv            # AlphaFold2 pLDDT ranking 
            - pairwise_ranking.tm              # Pairwise (APOLLO) ranking
            - pairwise_af_avg.ranking          # Average ranking of the two 

            # Rankings for the predictions generated by monomer structure prediction
            - alphafold_ranking_monomer.csv    # AlphaFold2 pLDDT ranking 
            - pairwise_af_avg_monomer.ranking  # Average ranking 

            # Rankings for the predictions extracted from multimer predictions
            - alphafold_ranking_multimer.csv   # AlphaFold2 pLDDT ranking 
            - pairwise_af_avg_multimer.ranking # Average ranking 

        - Subunit B
        - ...
    N8_multimer_structure_evaluation           # Working directory for evaluating multimer structural predictions
        - alphafold_ranking.csv                # AlphaFold2 pLDDT ranking
        - multieva.csv                         # Pairwise ranking using MMalign
        - pairwise_af_avg.ranking              # Average ranking of the two
    N9_multimer_structure_refinement           # Working directory for refining multimer structural predictions
    N9_multimer_structure_refinement_final     # Output directory for the refined multimer predictions
```

* The predictions and ranking files are saved in *N8_multimer_structure_evaluation*, similarly, you can check the AlphaFold-Multimer confidence score ranking file (alphafold_ranking.csv) to look for the structure with the highest predicted confidence score generated by AlphaFold-Multimer. The *multieva.csv* and *pairwise_af_avg.ranking* are the other two ranking files.

* The refined multimer predictions are saved in *N9_multimer_structure_refinement_final*.

* The monomer structures and ranking files are saved in *N7_monomer_structure_evaluation* if you want to check the predictions and rankings for the monomer structures.

# Some CASP15 Prediction Examples (MULTICOM versus the Standard AlphaFold method: NIBS-AF2-multimer)

![CASP15 pipeline](imgs/CASP15_good_examples1.png)
![CASP15 pipeline](imgs/CASP15_good_examples2.png)


# Citing this work

**If you use this package for tertiary or quaternary structure prediction, please cite:**

**Tertiary (monomer) structure prediction**

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583-589.

Liu, J., Guo, Z., Wu, T., Roy, R. S., Chen, C., & Cheng, J. (2023). Improving AlphaFold2-based Protein Tertiary Structure Prediction with MULTICOM in CASP15. bioRxiv, 2023-05-01. (https://www.biorxiv.org/content/10.1101/2023.05.01.538929v1) 

**Quaternary (multimer) structure prediction**

Evans, R., O’Neill, M., Pritzel, A., Antropova, N., Senior, A., Green, T., ... & Hassabis, D. (2021). Protein complex prediction with AlphaFold-Multimer. BioRxiv, 2021-10.
 
Liu, J., Guo, Z., Wu, T., Roy, R. S., Quadir, F., Chen, C., & Cheng, J. (2023). Enhancing AlphaFold-Multimer-based Protein Complex Structure Prediction with MULTICOM  in CASP15. bioRxiv, 2023-05-16. (https://www.biorxiv.org/content/10.1101/2023.05.16.541055v1)



