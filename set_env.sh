conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04
