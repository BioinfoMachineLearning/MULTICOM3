ARG CUDA=11.1.1
FROM nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu18.04
ARG CUDA

# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        cmake \
        cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
        git \
        # hmmer \
        # kalign \
        tzdata \
        wget \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

# Compile HHsuite from source.
# RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
#     && mkdir /tmp/hh-suite/build \
#     && pushd /tmp/hh-suite/build \
#     && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
#     && make -j 4 && make install \
#     && ln -s /opt/hhsuite/bin/* /usr/bin \
#     && popd \
#     && rm -rf /tmp/hh-suite

# Install Miniconda package manager.
RUN wget -q -P /tmp \
  https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh \
    && bash /tmp/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh

# Install conda packages.
ENV PATH="/opt/conda/bin:$PATH"
RUN conda install -qy conda==4.13.0 

RUN conda install -y -c conda-forge \
      openmm=7.5.1 \
      cudatoolkit==${CUDA_VERSION} \
      pdbfixer

RUN conda install tqdm && conda install -c conda-forge -c bioconda mmseqs2=14.7e284 -y

RUN conda install -c bioconda hmmer hhsuite==3.3.0 kalign2 -y

RUN conda clean --all --force-pkgs-dirs --yes

COPY . /app/MULTICOM3
RUN wget -q -P /app/MULTICOM3/tools/alphafold-v2.3.2/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Install pip packages.
RUN pip3 install --upgrade pip --no-cache-dir \
    && pip3 install -r /app/MULTICOM3/tools/alphafold-v2.3.2/requirements.txt --no-cache-dir \
    && pip3 install --upgrade --no-cache-dir \
      jax==0.3.25 \
      jaxlib==0.3.25+cuda11.cudnn805 \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Apply OpenMM patch.
WORKDIR /opt/conda/lib/python3.8/site-packages
RUN patch -p0 < /app/MULTICOM3/tools/alphafold-v2.3.2/docker/openmm.patch

# Add SETUID bit to the ldconfig binary so that non-root users can run it.
RUN chmod u+s /sbin/ldconfig.real

# We need to run `ldconfig` first to ensure GPUs are visible, due to some quirk
# with Debian. See https://github.com/NVIDIA/nvidia-docker/issues/1399 for
# details.
# ENTRYPOINT does not support easily running multiple commands, so instead we
# write a shell script to wrap them up.
WORKDIR /app/MULTICOM3
RUN chmod +x /app/MULTICOM3/docker/run_multicom3.sh
ENTRYPOINT ["/app/MULTICOM3/docker/run_multicom3.sh"]