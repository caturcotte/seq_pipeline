module purge
wget -O "${HOME}/Miniforge3.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash "${HOME}/Miniforge3.sh" -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
conda activate

mamba create -n seq_pipeline -c conda-forge -c bioconda snakemake
