#!/usr/bin/bash
#SBATCH --output=CB_Oppdal-Svalbard.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py CB_Oppdal-Svalbard.sfs
conda deactivate
exit 0
