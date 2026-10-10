#!/usr/bin/bash
#SBATCH --output=SAC_Svalbard_Lovisa-Alaska.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py SAC_Svalbard_Lovisa-Alaska.sfs
conda deactivate
exit 0
