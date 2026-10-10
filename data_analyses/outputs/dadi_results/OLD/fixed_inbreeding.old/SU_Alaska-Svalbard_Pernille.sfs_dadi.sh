#!/usr/bin/bash
#SBATCH --output=SU_Alaska-Svalbard_Pernille.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py SU_Alaska-Svalbard_Pernille.sfs
conda deactivate
exit 0
