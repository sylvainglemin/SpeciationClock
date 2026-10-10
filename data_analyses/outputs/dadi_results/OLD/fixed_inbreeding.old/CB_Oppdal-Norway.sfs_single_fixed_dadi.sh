#!/usr/bin/bash
#SBATCH --output=CB_Oppdal-Norway.sfs_single_fixed.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_single_fixed_inbreeding.py CB_Oppdal-Norway.sfs
conda deactivate
exit 0
