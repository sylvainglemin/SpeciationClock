#!/usr/bin/bash
#SBATCH --output=MT_Greece-Spain.sfs_single_fixed.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_single_fixed_inbreeding.py MT_Greece-Spain.sfs
conda deactivate
exit 0
