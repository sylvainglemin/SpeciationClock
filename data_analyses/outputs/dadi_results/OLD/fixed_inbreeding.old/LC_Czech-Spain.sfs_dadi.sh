#!/usr/bin/bash
#SBATCH --output=LC_Czech-Spain.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py LC_Czech-Spain.sfs
conda deactivate
exit 0
