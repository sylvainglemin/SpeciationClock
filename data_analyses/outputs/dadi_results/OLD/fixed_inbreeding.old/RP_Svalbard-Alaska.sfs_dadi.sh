#!/usr/bin/bash
#SBATCH --output=RP_Svalbard-Alaska.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py RP_Svalbard-Alaska.sfs
conda deactivate
exit 0
