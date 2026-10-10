#!/usr/bin/bash
#SBATCH --output=EV_Spain-Greece.sfs.log
#SBATCH --mem=20G
. /local/env/envconda.sh
conda activate ~/my_envs/env_dadi
python dadi_scenarios_with_fixed_inbreeding.py EV_Spain-Greece.sfs
conda deactivate
exit 0
