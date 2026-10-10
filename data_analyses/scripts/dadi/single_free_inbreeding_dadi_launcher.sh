#!/usr/bin/bash

for FILE in `ls *sfs`
    do
        echo "#!/usr/bin/bash" > ${FILE}_single_free_inbreeding_dadi.sh
        echo "#SBATCH --output=${FILE}_single_free_inbreeding_dadi.log" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo "#SBATCH --mem=20G" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo ". /local/env/envconda.sh" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo "conda activate ~/my_envs/env_dadi" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo "python dadi_scenarios_with_single_free_inbreeding.py $FILE" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo "conda deactivate" >> ${FILE}_single_free_inbreeding_dadi.sh
        echo "exit 0" >> ${FILE}_single_free_inbreeding_dadi.sh
        chmod u+x ${FILE}_single_free_inbreeding_dadi.sh
        sbatch -p ecobio ${FILE}_single_free_inbreeding_dadi.sh
    done

exit 0


