#!/usr/bin/bash

for FILE in `ls *sfs`
    do
        echo "#!/usr/bin/bash" > ${FILE}_fixed_dadi.sh
        echo "#SBATCH --output=${FILE}_fixed.log" >> ${FILE}_fixed_dadi.sh
        echo "#SBATCH --mem=20G" >> ${FILE}_fixed_dadi.sh
        echo ". /local/env/envconda.sh" >> ${FILE}_fixed_dadi.sh
        echo "conda activate ~/my_envs/env_dadi" >> ${FILE}_fixed_dadi.sh
        echo "python dadi_scenarios_with_fixed_inbreeding.py $FILE" >> ${FILE}_fixed_dadi.sh
        echo "conda deactivate" >> ${FILE}_fixed_dadi.sh
        echo "exit 0" >> ${FILE}_fixed_dadi.sh
        chmod u+x ${FILE}_fixed_dadi.sh
        sbatch -p ecobio ${FILE}_fixed_dadi.sh
    done

exit 0


