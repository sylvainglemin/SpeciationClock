#!/usr/bin/bash

for FILE in `ls *sfs`
    do
        echo "#!/usr/bin/bash" > ${FILE}_dadi.sh
        echo "#SBATCH --output=${FILE}.log" >> ${FILE}_dadi.sh
        echo ". /local/env/envconda.sh" >> ${FILE}_dadi.sh
        echo "conda activate ~/my_envs/env_dadi" >> ${FILE}_dadi.sh
        echo "python dadi_scenarios_with_inbreeding.py $FILE" >> ${FILE}_dadi.sh
        echo "conda deactivate" >> ${FILE}_dadi.sh
        echo "exit 0" >> ${FILE}_dadi.sh
        chmod u+x ${FILE}_dadi.sh
        sbatch -p ecobio ${FILE}_dadi.sh
    done

exit 0


