#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 8

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh

## Edit the line below as needed:
module load anaconda3

## substitute the command to run your code
## in the two lines below:

/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/salmon/salmon-latest_linux_x86_64/bin/salmon quant -i /u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/salmon/salmon-latest_linux_x86_64/Mus_musculus_c57bl6nj_index -l A -1 $1 -2 $2 -p 8 --gcBias --validateMappings -o ${1%.*}_quant


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####

