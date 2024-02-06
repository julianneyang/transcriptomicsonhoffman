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
#$ -pe shared 2

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load anaconda3
conda activate kneaddata

## substitute the command to run your code
## in the two lines below:



trimmomatic PE $1 $2 trim_${1%.*}_paired.fq.gz trim_${1%.*}_unpaired.fq.gz trim_${2%.*}_paired.fq.gz trim_${2%.*}_unpaired.fq.gz ILLUMINACLIP:/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/Trimmomatic/trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####

