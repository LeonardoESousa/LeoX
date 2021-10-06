#The script must include the specific command for job submission used in you batch system followed by $1.
#In this case here, batch.sh is another script that contains the details of submission and takes as 
#argument the input file

sbatch batch.sh $1