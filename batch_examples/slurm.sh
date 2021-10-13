#The script must include the specific command for job submission used in you batch system followed by $1.
#In this case here, batch.sh is another script that contains the details of submission and takes as 
#argument the input file. Both slurm.sh and batch.sh files must be present before running option 2.

sbatch batch.sh $1