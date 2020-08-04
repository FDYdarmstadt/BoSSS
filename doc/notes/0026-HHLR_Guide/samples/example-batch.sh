#!/bin/sh
# Job name
#SBATCH -J $job-name
#
# File / path where STDOUT will be written, the %J is the job id
#SBATCH -o /home/$TuID/$executionDir/$job-name.out%J
#
# File / path where STDERR will be written, the %J is the job id
#SBATCH -e /home/$TuID/$executionDir/$job-name.err%J
#
# Request the time you need for execution in minutes
# The format for the parameter is: [hour:]minute,
# that means for 80 minutes you could also use this: 1:20
#SBATCH -t 23:59:59
#
# Request virtual memory you need for your job in MB per CPU
# 1750 MB is the maximum per core to run CPU efficient. If you
# use more, then not all CPU on a node can run!
#SBATCH --mem-per-cpu=1750
#
# Request the number of compute slots you want to use
#SBATCH -n 16
#
# Specify your mail address - for activation replace < your email> and remove prepending "# "
#SBATCH --mail-user=$yourName@fdy.tu-darmstadt.de
#
# Send a mail when job is done - for activation remove prepending "# "
#SBATCH --mail-type=END
#
#SBATCH -C avx
#SBATCH --exclusive
#
# Unloading a predefined module
# module unload openmpi
# Loading the required module
  module load gcc
  module load openmpi/gcc/1.6.5
  module load acml
# Give an overview about all load modules
 module list

mpiexec mono /home/$TuID/$executionDir/$executable --control /home/$TuID/$executionDir/control-file.cs