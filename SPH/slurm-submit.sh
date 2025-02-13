#!/bin/bash -l

#SBATCH -J test
#SBATCH -p icelake
#SBATCH --nodes 1
#SBATCH --cpus-per-task=70
#SBATCH --tasks-per-node=1
#SBATCH --mail-user=pt426@cam.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --time=0-12:00:00

../SWIFT/swift -s -G -t 70 impact-parameters.yml 2>&1 | tee output.log