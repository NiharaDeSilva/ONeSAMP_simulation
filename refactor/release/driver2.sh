#!/bin/bash

###########################
###########################
## Simulation Parameters ##
###########################
###########################

# Number of trials to run simcoal2 (and simuPOP)
export NUMTRIALS=10

# List of different Ne values to be used in simulation. Must have 5 digits apiece.
export NeVals="00050 00100 00150 00200"
# 00128 00256 00512 01024
#export NeVals="00100"
export NeDigits=5

# Do not modify; used for cluster progress calculation
export NeNumberOfValues=`echo $NeVals | wc -w`


###########
# ONeSAMP #
###########

# Number of individuals simuPOP should input in to Ne estimators. (500)
#export outputSampleSize=100
export outputSampleSize=50

# Number of loci in simulated trials. Only needed if simulating own trials.
export loci=50

# Minimum allele frequency necessary in populations
ONESAMP2COAL_MINALLELEFREQUENCY=0.05

# Mutation rate specified for simuPOP and ONeSAMP
export mutationRate="0.000000012"

# We require a flat distribution of values to supply to the executable.
export rangeNe=100,500

# We require a specified range of theta values to supply to the executable.
export theta=0.000048,0.0048

# Number of populations that ONeSAMP generates
export numOneSampTrials=20000

# Duration of bottlenecks
export duration=2,8

# Unix niceness of processes used by ONeSAMP: 19 is low priority, 0 is normal priority
export processPriority=19

# Unix niceness of file I/O: 7 is low priority, 0 is normal
export filePriority=7

# Do not modify; used for checking length of input files.
export trialsplus1=$(($numOneSampTrials + 1))

# SNPs or microsats: s for SNPs, m for microsatellites
export microsatsOrSNPs=s

#########################
# Export some constants #
#########################
export suffix=".reduced"
export pop="pop"
export number="number"
export gen=".gen"
export PARAMETER=".par"
export ARPSUFFIX="_0.arp"
export DISTSUFFIX=".out"
export LDNeSUFFIX=".genLD.txt"
export myLogin=`whoami`
export hostx=`hostname`
#export numLoggedIn=`who | cut -f1 -d\  | uniq | wc -l`

# Do not modify
export popStrLen=${#pop}
export suffixLen=${#suffix}
export DISTSUFFIXLen=${#DISTSUFFIX}
export genLen=${#gen}

############################################################################
# EDITING STEP 1: enter the filename of this file
export DRIVERNAME=driver2.sh
############################################################################

# enter a python path for the simulator that contains the
# appropriate simuPOP library
#export PYTHONPATH=`echo ~`/simupop/lib64/python2.7/site-packages/simuPOP-1.1.6-py2.7-linux-x86_64.egg

# Do not modify
export CURRENTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

############################################################################
# EDITING STEP 2: if you are using a cluster setup, please supply local
# paths to be used for temporary file storage /home/sharique/Downloads/ONeSAMP2/refactor/release/tmp
#export in_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_in_filebuffer
#export out_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_out_filebuffer
#export singleton_out_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_filebuffer_safe
touch temp_in
touch temp_out
touch temp_safe
export in_bufferdir=temp_in
export out_bufferdir=temp_out
export singleton_out_bufferdir=temp_safe
############################################################################

# Do not use
#export SIMCOAL2EXECUTABLE=`echo ~`/simcoal2/simcoal2_1_2

# Do not modify 
export hostFile=$3

##################################
# Define a few programming tools #
##################################

############################################################################
# EDITING STEP 3: please supply the locations of the python, perl, and R
# interpreters, as well as the C and C++ compilers to be used
export PYTHONEXECUTABLE=python
# export PERLEXECUTABLE='perl -I/usr/lib64/perl5/vendor_perl -I/usr/share/perl5/vendor_perl'
# export RINTERPRETER=Rscript
# export CCCOMPILE=g++
# export CCOMPILE=gcc
############################################################################

############################################################################
# EDITING STEP 4: please supply the locations of the built ONeSAMP2 executable
# and the coalescent simulator, as well as the ONESAMP2 executable name
export ONESAMP2EXEC=refactor_main
export ONESAMP2=../release/${ONESAMP2EXEC}
export ONESAMP2COAL=./refactor_coalescent_simulator
#export RSCRIPT=rScript.r
############################################################################

#######################
# Usage error message #
#######################

if [ $# -eq 0 ]; then
  echo "Usage: Please use one of these sets of options"
  echo "./$DRIVERNAME generateSimulatedPopulationsOfKnownNe"
  echo "./$DRIVERNAME analyzeUnknownPopulation"
  echo "./$DRIVERNAME runOneSamp start"
  echo "./$DRIVERNAME runOneSamp progress"
  echo "./$DRIVERNAME estimateNe"
  echo "./$DRIVERNAME runLDNe"
  echo "./$DRIVERNAME onesampDataAnalytics outFile"
fi

export log_file="execution_stats_original"


############
############
## STEP 1 ##
############
############
if [ $1 ]; then
  if [ $1 == "generateSimulatedPopulationsOfKnownNe" ]; then
    log_file="simulation_log.txt"
    echo "Ne,Iteration,Time,Memory" > $log_file
    
    loci_values=(160) # Example loci values
    size_values=(100) # Example size values
    
    for numPOP in $(echo $NeVals); do
      for loci in "${loci_values[@]}"; do
        for outputSampleSize in "${size_values[@]}"; do
          for iteration in {1..10}; do
            #echo "Generating test cases for Ne = $numPOP, Loci = $loci, Size = $outputSampleSize"
            
            sleep 2

            #echo "Running ONESAMP2 coalescent simulator..."
            
            export genepopdata="ONeSAMP1Temp${numPOP}number"
            export simcoalout="ONeSAMP2Temp${numPOP}number"
            export simupopin="ONeSAMP3Temp${numPOP}number${gen}"
            export executablename="ONeSAMP4Temp${numPOP}number"
            export pythonscript="ONeSAMP5Temp${numPOP}number.py"
            export OUTPUT="${pop}${iteration}${numPOP}${number}${gen}"
  
            export expectedNe=$numPOP
  
            rm -fr $simcoalout
            rm -fr $simcoalout$PARAMETER
            rm -fr $genepopdata
            rm -fr $executablename
            rm -fr $pythonscript
            rm -fr $OUTPUT

            # Start timing in milliseconds
            start_time=$(date +%s%3N)
            
            pid=$(./ONESAMP2 -t1 -rC -b$numPOP -d1 -u$mutationRate -v${theta} -$microsatsOrSNPs -l$loci -i$outputSampleSize -o1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -p > $OUTPUT$suffix & echo $!)

            total_memory=0

            # Continuously check memory usage while the process is running
            while kill -0 $pid 2>/dev/null; do
              current_memory=$(ps -o rss= -p $pid)
              total_memory=$((total_memory + current_memory))
              sleep 0.1
            done

            # End timing in milliseconds
            end_time=$(date +%s%3N)
            execution_time=$((end_time - start_time))
            execution_time_sec=$(echo "scale=3; $execution_time_ms / 1000" | bc)
            
            # Convert total memory from KB to MB for better readability
            total_memory_mb=$((total_memory / 1024))

            # Write to log file
            echo "$numPOP,$iteration,$execution_time_sec,$total_memory_mb" >> $log_file

            echo "Trial complete"
          done # End of size loop
        done # End of loci loop
      done # End of numPOP loop
    done # End of iteration loop
  fi # End of inner if
fi # End of outer if
