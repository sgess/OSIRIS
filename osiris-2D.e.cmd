#!/bin/csh -f
#  osiris-2D.e.cmd
#
#  UGE job for osiris-2D.e built Wed Nov 20 02:37:59 PST 2013
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog
#$ -o /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/scratch/s/sgess/e_r3/osiris-2D.e
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Parallelism:  64-way parallel
#  Resources requested
#$ -pe dc* 64
#$ -l h_data=1024M,h_rt=2:00:00,highp
# # #
#
#  Name of application for log
#$ -v QQAPP=openmpi
#  Email address to notify
#$ -M sgess@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#  Uncomment the next line to have your environment variables used by SGE
#$ -V
#
# Initialization for mpi parallel execution
#
  unalias *
  set qqversion = 
  set qqapp     = "openmpi parallel"
  set qqptasks  = 64
  set qqidir    = /u/scratch/s/sgess/e_r3
  set qqjob     = osiris-2D.e
  set qqodir    = /u/scratch/s/sgess/e_r3
  cd     /u/scratch/s/sgess/e_r3
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for osiris-2D.e built Wed Nov 20 02:37:59 PST 2013"
  echo ""
  echo "  osiris-2D.e directory:"
  echo "    "/u/scratch/s/sgess/e_r3
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    $qqscratch"
  echo "  osiris-2D.e 64-way parallel job configuration:"
  echo "    $qqconfig" | tr "\\" "\n"
#
  echo ""
  echo "osiris-2D.e started on:   "` hostname -s `
  echo "osiris-2D.e started at:   "` date `
  echo ""
#
# Run the user program
#

  source /u/local/Modules/default/init/modules.csh
  module load intel/11.1
  module load openmpi/1.4

  echo osiris-2D.e -nt 64 "" \>\& osiris-2D.e.output.$JOB_ID
  echo ""

  setenv PATH /u/local/bin:$PATH

  time $MPI_BIN/mpiexec -n 64 -machinefile $QQ_NODEFILE  \
         /u/scratch/s/sgess/e_r3/osiris-2D.e  >& /u/scratch/s/sgess/e_r3/osiris-2D.e.output.$JOB_ID 

  echo ""
  echo "osiris-2D.e finished at:  "` date `

#
# Cleanup after mpi parallel execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/openmpi.log.parallel
 if (`wc -l /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
        echo " "  >> /u/local/apps/queue.logs/openmpi.log.parallel
        tail -10 /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  else
        cat /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  endif
#  cat            /u/scratch/s/sgess/e_r3/osiris-2D.e.joblog.$JOB_ID           >> /u/local/apps/queue.logs/openmpi.log.parallel
  exit (0)
