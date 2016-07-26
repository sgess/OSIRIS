#!/bin/bash

me=$USER
fl=`expr substr $USER 1 1`
d=$PWD

echo Hi $me!

day=`date +%d`
month=`date +%b`
year=`date +%Y`

echo Today is $day $month $year

echo Creating Simulation directory and copying files now...

cd /u/home/$fl/$me/d2scratch
mkdir $1
cd $1

cp $d/osinputs/$year/$month/$day/os-stdin_$1 os-stdin
cp /u/home/s/sgess/project/executables/OSIRIS/osiris-64-openmpi-production-develop-2.2.454M-2D.e osiris-2D.e

echo Writing .cmd file...

# file variables...
DAY=`date +%a`
month=`date +%b`
day=`date +%e`
TIME=`date +%T`
zone=`date +%Z`
year=`date +%Y`

me=$USER
fl=`expr substr $USER 1 1`
dirname=$1

nodes=$2
mem=$3
rt=$4

# write file using above variables
# NO TOUCHY!
cat > osiris-2D.e.cmd << EOF
#!/bin/csh -f
#  osiris-2D.e.cmd
#
#  UGE job for osiris-2D.e built $DAY $month $day $TIME $zone $year 
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog
#$ -o /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Parallelism:  $nodes-way parallel
#  Resources requested
#$ -pe dc* $nodes
#$ -l h_data=${mem}M,h_rt=$rt:00:00,highp
# # #
#
#  Name of application for log
#$ -v QQAPP=openmpi
#  Email address to notify
#$ -M $me@mail
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
  set qqptasks  = $nodes
  set qqidir    = /u/home/$fl/$me/d2scratch/$dirname
  set qqjob     = osiris-2D.e
  set qqodir    = /u/home/$fl/$me/d2scratch/$dirname
  cd     /u/home/$fl/$me/d2scratch/$dirname
  source /u/local/bin/qq.sge/qr.runtime
  if (\$status != 0) exit (1)
#
  echo "UGE job for osiris-2D.e built $DAY $month $day $TIME $zone $year"
  echo ""
  echo "  osiris-2D.e directory:"
  echo "    "/u/home/$fl/$me/d2scratch/$dirname
  echo "  Submitted to UGE:"
  echo "    "\$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    \$qqscratch"
  echo "  osiris-2D.e $nodes-way parallel job configuration:"
  echo "    \$qqconfig" | tr "\\\\" "\n"
#
  echo ""
  echo "osiris-2D.e started on:   "\` hostname -s \`
  echo "osiris-2D.e started at:   "\` date \`
  echo ""
#
# Run the user program
#

  source /u/local/Modules/default/init/modules.csh
  module load intel/11.1
  module load openmpi/1.4
  module load hdf5

  echo osiris-2D.e -nt $nodes "" \>\& osiris-2D.e.output.\$JOB_ID
  echo ""

  setenv PATH /u/local/bin:\$PATH

  time \$MPI_BIN/mpiexec -n $nodes -machinefile \$QQ_NODEFILE  \\
         /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e  >& /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.output.\$JOB_ID 

  echo ""
  echo "osiris-2D.e finished at:  "\` date \`

#
# Cleanup after mpi parallel execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID --------" >> /u/local/apps/queue.logs/openmpi.log.parallel
 if (\`wc -l /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID  | awk '{print \$1}'\` >= 1000) then
        head -50 /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
        echo " "  >> /u/local/apps/queue.logs/openmpi.log.parallel
        tail -10 /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  else
        cat /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  endif
#  cat            /u/home/$fl/$me/d2scratch/$dirname/osiris-2D.e.joblog.\$JOB_ID           >> /u/local/apps/queue.logs/openmpi.log.parallel
  exit (0)
EOF

echo Submitting Job...

#qsub osiris-2D.e.cmd

echo Job Submitted!
