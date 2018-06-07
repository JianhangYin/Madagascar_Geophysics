#!/bin/bash
#
#
# The SBATCH directive below says the name of the partition to be used.
# In most cases, you should use the queue named normal.
#
#SBATCH --partition=seiswave
#
# The SBATCH directive below says to use 1 CPU core of one CPU chip
# on one compute node, meaning that this batch jobs is non-parallel
# (serial). DON'T CHANGE THIS.
#
#SBATCH --ntasks=1
#
# The SBATCH directive below indicates to the scheduler the amount of
# memory your job will use in megabytes.  This is critical information
# for scheduling of nonexclusive jobs, since it prevents the scheduler
# from assigning more jobs to a given compute node than that node
# has memory for.
#
# Not supplying an accurate number for the required memory usage can
# have a dramatic negative impact on the speed of your code. Not
# supplying any value for the required memory usage will instruct
# the scheduler to use a modest value, currently 1024MB.
#
# If you need assistance determining the amount of memory your code
# requires, contact support@oscer.ou.edu
#
#SBATCH --mem=30000
#
# STDOUT is the output that normally would go to the terminal screen.
# STDERR is like STDOUT, except that it contains error messages instead
# of regular output.
#
# The SBATCH directive below says, send STDOUT and STDERR to the filenames
# listed below.
#
# Note that, in these filenames, %J will be replaced by the batch job ID
# number
#
# Everywhere throughout this batch script file:
# * CHANGE yourusername to your username;
# * CHANGE youremailaddress@yourinstitution.edu to your e-mail address;
# * CHANGE directory_to_run_in to the directory you want to run in;
# * CHANGE jobname to your preferred name for this batch job.
#
#SBATCH --output=./log/1.txt
#SBATCH --error=./log/2.txt
#
# The SBATCH directive below says to run for up to 12 hours (and zero
# minutes and zero seconds) of "wall clock" time (time experienced in
# real life).
# Currently, the maximum allowable wall clock time per batch job is
# 48:00 hours.
# Acceptable time formats include "minutes", "minutes:seconds",
# "hours:minutes:seconds", "days-hours", "days-hours:minutes" and
# "days-hours:minutes:seconds".
#
#SBATCH --time=20:00:00
#
# The SBATCH directive below says the name of the batch job, as it
# will appear in the batch partition when you do an squeue command.
#
# Everywhere throughout this batch script file:
# * CHANGE yourusername to your username;
# * CHANGE youremailaddress@yourinstitution.edu to your e-mail address;
# * CHANGE directory_to_run_in to the directory you want to run in;
# * CHANGE jobname to your preferred name for this batch job.
#
#SBATCH --job-name=GmRTM_sum_3_1
#
# The SBATCH directive below says the e-mail address to send
# notifications to, which should be changed to your e-mail address.
#
# Everywhere throughout this batch script file:
# * CHANGE yourusername to your username;
# * CHANGE youremailaddress@yourinstitution.edu to your e-mail address;
# * CHANGE directory_to_run_in to the directory you want to run in;
# * CHANGE jobname to your preferred name for this batch job.
#
#SBATCH --mail-user=yinjianhang@ou.edu
#
# The SBATCH directive below says to e-mail a notification when the
# batch job either completes or fails.
#
# If you only want e-mails when when a batch job fails, then delete
# the SBATCH directive below.
#
# Change to the directory that you want to run in.
#
# Everywhere throughout this batch script file:
# * CHANGE yourusername to your username;
# * CHANGE youremailaddress@yourinstitution.edu to your e-mail address;
# * CHANGE directory_to_run_in to the directory you want to run in;
# * CHANGE jobname to your preferred name for this batch job.
#
#SBATCH --workdir=./
#

#################################################

# If you like, you may uncomment the following line
# to display your working directory in your output,
# to make sure that you're in the correct directory.
# pwd

# If your job requires any "module load" commands,
# then please insert them here, without a pound sign (#).
# Typically, if your program required a
# module to compile it, then include the same module
# here to run your program.
# Run the executable, redirecting input from the given file.
# The date commands and the time command help track runtime
# and performance.
#
# NOTE: The line with a "<" sign, immediately before the
# second date command, redirects standard input
# (typically a user typing at the keyboard) from a file of
# human readable text.
# If your program uses no input, or if for all input your program
# explicitly opens its files inside the program (for example, via
# an OPEN statement in Fortran or an fopen command in C), then you
# don't need that line, so DELETE IT.
#
# Everywhere throughout this batch script file:
# * CHANGE yourusername to your username;
# * CHANGE youremailaddress@yourinstitution.edu to your e-mail address;
# * CHANGE directory_to_run_in to the directory you want to run in;
# * CHANGE jobname to your preferred name for this batch job.
#
# Here, also CHANGE executable_directory to the directory where your
# executable resides, CHANGE executable_filename to the name of the
# executable file, CHANGE input_directory to the directory where your
# standard input file resides, and CHANGE input_filename.txt to the
# name of your standard input file (that is, the file that you want
# to use for standard input, instead of a user typing live at the
# keyboard, which isn't possible for a batch job).
scons -f SConstruct_sum_1
