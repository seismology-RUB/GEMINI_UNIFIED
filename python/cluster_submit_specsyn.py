#!/usr/bin/env python3
#
#  Submit computeSyntheticsForSpecfemCoupling job to cluster machine
#------------------------------------------------------------------
import os
import os.path
import glob
import re
import shutil
import subprocess
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from modifyEventList import addCutoffToEventlist
#----------------------------------------------------------
def run_specsyn(machine,npe):
    """
    Submit job to grid engine using qsub command.
    Run directly from /data/Gemini using minos11-14.
    """
    rundir = os.getcwd()
    user = os.getlogin()
    print("You are logged in as: ",user)
#
#  parameter file and log file
#
    logfile = "specsyn.log"
    parfile = "parfile_gemini"
    subprocess.run(["readlink", "-f", parfile])
    print("Use parameter file: ",parfile)
#
#  some folder and path definitions
#
    executable = shutil.which("computeSyntheticsForSpecfemCoupling")
    print("Will run: ",executable)
#
#  job resources, options and command
#
    resources = ["-l","low","-l","hostname="+machine,"-l","memory=6G","-l","h_vmem=6G","-l","scratch_free=50G"]
    options = ["-cwd","-b","yes","-pe","mpi",npe,"-o",logfile,"-sync","yes", "-N", "coupling_synthetics"]
    command = ["mpirun","-n",npe,executable,parfile]
#
#  run job
#
    message = subprocess.check_output(["qsub"]+resources+options+command)
    print(message)
#-------------------------------------------------------------------------------------
#  read command line arguments
#
ap = ArgumentParser(description = "Submit computeSynthetcisForSpecfemCoupling job",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("machine", help = "Name of job machine")
ap.add_argument("npe", help = "Number of processes")
args = ap.parse_args()
#
# variable arguments
#
machine = args.machine
npe = args.npe

#
run_specsyn(machine,npe)
