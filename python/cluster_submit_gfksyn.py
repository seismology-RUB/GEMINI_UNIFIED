#!/usr/bin/env python3
#
#  Run computeGreenFKSpectraForSynthetics on cluster
#------------------------------------------------------------------
import os
import os.path
import glob
import re
import shutil
import subprocess
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
#----------------------------------------------------------
def run_gfksyn(machine,npe,gforversion):
    """
    Submit job to grid engine using qsub command.
    Assumes that parameters are specified in one available <parfile_gfk>.
    Assumes that there is one json file which contains earth model.
    Assumes that script is run from a folder whose name starts with
    2 standard path components such as /data/Gemini/subpath .
    A scratchdir with name /rscratch/machine/user/subpath is created.
    """
    rundir = os.getcwd()
    user = os.getlogin()
    print("You are logged in as: ",user)
#
#  split off /data/Gemini/ from rundir
#
    subdir = "/".join(re.split("/",rundir)[4:])
#
#  parameter file and log file
#
    logfile = "gfksyn.log"
    parfile = "parfile_gemini"
    print("Use parameter file: ",parfile)
#
#  earth model file, search for single available json file
#
    earthmodel = glob.glob("./*.json")[0]
    print("Use earth model: ",earthmodel)
#
#  some folder and path definitions
#
    executable = "/home/manuel/geminiUnified/bin/computeGreenFKSpectraForSynthetics"
    ldlibpath = "/rscratch/minos01/wolle/for-gfortran-"+gforversion+"/hdf5-1.10.5/hdf5/lib"
    print("Will run: ",executable)
#
#  create scratch folder
#
    scratchdir = "/".join(["/rscratch",machine,user,subdir,"/data"])
    if os.path.exists(scratchdir):
        print("ERROR: scratch directory exists. Please clean up before running job!")
        exit()
    os.makedirs(scratchdir)
#
#  copy parfile to scratch folder
#
    scrpar = shutil.copy(parfile,scratchdir)
    print(parfile+" copied to ",scrpar)
    scrpar = shutil.copy(earthmodel,scratchdir)
    print(earthmodel+" copied to ",scrpar)
#
#  cd to scratchdir
#
    os.chdir(scratchdir)
#
#  job resources, options and command
#
    resources = ["-l","low","-l","hostname="+machine,"-l","memory=500M","-l","h_vmem=500M","-l","scratch_free=5G"]
    options = ["-cwd","-b","yes","-pe","mpi",npe,"-o",logfile,"-v","LD_LIBRARY_PATH="+ldlibpath]
    command = ["mpirun","-n",npe,executable,"-s","6",parfile]
#
#  run job
#
    message = subprocess.check_output(["qsub"]+resources+options+command)
    print(message)
#
#  cd back to rundir
#
    os.chdir(rundir)
#-------------------------------------------------------------------------------------
#  read command line arguments
#
ap = ArgumentParser(description = "Submit computeGreenFKSpectraForSynthetics job",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("machine", help = "Name of job machine")
ap.add_argument("npe", help = "Number of processes")
ap.add_argument("--gforversion", help = "Version of gfortran, either 4.9 or 6.3", default = "4.9")
args = ap.parse_args()
#
# variable arguments
#
machine = args.machine
gforversion = args.gforversion
npe = args.npe
#
run_gfksyn(machine,npe,gforversion)
