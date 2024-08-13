#!/usr/bin/env python3
#
#  Submit computeGreenFKSpectraForASKI job to cluster machine
# ------------------------------------------------------------------
import os
import os.path
import glob
import re
import shutil
import subprocess
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


# ----------------------------------------------------------
def run_gfkask(machine, npe):
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
    print("You are logged in as: ", user)
    #
    #  split off /data/Gemini from rundir
    #
    subdir = "/".join(re.split("/", rundir)[3:])
    #
    #  parameter file and log file
    #
    logfile = "gfkask.log"
    parfile = "parfile_gemini"
    print("Use parameter file: ", parfile)
    #
    #  get information on earthmodel and boundary data file from parfile
    #
    with open(parfile, "r") as pf:
        for line in pf:
            if re.search("EARTH_MODEL", line): earthmodel = re.split("=", line)[1].strip()
            if re.search("EXTERNAL_NODES_HDF_FILE", line): boundaryData = re.split("=", line)[1].strip()
    print("Use earth model: ", earthmodel)
    print("Use boundary data: ", boundaryData)

    executable = shutil.which("computeGreenFKSpectraForASKI")
    print("Will run: ", executable)
    #
    #  create scratch folder
    #
    scratchdir = "/".join(["/rscratch", machine, user, subdir])
    if os.path.exists(scratchdir):
        print("ERROR: scratch directory exists. Please clean up before running job!")
        exit()
    os.makedirs(scratchdir)
    os.makedirs(scratchdir + "/gemini/")
    #
    #  copy parfile to scratch folder
    #
    scrpar = shutil.copy(parfile, scratchdir)
    print(parfile + " copied to ", scrpar)
    scrpar = shutil.copy(earthmodel, scratchdir + "/gemini/")
    print(earthmodel + " copied to ", scrpar)
    scrpar = shutil.copy(boundaryData, scratchdir)
    print(boundaryData + " copied to ", scrpar)
    #
    #  cd to scratchdir
    #
    os.chdir(scratchdir)
    #
    #  job resources, options and command
    #
    resources = ["-l", "low", "-l", "hostname=" + machine, "-l", "memory=500M", "-l", "h_vmem=500M", "-l",
                 "scratch_free=5G"]
    options = ["-cwd", "-b", "yes", "-pe", "mpi", npe, "-o", logfile, "-sync", "yes", "-N", "gfk-coupling"]
    command = ["mpirun", "-n", npe, executable, "-s", "6", parfile]
    #
    #  run job
    #
    message = subprocess.check_output(["qsub"] + resources + options + command)
    print(message)
    #
    #  cd back to rundir
    #
    os.chdir(rundir)


# -------------------------------------------------------------------------------------
#  read command line arguments
#
ap = ArgumentParser(description="Submit computeGreenFKSpectraForASKI job",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("machine", help="Name of job machine")
ap.add_argument("npe", help="Number of processes")
args = ap.parse_args()
#
# variable arguments
#
machine = args.machine
npe = args.npe
#
run_gfkask(machine, npe)
