#!/usr/bin/env python3
#
#  Copy output files of gfkask-run back from scratch to fileserver
#------------------------------------------------------------------
import os
import os.path
import glob
import re
import shutil
import subprocess
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
#-------------------------------------------------------------------
def run_copy(machine,cleanflag):
    """
    Copy HDF output files of computeGreenFKSpectraForASKI job
    from scratch back to file server /data/Gemini/subdir.
    Assumes that scratch is at /rscratch/machine/user/subdir.
    """
    rundir = os.getcwd()
    user = os.getlogin()
    print("You are logged in as: ",user)
#
#  split off /data/Gemini from rundir
#
    subdir = "/".join(re.split("/",rundir)[3:])
#
#  parameter file
#
    parfile = "parfile_gemini"
    print("Use parameter file: ",parfile)
#
#  get dsvbasename and source type and number of source nodes from parfile
#
    with open(parfile,"r") as pf:
        for line in pf:
            if re.search("DSVBASENAME",line): dsvbasename = re.split("=",line)[1].strip()
            if re.search("SOURCE_TYPE",line): sourcetype = re.split("=",line)[1].strip()
            if re.search("NSNOD",line): nsnod = int(re.split("=",line)[1].strip())
#
#  copy gfk meta file back to rundir
#
    print(dsvbasename)
    scratchdir = "/".join(["/rscratch",machine,user,subdir])
    gfkfile = shutil.copy(scratchdir+"/"+dsvbasename+".meta",rundir + "/gemini/")
    print("GFK meta file copied from scratch to: ",gfkfile)
#
#  copy gfk-files back to rundir
#
    for n in range(1,nsnod+1):
        gfkfile = shutil.copy(scratchdir+"/"+dsvbasename+"."+sourcetype+"."+str(n).zfill(3),rundir + "/gemini/")
        print("GFK-file for source node "+str(n)+" copied from scratch to: ",gfkfile)
#
#  copy gfk-file for FORCE back to rundir
#
    gfkfile = shutil.copy(scratchdir+"/"+dsvbasename+".FORCE"+"."+str(nsnod+1).zfill(3),rundir + "/gemini/")
    print("GFK-file for source node "+str(n)+" copied from scratch to: ",gfkfile)
#
#  delete last folder in scratchdir, if cleanflag is set
#
    if cleanflag:
        shutil.rmtree(scratchdir)
#-------------------------------------------------------------------------------------
#  read command line arguments
#
ap = ArgumentParser(description = "Copy output of computeGreenFKSpectraForSynthetics from scratch",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("machine", help = "Name of job machine")
ap.add_argument("--clean", action = "store_true", help = "Delete scratch folder")
args = ap.parse_args()
#
# variable arguments
#
machine = args.machine
cleanflag = args.clean
#
run_copy(machine,cleanflag)
