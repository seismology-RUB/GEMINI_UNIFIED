#--------------------------------------------------------------------------
#  Copyright 2017 Wolfgang Friederich (Ruhr-Universität Bochum)
#
#  This file is part of GeminiUnified.
#
#  GeminiUnified is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  any later version.
#
#  GeminiUnified is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with GeminiUnified.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------
#  Python module to run entire Gemini workflow from Green FK-Spectra
#  to synthetic seismograms
#
#  Usage: run_gemini_fk_to_synthetics setup_dir
#
#  It is assumed that parameter files, event and station file have already been created
#-----------------------------------------------------------
#  Parse commandline
#
import subprocess
import os
import argparse
#-----------------------------------------------------------
parser = argparse.ArgumentParser(description="Run full Gemini workflow to get synthetic seismograms")
parser.add_argument('setup_dir', help='Folder with setup files')
parser.add_argument("--gemini-bin", help="Folder with Gemini binaries (slash terminated)",default = "")
parser.add_argument("--skip_gfk", action="store_true", help="Skip computation of GFK spectra")
args = parser.parse_args()
#
bin = args.gemini_bin
setup = args.setup_dir
skipgfk = args.skip_gfk
#-----------------------------------------------------------
#  Change to setup folder, get original folder first
#
origdir = os.getcwd()
try:
    os.chdir(setup)
except OSError:
    print("Directory does not exist")
    exit()
print("Changed to folder: ",os.getcwd())
#-----------------------------------------------------------
# Run computeGreenFKSpectraForSynthetics
#
if (not skipgfk):
    if (not os.path.exists("parfile_gfk")):
        print("There is no parameter file called 'parfile_gfk'")
        exit()
    try:
        print("Running computeGreenFKSpectraForSynthetics")
        message = subprocess.check_output("mpirun -np 2 "+bin+"computeGreenFKSpectraForSynthetics "+"parfile_gfk",shell = True)
        print(message)
    except subprocess.CalledProcessError:
        print("execution of computeGreenFKSpectraForSynthetcs failed")
        exit()
#-----------------------------------------------------------
# Run computeEventFilter
#
if (not os.path.exists("parfile_syn")):
    print("There is no parameter file called 'parfile_syn'")
    exit()
if (not os.path.exists("event_filter")):
    os.mkdir("event_filter")
try:
    message = subprocess.check_output(bin+"computeEventFilter "+"parfile_syn",shell = True)
    print(message)
except subprocess.CalledProcessError:
    print("execution of computeEventFilter failed")
    exit()
#------------------------------------------------------------
# Run computeStationFilter
#
if (not os.path.exists("station_filter")):
    os.mkdir("station_filter")
try:
    message = subprocess.check_output(bin+"computeStationFilter "+"parfile_syn",shell = True)
    print(message)
except subprocess.CalledProcessError:
    print("execution of computeStationFilter failed")
    exit()
#------------------------------------------------------------
# Run computeSyntheticSeismograms
#
if (not os.path.exists("synseis")):
    os.mkdir("synseis")
try:
    message = subprocess.check_output(bin+"computeSyntheticSeismograms "+"parfile_syn",shell = True)
    print(message)
except subprocess.CalledProcessError:
    print("execution of computeSyntheticSeismograms failed")
    exit()
#------------------------------------------------------------
# Run buildSeismogramGather
#
try:
    message = subprocess.check_output(bin+"buildSeismogramGather "+"parfile_syn",shell = True)
    print(message)
except subprocess.CalledProcessError:
    print("execution of buildSeismogramGather failed")
    exit()

