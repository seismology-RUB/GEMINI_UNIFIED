#!/usr/bin/env python3
"""
Add earliest P-wave arrival time in the simulation Box (- Buffer) to the event list.
"""
from obspy.taup import TauPyModel
from seismicEvent import SeismicEvent
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
import os
from specfemBox import SpecfemBox


def addCutoffToEventlist(eventFile, buffer, folder = "/"):
    """
    Add a cutoff time to each event in the event file.

    :param eventFile: Name of the eventFile.
    :param buffer: Time ahead of the first arrival time of the P-wave.
    :param folder: Path to where the data is stored.
    :return:
    """
    box = SpecfemBox(file = folder + "specfem_boundary_data.hdf")
    model = TauPyModel(model="prem")
    tmpFile, tmpFile_path = mkstemp()
    with open(folder + eventFile, "r") as eventList, open(tmpFile, "w") as newFile:
        for line in eventList:
            event = SeismicEvent()
            eventData = line.split()
            if eventData[0] == "S" or eventData[0] == "C":
                event.readCoordSystem(eventData)
                newFile.write(line)
            else:
                event.readFromFile(eventData)
                event.calculateCutoffTime(box, model, buffer)
                if str(eventData[-1]) == str(event.cutoffTime):
                    print("Cutoff time already added.")
                    newFile.write(line)
                else:
                    newLine = line.rstrip("\n") + " " + str(event.cutoffTime) + "\n"
                    newFile.write(newLine)
                    print("Added {:6.2f} seconds as cutofftime.".format(event.cutoffTime))
    # Copy the file permissions from the old file to the new file
    copymode(folder + eventFile, tmpFile_path)
    # Remove original file
    remove(folder + eventFile)
    # Move new file
    move(tmpFile_path, folder + eventFile)


if __name__ == "__main__":
    ap = ArgumentParser(description = "Add earliest P-wave arrival time in the simulation Box (-Buffer) to the event list",
                    formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("buffer", help = "Padding for the cutoff time of the seismogram", type=float)
    ap.add_argument("eventFile", help = "File containing the events to be processed.")
    args = ap.parse_args()
    #
    # variable arguments
    #
    buffer = args.buffer
    eventFile = args.eventFile
    folder = os.getcwd() + "/"
    addCutoffToEventlist(eventFile, buffer, folder=folder)

