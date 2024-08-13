#!/usr/bin/env python3
#
#  Plot a Gemini synthetic seismogram gather
#-------------------------------------------------------------------
import altair as alt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from makeSynseisAltairChart import makeSynseisAltairChart 
#-------------------------------------------------------------------
#  Read commandline arguments
#
ap = ArgumentParser(description = "Plot Gemini synthetic seismogram gather",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("hdflist", help = "List of names of HDF files with synthetic seismograms")
ap.add_argument("eventid", help = "Event ID")
ap.add_argument("--comp", help = "choose one of Z, N or E component (default=Z)", default='Z')
args = ap.parse_args()
hdflist = args.hdflist.split(' ')
print(hdflist)
eventid = args.eventid
comp = args.comp

shift = [0.0,-0.05,+0.05,-0.1,0.1]                      # shift of traces
color = ['black','darkseagreen','blue','cyan','red']    # color of traces
dx = [0,30,60,90,120]                                   # station annotation dy

chart = makeSynseisAltairChart(hdflist[0],eventid,comp,shift[0],color[0],dx[0])
for j,f in enumerate(hdflist[1:]):
    chart = chart + makeSynseisAltairChart(f,eventid,comp,shift[j+1],color[j+1],dx[j+1])

chart.interactive().show()
