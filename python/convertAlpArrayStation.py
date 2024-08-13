#!/usr/bin/env python3
#--------------------------------------------------------------------------
#  Copyright 2021 Wolfgang Friederich (Ruhr-Universität Bochum)
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
#---------------------------------------------------------------
#  Tool to convert an AlpArray json station file to
#  one expected by Gemini
#---------------------------------------------------------------
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from operator import itemgetter
import json
import numpy as np
#
# read content of Json station file into dictionary
#
def readAlpArrayStation(stafile):
    with open(stafile,mode = 'r') as fp:
       stadict = json.load(fp)
    fp.close()
    return stadict
#
#  write content of dictionary into table
#
def writeGeminiStation(outfile,stadict):
    
    with open(outfile,mode = 'w') as fp:
       fp.write("{0:1}\n".format('S'))
       prevsta = 'None'
       for key,val in sorted(stadict.items()):
          [net,staname] = key.split('.')
          lat = val['latitude']
          alt = val['elevation']
          lon = val['longitude']
          fp.write("{1:8}{0:4}{2:12.5f}{3:12.5f}{4:12.1f}\n".format(net,staname,lat,lon,alt))
    fp.close()
#
#  main program
#
ap = ArgumentParser(description = "Convert AlpArray station Json file to Gemini type",
                    formatter_class = ArgumentDefaultsHelpFormatter)
ap.add_argument('aafile',help = 'Json AlpArray station file')
ap.add_argument('gemfile',help = 'Output station file for Gemini',default = 'aastations.out')
args = ap.parse_args()
aafile = args.aafile
gemfile = args.gemfile
#
stadict = readAlpArrayStation(aafile)
writeGeminiStation(gemfile,stadict)
