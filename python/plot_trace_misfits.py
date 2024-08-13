#!/usr/bin/env python3
#--------------------------------------------------------------------------
#  Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum)
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
#  Tool to convert an ASKI 1D background model to a Json model
#  file needed by Gemini
#---------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np

sms = []
stat = []
misfit = []
snr = []
rmsrat = []
stf = []
with open("trace_misfits_snr.txt","r") as fp:
   for line in fp:
      lsp = line.split()
      sms.append([lsp[0],float(lsp[1]),float(lsp[2]),float(lsp[3])])

smssort = sorted(sms, key=lambda station: station[1],reverse = True)

with open("trace_misfits_snr_sorted.txt","w") as fp:
   for station in smssort:
      stat.append(station[0])
      misfit.append(station[1])
      snr.append(station[2])
      rmsrat.append(station[3])
      fp.write("plotAsciiSeismogramGather -w 15 "+"'data_"+str(station[0])+"_UP "+"conv_"+str(station[0])+"_UP'   "
               +str(station[1])+'  '+str(station[2])+'  '+str(station[3])+'\n')

with open("stf.txt","r") as fp:
   line = fp.readline()
   line = fp.readline()
   lsp = fp.readline().split()
   nsamp = int(lsp[0])
   dt = float(lsp[1])
   tanf = float(lsp[2])
   for s in fp.readline().split():
      stf.append(float(s))

fig = plt.figure(figsize=(20,12))
ax = fig.add_subplot(2,2,1)
x = [tanf+i*dt for i in range(0,nsamp)]
ax.grid(b=True,which='major',axis='both',ls=':')
ax.plot(x,stf)
ax.set_xlabel("Time", fontsize=14)
ax.set_ylabel("Moment rate function", fontsize=14)

bx = fig.add_subplot(2,2,2)
bx.grid(b=True,which='major',axis='both',ls=':')
bx.hist(snr,bins=50)
bx.set_xlabel("SNR", fontsize=14)
bx.set_ylabel("Count", fontsize=14)

cx = fig.add_subplot(2,2,3)
x = range(1,len(stat)+1)
cx.grid(b=True,which='major',axis='both',ls=':')
cx.set_xlabel("Station number", fontsize=14)
cx.set_ylabel("Trace misfit (*100)",fontsize=14)
cx.plot(x,misfit,'.')

cx = fig.add_subplot(2,2,4)
x = range(1,len(stat)+1)
cx.grid(b=True,which='major',axis='both',ls=':')
cx.set_xlabel("Station number", fontsize=14)
cx.set_ylabel("Rms ratio Data/Syn", fontsize=14)
cx.plot(x,rmsrat,'.')
plt.savefig('plot_trace_misfits.png', dpi=300, bbox_inches='tight')
plt.show()
