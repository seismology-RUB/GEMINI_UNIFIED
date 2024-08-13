#!/usr/bin/env python3
#
#  Plot rays
#-------------------------------------------------------------------
import altair as alt
import h5py
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
#-------------------------------------------------------------------
#  Read commandline arguments
#
ap = ArgumentParser(description = "Plot rays",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("hdffile", help = "HDF file with ray table")
ap.add_argument("--pstep", help = "stepping of slowness index", default='10')
args = ap.parse_args()
hdffile = args.hdffile
step = int(args.pstep)

fid = h5py.File(hdffile,'r')                      # open file
slowness = fid['rayParameters']
rtp = fid['turningPointRadii']
rnod = fid['receiverRadii']
raytable = fid['rayTable']
ntp = np.size(rtp)

alt.data_transformers.disable_max_rows()
delta = raytable[0,:,0]*180./np.pi           # fixed slowness, varying receiver radius
p = slowness[0]
print('Slowness: ',p,' Turning radius: ',rtp[0])

data = pd.DataFrame({'delta': delta,'radius': rnod})  # create pandas data frame
chart = alt.Chart(data).mark_line(
    stroke='blue',strokeWidth=1.5
    ).encode(
        alt.X('radius',type='quantitative',title='Radius in km'),
        alt.Y('delta',type='quantitative',title='Distance in degrees'),
    ).properties(
        width=750, height=850
    )

for js in np.arange(step,ntp,step):
    print('Slowness: ',slowness[js],' Turning radius: ',rtp[js])
    delta = raytable[js,:,0]*180./np.pi           # fixed slowness, varying receiver radius
    data = pd.DataFrame({'delta': delta,'radius': rnod})  # create pandas data frame
    chart = chart + alt.Chart(data).mark_line(
        stroke='blue',strokeWidth=1.5
    ).encode(x='radius',y='delta'
    )
    
chart.interactive().show()
