#-------------------------------------------------------------------
#  produce an Altair chart from given hdf synseis file,
#  event an component
#
import h5py
import pandas as pd
import altair as alt
from altair import datum
import numpy as np

def makeSynseisAltairChart(hdffile,eventid,comp,shift,col,dx):
    """
    Produce an Altair chart from a given hdffile, event ID and component
    """
    fid = h5py.File(hdffile,'r')                      # open file
    if (eventid in fid.keys()):                       # check if event exists
        grp = fid[eventid]
        print("Group "+eventid+" was found in file!")
    else:
        print("Group "+eventid+" not found in file!")
        exit()
    stations = list(grp.keys())                       # List of datasets = station names
    
    ds = grp[stations[0]]                             # Take first station
    validcomps = str(ds.attrs["Components"])          # read components from dataset
    ic = validcomps.find(comp)-2                      # index of component (0=Z,1=N,2=E)
    print("Component index IC = ",ic)
    if (ic < 0):                                      # check validity of component
        print("Component "+comp+" is not available!")
        exit()
    nsamp = ds.shape[1]                               # get nsamp
    tanf,dt = ds.attrs["Time"]                        # get tanf and dt
    nstat = len(grp.keys())                           # number of stations
    print("Number of samples: ",nsamp)
    print("Start time and dt: ",tanf,dt)
    print("Number of stations found for event:",nstat)

    disp = np.zeros(nstat*nsamp)                      # space to concatenate all traces into one array
    ts = np.zeros(nstat*nsamp)                        # space for doing this for sample times too
    st = []                                           # empty list for station names
    ts1 = np.arange(0,nsamp)*dt                       # unique sample times
    for j,ds in enumerate(grp.values()):              # Loop over stations
        print("Name of station: ",ds.name.split('/')[2])
        ia = j*nsamp                                          # start in concatenated array
        ib = (j+1)*nsamp                                      # end in concantenated array
        disp[ia:ib] = 0.48*ds[ic,:]/np.amax(ds[ic,:])+float(j)+shift # reduce amplitudes to +/-0.48 and add 1 per trace
        ts[ia:ib] = ts1[:]                                    # insert sample times
        st1 = [ds.name.split('/')[2]] * nsamp                 # nsamp times the same station name
        st = st+st1
    
    data = pd.DataFrame({'station': st,'time': ts,'disp': disp})  # create pandas data frame

    alt.data_transformers.disable_max_rows()                  # allow more than 5000 rows
    lines = alt.Chart(data).mark_line(                        # chart with lines
        stroke=col,strokeWidth=1.5
    ).encode(
        alt.X('time',type='quantitative',title='Time in seconds'),
        alt.Y('disp',type='quantitative',title=comp+' displacement'),
        alt.Detail('station'),
    ).properties(
        width=750, height=850
    )
    
    text = lines.mark_text(                                   # add station names to start of traces
       align='left',baseline='bottom',dy=20,dx=dx,size=10,color=col
    ).transform_filter(
        datum.time == ts1[0]
    ).encode(
        text='station',y='disp:Q'
    )

    chart = lines+text
    return chart                                             # return Altair chart
