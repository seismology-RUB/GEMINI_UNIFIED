#--------------------------------------------------------------------------
#  Copyright 2021 Wolfgang Friederich (Ruhr-Universitaet Bochum)
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
#  Tool to convert the polynomial standard earth model PREM into
#  an ASKI 1D background model
#---------------------------------------------------------------
import sys
import argparse
import json
import numpy as np
#
def convertPolypremToASKIBackground(argv):
    parser = argparse.ArgumentParser(description="Convert polynomial PREM to an ASKI background model")
    parser.add_argument('polyprem', help='Polynomial PREM JSON file')
    parser.add_argument('--outfile', help='List output according to ASKI specs', default = 'prem_1d_aski.out')
    parser.add_argument('--rbot', help='radius of deepest node', default = '100.0')
    args = parser.parse_args()
    askifile = args.outfile
    polyfile = args.polyprem
    rbot = float(args.rbot)
#------------------------------------------------------------------
#  Read polynomial model
#
    fp = open(polyfile,mode="r")
    ct = json.load(fp)
    fp.close()
#-------------------------------------------------------------------
#  Extract some variables from file content
#
    units = ct["units"]
    units["depth"] = "m"                             # add depth:m pair
    c0 = np.array(ct["layprop_order_zero"])          # polynomial coefficients
    c1 = np.array(ct["layprop_order_one"])
    c2 = np.array(ct["layprop_order_two"])
    c3 = np.array(ct["layprop_order_three"])
    nlay = ct["num_layers"]
    nic = ct["layer_inner_core"]
    noc = ct["layer_outer_core"]
    ric = ct["radius_inner_core"]
    roc = ct["radius_outer_core"]
    numival = np.array(ct["intervals_per_layer"],dtype=int)
    rb = np.array(ct["radius_upper_layer_boundary_km"])
    isfluid = ct["isfluid"]
    re = ct["earth_radius"]
#-------------------------------------------
    nlbot = np.amin(np.where(rb > rbot))             # Layer of deepest node
    if nlbot > noc:                                  # Model wthout core
        noc = -1
        roc = -9999.
        ric = -9999.
    elif nlbot == noc:                               # inner sphere in outer core
        nic = -1                                     # model without inner core
        noc = 2
        ric = -9999.
    elif nlbot < noc:                                # inner sphere in inner core
        nic = 2                                      # model with inner and outer core
        noc = 3
    #   
    print("Original layer index of rbot: ",nlbot)
    if nlbot == 0:                              
        dr = rb[nlbot]/numival[nlbot]                # dr in layer
        nval_skipped = np.floor(rbot/dr)             # number of skipped intervals in layer nlbot
        rstart = dr*nval_skipped                     # radius of first node
        print("deepest node in inner core, start at :",rstart)
    else:
        dr = (rb[nlbot]-rb[nlbot-1])/numival[nlbot]     # dr in layer
        nval_skipped = np.floor((rbot-rb[nlbot-1])/dr)
        rstart = rb[nlbot-1]+dr*nval_skipped         # radius of first node
        print("layer bottom at ",rb[nlbot-1]," layer top at ",rb[nlbot]," start at :",rstart)
    #
    if nval_skipped == 0 and nlbot == 0:
        print("Error: ","dr in inner core is greater than rbot")
        print("Error: ","Choose more intervals in the inner core")
        exit()
#
#  calculate new values according to bottom of model
#
    nlay_new = nlay-nlbot                            # new number of layers
    isfluid_new = [None]*nlay_new                    # allocate fluid flags
    isfluid_new[0] = isfluid[nlbot]
    isfluid_new[1:nlay_new] = isfluid[nlbot+1:nlay]
    numival_new = np.zeros(nlay_new,dtype=int)       # allocate new interval counts
    numival_new[0] = numival[nlbot]-nval_skipped     # bottom layer has fewer intervals now
    numival_new[1:nlay_new] = numival[nlbot+1:nlay]  # remaining layers are unchanged
    rb_new = np.zeros(nlay_new)                      # allocate new radii of upper layer boundaries
    rb_new[0:nlay_new] = rb[nlbot:nlay]              # remaining ones
    print("Top radii of remaining layers: ",rb_new)
    print("Intervals in remaining layers: ",numival_new)
#
    nk = np.sum(numival_new)+nlay_new                # number of intervals plus nodes at layer top
    print("Expected number of nodes ",nk)
#
#  open ASKI background model output file
#
    fp = open(askifile, mode = "w")
    fp.write("PREM as ASKI 1D background model\n")
    fp.write(str(0.0)+"\n")                                        # zmax
    fp.write(str(nlay_new)+"\n")                                   # number of layers
    fp.write(" ".join(map(str, reversed(numival_new+1)))+'\n')     # number of nodes per layer
#
#  initialize arrays
#
    depth = np.zeros(nk)                             # allocate depth of nodes in m
    iktop = [0]*nlay_new                             # allocate node indices upper layer boundaries
    layind = []                                      # list with layer indices
    dens = []                                        # list with densities
    vpv = []                                         # list with vpv-values
    vph = []                                         # list with vph-values
    vsv = []                                         # list with vsv-values
    vsh = []                                         # list with vsh-values
    eta = []                                         # list with eta-values
    qkap = []                                        # list with qkap-values
    qmu = []                                         # list with qmu-values
    jn = 0                                           # next node index to be used
#-----------------------------------------------------------------------------------------------
#  Start calculation of nodes and properties in all layers
#
    for nl in range(nlay_new-1,-1,-1):                        # Loop over layers (top down)
        n = numival_new[nl]                                   # number of intervals in current layer
        if nl > 0:
            dr = (rb_new[nl]-rb_new[nl-1])/n                  # current dr
        else:
            dr = (rb_new[nl]-rstart)/n
        nlo = nl+nlbot                                        # original layer index
        r2 = rb_new[nl]                                       # upper layer boundary
        for j in range(jn,jn+n):                              # nodes in layer
            depth[j] = re-r2+(j-jn)*dr
        if nl > 0:
            depth[jn+n] = re-rb_new[nl-1]                     # add bottom node of layer
        else:
            depth[jn+n] = re-rstart
    #
        for j in range(jn,jn+n+1):                            # loop over nodes in this layer
            x = (re-depth[j])/re
            layind.append(nl+1)
            dens.append(c0[nlo,1]+c1[nlo,1]*x+c2[nlo,1]*x**2+c3[nlo,1]*x**3)
            vpv.append(c0[nlo,2]+c1[nlo,2]*x+c2[nlo,2]*x**2+c3[nlo,2]*x**3)
            vph.append(c0[nlo,3]+c1[nlo,3]*x+c2[nlo,3]*x**2+c3[nlo,3]*x**3)
            vsv.append(c0[nlo,4]+c1[nlo,4]*x+c2[nlo,4]*x**2+c3[nlo,4]*x**3)
            vsh.append(c0[nlo,5]+c1[nlo,5]*x+c2[nlo,5]*x**2+c3[nlo,5]*x**3)
            qmu.append(c0[nlo,6]+c1[nlo,6]*x+c2[nlo,6]*x**2+c3[nlo,6]*x**3)
            qkap.append(c0[nlo,7]+c1[nlo,7]*x+c2[nlo,7]*x**2+c3[nlo,7]*x**3)
            eta.append(c0[nlo,8]+c1[nlo,8]*x+c2[nlo,8]*x**2+c3[nlo,8]*x**3)
    #
    #  write to file 
    #
            fp.write("{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12.3f}".format(
                   depth[j]*1.e3,dens[j]*1.e3,vpv[j]*1.e3,vsv[j]*1.e3,qmu[j],qkap[j])+'\n')
    #
        jn = jn+n+1                                           # next node index to be used
#
#  close file
#
    fp.close()
#
if __name__ == "__main__":
    convertPolypremToASKIBackground(sys.argv)
