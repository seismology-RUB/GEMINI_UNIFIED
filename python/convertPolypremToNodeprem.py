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
#---------------------------------------------------------------
#  Tool to convert the polynomial standard earth model PREM into
#  a node model (written as JSON file)
#---------------------------------------------------------------
import sys
import argparse
import json
import numpy as np
#
def convertPolypremToNodeprem(argv):
    parser = argparse.ArgumentParser(description="Convert polynomial PREM to a node model")
    parser.add_argument('polyprem', help='Polynomial PREM JSON file')
    parser.add_argument('--outfile', help='JSON output with node PREM model', default = 'nodeprem.out')
    parser.add_argument('--rhs', help='Radius of top of halfspace or inner sphere', default = '100.0')
    args = parser.parse_args()
    nodefile = args.outfile
    polyfile = args.polyprem
    rhs = float(args.rhs)
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
    nlhs = np.amin(np.where(rb > rhs))               # Layer of rhs
    if nlhs > noc:                                   # Model wthout core
        noc = -1
        roc = -9999.
        ric = -9999.
    elif nlhs == noc:                                # inner sphere in outer core
        nic = -1                                     # model without inner core
        noc = 2
        ric = -9999.
    elif nlhs < noc:                                 # inner sphere in inner core
        nic = 2                                      # model with inner and outer core
        noc = 3
    #   
    print("Original layer index of rhs: ",nlhs+1)
    if nlhs == 0:                              
        rbot = 0.0                                   # bottom of layer
        dr = rb[nlhs]/numival[nlhs]                  # dr in layer
    else:
        rbot = rb[nlhs-1]                            # bottom of layer
        dr = (rb[nlhs]-rbot)/numival[nlhs]           # dr in layer
    #
    nval_skipped = np.floor((rhs-rbot)/dr)           # number of skipped intervals in layer nlhs
    if nval_skipped == 0 & nlhs == 0:
        print("Error: ","dr in inner core is greater than rhs")
        print("Error: ","Choose more intervals in the inner core")
        exit()
    rstart = rbot+dr*nval_skipped                    # radius of first node
#
#  calculate new values according to bottom of model
#
    nlay_new = nlay-nlhs+1                           # new number of layers: remove nlhs layers but add halfspace
    isfluid_new = [None]*nlay_new                    # allocate fluid flags
    isfluid_new[0] = isfluid[nlhs]
    isfluid_new[1] = isfluid[nlhs]
    isfluid_new[2:nlay_new] = isfluid[nlhs+1:nlay]
    numival_new = np.zeros(nlay_new,dtype=int)       # allocate new interval counts
    numival_new[0] = 0                               # zero intervals in halfspace (only top node)
    numival_new[1] = numival[nlhs]-nval_skipped      # layer above halfspace has fewer intervals now
    numival_new[2:nlay_new] = numival[nlhs+1:nlay]   # remaining layers are unchanged
    rb_new = np.zeros(nlay_new)                      # allocate new radii of upper layer boundaries
    rb_new[0] = rstart                               # top of halfspace is at rstart
    rb_new[1:nlay_new] = rb[nlhs:nlay]               # remaining ones
#
    nk = np.sum(numival_new)+nlay_new                # number of intervals plus nodes at layer top
    print("Expected number of nodes ",nk)
    rk = np.zeros(nk)                                # allocate node radii
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
    print(rb_new)
#-----------------------------------------------------------------------------------------------
#  Start calculation of nodes and properties in all layers
#
    for nl in range(0,nlay_new):                              # Loop over layers
        n = numival_new[nl]                                   # number of intervals in current layer
        if nl > 0:
            dr = (rb_new[nl]-rb_new[nl-1])/n                  # current dr
            nlo = nl+nlhs-1                                   # original layer index
            rbot = rb_new[nl-1]
        else:
            nlo = nlhs                                        # halfspace belongs to nlhs
            rbot = rstart
        for j in range(jn,jn+n):                              # nodes in layer
            rk[j] = rbot+(j-jn)*dr
        iktop[nl] = int(jn+n+1)                               # index of layer top node (add one for Fortran use later)
        rk[jn+n] = rb_new[nl]                                 # add top node of layer
#
        for j in range(jn,jn+n+1):                            # loop over nodes in this layer
            x = rk[j]/re
            layind.append(nl+1)
            dens.append(c0[nlo,1]+c1[nlo,1]*x+c2[nlo,1]*x**2+c3[nlo,1]*x**3)
            vpv.append(c0[nlo,2]+c1[nlo,2]*x+c2[nlo,2]*x**2+c3[nlo,2]*x**3)
            vph.append(c0[nlo,3]+c1[nlo,3]*x+c2[nlo,3]*x**2+c3[nlo,3]*x**3)
            vsv.append(c0[nlo,4]+c1[nlo,4]*x+c2[nlo,4]*x**2+c3[nlo,4]*x**3)
            vsh.append(c0[nlo,5]+c1[nlo,5]*x+c2[nlo,5]*x**2+c3[nlo,5]*x**3)
            qmu.append(c0[nlo,6]+c1[nlo,6]*x+c2[nlo,6]*x**2+c3[nlo,6]*x**3)
            qkap.append(c0[nlo,7]+c1[nlo,7]*x+c2[nlo,7]*x**2+c3[nlo,7]*x**3)
            eta.append(c0[nlo,8]+c1[nlo,8]*x+c2[nlo,8]*x**2+c3[nlo,8]*x**3)
        jn = jn+n+1                                           # next node index to be used
#------------------------------------------------------------------------
#  Write result to json node model output file
#  Pack everything into a dictionary and the use the json.encoder
#
    depth = (re-rk)*1.e3
    dict = {
        "description" : "Standard earth model PREM converted to node model",
        "type" : "nm",
        "reference_frequency_hz" : ct["reference_frequency_hz"],
        "anisotropic" : ct["anisotropic"],
        "earth_radius" : ct["earth_radius"],
        "full_space" : ct["full_space"],
        "equidistant_nodes" : ct["equidistant_nodes"],
        "num_layers" : int(nlay_new),
        "num_nodes" : int(nk),
        "layer_inner_core" : nic,
        "layer_outer_core" : noc,
        "radius_inner_core" : ric,
        "radius_outer_core" : roc,
        "radius_upper_layer_boundary_km" : list(rb_new),
        "isfluid" : isfluid_new,
        "index_upper_boundary" : iktop,
        "units" : units,
        "node_properties" : {
            "layer_index" : layind,
            "depth" : list(depth),
            "density" : dens,
            "PV-velocity" : vpv,
            "PH-velocity" : vph,
            "SV-velocity" : vsv,
            "SH-velocity" : vsh,
            "Qkappa" : qkap,
            "Qmu" : qmu,
            "eta" : eta,
        }
    }
    #----------------------------------------------------------------------
    #  Write to output file
    #
    fp = open(nodefile, mode = "w")
    print(json.JSONEncoder(indent=4).encode(dict),file = fp)
    fp.close()
#
if __name__ == "__main__":
    convertPolypremToNodeprem(sys.argv)
