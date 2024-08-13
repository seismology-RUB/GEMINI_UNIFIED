#!/usr/bin/env python3
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
#  Tool to convert an ASKI 1D background model to a Json model
#  file needed by Gemini
#---------------------------------------------------------------
import sys
import argparse
import json
import numpy as np
#
def convertASKIBackgroundToJson(argv):
    parser = argparse.ArgumentParser(description="Convert ASKI background model to GEMINI Json file")
    parser.add_argument('askibg', help='ASKI background model file')
    parser.add_argument('--outfile', help='Gemini Json model file output', default = 'askibg.json')
    parser.add_argument('--rbot', help='radius of homogeneous inner sphere in m', default = '100000.0')
    args = parser.parse_args()
    jsonfile = args.outfile
    askibg = args.askibg
    rbot = float(args.rbot)
#------------------------------------------------------------------
#  Read ASKI background model
#
    fp = open(askibg,mode="r")
    content = fp.readlines()
    fp.close()
#
#  extract file content
#
    description = content[0]
    nlay = int(content[2])
    nodes_per_layer = list(map(int,content[3].split()))
    nnodes = sum(nodes_per_layer)
    prop = np.zeros((nnodes,6))
    for j,line in enumerate(content[4:]):
       prop[j,:] = np.array(list(map(float,line.split())))
#
#  start to fill dictionary
#
    Type = "nm"
    reference_frequency_hz = 1.0
    units = {"density": "kg/m3","velocity": "m/s","depth": "m"}
    anisotropic = False
    earth_radius = 6371000.0
    full_space = 0
    equidistant_nodes = 0
    layer_inner_core = 2
    layer_outer_core = 3
    radius_inner_core = earth_radius-prop[sum(nodes_per_layer[0:nlay-1]),0]
    radius_outer_core = earth_radius-prop[sum(nodes_per_layer[0:nlay-2]),0]
#
#  go from bottom to top through ASKI background model 
#  and assign values to Gemini earth model
#
    depth = []
    rho =[]
    vp = []
    vs = []
    qmu = []
    qkap = []
    layer_index = []
    isfluid = []
    radius_upper_layer_boundary = []
    index_upper_boundary = []
    nl = 1
    kg = 0
    for k in range(nnodes,0,-1): 
        if prop[k-1,0] > earth_radius-rbot:   # skip nodes below rbot
            continue
        depth.append(prop[k-1,0])
        rho.append(prop[k-1,1])
        vp.append(prop[k-1,2])
        vs.append(prop[k-1,3])
        qkap.append(prop[k-1,5])
        if prop[k-1,3] > 0.0:
            qmu.append(prop[k-1,4])
        else:
            qmu.append(-1.0)
        layer_index.append(nl)
        kg = kg+1
        if nl == 1:                      # add double node at rbot
            depth.append(prop[k-1,0])
            rho.append(prop[k-1,1])
            vp.append(prop[k-1,2])
            vs.append(prop[k-1,3])
            qkap.append(prop[k-1,5])
            layer_index.append(nl+1)
            if prop[k-1,3] > 0.0:
                qmu.append(prop[k-1,4])
            else:
                qmu.append(0.0)
            kg = kg+1
            radius_upper_layer_boundary.append(earth_radius-depth[-2])
            index_upper_boundary.append(kg-1)
            isfluid.append(int(prop[k-1,3] == 0.0)) 
            nl = nl+1
            continue
        if abs(depth[-1]-depth[-2]) < 1.e-3:  # new layer if double node
            nl = nl+1
            isfluid.append(int(prop[k,3] == 0.0))
            radius_upper_layer_boundary.append(earth_radius-depth[-2])
            index_upper_boundary.append(kg-1)
        if prop[k-1,0] <= 0.00:   # surface reached
            radius_upper_layer_boundary.append(earth_radius-depth[-1])
            index_upper_boundary.append(kg)
            isfluid.append(int(prop[k-1,3] == 0.0))
            break
    ngnodes = kg
    nlgem = nl
#
    print("Core radii: ",radius_inner_core,radius_outer_core)
    print("Fluid indices: ",isfluid)
    print("Upper layer boundaries: ",radius_upper_layer_boundary)
    print("Node indices of upper boundaries: ",index_upper_boundary)
#
    node_properties ={"layer_index":layer_index, "depth":depth, "density":rho, 
                      "P-velocity":vp, "S-velocity":vs, "Qkappa":qkap, "Qmu":qmu}
    ak135 = {"description":description,
             "type":Type,
             "reference_frequency_hz":reference_frequency_hz,
             "anisotropic":anisotropic,
             "earth_radius":earth_radius,
             "full_space": full_space,
             "equidistant_nodes":equidistant_nodes,
             "num_layers":nlgem,
             "num_nodes":ngnodes,
             "layer_inner_core":layer_inner_core,
             "layer_outer_core":layer_outer_core,
             "radius_inner_core":radius_inner_core,
             "radius_outer_core":radius_outer_core,
             "radius_upper_layer_boundary_m":radius_upper_layer_boundary,
             "isfluid":isfluid,
             "index_upper_boundary":index_upper_boundary,
             "units":units,
             "node_properties":node_properties
             }
    with open(jsonfile, "w") as fout:
        json.dump(ak135, fout, indent=4)

if __name__ == "__main__":
    convertASKIBackgroundToJson(sys.argv)
