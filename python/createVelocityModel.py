import json
from prem import prem
import numpy as np

model = list(np.loadtxt("/Users/tommy/git/Gemini/earthModels/ak135_diehl_v2.vel", unpack=True))

description = "Standard AK135 Earth model, modified with a crustal model for the Alps by Diehl et. al."
Type = "nm"
reference_frequency_hz = 1.0
anisotropic = False
earth_radius = 6371000.0
full_space = 0
equidistant_nodes = 0
num_layers = 0
num_nodes = len(model[0])
layer_inner_core = 2
layer_outer_core = 3
radius_inner_core = 1217500.0
radius_outer_core = 3479500.0
isfluid=[0,0,1,0,0,0,0,0,0,0]
radius_upper_layer_boundary =[]
index_upper_boundary = []
units = {"density": "kg/m3","velocity": "m/s","depth": "m"}
layerIndex = []

depthModel = list(model[0])
p_velocity = list(model[1])
#s_velocity = list(model[2])
# Sort velocity model such that it starts at Radius = 0 /maximum depth
depthModel = depthModel[::-1]
p_velocity = p_velocity[::-1]
#s_velocity = s_velocity[::-1]
# Convert depth from km to m
depthModel = [d*1000 for d in depthModel]
p_velocity = [d*1000 for d in p_velocity]
# Same aproximations on the S-wave velocity and densitiy as in the tomo model
s_velocity = [d/np.sqrt(3.) for d in p_velocity]
density = [310.*(d**0.25) for d in p_velocity]
q_kappa = []
q_mu = []

idx = 1
prevValue = 0
for depth, number in zip(depthModel, range(1, len(depthModel)+1)):
    cnt = list(depthModel).count(depth)
    tmprho, tmpvp, tmpvs, qMu, qKappa = prem(depth)
    if (depth == prevValue and cnt > 1) or depth == 0.:
        num_layers += 1
        index_upper_boundary.append(depthModel.index(depth)+1)
        radius_upper_layer_boundary.append(earth_radius - depthModel[depthModel.index(depth)])
        tmprho, tmpvp, tmpvs, qMu, qKappa = prem(depth-0.00000001)
    if depth == prevValue:
        idx += 1
    if qMu == 0:
        qMu = -1
    #print(number, earth_radius-depth, rho, p_velocity[depthModel.index(depth)], s_velocity[depthModel.index(depth)], qKappa, qMu)
    #density.append(rho*1000)
    q_kappa.append(qKappa)
    q_mu.append(qMu)
    layerIndex.append(idx)
    prevValue = depth


node_properties ={"layer_index":layerIndex, "depth":depthModel, "density":density, "P-velocity":p_velocity, "S-velocity":s_velocity, "Qkappa":q_kappa, "Qmu":q_mu}
ak135 = {"description":description,
         "type":Type,
         "reference_frequency_hz":reference_frequency_hz,
         "anisotropic":anisotropic,
         "earth_radius":earth_radius,
         "full_space": full_space,
         "equidistant_nodes":equidistant_nodes,
         "num_layers":num_layers,
         "num_nodes":num_nodes,
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
with open("/Users/tommy/git/Gemini/earthModels/ak135_diehl_v2.json", "w") as outfile:
    json.dump(ak135, outfile, indent=4)
