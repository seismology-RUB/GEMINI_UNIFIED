# -*- coding: utf-8 -*-
#-------------------------------------------------------------
#   A template for the GFKASKI parameter file used to generate
#   the parameter file using a python notebook.
#-------------------------------------------------------------
def create_gfkaski_parfile(tlen,xlen,fmax,earth_model,dsvbasename,
                   source_type,attn_mode,eps,
                   slowness_limit_1,slowness_limit_2,
                   wavenumber_margin_fraction,
                   rearth,nblocks,nnod,dr,shift,min_source_node,rec_node_index):
    """creates string containing parfile for gfdsv-computations"""
    parfile_string = """#
#------------------------------------------------------------------------
#  GFKASKI parameter file
#------------------------------------------------------------------------------
#  Maximum length of synthetics to be calculated from FK-spectra in seconds
#  Determines the frequency stepping: df = 1/tlen
#
  TLEN = {0:f}
#------------------------------------------------------------------------
#  Distance from which wavenumber stepping is calculated in km
#  (should be about 10 times the maximum source receiver distance)
#
  XLEN = {1:f}
#------------------------------------------------------------------------
#  Upper frequency limit of FK-spectra in Hz 
#  (should be twice the Nyquist frequency desired for synthetics)
#
  FMAX = {2:f}
#------------------------------------------------------------------------
#  1D earth model in nm-format
#
  EARTH_MODEL = {3:s}
#------------------------------------------------------------------------
#  Basename of output files for GFK-spectra
#  Program appends '.sourcetype.crank' in parallel version
#
  DSVBASENAME = {4:s}
#------------------------------------------------------------------------
#  source type: set 'FORCE' for single force or 'MOMENT' for moment tensor source
#
  SOURCE_TYPE = {5:s}
#------------------------------------------------------------------------
#  Mode of attenuation 
#  allowed values are: ELASTIC,ATTENUATION_ONLY,DISPERSION_ONLY,ATTENUATION_AND_DISPERSION
#
  ATTENUATION_MODE = {6:s}
#------------------------------------------------------------------------
#  Accuracy of computation
#
  ACCURACY = {7:f}
#------------------------------------------------------------------------
#  Slowness limit for FK spectra used from 0 to 0.5*fmax in s/km
#  determines the maximum wavenumber according to kmax = om*plim
#
  SLOWNESS_LIMIT_1 = {8:f}
#------------------------------------------------------------------------
#  Slowness limit for FK-spectra used from 0.5*fmax to fmax
#
  SLOWNESS_LIMIT_2 = {9:f}
#------------------------------------------------------------------------
#  Wavenumber margin beyond slowness limit in fraction of max wavenumber
#
  WAVENUMBER_MARGIN_FRACTION = {10:f}
#------------------------------------------------------------------------
  EXTERNAL_NODES_REARTH = {11:f}
  EXTERNAL_NODES_NBLOCKS = {12:d}
  EXTERNAL_NODES_NNOD =  {13:s}
  EXTERNAL_NODES_DR =  {14:s}
#
#  Shift the nodes to greater depth by value EXTERNAL_NODES_SHIFT (is expected to be > 0.)
#  Value EXTERNAL_NODES_SHIFT is the depth of the top node.
#
  EXTERNAL_NODES_SHIFT = {15:f}
#
#  lower limit to source node range
#
  MIN_SOURCE_NODE = {16:d}
#
#  index of receiver node 
#
  REC_NODE_IDX = {17:d}
#--------------------------------------------------------------------------
""".format(tlen,xlen,fmax,earth_model,dsvbasename,
           source_type,attn_mode,
           eps,slowness_limit_1,slowness_limit_2,
           wavenumber_margin_fraction,
           rearth,nblocks,nnod,dr,shift,
           min_source_node,rec_node_index)
    return parfile_string
