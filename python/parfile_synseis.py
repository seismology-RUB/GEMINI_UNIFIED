# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#   A template for the computeSyntheticSeismograms parameter file used to generate
#   the parameter file using a python notebook.
#-------------------------------------------------------------
def create_synseis_parfile(basename,path_synthetic_seis,file_event_list,
                       file_station_list,path_event_filter,path_station_filter,
                       components,event_filter_type,event_filter_specs,
                       station_filter_type):
    """creates string containing parfile for gfdsvseis-computations"""
    parfile_string = """#
#------------------------------------------------------------------------
#  computeSyntheticSeismograms parameter file:
#------------------------------------------------------------------------
#  Base name for output file of computeGreenFKSpectraForSynthetics 
#  with Green functions in the frequency-wavenumber domain. 
#  Provide full path and first part
#  of file name including source type specification which is either
#  "FORCE" for a single force or "MOMENT" for a moment tensor source.
#  Further extensions such as ".000" are automatically provided by
#  computeSyntheticSeismograms.
#
  GFDSV_BASENAME = {0:s}
#------------------------------------------------------------------------
#  Path to folder where synthetic seismograms are stored.
#
  PATH_SYNTHETIC_SEIS =  {1:s}
#------------------------------------------------------------------------
#  Name of event list file including full path
#
  FILE_EVENT_LIST =  {2:s}
#------------------------------------------------------------------------
#  Name of station list file including full path
#
  FILE_STATION_LIST =  {3:s}
#------------------------------------------------------------------------
#  Path to folder where event filters are stored.
#  Filter file names follow the convention "filter_eventid".
#
  PATH_EVENT_FILTER =  {4:s}
#-------------------------------------------------------------------------
#  Path to folder where station filters are stored.
#  Filter files follow the convention "filter_staname_comp".
#
  PATH_STATION_FILTER =  {5:s}
#-------------------------------------------------------------------------
#  String containing sequence of components to be calculated.
#  Example: ZNEHRS
#
  COMPONENTS =  {6:s}
#-------------------------------------------------------------------------
#  Information specifying the way the filter is computed.
#  "EVENT_FILTER_TYPE" contains a keyword for the algorithm.
#  "EVENT_FILTER_SPECS" contains a real array with necessary data.
#  Here, it is order and corner frequency (Hz) for highpass and lowpass,
#  respectively.
#
  EVENT_FILTER_TYPE =  {7:s}
  EVENT_FILTER_SPECS =  {8:s}
#-------------------------------------------------------------------------
#  Information specifying the way station filters are computed.
#  "STATION_FILTER_TYPE" contains a keyword for the algorithm.
#  More should be implemented later
#
  STATION_FILTER_TYPE =  {9:s}
#-------------------------------------------------------------------------
""".format(basename,path_synthetic_seis,file_event_list,
           file_station_list,path_event_filter,path_station_filter,
           components,event_filter_type,event_filter_specs,
           station_filter_type)
    return parfile_string
