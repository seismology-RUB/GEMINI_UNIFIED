from seismicEvent import SeismicEvent
from specfemBox import SpecfemBox

basefolder = "/Users/tommy/git/Gemini/shyam_test/comparisonGeminiHybrid/"
geminiFolder = basefolder + "gemini/"
specfemFolder = basefolder + "OUTPUT_FILES/"

gemini_synseis_file = "seismograms_at_stations.hdf"

# Only one station file is needed to extract the data from Gemini/Specfem output.
# Since some calculations require the station-coordinates in geographical coordinates,
# the gemini station file is used here.
stationFile = "p11.sta"


box = SpecfemBox(file = geminiFolder + "/specfem_boundary_data.hdf", readLatLonOnly=True)
with open (basefolder + "demoEvent.ev") as eventFile:
    for eventLine in eventFile:
        eventData = eventLine.split()
        if eventData[0] == "S" or eventData[0] == "C":
            continue
        else:
            print(box.centerLat, box.centerLon)
            event = SeismicEvent()
            event.readFromFile(eventData)
            print("Event located at: ", event.lat, event.lon, event.depth)
            event.createStreamData(box.centerLat, box.centerLon, geminiFolder + stationFile,
                                   dataFolder=specfemFolder,
                                   rotate=True, datatype="velocity", typ="SPECFEM")
            event.createStreamData(box.centerLat, box.centerLon, geminiFolder + stationFile,
                                   geminiFile=geminiFolder + gemini_synseis_file,
                                   dataFolder=geminiFolder, rotate=False, datatype="velocity",
                                   typ="GEMINI")





