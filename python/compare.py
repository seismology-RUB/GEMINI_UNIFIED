from obspy import read
import matplotlib.pyplot as plt
from scipy import interpolate
from obspy.signal.tf_misfit import plot_tf_misfits
from sklearn.metrics import mean_squared_error
import numpy as np

def compareGeminiUnifiedToSpecfemCoupling(geminiFile, specfemFile):
    geminiStream = read(geminiFile)
    specfemStream = read(specfemFile)
    component = ["Z", "N", "E"]
    #stations = ["ST001", "ST005", "ST010", "ST015", "ST020", "ST025", "ST030", "ST035", "ST040"]
    stations = ["ST001", "ST050", "ST100",  "ST150",  "ST200"]

    fig, ax = plt.subplots(1, 3, sharey=True, sharex=True, constrained_layout=True, figsize=(24, 6))
    for nc, comp in enumerate(component):
        #ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])
        rmsax = ax[nc].secondary_yaxis("right")
        rmsTickLabels = []
        yticklabels = []
        yticks = []
        i = 0
        gemMax = 0
        specmax = 0
        tbuf = 500.0
        for gemTrace, specfemTrace in zip(geminiStream, specfemStream):
            if gemTrace.stats.station in stations:
                if max(gemTrace.data) > gemMax:
                    gemMax = max(gemTrace.data)
                if max(specfemTrace.data) > specmax:
                    specmax = max(specfemTrace.data)

        offset = max(gemMax, specmax)

        for gemTrace, specfemTrace in zip(geminiStream, specfemStream):
            if gemTrace.stats.station in stations:
                if gemTrace.stats.channel == comp and specfemTrace.stats.channel == comp:
                    i += 0.5
                    yticklabels.append(gemTrace.stats.station)
                    yticks.append(i*offset)
                    # The 10 is added as the specfemseismogram starts at -10 seconds.
                    geminiTime = [t * gemTrace.stats.delta + tbuf + 10 for t in range(len(gemTrace.data))]
                    specfemTime = [t * specfemTrace.stats.delta + tbuf + 10 for t in range(len(specfemTrace.data))]
                    specfemData = specfemTrace.data

                    inter = interpolate.InterpolatedUnivariateSpline(specfemTime, specfemData)
                    #geminiTime = geminiTime[:int(len(gemTrace.data) / 4.4)]
                    specfemData = inter(geminiTime)
                    #geminiData = gemTrace.data[:int(len(gemTrace.data)/4.4)]*1e6

                    geminiData = gemTrace.data

                    #print(max(geminiData))
                    #geminiTime = geminiTime[150:449]
                    #geminiData = geminiData[150:449]
                    #specfemData = specfemData[150:449]

                    dt = specfemTrace.stats.delta
                    #plot_tf_misfits(specfemData, geminiData, dt=dt)

                    normGeminiData = geminiData/max(geminiData)
                    normSpecfemData = specfemData/max(specfemData)

                    rms = mean_squared_error(geminiData, specfemData, squared=False)
                    nrms = rms/(max(specfemData) - min(specfemData))
                    rmsTickLabels.append(format(round(nrms, 3), '.3f'))
                    diff = geminiData - specfemData
                    factor = 1.2
                    if i <= 0.5:
                        ax[nc].plot(geminiTime, geminiData*factor + (i*offset), c="red", lw=1, label="Gemini (1D)")
                        ax[nc].plot(geminiTime, specfemData*factor + (i*offset), c="blue", lw=1, label="Hybrid (3D)", ls="-.")
                        ax[nc].plot(geminiTime, diff*factor + (i*offset), c="black", lw=1, label="Residual", ls="--")
                    else:
                        ax[nc].plot(geminiTime, geminiData*factor + (i*offset), c="red", lw=1)
                        ax[nc].plot(geminiTime, specfemData*factor + (i*offset), c="blue", lw=1, ls="-.")
                        ax[nc].plot(geminiTime, diff*factor + (i*offset), c="black", lw=1, ls="--")
        time = geminiStream[0].stats.starttime.strftime("%d. %B %Y %H:%M:%S")
        network = geminiStream[0].stats.network
        ax[nc].set_yticks(yticks)
        ax[nc].set_yticklabels(yticklabels)

        ax[nc].set_xlabel("Time in seconds after " + time)
        rmsax.set_yticks(yticks)
        rmsax.set_yticklabels(rmsTickLabels)
        ax[nc].set_title(comp + "-Component")
        #ax[nc].annotate('Stations',
        #            xy=(.018, .955), xycoords='figure fraction',
        #            horizontalalignment='left', verticalalignment='top',
        #            fontsize=12)
        #ax[nc].annotate('RMSD',
        #            xy=(.32*(nc+1), .9), xycoords='figure fraction',
        #            horizontalalignment='left', verticalalignment='top',
        #            fontsize=12)
        leg = ax[nc].legend(loc="lower left")
        #leg2 = ax[1].legend(loc="lower left")

    ax[0].set_ylabel("Particle velocity")
    ax[0].annotate('Stations',
                   xy=(.025, .89), xycoords='figure fraction',
                   horizontalalignment='left', verticalalignment='top',
                   fontsize=12)
    ax[0].annotate('NRMSD',
                    xy=(.334, .89), xycoords='figure fraction',
                    horizontalalignment='left', verticalalignment='top',
                    fontsize=12)
    ax[1].annotate('NRMSD',
                  xy=(.636, .89), xycoords='figure fraction',
                  horizontalalignment='left', verticalalignment='top',
                  fontsize=12)
    ax[2].annotate('NRMSD',
                  xy=(.955, .89), xycoords='figure fraction',
                  horizontalalignment='left', verticalalignment='top',
                  fontsize=12)

    ax[1].annotate("Profile " + network + " - Application",
                   xy=(.505, .985), xycoords='figure fraction',
                   horizontalalignment='center', verticalalignment='top',
                   fontsize=18)
    plt.show()
    #fig.savefig("/Users/tommy/Desktop/plot_comparison_P10_ak135.png", dpi=300, bbox_inches="tight")


#geminiFile = "/Users/tommy/git/Gemini/syntheticSeismograms/onEquator-prem-no-ocean-0.2Hz/P10/AK135/geminiData171117_223425.dat"
#specfemFile = "/Users/tommy/git/Gemini/specfemCoupling/box-80-60-30-prem-no-ocean-0.2Hz/P10/AK135/specfemData171117_223425.dat"


#specfemFile = "/Users/tommy/git/Gemini/shyam_test/comparisonGeminiHybrid/OUTPUT_FILES/specfemData171117_223425.dat"

specfemFile = "/Users/tommy/git/Gemini/shyam_test/comparisonGeminiHybrid/gemini/horizontal-profile_side_right_radius_0_velocity_si.dat"
geminiFile = "/Users/tommy/git/Gemini/shyam_test/comparisonGeminiHybrid/gemini/geminiData171117_223425.dat"
compareGeminiUnifiedToSpecfemCoupling(geminiFile, specfemFile)