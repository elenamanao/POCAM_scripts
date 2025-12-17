#!/usr/bin/env python3

import argparse
import numpy as np
from os.path import expandvars

parser = argparse.ArgumentParser(description="Detector simulator")
parser.add_argument("-i", "--infile", nargs="+", default=["saturation1.i3.bz2"],
                    dest="INFILE", help="Input file(s) (.i3{.gz} format)")
parser.add_argument("-o", "--outfile", default="test_detsim.i3",
                    dest="OUTFILE", help="Output file (.i3{.gz} format)")
parser.add_argument("-g", "--gcd", default=expandvars("$I3_DATA/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz"),
                    dest="GCDFILE", help="GCD file (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", type=int, default=1,
                    dest="RUNNUMBER", help="Run number")
parser.add_argument("--subrunnumber", type=int, default=1,
                    dest="SUBRUNNUMBER", help="Subrun number")
parser.add_argument("-l", "--holeice", default="angsens/as.flasher_p1_0.25_p2_0",
                    dest="HOLEICE", help="Hole ice model")
parser.add_argument("-s", "--storeintermediate", action="store_true", default=False,
                    dest="STOREINTERMEDIATE", help="Store intermediate MCPEs")

options = parser.parse_args()

# parse cmd line args, bail out if anything is not understood

# (options,args) = parser.parse_args()
# if len(args) != 0:
#     errmsg = "Got undefined options:"
#     for a in args:
#             errmsg += a
#             errmsg += " "
#     parser.error(errmsg)

from icecube.icetray import I3Tray, I3Units
tray = I3Tray()

from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services, simclasses
from icecube import vuvuzela
from icecube import DOMLauncher, DomTools, trigger_sim
from icecube import icetray, DomTools, WaveCalibrator, topeventcleaning, tpx, wavedeform
from icecube.online_filterscripts.base_segments.superdst import Unify
tray.AddService("I3SPRNGRandomServiceFactory","random",
                Seed = options.RUNNUMBER,
                StreamNum = options.SUBRUNNUMBER,
                NStreams = 10000,
                )
randomService = phys_services.I3SPRNGRandomService(
    seed = options.RUNNUMBER,
    nstreams = 10000,
    streamnum = options.SUBRUNNUMBER,
    )

infiles = options.INFILE
filenamelist = [options.GCDFILE]
for file in infiles:
    filenamelist = np.append(filenamelist,file)
print("Using hole ice model: ", options.HOLEICE)
print("Using random number seed: ", options.RUNNUMBER)
print("Using random number stream: ", options.SUBRUNNUMBER)
print("Using input file: ", infiles)
print("Using GCD file: ", options.GCDFILE)
print("Storing intermediate MCPEs: ", options.STOREINTERMEDIATE)
tray.AddModule("I3Reader","reader",
                FilenameList = filenamelist)

mcpe_to_pmt = "MCPESeriesMap"
tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons,"makeHitsFromPhotons",
                PhotonSeriesName = "I3PhotonSeriesMap",
                MCPESeriesName = "MCPESeriesMap",
                RandomService = randomService,
                DOMOversizeFactor=1.,
                DOMEfficiency=1.,
                IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_bfr-v2"),
                GCDFile=options.GCDFILE,
                HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/ANGSENS/%s"%options.HOLEICE),
)

mcpe_with_noise = mcpe_to_pmt + "_withNoise"
tray.AddModule("Vuvuzela", "vuvuzela_noise" ,
    InputHitSeriesMapName  = mcpe_to_pmt,
    OutputHitSeriesMapName = mcpe_with_noise,
    StartWindow            = -11*I3Units.microsecond,
    EndWindow              = 11*I3Units.microsecond,
    # IceTop                 = False, # arg no longer exists
    # InIce                  = True, # arg no longer exists
    ScaleFactor            = 1.0,
    DeepCoreScaleFactor    = 1,
    DOMsToExclude          = [], # This will be cleaned later by DOM launch cleaner
    RandomService          = "I3RandomService",
    SimulateNewDOMs        = True,
    DisableLowDTCutoff     = True,
    UseIndividual          = True
)

tray.AddModule("PMTResponseSimulator","rosencrantz",
    Input=mcpe_with_noise,  
    Output="I3MCPulses",
    MergeHits=True,
    )

tray.AddModule("DOMLauncher", "guildenstern",
    Input= "I3MCPulses",
    Output="InIceRawData_unclean",
    UseTabulatedPT=True,
    )

tray.AddModule("I3DOMLaunchCleaning","launchcleaning",
    InIceInput="InIceRawData_unclean",
    InIceOutput="InIceRawData",
    FirstLaunchCleaning=False
    )

#If = lambda f: True
tray.AddModule("I3WaveCalibrator", "_wavecal",
    Launches="InIceRawData")#,
    #If=If)
tray.AddModule("I3PMTSaturationFlagger", "_saturationflagger",
    Output="SaturationWindows")#,
    #If=If)

tray.AddModule("I3Wavedeform", "_wavedeform_without_reduce",
    Output="InIcePulses",
    Reduce=False)#,
    #If = lambda frame: If(frame) and  wvdfrm_flag in frame and frame[wvdfrm_flag].value)


###### triggering 

dellist = ['I3TriggerHierarchy', 'TimeShift', 'CleanIceTopRawData', 
                         'I3MCPulsesParticleIDMap', 'InIceRawData_unclean',
                        'I3PhotonSeriesMap', mcpe_with_noise]
if not options.STOREINTERMEDIATE:
    dellist.extend([mcpe_to_pmt, "I3MCPulses"])
tray.AddModule('Delete', 'delete_triggerHierarchy',
                Keys = dellist)

tray.AddSegment(trigger_sim.TriggerSim, 'trig', 
                gcd_file = options.GCDFILE,
                time_shift_args = {"SkipKeys" : ["BundleGen"]},
                run_id=1)

tray.AddModule("I3Writer","writer",
                Filename = options.OUTFILE,
                Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
tray.Execute()
