#!/usr/bin/env python3
from optparse import OptionParser
from os.path import expandvars
import numpy as np
#
usage = "usage: %prog [options] inputfile"

#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    description="Generate POCAM flash simulation input."
)

parser.add_argument("-o", "--outfile",
                    default="test_flashes.i3",
                    help="Write output to OUTFILE (.i3{.gz} format)")

parser.add_argument("-s", "--seed",
                    type=int, default=12344,
                    help="Initial seed for the random number generator")

parser.add_argument("-g", "--gcd",
                    type=str,
                    default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
                    help="Read geometry from GCDFILE (.i3{.gz} format)")

parser.add_argument("-r", "--runnumber",
                    type=int, default=1,
                    help="The run number for this simulation")

parser.add_argument("-n", "--numevents",
                    type=int, default=100,
                    help="Number of events per run")

parser.add_argument("-p", "--position", type=str,
                    help="Flasher position in detector (x y z)")

parser.add_argument("-f", "--flashertype",
                    type=str,
                    default="POCAMKapu405nmIsotropic",
                    help="Flasher type (see I3FlasherPulse.h)")

parser.add_argument("--nphotons",
                    type=int, default=1_000_000_000,
                    help="Number of photons to generate per event")

parser.add_argument("--pulsewidth",
                    type=float, default=1.0,
                    help="Pulse width of a single flash")

# weighted/unweighted toggle — argparse handles this cleanly
weighted_group = parser.add_mutually_exclusive_group()
weighted_group.add_argument("--weighted",
                            action="store_true", default=True,
                            help="Use weighted photons (default)")
weighted_group.add_argument("--unweighted",
                            action="store_false",
                            dest="weighted",
                            help="Use unweighted photons")

args = parser.parse_args()

print(args)

from icecube.icetray import I3Tray, I3Units
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services, clsim, sim_services, simclasses

tray = I3Tray()
tray.AddModule("I3InfiniteSource","streams",
               Prefix=args.gcd,
               Stream=icetray.I3Frame.DAQ)
tray.AddModule("I3MCEventHeaderGenerator","gen_header",
               Year=2012,
               DAQTime=7968509615844458,
               RunNumber=1,
               EventID=1,
               IncrementEventID=True)

flasher_pos = dataclasses.I3Position(*(float(v) for v in args.position.split(",")))
print("Source position: ", flasher_pos)
flasher_type = simclasses.I3FlasherPulse.FlasherPulseType.names[args.flashertype]
print("Flasher type: ", flasher_type)
print("Number of flash events: ", args.numevents)
print("Number of photons per event: ", args.nphotons)
print("Weighted photons: ", args.weighted)
# a random number generator
try:
    randomService = phys_services.I3SPRNGRandomService(
        seed = args.seed,
        nstreams = 10000,
        streamnum = args.runnumber)
except AttributeError:
    randomService = phys_services.I3GSLRandomService(
        seed = args.seed*10000 + args.runnumber,
    )

from icecube.clsim import POCAMFlashInfoIsotropic
tray.AddModule(POCAMFlashInfoIsotropic, "POCAMFlasher", 
               FlasherPulseType=flasher_type,
               PhotonPulseSeriesName="I3FlasherPulseSeries",
               NumberOfPhotons=args.nphotons,
               FlasherPosition=flasher_pos,
               PulseWidth= args.pulsewidth*I3Units.ns, 
               NEvents=args.numevents
               )
print("Uniform AS")
# start photon propagation with CLSim
tray.AddSegment(clsim.I3CLSimMakePhotons, "goCLSIM",
    UseGPUs=True,
    UseCPUs=False,
    UseGeant4=False,
    UseI3PropagatorService=False,
    RandomService=randomService,
    DoNotParallelize=False,
    UnweightedPhotons=(not args.weighted),
    StopDetectedPhotons=True,
    #OMKeyMaskName=I3Vector,
    PhotonSeriesName='I3PhotonSeriesMap',
    MCPESeriesName='',
    FlasherPulseSeriesName="I3FlasherPulseSeries",
    GCDFile=args.gcd,
    IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/ICEMODEL/spice_ftp-v3"),
    # HoleIceParameterization=expandvars("$I3_BUILD/ice-models/resources/models/ANGSENS/angsens/as.uniform"),
    DOMOversizeFactor=1.,
    DOMEfficiency=1.5,
    IgnoreSubdetectors=['IceTop', 'NotOpticalSensor'],
    MCTreeName= None,
   )
def EmptyI3MCTree(frame):
    if not frame.Has('I3MCTree'):
        frame.Put('I3MCTree', dataclasses.I3MCTree())
    return
tray.AddModule(EmptyI3MCTree, "EmptyI3MCTree", streams=[icetray.I3Frame.DAQ])
tray.AddModule("I3Writer","writer",
    Filename = args.outfile,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo])

tray.Execute()
