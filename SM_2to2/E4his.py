from __future__ import print_function
import sys
import os
import random

import ROOT as r
r.gROOT.ProcessLine('.include /Users/isaac/Work/MG5_aMC_v2_6_7')
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/Delphes-3.4.2/external/ExRootAnalysis/ExRootTreeReader.h"'
)
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v2_6_7/ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"')
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v2_6_7/ExRootAnalysis/ExRootAnalysis/ExRootLHEFReader.h"'
)
r.gSystem.Load(
    "~/Work/MG5_aMC_v2_6_7/ExRootAnalysis/libExRootAnalysis.so")


os.environ['TERM'] = 'linux'
random.seed(1)

rootfilepath = sys.argv[1]

newrootfilepath = sys.argv[2]
newrootfile = r.TFile(newrootfilepath, 'recreate')

chain = r.TChain("LHEF")
chain.Add(rootfilepath)

E4_Hisall = r.TH1F("E4 distribution", "E4 distribution", 45, 5, 50)

for event in chain:
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 11):
            E4_Hisall.Fill(particle.E)

newrootfile.Write()
newrootfile.Close()
