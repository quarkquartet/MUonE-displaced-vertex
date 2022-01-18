import csv
import os
import random
import sys
import numpy as np
import ROOT as r
r.gROOT.ProcessLine(".include /Users/isaac/Work/MG5_aMC_v3_1_0")
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v3_1_0/ExRootAnalysis/\
ExRootAnalysis/ExRootTreeReader.h"'
)
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v3_1_0/\
ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"'
)
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v3_1_0/ExRootAnalysis/\
ExRootAnalysis/ExRootLHEFReader.h"'
)
r.gSystem.Load(
    "/Users/isaac/Work/MG5_aMC_v3_1_0/\
ExRootAnalysis/libExRootAnalysis.so"
)
os.environ["TERM"] = "linux"
random.seed(1)

mass_list = [2, 20, 200]
coupling = 1.0e-04
event_path = "/Users/isaac/Work/MUonE/MG5_events/check_iftah/Events/"

for mass in mass_list:
    chain = r.TChain("LHEF")
    chain.Reset()
    root_file_subfolder = "X_" + str(mass) + "mev"
    root_file_name = "unweighted_events.root"
    root_file_path = os.path.join(event_path, root_file_subfolder, root_file_name)
    chain.Add(root_file_path)
    gamma_file_path = "./" + str(mass) + "gamma.csv"
    gamma_list = []
    for event in chain:
        for particle in event.Particle:
            if particle.PID == 103:
                E_X = particle.E
                gamma = 1000*E_X/mass
        gamma_list.append([gamma])
    with open(gamma_file_path, "w") as f:
        data_writer = csv.writer(f)
        data_writer.writerows(gamma_list)