import csv
import os
import random
import sys

import numpy as np
import ROOT as r

# -----------------------------------------------------------------
# Import ROOT libraries for MadGraph generated rootfiles.
# -----------------------------------------------------------------
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

# -----------------------------------------------------------------
# set up environment
# -----------------------------------------------------------------
os.environ["TERM"] = "linux"
random.seed(1)

# -----------------------------------------------------------------
# Basic setup
# -----------------------------------------------------------------
event_path = "/Users/isaac/Work/MUonE/MG5_events/2to3vec/Events/"


mass_list = np.logspace(np.log10(0.002), np.log10(0.29), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]

for mass in mass_list:
    csv_folder = "/Users/isaac/Work/MUonE/data/"
    electron_file = os.path.join(csv_folder, str(mass) + "mev_e.csv")
    muon_file = os.path.join(csv_folder, str(mass) + "mev_mu.csv")
    X_file = os.path.join(csv_folder, str(mass) + "mev_X.csv")
    chain = r.TChain("LHEF")
    chain.Reset()
    root_file_subfolder = "X_mass_" + str(mass)
    root_file_name = "unweighted_events.root"
    root_file_path = os.path.join(event_path, root_file_subfolder, root_file_name)
    chain.Add(root_file_path)
    e_list = []
    mu_list = []
    x_list = []
    for event in chain:
        for particle in event.Particle:
            if particle.PID == 103:
                x_list.append([particle.E, particle.PT, particle.Phi, particle.Eta])
            elif particle.PID == 13 and particle.Status == 1:
                mu_list.append([particle.E, particle.PT, particle.Phi, particle.Eta])
            elif particle.PID == 11 and particle.Status == 1:
                e_list.append([particle.E, particle.PT, particle.Phi, particle.Eta])

    with open(electron_file, "w") as f:
        header = ["E", "PT", "Phi", "Eta"]
        data_writer = csv.writer(f)
        data_writer.writerow(header)
        data_writer.writerows(e_list)

    with open(muon_file, "w") as f:
        header = ["E", "PT", "Phi", "Eta"]
        data_writer = csv.writer(f)
        data_writer.writerow(header)
        data_writer.writerows(mu_list)

    with open(X_file, "w") as f:
        header = ["E", "PT", "Phi", "Eta"]
        data_writer = csv.writer(f)
        data_writer.writerow(header)
        data_writer.writerows(x_list)
