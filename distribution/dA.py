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

luminosity = 1.5e04
event_path = "/Users/isaac/Work/MUonE/MG5_events/2to3vec/Events/"
m_e = 5.11e-04
meter = 1
GeV = 1 / (1.97 * meter * 1.0e-16)
Radius = 0.05

lowest_angle_cut = 1.0 / 1000.0  # This requires that the opening angle should be within the resolution
L_min = 0.02
L_max = 14.5 / 100.0  # input cm, change to m
E_cut = 1.0  # energy cuts in unit of GeV

def decay_rate(mass, g_e):
    # the \Gamma, eq 10
    alpha_e = (g_e ** 2) / (4 * np.pi)
    return (
        alpha_e
        * mass
        * (1 + 2 * ((m_e ** 2) / (mass ** 2)))
        * np.sqrt(1 - 4 * ((m_e ** 2) / (mass ** 2)))
        * GeV
        / 3
    )


def distance(momentum, mass, g_e):
    # the d_{A^\prime}, eq 11
    return momentum / (decay_rate(mass, g_e) * mass)

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
    dApath = "./" + str(mass) + "dA.csv"
    X_list = []
    for event in chain:
        for particle in event.Particle:
            if particle.PID == 103:
                E_X = particle.E
                p_X = np.sqrt(particle.E**2 - particle.M**2)
                X_vector = r.TLorentzVector()
                X_vector.SetPtEtaPhiM(particle.PT, particle.Eta, particle.Phi, particle.M)
                theta_X = X_vector.Theta()
        X_list.append([E_X, p_X, theta_X])

    dAlist = []
    for j in range(0, 40):
        axis = j+1
        for i in range(0, len(X_list)):
            E_X = X_list[i][0]
            p_X = X_list[i][1]
            theta_X = X_list[i][2]
            d_A = np.random.exponential(scale=distance(p_X, mass/1000.0, 1.0e-04))
            if (L_min <= d_A <= L_max and axis*theta_X<=Radius):
                BD = axis*theta_X
                dAlist.append([BD])

    with open(dApath, "w") as f:
        data_writer = csv.writer(f)
        data_writer.writerows(dAlist)
    
