"""This code calculates the number of events that passes the geometrical
cuts for the MuonE system for all the masses and coupling
parameters. For each parameter points, the expected number of events
are calculated and write into the .csv file.

The geometrical cuts are:

- The opening angle of the e+e- should be large enough for resolution.
Set as argument 1 of this code. Suggested value: 0.1 mrad, 0.5 mrad, 1 mrad.

- The decay volume is between the target and the first tracker.
That is, L_min = 1.6cm.
The L_max depends on what is the real distance of the tracker: 10cm or 15 cm.
This value is set as the second argument. Possible choices: 10 or 15. 0.5cm of
uncertainty is applied. That is, L_min = 9.5 or 14.5.

- Lowest energy cuts on X.

The file is written in .csv file. Path is argument 4.
"""
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
luminosity = 1.5e04
event_path = "/Users/isaac/Work/MUonE/MG5_events/2to3vec_v3/Events/"
m_e = 5.11e-04
meter = 1
GeV = 1 / (1.97 * meter * 1.0e-16)
Radius = 0.05

angle_cut = 1.0 / 1000.0  # input mrad, change to rad
L_min = 0.02
L_max = 14.5 / 100.0  # input cm, change to m
E_cut = 1.0  # energy cuts in unit of GeV
csv_file_path = sys.argv[1]
print("Calculating sensitivities!")


# -----------------------------------------------------------------
# Define useful functions
# -----------------------------------------------------------------


def decay_rate(mass, g_e):
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
    return momentum / (decay_rate(mass, g_e) * mass)


def probability(momentum, theta, mass, g_e):
    """This function calculates the probability for one
    single X particle to travel into the detectable region.
    open angle cuts is not applied here."""
    travel_distance = distance(momentum, mass, g_e)
    # Here loop over all the targets
    if L_max * np.tan(theta) <= Radius:
        L_prob = np.exp(-L_min / travel_distance) * (
            1.0 - np.exp((L_min - L_max) / travel_distance)
        )
    else:
        L_prob = 0.0
    return L_prob


def read_x_section(mass):
    """Read the cross section from madgraph banner file"""
    banner_file_subfolder = "X_mass_" + str(mass)
    banner_file_name = "X_mass_" + str(mass) + "_tag_1_banner.txt"
    banner_file_path = os.path.join(event_path, banner_file_subfolder, banner_file_name)
    with open(banner_file_path) as banner_file:
        for line in banner_file:
            if line.startswith("#  Integrated weight (pb)"):
                for x in line:
                    if x in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                        index = line.index(x)
                        break
                x_section_value = float(line[index:])
    return x_section_value


def get_sensitivity(mass, g_e, x_section):
    """Calculate the total number of signal events, given the probability for
    a single X particle to decay inside the volume.
    Open angle cuts is applied here."""
    rescaled_number = x_section * luminosity * (g_e ** 2)
    avg_probability = 0
    for i in range(len(X_list)):
        E = X_list[i][0]
        p = X_list[i][1]
        theta = X_list[i][2]
        open_angle = 2.0 * mass / E
        if open_angle >= angle_cut and E >= E_cut:
            avg_probability += probability(p, theta, mass, g_e)
    sensitivity = (avg_probability / 10000.0) * rescaled_number
    return sensitivity


# -----------------------------------------------------------------
# Initiate csv file writing, then loop over parameters
# -----------------------------------------------------------------
mass_list = np.logspace(np.log10(0.002), np.log10(0.29), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]
coupling_list = [10.0 ** power for power in np.arange(-8.0, 0.0, 0.2)]
csv_file = open(csv_file_path, "w")
with csv_file:
    header = ["Mass", "Coupling", "Number of Events"]
    data_writer = csv.writer(csv_file)
    data_writer.writerow(header)
    for mass in mass_list:
        # read cross section and rootfile
        x_section_mass = read_x_section(mass)
        chain = r.TChain("LHEF")
        chain.Reset()
        root_file_subfolder = "X_mass_" + str(mass)
        root_file_name = "unweighted_events.root"
        root_file_path = os.path.join(event_path, root_file_subfolder, root_file_name)
        chain.Add(root_file_path)
        print("search rootfile " + root_file_path)
        if float(chain.GetEntries()) == 10000.0:
            print(
                "Read rootfile for mass "
                + str(mass)
                + " successful! "
                + "Cross section: "
                + str(x_section_mass)
            )
        elif float(chain.GetEntries()) > 10000.0:
            sys.exit("Error occurred: rootfile combined by mistake!")
        else:
            sys.exit("Error occurred: did not find rootfile!")

        # Extract momentum and theta for all X particles
        X_list = []
        for event in chain:
            for particle in event.Particle:
                if particle.PID == 103:
                    E_X = particle.E
                    p_X = np.sqrt(particle.E ** 2 - particle.M ** 2)
                    X_vector = r.TLorentzVector()
                    X_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M
                    )
                    theta_X = X_vector.Theta()
                    X_list.append([E_X, p_X, theta_X])

        for coupling in coupling_list:
            sensitivity = get_sensitivity(mass, g_e=coupling, x_section=x_section_mass)
            print(
                "Coupling: " + str(coupling) + " number of events: " + str(sensitivity)
            )
            data_writer.writerow([mass, coupling, sensitivity])

print("Done!")
