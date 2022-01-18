"""This code calculates the number of events that passes the geometrical cuts of the MUonE system for all the masses and couplings.

The angular acceptance requires:
- Case 1: require the decay products to pass all three trackers in the corresponding module.
- Case 2: require the all the decay products to enter the final ECAL.

Other geometrical cuts are:
- The opening angle of the e+e- should be large enough for resolution.
Set as argument 1 of this code. Suggested value: 1 mrad.

- The decay volume is between the target and the first tracker.
That is, L_min = 2.0 cm.
The L_max depends on what is the real distance of the tracker: 15 cm.
This value is set as the second argument. 0.5cm of
uncertainty is applied. That is, L_max = 14.5.

- Lowest energy cuts on X: this offers enough resolution for that.

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
event_path = "/Users/isaac/Work/MUonE/MG5_events/20211215/Events/"
m_e = 5.11e-04
meter = 1
GeV = 1 / (1.97 * meter * 1.0e-16)
Radius = 0.05

angle_cut = 1.0 / 1000.0  # input mrad, change to rad
L_min = 0.02
L_max = 14.5 / 100.0  # input cm, change to m
E_cut = round(float(sys.argv[1]), 1)  # energy cuts in unit of GeV
nlayers = int(sys.argv[2])
ECAL_result_path = "./reach_tight_" + sys.argv[1] + "GeV_" + sys.argv[2] + "layers.csv"

print("Calculating sensitivities, file are written in " + ECAL_result_path)
print("Energy cut: " + str(E_cut) + " GeV")
print("Number of layers: " + str(nlayers))


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


def boost(theta, mass, energy):
    beta = np.sqrt(energy ** 2.0 - mass ** 2.0) / energy
    gamma = energy / mass
    pcom = np.sqrt(mass ** 2.0 - 4.0 * (m_e ** 2.0)) / 2.0
    nominator = gamma * (pcom * np.cos(theta) + beta * mass / 2.0)
    denominator = np.sqrt(
        (pcom ** 2.0) * (np.sin(theta) ** 2.0)
        + (gamma ** 2.0) * ((pcom * np.cos(theta) + beta * mass * 0.5) ** 2.0)
    )
    return nominator / denominator


def energy_boost(theta, mass, energy):
    beta = np.sqrt(energy ** 2.0 - mass ** 2.0) / energy
    gamma = energy / mass
    pcom = np.sqrt(mass ** 2.0 - 4.0 * (m_e ** 2.0)) / 2.0
    return gamma * (0.5 * mass + beta * pcom * np.cos(theta))


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


mass_list = np.logspace(np.log10(0.002), np.log10(0.2), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]
coupling_list = np.logspace(-6.0, -2.0, 30).tolist()


with open(ECAL_result_path, "w") as f:
    header = ["mass", "coupling", "number of events"]
    data_writer = csv.writer(f)
    data_writer.writerow(header)
    for mass in mass_list:
        x_section_mass = read_x_section(mass)
        chain = r.TChain("LHEF")
        chain.Reset()
        root_file_subfolder = "X_mass_" + str(mass)
        root_file_name = "unweighted_events.root"
        root_file_path = os.path.join(event_path, root_file_subfolder, root_file_name)
        chain.Add(root_file_path)
        print("search rootfile " + root_file_path)
        if float(chain.GetEntries()) == 100000.0:
            print(
                "Read rootfile for mass "
                + str(mass)
                + " successful! "
                + "Cross section: "
                + str(x_section_mass)
            )
        elif float(chain.GetEntries()) > 100000.0:
            sys.exit("Error occurred: rootfile combined by mistake!")
        else:
            sys.exit("Error occurred: did not find rootfile!")

        event_list = []
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
                elif particle.PID == 13 and particle.Status == 1:
                    mu_vector = r.TLorentzVector()
                    mu_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M
                    )
                    theta_mu = mu_vector.Theta()
                    e_mu = particle.E
                elif particle.PID == 11 and particle.Status == 1:
                    e_vector = r.TLorentzVector()
                    e_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M
                    )
                    theta_e = e_vector.Theta()
                    e_e = particle.E
            event_list.append([E_X, p_X, theta_X, theta_mu, theta_e, e_mu, e_e])

        for coupling in coupling_list:
            event_number = 0
            for j in range(0, nlayers):
                axis = j + 1
                for i in range(0, len(event_list)):
                    E_X = event_list[i][0]
                    p_X = event_list[i][1]
                    theta_X = event_list[i][2]
                    theta_mu = event_list[i][3]
                    theta_e = event_list[i][4]
                    e_mu = event_list[i][5]
                    e_e = event_list[i][6]
                    d_A = np.random.exponential(scale=distance(p_X, mass, coupling))
                    if (
                        axis * theta_X <= Radius
                        and L_min <= d_A <= L_max
                        and e_mu >= E_cut
                        and e_e >= E_cut
                        and axis * theta_e <= Radius
                        and axis * theta_mu <= Radius
                    ):
                        theta_e = np.arccos(np.random.uniform(-1.0, 1.0))
                        phi_e = np.random.uniform(-np.pi * 0.5, np.pi * 0.5)
                        theta_elab1 = np.arccos(boost(theta_e, mass, E_X))
                        theta_elab2 = np.arccos(boost(theta_e + np.pi, mass, E_X))
                        open_angle = theta_elab1 + theta_elab2
                        energy_elab1 = energy_boost(theta_e, mass, E_X)
                        energy_elab2 = energy_boost(theta_e + np.pi, mass, E_X)
                        if (
                            open_angle >= angle_cut
                            and energy_elab1 >= E_cut
                            and energy_elab2 >= E_cut
                        ):
                            AB = axis - d_A
                            BC = AB * theta_elab1
                            BD = axis * theta_X
                            BF = AB * theta_elab2
                            sin_alpha = np.sin(phi_e) * BD / Radius
                            alpha = np.arcsin(sin_alpha)
                            BCmax = Radius * (np.sin(phi_e - alpha) / np.sin(phi_e))
                            BFmax = Radius * (np.sin(phi_e + alpha) / np.sin(phi_e))
                            if BC <= BCmax and BF <= BFmax:
                                event_number += 1
            event_number = event_number / float(nlayers)
            sensitivity = (
                x_section_mass
                * luminosity
                * (float(nlayers) / 40.0)
                * (coupling ** 2)
                * (event_number / 100000.0)
            )
            data_writer.writerow([mass, coupling, sensitivity])
            print(
                "mass: "
                + str(mass)
                + ", coupling: "
                + str(coupling)
                + ", sensitivity: "
                + str(sensitivity)
            )

print("Done!")
