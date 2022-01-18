"""This code calculate how much decay product can finally enters the required augular
acceptance. The angular acceptance depends on what we want: require the electron-positron pair to enter the 3 trackers, or enter the final ECAL."""

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

lowest_angle_cut = (
    1.0 / 1000.0
)  # This requires that the opening angle should be within the resolution
L_min = 0.02
L_max = 14.5 / 100.0  # input cm, change to m
E_cut = 1.0  # energy cuts in unit of GeV

tracker_csv_path = "./tracker.csv"
ECAL_csv_path = "./ECAL.csv"

FLAG = sys.argv[1]


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


def boost(theta, mass, energy):
    beta = np.sqrt(1.0 - mass / energy)
    gamma = energy / mass
    pcom = np.sqrt(mass ** 2.0 - 4.0 * (m_e ** 2.0)) / 2.0
    nominator = gamma * (pcom * np.cos(theta) + beta * mass / 2.0)
    denominator = np.sqrt(
        (pcom ** 2.0) * (np.sin(theta) ** 2.0)
        + (gamma ** 2.0) * ((pcom * np.cos(theta) + beta * mass * 0.5) ** 2.0)
    )
    return nominator / denominator


mass_list = np.logspace(np.log10(0.002), np.log10(0.29), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]
coupling = 1.0e-04

if FLAG == "tracker":
    # Require pair to enter trackers
    tracker_file = open(tracker_csv_path, "w")
    with tracker_file:
        header = ["mass", "efficiency rate"]
        data_writer = csv.writer(tracker_file)
        data_writer.writerow(header)
        for mass in mass_list:
            chain = r.TChain("LHEF")
            chain.Reset()
            root_file_subfolder = "X_mass_" + str(mass)
            root_file_name = "unweighted_events.root"
            root_file_path = os.path.join(
                event_path, root_file_subfolder, root_file_name
            )
            chain.Add(root_file_path)
            print("search rootfile " + root_file_path)
            if float(chain.GetEntries()) == 100000.0:
                print("Read rootfile for mass " + str(mass) + " successful! ")
            elif float(chain.GetEntries()) > 100000.0:
                sys.exit("Error occurred: rootfile combined by mistake!")
            else:
                sys.exit("Error occurred: did not find rootfile!")

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

            good_photon_number = 0
            good_decay_number = 0
            for i in range(0, len(X_list)):
                E_X = X_list[i][0]
                p_X = X_list[i][1]
                theta_X = X_list[i][2]
                d_A = np.random.exponential(scale=distance(p_X, mass, coupling))
                if L_min <= d_A <= L_max:
                    good_photon_number += 1
                    theta_e = np.arccos(np.random.uniform(-1.0, 1.0))
                    phi_e = np.random.uniform(-np.pi * 0.5, np.pi * 0.5)
                    theta_elab1 = np.arccos(boost(theta_e, mass, E_X))
                    theta_elab2 = np.arccos(boost(theta_e + np.pi, mass, E_X))
                    AB = 1.0 - d_A
                    BC = AB * theta_elab1
                    BD = 1.0 * theta_X
                    BF = AB * theta_elab2
                    if BD <= Radius:
                        sin_alpha = np.sin(phi_e) * BD / Radius
                        alpha = np.arcsin(sin_alpha)
                        BCmax = Radius * (np.sin(phi_e - alpha) / np.sin(phi_e))
                        BFmax = Radius * (np.sin(phi_e + alpha) / np.sin(phi_e))
                        if BC <= BCmax and BF <= BFmax:
                            good_decay_number += 1
            if good_photon_number > 0:
                efficiency_rate = float(good_decay_number) / float(good_photon_number)
                print(
                    "mass = " + str(mass) + ", efficiency rate: " + str(efficiency_rate)
                )
                data_writer.writerow([mass, efficiency_rate])
            else:
                print("mass = " + str(mass) + ", no dark photon decays in volume.")

elif FLAG == "ECAL" or "ecal":
    # Require pair to enter ECAL
    ECAL_file = open(ECAL_csv_path, "w")
    with ECAL_file:
        header = [
            "mass",
            "small angle photon ratio",
            "small open angle rate",
            "efficiency rate",
        ]
        data_writer = csv.writer(ECAL_file)
        data_writer.writerow(header)
        for mass in mass_list:
            chain = r.TChain("LHEF")
            chain.Reset()
            root_file_subfolder = "X_mass_" + str(mass)
            root_file_name = "unweighted_events.root"
            root_file_path = os.path.join(
                event_path, root_file_subfolder, root_file_name
            )
            chain.Add(root_file_path)
            print("search rootfile " + root_file_path)
            if float(chain.GetEntries()) == 100000.0:
                print("Read rootfile for mass " + str(mass) + " successful! ")
            elif float(chain.GetEntries()) > 100000.0:
                sys.exit("Error occurred: rootfile combined by mistake!")
            else:
                sys.exit("Error occurred: did not find rootfile!")

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

            good_photon_number = 0
            small_angle_photon_number = 0
            good_decay_number = 0
            for j in range(0, 40):
                axis = (
                    j + 1
                )  # The distance from the j-th target to the ECAL. j is counted from the one that is closest to the ECAL
                for i in range(0, len(X_list)):
                    E_X = X_list[i][0]
                    p_X = X_list[i][1]
                    theta_X = X_list[i][2]
                    d_A = np.random.exponential(scale=distance(p_X, mass, 1.0e-04))
                    if L_min <= d_A <= L_max:
                        good_photon_number += 1
                        theta_e = np.arccos(np.random.uniform(-1.0, 1.0))
                        phi_e = np.random.uniform(-np.pi * 0.5, np.pi * 0.5)
                        theta_elab1 = np.arccos(boost(theta_e, mass, E_X))
                        theta_elab2 = np.arccos(boost(theta_e + np.pi, mass, E_X))
                        AB = axis - d_A
                        BC = AB * theta_elab1
                        BD = axis * theta_X
                        BF = AB * theta_elab2
                        if BD <= Radius:
                            small_angle_photon_number += 1
                            sin_alpha = np.sin(phi_e) * BD / Radius
                            alpha = np.arcsin(sin_alpha)
                            BCmax = Radius * (np.sin(phi_e - alpha) / np.sin(phi_e))
                            BFmax = Radius * (np.sin(phi_e + alpha) / np.sin(phi_e))
                            if BC <= BCmax and BF <= BFmax:
                                good_decay_number += 1
            good_photon_number = good_photon_number / 40
            small_angle_photon_number = small_angle_photon_number / 40
            good_decay_number = good_decay_number / 40
            if good_photon_number > 0:
                small_angle_photon_rate = float(small_angle_photon_number) / float(
                    good_photon_number
                )
                small_open_angle_rate = float(good_decay_number) / float(
                    small_angle_photon_number
                )
                efficiency_rate = float(good_decay_number) / float(good_photon_number)
                data_writer.writerow(
                    [
                        mass,
                        small_angle_photon_rate,
                        small_open_angle_rate,
                        efficiency_rate,
                    ]
                )
                print(
                    "mass = " + str(mass) + ", efficiency rate: " + str(efficiency_rate)
                )
            else:
                print("mass = " + str(mass) + ", no dark photon decays in volume.")

else:
    print("Please input whether require the decay products to enter tracker or ECAL.")
