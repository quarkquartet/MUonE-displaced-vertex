"""Compare the background process 'mu- e- > mu- e- gamma' to
the signal process 'mu- e- > mu- e- sc'.
Basic cuts: cuts on pt4 and pta. Different value of this cuts are applied
and looped over.
Further cuts: require theta_e smaller than 30 mrad.
This code does the following to the events that pass the cuts:
1. Plot the events in the theta_e theta_mu plane.
2. Count the number of events that above the elastic curve, rescale
it according to the cross section and the luminosity to get the
real world number of events.
3. Calculate the sensitivities for the signal events vs the background
events."""

import os
import random

import numpy as np
import ROOT as r
from matplotlib import pyplot as plt

# ---------------------------------------------------------------------------
# Import ROOT libraries for MadGraph generated rootfiles.
# ---------------------------------------------------------------------------
r.gROOT.ProcessLine('.include /Users/isaac/Work/MG5_aMC_v2_6_7')
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/Delphes-3.4.2/external/ExRootAnalysis/ExRootTreeReader.h"')
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v2_6_7/ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"'
)
r.gInterpreter.Declare(
    '#include "/Users/isaac/Work/MG5_aMC_v2_6_7/ExRootAnalysis/ExRootAnalysis/ExRootLHEFReader.h"'
)
r.gSystem.Load(
    "/Users/isaac/Work/MG5_aMC_v2_6_7/ExRootAnalysis/libExRootAnalysis.so")

# ---------------------------------------------------------------------------
# Set up environment
# ---------------------------------------------------------------------------
os.environ['TERM'] = 'linux'
random.seed(1)

# ---------------------------------------------------------------------------
# Basic setups and load files
# ---------------------------------------------------------------------------
ptcuts = [0.001, 0.01, 0.05]  # Input pt cuts
luminosity = 1.5e+04  # Luminosity, which is 1.5 \times 10^4 pb^-1
signal_x_section = 248180.0  # Signal event cross section

# Create TChain to read the rootfile. We only use one signal rootfile, which
# is the 100 MeV one.
signal_chain = r.TChain("LHEF")
bg_chain = r.TChain("LHEF")
signal_chain.Add(
    "~/Work/MuonE/2to3sc/Events/nocut_0.1GeV/unweighted_events.root")
event_path = os.path.expanduser('~/Work/MuonE/SMmuegamma/Events/')
efficiency_file_path = "/Users/isaac/Work/MuonE/sensitivity.dat"
efficiency_file = open(efficiency_file_path, "w")
efficiency_file.write("pt4 cut    ")
efficiency_file.write("pta cut    ")
efficiency_file.write("signal number  ")
efficiency_file.write("signal efficiency  ")
efficiency_file.write("signal number rescaled  ")
efficiency_file.write("background number  ")
efficiency_file.write("background efficiency  ")
efficiency_file.write("background cross section  ")
efficiency_file.write("background number rescaled   ")
efficiency_file.write("sensitivity\n")

# ---------------------------------------------------------------------------
# Define necessary functions
# ---------------------------------------------------------------------------


def x_section(pt4cut, ptacut):
    """Read the background cross section from the madgraph"""
    banner_file_subfolder = 'pt4_' + str(pt4cut) + '_pta_' + str(ptacut)
    banner_file_name = 'pt4_' + \
        str(pt4cut) + '_pta_' + str(ptacut) + '_tag_1_banner.txt'
    banner_file_path = os.path.join(
        event_path, banner_file_subfolder, banner_file_name)
    with open(banner_file_path) as banner_file:
        for line in banner_file:
            if line.startswith('#  Integrated weight (pb)'):
                for x in line:
                    if x in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                        index = line.index(x)
                        break
                x_section_value = float(line[index:])
    return x_section_value


def T(theta_e):
    """T: defined as $E_4 - m_e$."""
    t = 2 * 0.000511 / (((150 + 0.000511)**2) /
                        ((np.cos(theta_e / 1000)**2) * (150**2 - 0.1057**2))
                        - 1)
    return t


def theta_mu_theory(theta_e):
    """The expression for the theta_mu vs theta_e relation in SM 2 to 2
    QED process"""
    theta = 1000 * np.arctan((np.sqrt(T(theta_e) * (2 * 0.000511 + T(theta_e)))
                              * np.sin(theta_e / 1000)) /
                             (np.sqrt(150**2 - 0.1057**2)
                              - np.sqrt(T(theta_e)
                                        * (2 * 0.000511 + T(theta_e)))
                              * np.cos(theta_e / 1000)))
    return theta


# Generate variable and prepare for the plots.
theta_mu2 = np.vectorize(theta_mu_theory)
theta_e_variable = np.linspace(0, 30, 500)
mu_variable = theta_mu2(theta_e_variable)


# ---------------------------------------------------------------------------
# Create figure panel
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=[12.8, 9.6])
fig.suptitle(r'$\theta_\mu - \theta_e$')


# ---------------------------------------------------------------------------
# Loop over the cuts. First loop over the pt4 cuts. For each pt4 cuts, apply
# it to the signal events.
# ---------------------------------------------------------------------------
plot_index = 0  # Count the number of plot, start with 0
for pt4_cuts in ptcuts:
    signal_theta_e_list = []
    signal_theta_mu_list = []
    signal_above_curve_number = 0

    # For each pt4_cuts, do analysis for the signal process
    for event in signal_chain:
        for particle in event.Particle:
            # Pick out final state electron
            if(particle.Status == 1 and particle.PID == 11):
                electron_vector = r.TLorentzVector()  # electron 4-vector
                pt4 = particle.PT  # electron pt
                electron_vector.SetPtEtaPhiM(
                    particle.PT, particle.Eta, particle.Phi, particle.M)
                theta_e = electron_vector.Theta() * 1000  # theta_e
            # Pick out final state muon
            elif(particle.Status == 1 and particle.PID == 13):
                muon_vector = r.TLorentzVector()  # muon 4-vector
                muon_vector.SetPtEtaPhiM(
                    particle.PT, particle.Eta, particle.Phi, particle.M)
                theta_mu = muon_vector.Theta() * 1000  # theta_mu

        # Apply cuts on pt
        if(pt4 > pt4_cuts and theta_e <= 30):
            signal_theta_e_list.append(theta_e)
            signal_theta_mu_list.append(theta_mu)
            if(theta_mu > theta_mu_theory(theta_e)):
                signal_above_curve_number += 1
    # Calculate the signal efficiency and rescaled number of events
    signal_efficiency = float(signal_above_curve_number) / 10000.0
    signal_rescaled_number = luminosity * signal_x_section * signal_efficiency
    # Then loop over the gamma cuts. For each cuts, we have a background file.
    # Create one plot for each cut.
    for pta_cuts in ptcuts:
        bg_above_curve_number = 0  # counting background events above the curve
        # Create plots.
        plot_index += 1  # Count the number of plots
        subplot_title = "pt4=" + str(pt4_cuts) + ", pta=" + str(pta_cuts)
        ax = fig.add_subplot(len(ptcuts), len(
            ptcuts), plot_index)  # Create subplot
        ax.set_xlim(0, 30)  # set limit on the theta_e
        ax.set_title(subplot_title)

        # Add rootfile of background events.
        event_folder_name = "pt4_" + str(pt4_cuts) + "_pta_" + str(pta_cuts)
        bg_root_file = os.path.join(
            event_path,
            event_folder_name, "unweighted_events.root")
        bg_chain.Add(bg_root_file)
        bg_theta_e_list = []
        bg_theta_mu_list = []
        bg_x_section = x_section(pt4_cuts, pta_cuts)
        for event in bg_chain:
            for particle in event.Particle:
                if(particle.Status == 1 and particle.PID == 11):
                    electron_vector = r.TLorentzVector()
                    phi_e = particle.Phi * 1000
                    electron_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M)
                    theta_e = electron_vector.Theta() * 1000
                elif(particle.Status == 1 and particle.PID == 13):
                    muon_vector = r.TLorentzVector()
                    muon_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M)
                    phi_mu = particle.Phi * 1000
                    theta_mu = muon_vector.Theta() * 1000
            if(theta_e <= 30):
                bg_theta_e_list.append(theta_e)
                bg_theta_mu_list.append(theta_mu)
                if(theta_mu > theta_mu_theory(theta_e)):
                    bg_above_curve_number += 1
        # This step must be applied. Because we need a new TChain to read
        # the new rootfile.
        bg_chain.Reset()
        # Calculate the background efficiency
        bg_efficiency = float(bg_above_curve_number) / 10000.0
        bg_rescaled_number = luminosity * bg_x_section * bg_efficiency
        sensitivity = signal_rescaled_number / \
            (np.sqrt(bg_rescaled_number + signal_rescaled_number))
        bg_x_section = "{:.2e}".format(bg_x_section)
        signal_rescaled_number_format = "{:.2e}".format(signal_rescaled_number)
        bg_rescaled_number = "{:.2e}".format(bg_rescaled_number)
        efficiency_file.write(str(pt4_cuts) + '       ')
        efficiency_file.write(str(pta_cuts) + '          ')
        efficiency_file.write(str(signal_above_curve_number) + '             ')
        efficiency_file.write(str(signal_efficiency) + '                ')
        efficiency_file.write(signal_rescaled_number_format + '              ')
        efficiency_file.write(str(bg_above_curve_number) + '               ')
        efficiency_file.write(str(bg_efficiency) + '                ')
        efficiency_file.write(bg_x_section + '                ')
        efficiency_file.write(bg_rescaled_number + '               ')
        efficiency_file.write(str(sensitivity) + '\n')
        print("pt4 cut: " + str(pt4_cuts))
        print("pta cut: " + str(pta_cuts))
        print("signal efficiency: " + str(signal_efficiency))
        print("background efficiency: " + str(bg_efficiency))
        print("sensitivity: " + str(sensitivity) + '\n')

        # Plot the points.
        ax.scatter(signal_theta_e_list, signal_theta_mu_list, color="blue")
        ax.scatter(bg_theta_e_list, bg_theta_mu_list, color="orange")
        ax.plot(theta_e_variable, mu_variable, color="r")

# ---------------------------------------------------------------------------
# Finish
# ---------------------------------------------------------------------------
fig.subplots_adjust(hspace=0.3, wspace=0.3)
fig.savefig("compare_all_ptcuts", dpi=300)
efficiency_file.close()
