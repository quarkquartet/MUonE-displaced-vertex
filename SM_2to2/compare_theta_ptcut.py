"""This code compares the signal events 'mu- e- > mu- e- sc' with the
background events 'mu- e- > mu- e- gamma' in the theta_e vs theta_mu plane.
The efficiencies that events who is above the elastic band is calculated
for both the signal events and background events.
Different cuts are applied. The background events receive cuts during event
generation.  The signal events receive cuts in the analysis.
Different cuts on pt_electron and pt_gamma are applied and looped over."""
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
    "~/Work/MG5_aMC_v2_6_7/ExRootAnalysis/libExRootAnalysis.so")

# ---------------------------------------------------------------------------
# Set up environment
# ---------------------------------------------------------------------------
os.environ['TERM'] = 'linux'
random.seed(1)


# ---------------------------------------------------------------------------
# Basic setups and load files
# ---------------------------------------------------------------------------
ptcuts = [0.001, 0.01, 0.05]  # Input pt cuts

# Create TChain to read the rootfile. We only use one signal rootfile, which
# is the 100 MeV one.
signal_chain = r.TChain("LHEF")
bg_chain = r.TChain("LHEF")
signal_chain.Add(
    "~/Work/MuonE/2to3sc/Events/nocut_0.1GeV/unweighted_events.root")
efficiency_file_path = "/Users/isaac/Work/MuonE/ptcut_efficiency.dat"
efficiency_file = open(efficiency_file_path, "w")
efficiency_file.write("pt4 cut    ")
efficiency_file.write("pta cut    ")
efficiency_file.write("signal_pass_cut  ")
efficiency_file.write("signal_above_curve  ")
efficiency_file.write("background_pass_cut  ")
efficiency_file.write("background_above_curve\n")


# ---------------------------------------------------------------------------
# Define function to theoretically calculate the elastic band (theta_e and
# theta_mu relation in the SM process)
# ---------------------------------------------------------------------------

# T: defined as $E_4 - m_e$.


def T(theta_e):
    t = 2 * 0.000511 / (((150 + 0.000511)**2) /
                        ((np.cos(theta_e / 1000)**2) * (150**2 - 0.1057**2))
                        - 1)
    return t


def theta_mu_theory(theta_e):
    theta = 1000 * np.arctan((np.sqrt(T(theta_e) * (2 * 0.000511 + T(theta_e)))
                              * np.sin(theta_e / 1000)) /
                             (np.sqrt(150**2 - 0.1057**2)
                              - np.sqrt(T(theta_e) * (2 * 0.000511 + T(theta_e))) * np.cos(theta_e / 1000)))
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

    signal_pass_cut_number = 0  # counting number of signal events that pass the cut
    # counting number of singal events that above the event
    signal_above_curve_number = 0
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
            signal_pass_cut_number += 1
            if(theta_mu > theta_mu_theory(theta_e)):
                signal_above_curve_number += 1

    # Then loop over the gamma cuts. For each cuts, we have a background file.
    # Create one plot for each cut.
    for pta_cuts in ptcuts:

        bg_pass_cut_number = 0  # counting background event that pass the cuts
        bg_above_curve_number = 0  # counting background events that above the curve
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
            "/Users/isaac/Work/MuonE/SMmuegamma/Events/",
            event_folder_name, "unweighted_events.root")
        bg_chain.Add(bg_root_file)
        bg_theta_e_list = []
        bg_theta_mu_list = []
        for event in bg_chain:
            for particle in event.Particle:
                if(particle.Status == 1 and particle.PID == 11):
                    electron_vector = r.TLorentzVector()
                    electron_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M)
                    theta_e = electron_vector.Theta() * 1000
                elif(particle.Status == 1 and particle.PID == 13):
                    muon_vector = r.TLorentzVector()
                    muon_vector.SetPtEtaPhiM(
                        particle.PT, particle.Eta, particle.Phi, particle.M)
                    theta_mu = muon_vector.Theta() * 1000
            if(theta_e <= 30):
                bg_theta_e_list.append(theta_e)
                bg_theta_mu_list.append(theta_mu)
                bg_pass_cut_number += 1
                if(theta_mu > theta_mu_theory(theta_e)):
                    bg_above_curve_number += 1
        # This step must be applied. Because we need a new TChain to read
        # the new rootfile.
        bg_chain.Reset()
        efficiency_file.write(str(pt4_cuts) + '       ')
        efficiency_file.write(str(pta_cuts) + '       ')
        efficiency_file.write(
            str(signal_pass_cut_number) + '                    ')
        efficiency_file.write(
            str(signal_above_curve_number) + '                 ')
        efficiency_file.write(str(bg_pass_cut_number) + '                    ')
        efficiency_file.write(str(bg_above_curve_number) + '\n')

        # Plot the points.
        ax.scatter(signal_theta_e_list, signal_theta_mu_list, color="blue")
        ax.scatter(bg_theta_e_list, bg_theta_mu_list, color="orange")
        ax.plot(theta_e_variable, mu_variable, color="r")

# ---------------------------------------------------------------------------
# Finish
# ---------------------------------------------------------------------------
fig.subplots_adjust(hspace=0.3, wspace=0.3)
fig.savefig("compare_all_cuts.png", dpi=300)
efficiency_file.close()
