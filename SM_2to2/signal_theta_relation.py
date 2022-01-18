"""Plot out the events in the \theta_\mu \theta_e plane for the signal events"""
from __future__ import print_function
import os
import random
from matplotlib import pyplot as plt
import ROOT as r

# Import PyROOT classes from madgraph.
# Madgraph uses its own defined ROOT classes to write ROOT files. Import these
# libraries so that PyROOT can recognize the file.
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

CHAIN1 = r.TChain("LHEF")
CHAIN1.Add("~/Work/MuonE/2to3sc/Events/nocut_0.001GeV/unweighted_events.root")
CHAIN2 = r.TChain("LHEF")
CHAIN2.Add("~/Work/MuonE/2to3sc/Events/nocut_0.01GeV/unweighted_events.root")
CHAIN3 = r.TChain("LHEF")
CHAIN3.Add("~/Work/MuonE/2to3sc/Events/nocut_0.1GeV/unweighted_events.root")

theta_e_list1 = []
theta_mu_list1 = []
theta_e_list2 = []
theta_mu_list2 = []
theta_e_list3 = []
theta_mu_list3 = []

for event in CHAIN1:
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 11):
            electron_vector = r.TLorentzVector()
            electron_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_e = electron_vector.Theta() * 1000
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 13):
            muon_vector = r.TLorentzVector()
            muon_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_mu = muon_vector.Theta() * 1000
    theta_e_list1.append(theta_e)
    theta_mu_list1.append(theta_mu)

for event in CHAIN2:
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 11):
            electron_vector = r.TLorentzVector()
            electron_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_e = electron_vector.Theta() * 1000
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 13):
            muon_vector = r.TLorentzVector()
            muon_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_mu = muon_vector.Theta() * 1000
    theta_e_list2.append(theta_e)
    theta_mu_list2.append(theta_mu)

for event in CHAIN3:
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 11):
            electron_vector = r.TLorentzVector()
            electron_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_e = electron_vector.Theta() * 1000
    for particle in event.Particle:
        if(particle.Status == 1 and particle.PID == 13):
            muon_vector = r.TLorentzVector()
            muon_vector.SetPtEtaPhiM(
                particle.PT, particle.Eta, particle.Phi, particle.M)
            theta_mu = muon_vector.Theta() * 1000
    theta_e_list3.append(theta_e)
    theta_mu_list3.append(theta_mu)

fig = plt.figure()
fig.suptitle(r'$\theta_\mu - \theta_e$')
ax1 = fig.add_subplot(131)
ax1.scatter(theta_e_list1, theta_mu_list1, color="blue")
ax1.set_xlim(0, 30)
ax1.set_title("mX=1 MeV")
ax2 = fig.add_subplot(132)
ax2.scatter(theta_e_list2, theta_mu_list2, color="orange")
ax2.set_xlim(0, 30)
ax2.set_title("mX = 10 MeV")
ax3 = fig.add_subplot(133)
ax3.scatter(theta_e_list3, theta_mu_list3, color="green")
ax3.set_xlim(0, 30)
ax3.set_title("mX = 100 MeV")
fig.savefig("signal_theta_relation.png", dpi=300)
