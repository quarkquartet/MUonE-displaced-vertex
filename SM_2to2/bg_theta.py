"""Plot out the events in the \theta_\mu \theta_e plane."""
from __future__ import print_function
import sys
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
CHAIN1.Add(
    "~/Work/MuonE/SMmuegamma/Events/pt4_0.01_pta_0.01/unweighted_events.root")
CHAIN2 = r.TChain("LHEF")
CHAIN2.Add(
    "~/Work/MuonE/SMmuegamma/Events/pt4_0.05_pta_0.01/unweighted_events.root")
CHAIN3 = r.TChain("LHEF")
CHAIN3.Add(
    "~/Work/MuonE/SMmuegamma/Events/pt4_0.01_pta_0.05/unweighted_events.root")
CHAIN4 = r.TChain("LHEF")
CHAIN4.Add(
    "~/Work/MuonE/SMmuegamma/Events/pt4_0.05_pta_0.05/unweighted_events.root")

theta_e_list1 = []
theta_mu_list1 = []
theta_e_list2 = []
theta_mu_list2 = []
theta_e_list3 = []
theta_mu_list3 = []
theta_e_list4 = []
theta_mu_list4 = []

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

for event in CHAIN4:
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
    theta_e_list4.append(theta_e)
    theta_mu_list4.append(theta_mu)

fig = plt.figure()
fig.suptitle(r'$\theta_\mu - \theta_e$')
ax1 = fig.add_subplot(221)
ax1.scatter(theta_e_list1, theta_mu_list1, color="blue")
ax1.set_xlim(0, 30)
ax1.set_title("pt4=0.01, pta=0.01")
ax2 = fig.add_subplot(222)
ax2.scatter(theta_e_list2, theta_mu_list2, color="magenta")
ax2.set_xlim(0, 30)
ax2.set_title("pt4=0.05, pta=0.01")
ax3 = fig.add_subplot(223)
ax3.scatter(theta_e_list3, theta_mu_list3, color="orange")
ax3.set_xlim(0, 30)
ax3.set_title("pt4=0.01, pta=0.05")
ax4 = fig.add_subplot(224)
ax4.scatter(theta_e_list4, theta_mu_list4, color="green")
ax4.set_xlim(0, 30)
ax4.set_title("pt4=0.05, pta=0.05")
fig.subplots_adjust(hspace=0.6)

fig.savefig("bg_theta_relation.png", dpi=300)
