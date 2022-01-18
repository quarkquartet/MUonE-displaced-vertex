"""
Simulate sample events for different masses. Generation commands
and "output" are already done in MG5 interface. This part clears all
cuts and adjusts the mass parameters.
"""
import os

import numpy as np

mg5_command_file_path = "/Users/isaac/Work/MuonE/code/2to3vec/vec_mass.mg5"
output_path = "/Users/isaac/Work/MUonE/MG5_events/2to3vec10k"

mass_spectrum = np.logspace(np.log10(0.002), np.log10(0.29), 30).tolist()
mass_spectrum = [round(x, 4) for x in mass_spectrum]


def run_madgraph(X_mass):
    mg5_command_file = open(mg5_command_file_path, "w")
    run_name = "X_mass_" + str(X_mass)
    mg5_command = "launch -n " + run_name + " " + output_path + " \n"
    # mg5_command = mg5_command + "set nevents = 100000 \n"
    mg5_command = mg5_command + "set ebeam1 = 150.0 \n"
    mg5_command = mg5_command + "set ebeam2 = 5.11e-04 \n"
    mg5_command = mg5_command + "set Mvec = " + str(X_mass) + " \n"
    mg5_command = mg5_command + "set ptl = 0.0 \n"
    mg5_command = mg5_command + "set etal = -1.0 \n"
    mg5_command = mg5_command + "set drll = 0.0 \n"
    mg5_command = mg5_command + "done"
    mg5_command_file.write(mg5_command)
    mg5_command_file.close()
    run_command = "~/Work/MG5_aMC_v3_1_0/bin/mg5_aMC " + mg5_command_file_path
    os.system(run_command)


for X_mass in mass_spectrum:
    run_madgraph(X_mass)
