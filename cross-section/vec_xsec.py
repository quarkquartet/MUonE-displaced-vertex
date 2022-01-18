"""Read the cross sections for 2-3 signal process"""
import os

import numpy as np

event_path = "/Users/isaac/Work/MUonE/MG5_events/2to3vec/Events"


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


xsc_file = open(
    "/Users/isaac/Work/MUonE/share/2to3vec/vec_displaced_vertex/vec_xsec.txt", "w"
)
mass_list = np.logspace(np.log10(0.002), np.log10(0.29), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]
with xsc_file:
    xsc_file.write("X_mass (GeV)        x_sec (pb)\n")
    for mass in mass_list:
        x_sec = read_x_section(mass)
        x_sec = "{:.2e}".format(x_sec)
        xsc_file.write(str(mass) + "           " + x_sec + "\n")

xsc_file.close()
