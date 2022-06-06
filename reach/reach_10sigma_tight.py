import pylhe
import numpy as np
import os
import csv
import sys

l1=1.5*0.01 + 10*0.001
l2 = 15.0*0.01 - 10*0.001
m_e = 0.000511
m_mu = 0.105
angle_cut = 0.001
meter = 1
GeV = 1 / (1.97 * meter * 1.0e-16)
Radius = 0.05
E_cut = 5.0
total_fixed_target_luminosity = 1.5*10**4
ee = np.sqrt(4*np.pi/137)
event_path = "/Users/isaac/Work/MUonE-displaced-vertex/Signal/Events"
result_path = "./reach_10sigma_tight.csv"
nlayers = 5

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


def read_x_section(mass):
    """Read the cross section from madgraph banner file"""
    event_file_folder = "/Users/isaac/Work/MUonE-displaced-vertex/Signal/Events"
    banner_file_subfolder = "X_mass_"+str(mass)
    banner_file_name = "X_mass_"+ str(mass) + '_tag_1_banner.txt'
    banner_file_path = os.path.join(event_file_folder,banner_file_subfolder, banner_file_name)
    with open(banner_file_path) as banner_file:
        for line in banner_file:
            if line.startswith("#  Integrated weight (pb)"):
                for x in line:
                    if x in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                        index = line.index(x)
                        break
                x_section_value = float(line[index:])
    return x_section_value

def angle_boost(theta, mass, energy):
    """theta angle of the decay products in the lab frame. Input the rest frame theta, dark photon masss and dark photon energy, we boost the angle theta angle of the decay product into the lab frame."""
    beta = np.sqrt(energy ** 2.0 - mass ** 2.0) / energy
    gamma = energy / mass
    pcom = np.sqrt(mass ** 2.0 - 4.0 * (m_e ** 2.0)) / 2.0
    nominator = gamma * (pcom * np.cos(theta) + beta * mass / 2.0)
    denominator = np.sqrt(
        (pcom ** 2.0) * (np.sin(theta) ** 2.0)
        + (gamma ** 2.0) * ((pcom * np.cos(theta) + beta * mass * 0.5) ** 2.0)
    )
    return np.arccos(nominator / denominator)


def energy_boost(theta, mass, energy):
    """energy of the decay product in the lab frame. Input the rest frame theta, dark photon masss and dark photon energy, we boost the energy of the decay product into the lab frame. """
    beta = np.sqrt(energy ** 2.0 - mass ** 2.0) / energy
    gamma = energy / mass
    pcom = np.sqrt(mass ** 2.0 - 4.0 * (m_e ** 2.0)) / 2.0
    return gamma * (0.5 * mass + beta * pcom * np.cos(theta))

class LHEdataset:

    def __init__(self, mA):
        self.mA = mA
        self.xsec = read_x_section(self.mA)
        self.four_momentum = np.zeros((100000,3,4), dtype="float64")
        self.conv_data = np.zeros((100000,17), dtype = "float64")
        self.event_file_folder = "/Users/isaac/Work/MUonE-displaced-vertex/Signal/Events"
        self.lhe_file_subfolder = "X_mass_"+str(self.mA)
        self.lhe_file_name = "unweighted_events.lhe"
        self.lhereader = pylhe.readLHE(os.path.join(self.event_file_folder, self.lhe_file_subfolder, self.lhe_file_name))
        self.read_lhe()
        self.convert_4_momenta()

    def read_lhe(self):
        for ievent, event in enumerate(self.lhereader):
            for particle in event.particles:
                if particle.id == 103:
                    self.four_momentum[ievent,0] = np.array([particle.e, particle.px, particle.py, particle.pz])
                if particle.id == 13 and particle.status == 1:
                    self.four_momentum[ievent,1] = np.array([particle.e, particle.px, particle.py, particle.pz])
                if particle.id == 11 and particle.status == 1:
                    self.four_momentum[ievent,2] = np.array([particle.e, particle.px, particle.py, particle.pz])
                
    def convert_4_momenta(self):
        # dark photon e, p, theta, phi
        self.conv_data[:,0] = self.four_momentum[:,0,0] # dark photon energy
        self.conv_data[:,1] = np.linalg.norm(self.four_momentum[:,0,1:4], axis = 1) # dark photon p
        self.conv_data[:,2] = np.arccos(self.four_momentum[:,0,3]/self.conv_data[:,1]) # dark photon theta
        self.conv_data[:,3] = np.arctan2(self.four_momentum[:,0,2], self.four_momentum[:,0,1]) # drk photon phi

        # muon e, p, theta, phi
        self.conv_data[:,4] = self.four_momentum[:,1,0] # muon energy
        self.conv_data[:,5] = np.linalg.norm(self.four_momentum[:,1,1:4], axis=1) # muon p
        self.conv_data[:,6] = np.arccos(self.four_momentum[:,1,3]/self.conv_data[:,5]) # muon theta
        self.conv_data[:,7] = np.arctan2(self.four_momentum[:,1,2], self.four_momentum[:,1,1]) # muon phi

        # electron e, p, theta, phi
        self.conv_data[:,8] = self.four_momentum[:,2,0] # electron energy
        self.conv_data[:,9] = np.linalg.norm(self.four_momentum[:,2,1:4], axis=1) # electron p
        self.conv_data[:,10] = np.arccos(self.four_momentum[:,2,3]/self.conv_data[:,9]) # electron theta
        self.conv_data[:,11] = np.arctan2(self.four_momentum[:,2,2], self.four_momentum[:,2,1]) # electron phi

        # simulation of decay products
        theta_restframe = np.arccos(np.random.uniform(-1.0,1.0,size=(100000,))) # theta in rest frame

        # lab frame boost
        self.conv_data[:,12] = angle_boost(theta_restframe, self.mA, self.conv_data[:,0]) # theta angle of decay product 1 in lab frame
        self.conv_data[:,13] = energy_boost(theta_restframe, self.mA,self.conv_data[:,0]) # energy of decay product 1 in lab frame
        self.conv_data[:,14] = angle_boost(theta_restframe + np.pi, self.mA, self.conv_data[:,0]) # theta angle of decay product 2 in lab frame
        self.conv_data[:,15] = energy_boost(theta_restframe + np.pi, self.mA, self.conv_data[:,0]) # energy of decay product 2 in lab frame

        # opening angle
        self.conv_data[:,16] = self.conv_data[:,12] + self.conv_data[:,14] # opening angle between decay products

couplings = np.logspace(-6.0, -2.0, 30).tolist()
mass_list = mass_list = np.logspace(np.log10(0.002), np.log10(0.2), 30).tolist()
mass_list = [round(x, 4) for x in mass_list]

with open(result_path, "w") as f:
    header = ['mass', 'coupling', 'number of events']
    data_writer = csv.writer(f)
    data_writer.writerow(header)
    for mA in mass_list:
        LHE = LHEdataset(mA)
        x_section_mass = LHE.xsec
        if float(len(LHE.four_momentum)) == 100000.0:
            print("LHE file of mA = " + str(mA) + " is loaded.")
        else: sys.exit("Error occured: did not load the data file successfully!")

        for coupling in couplings:
            efficiency_list = []
            for j in range(0,nlayers):
                axis = j+1
                selection = (((LHE.conv_data[:,2] <= (axis*Radius))&(LHE.conv_data[:,6] <= (axis*Radius)) & (LHE.conv_data[:,10] <= (axis*Radius)) & (LHE.conv_data[:,13] >= E_cut) & (LHE.conv_data[:,15] >= E_cut)&(LHE.conv_data[:,4] >= E_cut)&(LHE.conv_data[:,8] >= E_cut) & (LHE.conv_data[:,16] >= angle_cut)))

                LHE_pass = LHE.conv_data[selection]
                pass_data = np.zeros((len(LHE_pass), 7), dtype="float64")
                pass_data[:,0] = LHE_pass[:,1] # dark photon momentum
                pass_data[:,1] = LHE_pass[:,2] # dark photon theta
                pass_data[:,2] = distance(pass_data[:,0], mA, coupling) # dark photon travel distance
                AB = axis - pass_data[:,2]
                BC = AB * LHE_pass[:,12]
                pass_data[:,3] = BC
                BD = axis * pass_data[:,1]
                BF = AB * LHE_pass[:,14]
                pass_data[:,4] = BF
                phi_e = np.random.uniform(-np.pi*0.5, np.pi*0.5, size=(len(LHE_pass),))
                sin_alpha = np.sin(phi_e) * BD/Radius
                alpha = np.arcsin(sin_alpha)
                BCmax = Radius * np.sin(phi_e - alpha)/np.sin(phi_e)
                BFmax = Radius * np.sin(phi_e + alpha)/np.sin(phi_e)
                pass_data[:,5] = BCmax
                pass_data[:,6] = BFmax
                selection_2 = ((pass_data[:,3] <= pass_data[:,5]) & (pass_data[:,4] <= pass_data[:,6]) & (l1 <= pass_data[:,2]) & (pass_data[:,2] <= l2))
                efficiency = float(len(pass_data[selection_2]))/100000.0
                efficiency_list.append(efficiency)
            total_efficiency = sum(efficiency_list)/float(nlayers)
            sensitivity = x_section_mass * total_fixed_target_luminosity * coupling**2 * total_efficiency * float(nlayers)/40.0
            data_writer.writerow([mA, coupling, sensitivity])
            print("mass: " + str(mA) + ", coupling: " + str(coupling) + ", sensitivity: " + str(sensitivity))

print("Done!")
