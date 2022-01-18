import os

mg5_command_file_path = "/Users/isaac/Work/MuonE/bg.mg5"
output_path = "/Users/isaac/Work/MuonE/SMmuegamma/"


ptcuts = [0.001, 0.01, 0.05]


def run_madgraph_with_ptcuts(pt4cut, ptacut):
    mg5_command_file = open(mg5_command_file_path, "w")
    run_name = "pt4_" + str(pt4cut) + "_pta_" + str(ptacut)
    mg5_command = "launch -n " + run_name + " " + output_path + " \n"
    mg5_command = mg5_command + "set ebeam1 = 150.0 \n"
    mg5_command = mg5_command + "set ebeam2 = 5.11e-04 \n"
    mg5_command = mg5_command + "set ptj = 0.0 \n"
    mg5_command = mg5_command + "set ptl = 0.0 \n"
    mg5_command = mg5_command + "set etaj = -1.0 \n"
    mg5_command = mg5_command + "set etaa = -1.0 \n"
    mg5_command = mg5_command + "set etal = -1.0 \n"
    mg5_command = mg5_command + "set drjj = 0.0 \n"
    mg5_command = mg5_command + "set drll = 0.0 \n"
    mg5_command = mg5_command + "set draa = 0.0 \n"
    mg5_command = mg5_command + "set draj = 0.0 \n"
    mg5_command = mg5_command + "set drjl = 0.0 \n"
    mg5_command = mg5_command + "set dral = 0.0 \n"
    mg5_command = mg5_command + "set pt_min_pdg {11: " + str(pt4cut) + "} \n"
    mg5_command = mg5_command + "set pta " + str(ptacut) + " \n"
    mg5_command = mg5_command + "done"
    mg5_command_file.write(mg5_command)
    mg5_command_file.close()
    run_command = '~/Work/MG5_aMC_v2_6_7/bin/mg5_aMC ' + mg5_command_file_path
    os.system(run_command)


for pt4cuts in ptcuts:
    for ptacuts in ptcuts:
        # run madgraph with such parameter
        run_madgraph_with_ptcuts(pt4cuts, ptacuts)
