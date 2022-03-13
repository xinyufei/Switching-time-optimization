import os
import numpy as np
import argparse

import sys

sys.path.append("../")
from switching.switch_time import *

parser = argparse.ArgumentParser()
# name of example
parser.add_argument('--name', help='example name', type=str, default='MoleculeST')
# name of molecule
parser.add_argument('--molecule', help='molecule name', type=str, default='H2')
# number of quantum bits
parser.add_argument('--qubit_num', help='number of quantum bits', type=int, default=2)
# evolution time
parser.add_argument('--evo_time', help='evolution time', type=float, default=4)
# time steps
parser.add_argument('--n_ts', help='time steps', type=int, default=80)
# initial type: ave, warm, rnd
parser.add_argument('--initial_type', help='initial type of control variables (rnd, ave, warm)', type=str,
                    default="warm")
# initial control file obtained from ADMM algorithm
parser.add_argument('--c_control', help='file name of initial control', type=str,
                    default="../result/control/ADMM/MoleculeNewH2_evotime4_n_ts80_ptypeWARM_offset0_objUNIT_penalty0.001_sum_penalty0.01.csv")
# method to extract control sequence: obj, la, sur
parser.add_argument('--extract', help='method to extract control sequence', type=str, default='obj')
# threshold for choosing the maximum controller
parser.add_argument('--thre_ratio', help='threshold for choosing the maximum controller', type=float, default=0)
# file store the target circuit
parser.add_argument('--target', help='unitary matrix of target circuit', type=str, default=None)

args = parser.parse_args()

# args.molecule = "LiH"
# args.qubit_num = 4

d = 2
Hops, H0, U0, U = generate_molecule_func(args.qubit_num, d, args.molecule)

if args.target is not None:
    U = np.loadtxt(args.target, dtype=np.complex_, delimiter=',')
else:
    print("Please provide the target file!")
    exit()

# n_ts = 80
# evo_time = 4

min_up_time = 0

# The control Hamiltonians (Qobj classes)
H_c = [Qobj(hops) for hops in Hops]
# Drift Hamiltonian
H_d = Qobj(H0)
# start point for the gate evolution
X_0 = Qobj(U0)
# Target for the gate evolution
X_targ = Qobj(U)

if args.c_control is None:
    print("Must provide continuous control results!")
    exit()

switches = Switches(args.c_control, delta_t=args.evo_time / args.n_ts)
start1 = time.time()
switches.init_gradient_computer(H0, Hops, U0, U, args.n_ts, args.evo_time, 'fid')
warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches(args.extract, thre_ratio=args.thre_ratio)
end1 = time.time()

# sequence of control hamiltonians
ctrl_hamil = [(H_d + H_c[j]).full() for j in range(len(H_c))]

# initial control
if args.initial_type == "ave":
    initial = np.ones(num_switch + 1) * args.evo_time / (num_switch + 1)
if args.initial_type == "rnd":
    initial_pre = np.random.random(num_switch + 1)
    initial = initial_pre.copy() / sum(initial_pre) * args.evo_time
if args.initial_type == "warm":
    initial = warm_start_length

# build optimizer
molecule_opt = SwitchTimeOpt()
min_time = 1
molecule_opt.build_optimizer(
    ctrl_hamil, ctrl_hamil_idx, initial, X_0.full(), X_targ.full(), args.evo_time, num_switch, min_up_time, None)
start2 = time.time()
res = molecule_opt.optimize()
end2 = time.time()

if not os.path.exists("../result/output/SwitchTime/"):
    os.makedirs("../result/output/SwitchTime/")
if not os.path.exists("../result/control/SwitchTime/"):
    os.makedirs("../result/control/SwitchTime/")
if not os.path.exists("../result/figure/SwitchTime/"):
    os.makedirs("../result/figure/SwitchTime/")

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}_extraction".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, str(min_up_time), args.thre_ratio) + ".png"
switches.draw_extracted_control(fig_name)

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}_metric".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, str(min_up_time), args.thre_ratio) + ".png"
switches.draw_metric(fig_name, args.extract)
# exit()

# output file
output_name = "../result/output/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, str(min_up_time), args.thre_ratio) + ".log"
output_file = open(output_name, "a+")
if args.extract == "obj":
    print("difference", sum(min(switches.cost[k, :]) for k in range(args.n_ts)), file=output_file)
    print("summation of all the changes of all the controllers", sum([sum(switches.cost[k, j] for j in range(len(H_c)))
                                                                      for k in range(args.n_ts)]), file=output_file)
    print("maximum summation change of all the controllers", max([sum(switches.cost[k, j] for j in range(len(H_c)))
                                                                  for k in range(args.n_ts)]), file=output_file)
print("objective function before optimization", compute_obj_by_switch(
    Hops, warm_start_length, ctrl_hamil_idx, X_0.full(), X_targ.full(), 'fid'), file=output_file)
print("TV regularizer before optimization", 2 * len(warm_start_length) - 2, file=output_file)
print(res, file=output_file)
print("switching time points", molecule_opt.switch_time, file=output_file)
print("computational time of retrieving switches", end1 - start1, file=output_file)
print("computational time of optimization", end2 - start2, file=output_file)
print("total computational time", end2 - start1, file=output_file)

# retrieve control
control_name = "../result/control/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, str(min_up_time), args.thre_ratio) + ".csv"
control = molecule_opt.retrieve_control(args.n_ts)
np.savetxt(control_name, control, delimiter=",")

tv_norm = molecule_opt.tv_norm()
print("tv norm", tv_norm, file=output_file)

# figure file
figure_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, str(min_up_time), args.thre_ratio) + ".png"
molecule_opt.draw_control(figure_name)

b_bin = np.loadtxt(control_name, delimiter=",")
f = open(output_name, "a+")
print("total tv norm", compute_TV_norm(b_bin), file=f)
print("initial file", args.c_control, file=f)
f.close()