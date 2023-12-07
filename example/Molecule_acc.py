import os
import numpy as np
import argparse

import sys
import time
from qutip import Qobj

sys.path.append("../")
from utils.auxiliary_molecule import *
from utils.evolution import compute_TV_norm, compute_obj_by_switch, time_evolution, compute_obj_fid, test_trace
from switching.switch_time_acc import SwitchTimeOptAcc
from switching.obtain_switches_acc_sos1 import SwitchesAccSOS1

parser = argparse.ArgumentParser()
# name of example
parser.add_argument('--name', help='example name', type=str, default='MoleculeNew')
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
                    default=None)
# method to extract control sequence: binary, obj, sur
parser.add_argument('--extract', help='method to extract control sequence', type=str, default='obj')
# # threshold for choosing the maximum controller
# parser.add_argument('--thre_ratio', help='threshold for choosing the maximum controller', type=float, default=0)
# weight parameter for TV regularizer
parser.add_argument('--alpha', help='weight parameter for TV regularizer', type=float, default=0)
# if extracted from binary control, give the parameter
parser.add_argument('--b_control', help='binary control for extraction', type=str, default=None)
# parser.add_argument('--warm_start', help='warm start for TV extraction', type=int, default=0)
# file store the target circuit
parser.add_argument('--target', help='unitary matrix of target circuit', type=str,
                    default="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv")
# available option: SQ, SU, PSU
parser.add_argument('--phase', help='computational way for infidelity function', type=str, default="PSU")
args = parser.parse_args()

# args.target = "../result/control/Target/MoleculeVQE_BeH2_evotime20.0_n_ts200_target.csv"
# args.molecule = "BeH2"
# args.qubit_num = 6

# args.n_ts = 64
# args.alpha = 0.001
# args.extract = "tv"
# args.warm_start = 1

d = 2
if args.target is not None:
    Hops, H0, U0, U = generate_molecule_func(args.qubit_num, d, args.molecule, optimize=True, target=args.target)
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

# for hops in Hops:
#     print(test_trace(hops + H0, U0, U, 20))
# print(test_trace(sum(hops for hops in Hops) / len(Hops) + H0, U0, U, 100))
# exit()

if args.c_control is None:
    print("Must provide continuous control results!")
    exit()

if args.extract == 'binary':
    parameter = 'alpha_' + str(args.alpha) + '_mintv'
    param = {'binary_file': args.b_control}
if args.extract in ['obj', 'sur']:
    parameter = 'alpha_' + str(args.alpha)
    param = {'alpha': args.alpha}

switches = SwitchesAccSOS1(args.c_control, delta_t=args.evo_time / args.n_ts)
start1 = time.time()
switches.init_gradient_computer(H0, Hops, U0, U, args.n_ts, args.evo_time, 'fid', phase_option=args.phase)
warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches(args.extract, param=param)
diag, matrix = switches.get_decomposition()
end1 = time.time()

# print("extracted!")
# sequence of control hamiltonians
# ctrl_hamil = [(H_d + H_c[j]).full() for j in range(len(H_c))]

# initial control
if args.initial_type == "ave":
    initial = np.ones(num_switch + 1) * args.evo_time / (num_switch + 1)
if args.initial_type == "rnd":
    initial_pre = np.random.random(num_switch + 1)
    initial = initial_pre.copy() / sum(initial_pre) * args.evo_time
if args.initial_type == "warm":
    initial = warm_start_length

# build optimizer
molecule_opt = SwitchTimeOptAcc()
min_time = 1
molecule_opt.build_optimizer(
    Hops, H0, ctrl_hamil_idx, initial, X_0.full(), X_targ.full(), args.evo_time, num_switch, min_up_time, None,
    phase=args.phase, decompose=None)
# molecule_opt.build_optimizer(
#     ctrl_hamil, ctrl_hamil_idx, initial, X_0.full(), X_targ.full(), args.evo_time, num_switch, min_up_time, None,
#     phase=args.phase, decompose={'diag':diag, 'matrix':matrix})

start2 = time.time()
res = molecule_opt.optimize()
end2 = time.time()

if not os.path.exists("../result/output/SwitchTime/"):
    os.makedirs("../result/output/SwitchTime/")
if not os.path.exists("../result/control/SwitchTime/"):
    os.makedirs("../result/control/SwitchTime/")
if not os.path.exists("../result/figure/SwitchTime/"):
    os.makedirs("../result/figure/SwitchTime/")

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}_extraction".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, parameter) + ".png"
# if args.extract == 'tv':
#     fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_alpha{}_extraction".format(
#         args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
#         args.initial_type, str(min_up_time), args.alpha) + ".png"
switches.draw_extracted_control(fig_name)

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}_metric".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, parameter) + ".png"

if args.extract in ['obj', 'sur']:
    switches.draw_metric(fig_name, args.extract)

output_name = "../result/output/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, parameter) + ".log"

output_file = open(output_name, "a+")
# if args.extract == "la":
#     print("difference", switches.binary_obj - switches.continuous_obj, file=output_file)
if args.extract == "obj":
    print("difference with alpha=0", sum(min(switches.cost[k, :]) for k in range(args.n_ts)), file=output_file)
    print("summation of all the changes of all the controllers", sum([sum(switches.cost[k, j] for j in range(len(H_c)))
                                                                      for k in range(args.n_ts)]), file=output_file)
    print("maximum summation change of all the controllers", max([sum(switches.cost[k, j] for j in range(len(H_c)))
                                                                  for k in range(args.n_ts)]), file=output_file)

if args.extract == 'sur':
    # print("difference of obj", switches.binary_obj - switches.continuous_obj, file=output_file)
    print("all objective value of milp", switches.overall_obj, file=output_file)
    print("difference of control", switches.diff_obj, file=output_file)
    print("tv norm obj", switches.tv_obj, file=output_file)

print("controller sequence is", ctrl_hamil_idx, file=output_file)
print("warm start length is", warm_start_length, file=output_file)
print("objective of continuous control", switches.continuous_obj, file=output_file)
print("objective function before optimization", compute_obj_by_switch(
    Hops, H0, warm_start_length, ctrl_hamil_idx, X_0.full(), X_targ.full(), 'fid', phase=args.phase), file=output_file)
print("TV regularizer before optimization", 2 * len(warm_start_length) - 2, file=output_file)
print(res, file=output_file)
print("switching time points", molecule_opt.switch_time, file=output_file)
print("computational time of retrieving switches", end1 - start1, file=output_file)
print("computational time of optimization", end2 - start2, file=output_file)
print("pre-computing time, objective function evaluation time, gradient time", molecule_opt.pre_compute_time,
      molecule_opt.evaluate_time, molecule_opt.gradient_time, file=output_file)
print("total computational time", end2 - start1, file=output_file)

# retrieve control
control_name = "../result/control/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, parameter) + ".csv"
control = molecule_opt.retrieve_control(args.n_ts)
np.savetxt(control_name, control, delimiter=",")

convert_switch_control, filter_obj = molecule_opt.filter_control()
tv_norm = molecule_opt.tv_norm()
print("tv norm", tv_norm, file=output_file)

# figure file
# figure_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}".format(
#     args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
#     args.initial_type, str(min_up_time), args.thre_ratio) + ".png"
figure_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}".format(
    args.name + args.extract + "_" + args.molecule, str(args.evo_time), str(args.n_ts), str(num_switch),
    args.initial_type, parameter) + ".png"
molecule_opt.draw_control(figure_name)

b_bin = np.loadtxt(control_name, delimiter=",")
bin_result = time_evolution(H0, Hops, args.n_ts, args.evo_time, b_bin, X_0.full(), False, 1)
f = open(output_name, "a+")
print("filtered objective value", filter_obj, file=f)
print("binary objective value", compute_obj_fid(X_targ, bin_result, phase=args.phase), file=f)
print("total binary tv norm", compute_TV_norm(b_bin), file=f)
print("initial file", args.c_control, file=f)
f.close()