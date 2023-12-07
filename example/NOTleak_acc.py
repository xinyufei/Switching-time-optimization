import argparse
import os
import sys
import math
import time
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np

from qutip import identity, Qobj
from qutip.qip.operations.gates import cnot

sys.path.append("../")
from utils.auxiliary_energy import *
from utils.evolution import compute_TV_norm, compute_obj_by_switch, time_evolution, compute_obj_fid, test_trace
from utils.modify_control import remove_sos1
from switching.switch_time_acc import SwitchTimeOptAcc
from switching.obtain_switches_acc import SwitchesAcc

parser = argparse.ArgumentParser()
# name of example
parser.add_argument('--name', help='example name', type=str, default='NOTleak')
# evolution time
parser.add_argument('--evo_time', help='evolution time', type=float, default=5)
# time steps
parser.add_argument('--n_ts', help='time steps', type=int, default=20)
# initial type
parser.add_argument('--initial_type', help='initial controls type', type=str, default='CONSTANT')
# initial control file obtained from ADMM algorithm
parser.add_argument('--c_control', help='file name of initial control', type=str, default=None)
# method to extract control sequence: binary, obj, sur
parser.add_argument('--extract', help='method to extract control sequence', type=str, default='obj')
# # threshold for choosing the maximum controller
# parser.add_argument('--thre_ratio', help='threshold for choosing the maximum controller', type=float, default=0)
# if extracted from binary control, give the parameter
parser.add_argument('--b_control', help='binary control for extraction', type=str, default=None)
# weight parameter for TV regularizer
parser.add_argument('--alpha', help='weight parameter for TV regularizer', type=float, default=0)
# available option: SQ, SU, PSU
parser.add_argument('--phase', help='computational way for infidelity function', type=str, default="PSU")

args = parser.parse_args()
# QuTiP control modules
# a two-qubit system with target control mode as CNOT gate with summation one constraint
# The control Hamiltonians (Qobj classes)
# Drift Hamiltonian
H_d = 0 * np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]], dtype=complex) + 2 * math.pi * np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=complex)
h_c = [1 / 2 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
       + np.sqrt(2) / 2 * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
       1 / 2 * np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
       + np.sqrt(2) / 2 * np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)]
H_c = [np.zeros((3, 3), dtype=complex)]
s = [j for j in range(2)]
for n in range(1, 3):
    for i in combinations(s, n):
        cur_list = list(i)
        H_c.append(sum(h_c[j] for j in cur_list))
# print(H_d, H_c)
H_d_qobj = Qobj(H_d)
H_c_qobj = [Qobj(hc) for hc in H_c]
# start point for the gate evolution
X_0_qobj = identity(3)
X_0 = X_0_qobj.full()
# Target for the gate evolution
X_targ = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
X_targ_qobj = Qobj(X_targ)

if args.c_control is None:
    print("Must provide continuous control results!")
    exit()

if args.extract == 'binary':
    parameter = 'alpha_' + str(args.alpha) + '_mintv'
    param = {'binary_file': args.b_control}
if args.extract in ['obj', 'sur']:
    parameter = 'alpha_' + str(args.alpha)
    param = {'alpha': args.alpha}

switches = SwitchesAcc(args.c_control, delta_t=args.evo_time / args.n_ts)
start1 = time.time()
switches.init_gradient_computer(H_d, H_c, X_0, X_targ, args.n_ts, args.evo_time, 'fid', phase_option=args.phase)
warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches(args.extract, param=param)
end1 = time.time()

# ctrl_hamil = [H_d + H_c[j] for j in range(len(H_c))]

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
    H_c, H_d, ctrl_hamil_idx, initial, X_0, X_targ, args.evo_time, num_switch, 0, None,
    phase=args.phase, decompose=None)

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
    args.name + args.extract, args.evo_time, args.n_ts, str(num_switch), args.initial_type, parameter) + ".png"
switches.draw_extracted_control(fig_name)

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}_metric".format(
    args.name + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type, parameter) + ".png"

if args.extract in ['obj', 'sur']:
    switches.draw_metric(fig_name, args.extract)

output_num = "../result/output/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}".format(
    args.name + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type, parameter) + ".log"

output_file = open(output_num, "a+")
# solve the optimization model
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
    H_c, H_d, warm_start_length, ctrl_hamil_idx, X_0, X_targ, 'fid', phase=args.phase), file=output_file)
print("TV regularizer before optimization", 2 * len(warm_start_length) - 2, file=output_file)
print(res, file=output_file)
print("switching time points", molecule_opt.switch_time, file=output_file)
print("computational time of retrieving switches", end1 - start1, file=output_file)
print("computational time of optimization", end2 - start2, file=output_file)
print("pre-computing time, objective function evaluation time, gradient time", molecule_opt.pre_compute_time,
      molecule_opt.evaluate_time, molecule_opt.gradient_time, file=output_file)
print("total computational time", end2 - start1, file=output_file)

# retrieve control
control_name = "../result/control/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}_metric".format(
    args.name + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type, parameter) + ".csv"
control = molecule_opt.retrieve_control(args.n_ts)
np.savetxt(control_name, control, delimiter=",")

convert_switch_control, filter_obj = molecule_opt.filter_control()
tv_norm = molecule_opt.tv_norm()
print("tv norm", tv_norm, file=output_file)

# figure file
figure_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_{}_metric".format(
    args.name + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type, parameter) + ".png"
molecule_opt.draw_control(figure_name)

b_bin = np.loadtxt(control_name, delimiter=",")
bin_result = time_evolution(H_d, H_c, args.n_ts, args.evo_time, b_bin, X_0, False, 1)
f = open(output_num, "a+")
print("binary objective value", compute_obj_fid(X_targ_qobj, bin_result, phase=args.phase), file=f)
print("total binary tv norm", compute_TV_norm(b_bin), file=f)
print("initial file", args.c_control, file=f)

# convert new control
no_sos1_control = remove_sos1(convert_switch_control)
print("filtered objective value", filter_obj, file=f)
print("converted tv norm", compute_TV_norm(no_sos1_control), file=f)
f.close()