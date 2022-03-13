import argparse
import os
import sys
import time


sys.path.append("../")
from tools.auxiliary_energy import *
from tools.evolution import compute_TV_norm, compute_obj_by_switch
from switching.switch_time import SwitchTimeOpt
from switching.obtain_switches import Switches

parser = argparse.ArgumentParser()
# name of example
parser.add_argument('--name', help='example name', type=str, default='EnergyST')
# number of qubits
parser.add_argument('--n', help='number of qubits', type=int, default=2)
# number of edges for generating regular graph
parser.add_argument('--num_edges', help='number of edges for generating regular graph', type=int, default=1)
# if generate the graph randomly
parser.add_argument('--rgraph', help='if generate the graph randomly', type=int, default=0)
# number of instances
parser.add_argument('--seed', help='random seed', type=int, default=0)
# evolution time
parser.add_argument('--evo_time', help='evolution time', type=float, default=2)
# time steps
parser.add_argument('--n_ts', help='time steps', type=int, default=40)
# initial type: ave, warm, rnd
parser.add_argument('--initial_type', help='type of initialized variables', type=str, default='warm')
# initial control file obtained from ADMM algorithm
parser.add_argument('--c_control', help='file name of initial control', type=str,
                    default="../result/control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv")
# method to extract control sequence: obj, la, sur
parser.add_argument('--extract', help='method to extract control sequence', type=str, default='obj')
# threshold for choosing the maximum controller
parser.add_argument('--thre_ratio', help='threshold for choosing the maximum controller', type=float, default=0)
# # minimum up time constraint
# parser.add_argument('--min_up_time', help='minimum up time', type=float, default=0)
# # tv regularizer parameter
# parser.add_argument('--alpha', help='tv regularizer parameter', type=float, default=0.05)

args = parser.parse_args()

if args.rgraph == 0:
    Jij, edges = generate_Jij_MC(args.n, args.num_edges, 100)

    C = get_ham(args.n, True, Jij)
    B = get_ham(args.n, False, Jij)

    args.seed = 0

if args.rgraph == 1:
    Jij = generate_Jij(args.n, args.seed)
    C = get_ham(args.n, True, Jij)
    B = get_ham(args.n, False, Jij)

y0 = uniform(args.n)

if args.c_control is None:
    print("Must provide continuous control results!")
    exit()

switches = Switches(args.c_control, delta_t=args.evo_time / args.n_ts)
start1 = time.time()
switches.init_gradient_computer(None, [B, C], y0[0:2 ** args.n], None, args.n_ts, args.evo_time, 'energy')
warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches(args.extract, thre_ratio=args.thre_ratio)
end1 = time.time()

# sequence of control hamiltonians
ctrl_hamil = [B, C]

# X_0 = np.expand_dims(y0[0:2**n], 1)
X_0 = y0[0:2 ** args.n]

if args.initial_type == "ave":
    initial = np.ones(num_switch + 1) * args.evo_time / (num_switch + 1)
if args.initial_type == "rnd":
    initial_pre = np.random.random(num_switch + 1)
    initial = initial_pre.copy() / sum(initial_pre) * args.evo_time
if args.initial_type == "warm":
    initial = warm_start_length

# build optimizer
min_up_time = 0
energy_opt = SwitchTimeOpt()
start2 = time.time()
energy_opt.build_optimizer(
    ctrl_hamil, ctrl_hamil_idx, initial, X_0, None, args.evo_time, num_switch, min_up_time, None,
    obj_type='energy')
res = energy_opt.optimize()
end2 = time.time()

if not os.path.exists("../result/output/SwitchTime/"):
    os.makedirs("../result/output/SwitchTime/")
if not os.path.exists("../result/control/SwitchTime/"):
    os.makedirs("../result/control/SwitchTime/")
if not os.path.exists("../result/figure/SwitchTime/"):
    os.makedirs("../result/figure/SwitchTime/")

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_instance{}_extraction_thre{}".format(
    args.name + str(args.n) + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type,
    str(min_up_time), args.seed, args.thre_ratio) + ".png"
switches.draw_extracted_control(fig_name)

fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_instance{}_extraction_thre{}".format(
    args.name + str(args.n) + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type,
    str(min_up_time), args.seed, args.thre_ratio) + ".png"
switches.draw_metric(fig_name, args.extract)

# output file
output_name = "../result/output/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_instance{}_extraction_thre{}".format(
    args.name + str(args.n) + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type,
    str(min_up_time), args.seed, args.thre_ratio) + ".log"
output_file = open(output_name, "a+")

if args.extract == "obj":
    print("difference", sum(min(switches.cost[k, :]) for k in range(args.n_ts)), file=output_file)
    print("summation of all the changes of all the controllers", sum([sum(switches.cost[k, j] for j in range(2))
                                                                      for k in range(args.n_ts)]), file=output_file)
    print("maximum summation change of all the controllers", max([sum(switches.cost[k, j] for j in range(2))
                                                                  for k in range(args.n_ts)]), file=output_file)
# exit()

print("objective function before optimization", compute_obj_by_switch(
    [B, C], warm_start_length, ctrl_hamil_idx, X_0, None, 'energy'), file=output_file)
print("TV regularizer before optimization", 2 * len(warm_start_length) - 2, file=output_file)
print(res, file=output_file)
print("switching time points", energy_opt.switch_time, file=output_file)
print("computational time of retrieving switches", end1 - start1, file=output_file)
print("computational time of optimization", end2 - start2, file=output_file)
print("total computational time", end2 - start1, file=output_file)

# retrieve control
control_name = "../result/control/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_instance{}_extraction_thre{}".format(
    args.name + str(args.n) + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type,
    str(min_up_time), args.seed, args.thre_ratio) + ".csv"
control = energy_opt.retrieve_control(args.n_ts)
np.savetxt(control_name, control, delimiter=",")

# print("alpha", alpha, file=output_file)
tv_norm = energy_opt.tv_norm()
print("tv norm", tv_norm, file=output_file)
# print("objective with tv norm", energy_opt.obj + alpha * tv_norm, file=output_file)

# figure file
figure_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_instance{}_extraction_thre{}".format(
    args.name + str(args.n) + args.extract, str(args.evo_time), str(args.n_ts), str(num_switch), args.initial_type,
    str(min_up_time), args.seed, args.thre_ratio) + ".png"
energy_opt.draw_control(figure_name)

b_bin = np.loadtxt(control_name, delimiter=",")
f = open(output_name, "a+")
print("total tv norm", compute_TV_norm(b_bin), file=f)
print("initial file", args.c_control, file=f)
f.close()
