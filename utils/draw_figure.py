import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple, HandlerBase
from utils.auxiliary_molecule import generate_molecule_func
from utils.auxiliary_energy import *
from utils.evolution import compute_obj_fid, time_evolution

sys.path.append("../")
from wirte_energy_csv import write_csv, get_all_results
from qutip import Qobj


def draw_control(evo_time, n_ts, control, output_fig):
    plt.figure(dpi=300)
    plt.xlabel("Time")
    plt.ylabel("Control amplitude")
    # plt.ylim([-0.04, 1.04])
    marker_list = ['o', '^', '*', 's', 'P']
    line = ['-', 'dotted']
    marker_size_list = [5, 5, 8, 5, 8]
    for j in range(control.shape[1]):
        # plt.step(np.linspace(0, evo_time, n_ts + 1), np.hstack((control[:, j], control[-1, j])), marker_list[j % 5],
        #          where='post', linewidth=2, label='controller ' + str(j + 1), markevery=5,
        #          markersize=marker_size_list[j % 5], linestyle=line[j % 10])
        plt.step(np.linspace(0, evo_time, n_ts + 1), np.hstack((control[:, j], control[-1, j])),
                 where='post', linewidth=2, label='controller ' + str(j + 1), linestyle=line[j // 10])

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.legend()
    plt.tight_layout()
    plt.savefig(output_fig)


def draw_switchtime(evo_time, num_steps, num_control, control_sequence, length, output_fig, extraction=False):
    # plt.figure(dpi=300)
    # plt.xlabel("Time")
    # plt.ylabel("Control amplitude")
    # plt.ylim([-0.04, 1.04])
    marker_list = ['o', '^', '*', 's', 'P']
    line = ['-', 'dotted']
    marker_size_list = [5, 5, 8, 5, 8]
    points = int(1000 * evo_time)
    x = np.linspace(0, evo_time, points)
    y = np.zeros((num_control, points))
    cur_point = 0
    mark_step = num_steps / evo_time * 5 * 1000
    # for i in range(len(control_sequence)):
    #     for cur_control in control_sequence[i]:
    #         y[cur_control, cur_point:cur_point + int(length[i] * 1000)] = 1
    #     cur_point += int(length[i] * 1000)
    for i in range(len(control_sequence)):
        cur_control = control_sequence[i]
        y[cur_control, cur_point:cur_point + int(length[i] * 1000)] = 1
        cur_point += int(length[i] * 1000)
    # if cur_point != 1000 * evo_time:
    #     for cur_control in control_sequence[len(control_sequence) - 1]:
    #         y[cur_control, cur_point:points] = 1
    if cur_point != 1000 * evo_time:
        cur_control = control_sequence[len(control_sequence) - 1]
        y[cur_control, cur_point:points] = 1
    return x, y
    # for j in range(num_control):
    #     # plt.plot(x, y[j, :], linewidth=2, label='controller ' + str(j + 1), marker=marker_list[j % 5],
    #     #          markevery=round(points/10), markersize=marker_size_list[j % 5], linestyle=line[j // 10])
    #     if extraction:
    #         # plt.plot(x, y[j, :], linewidth=2, label='controller ' + str(j + 1), marker=marker_list[j % 5],
    #         #          markevery=mark_step, markersize=marker_size_list[j % 5], linestyle=line[j % 10])
    #         plt.plot(x, y[j, :], linewidth=2, label='controller ' + str(j + 1), linestyle=line[j // 10])
    #     else:
    #         plt.plot(x, y[j, :], linewidth=2, label='controller ' + str(j + 1), linestyle=line[j // 10])
    # # plt.show()
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # # plt.legend()
    # plt.tight_layout()
    # plt.savefig(output_fig)


def convert_control_idx(idx):
    cur_control_list = []
    if idx == 0:
        cur_control_list = []
    if idx in [1, 2]:
        cur_control_list = [idx - 1]
    if idx == 3:
        cur_control_list = [0, 1]
    return cur_control_list


def draw_switchtime_cnot(evo_time, num_control, control_sequence, length, output_fig):
    # plt.figure(dpi=300)
    # plt.xlabel("Time")
    # plt.ylabel("Control amplitude")
    marker_list = ['o', '^', '*', 's', 'P']
    line = ['-', 'dotted']
    marker_size_list = [5, 5, 8, 5, 8]
    points = int(1000 * evo_time)
    x = np.linspace(0, evo_time, points)
    y = np.zeros((num_control, points))
    cur_point = 0
    for i in range(len(control_sequence)):
        cur_control_list = convert_control_idx(control_sequence[i])
        # cur_control_list = control_sequence
        for cur_control in cur_control_list:
            y[cur_control, cur_point:cur_point + int(length[i] * 1000)] = 1
        cur_point += int(length[i] * 1000)
    if cur_point != 1000 * evo_time:
        cur_control_list = convert_control_idx(control_sequence[len(control_sequence) - 1])
        for cur_control in cur_control_list:
            y[cur_control, cur_point:points] = 1
    return x, y
    # for j in range(num_control):
    #     plt.plot(x, y[j, :], linewidth=2, label='controller ' + str(j + 1), linestyle=line[j % 10])
    #     # marker=marker_list[j % 5],
    #     #          markevery=round(points / 10), markersize=marker_size_list[j % 5], linestyle=line[j // 10])
    # # plt.show()
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # # plt.legend()
    # plt.tight_layout()
    # plt.savefig(output_fig)
