import sys
import time
import numpy as np
from scipy.linalg import expm
from qutip import Qobj
import qutip.control.pulseoptim as cpo
import matplotlib.pyplot as plt
import gurobipy as gb
from itertools import combinations

sys.path.append("..")
from utils.draw_figure import draw_switchtime
from utils.auxiliary_energy import *

from utils.auxiliary_molecule import generate_molecule_func


def obtain_controller(control_arr, thre_max=0.65, thre_min=0.1, seed=None):
    if np.max(control_arr) >= thre_max:
        return np.argmax(control_arr)
    else:
        pool = []
        prob = []
        for j in range(len(control_arr)):
            if control_arr[j] >= thre_min:
                pool.append(j)
                prob.append(control_arr[j])
        prob = np.array(prob) / sum(prob)
        if seed:
            np.random.seed(seed)
        return np.random.choice(pool, 1, p=prob)


class SwitchesAccSOS1:
    def __init__(self, initial_control_file=None, delta_t=0.05):
        self.initial_control_file = initial_control_file
        self.delta_t = delta_t

        self.j_hat = None
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []

        self.obj_type = None
        self.H_d = None
        self.H_c = None
        self.X_0 = None
        self.X_targ = None
        self.n_ts = None
        self.evo_time = None
        self.optim = None
        self._into = None
        self._onto = None

        self.n_ctrl = None
        self.phase_option = None
        self.time_limit = 3600

    def init_gradient_computer(self, H_d, H_c, X_0, X_targ, n_ts, evo_time, obj_type='fid', phase_option="SU",
                               min_energy=-1):
        self.H_d = H_d
        self.H_c = H_c
        self.X_0 = X_0
        self.X_targ = X_targ
        self.n_ts = n_ts
        self.j_hat = np.zeros(n_ts)
        self.evo_time = evo_time
        self.delta_t = self.evo_time / self.n_ts
        self.obj_type = obj_type
        self.phase_option = phase_option
        if self.obj_type == 'fid':
            # self.optim = cpo.create_pulse_optimizer(Qobj(self.H_d), [Qobj(hc) for hc in self.H_c],
            #                                         Qobj(self.X_0), Qobj(self.X_targ), self.n_ts, self.evo_time,
            #                                         amp_lbound=0, amp_ubound=1, dyn_type='UNIT',
            #                                         phase_option=self.phase_option,
            #                                         init_pulse_params={"offset": 0}, gen_stats=True,
            #                                         init_pulse_type="ZERO")
            self.n_ctrl = len(self.H_c)

        if self.obj_type in ['energyratio', 'energy']:
            self.n_ctrl = len(self.H_c) - 1
            self._into = None
            self._onto = None
            self.min_energy = min_energy

    def obtain_switches(self, mode='random', param=None):
        if mode == 'obj':
            return self._obtain_switches_reduction(param)
        if mode == 'sur':
            return self._obtain_switches_sur(param)
        if mode == 'binary':
            return self._obtain_switches_binary(param)

    def _obtain_switches_binary(self, param):
        binary_control = np.loadtxt(param['binary_file'], delimiter=',')
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(binary_control.shape[0])
        self.j_hat[0] = np.argmax(binary_control[0, :])
        self.hsequence.append(int(self.j_hat[0]))
        # if self.obj_type == 'fid':
        #     dyn = self.optim.dynamics
        #     dyn.initialize_controls(initial_control)
        self._initialize_matrix_exp(initial_control)
        self.continuous_obj = self._initialize_continuous_state()
        prev_switch = 0
        for k in range(1, binary_control.shape[0]):
            self.j_hat[k] = np.argmax(binary_control[k, :])
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(int(self.j_hat[k]))
        self.tau0.append(binary_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _choose_control_obj(self, cost):
        s = [j for j in range(len(H_c))]
        for n in range(len(H_c) + 1):
            for i in combinations(s, n):
                cur_list = list(i)
                cur_list.sort()
        if self.min_control == 1 and self.max_control == 1:
            cur_j_hat = [np.argmin(cost)]
        else:
            cur_j_hat = []
            sorted_idx = np.argsort(cost)
            sorted_cost = cost[sorted_idx]
            num_select = 0
            for sortj in range(len(cost)):
                if num_select >= self.max_control:
                    break
                elif sorted_cost[sortj] < 0:
                    cur_j_hat.append(sorted_idx[sortj])
                    num_select += 1
                elif num_select < self.min_control:
                    cur_j_hat.append(sorted_idx[sortj])
                    num_select += 1
                else:
                    break
        return cur_j_hat


    def _obtain_switches_reduction(self, param):
        alpha = param['alpha']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        num_ctrl = initial_control.shape[1]
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = [0] * initial_control.shape[0]
        self._initialize_matrix_exp(initial_control)
        self.continuous_obj = self._initialize_continuous_state()
        updated_obj = self.continuous_obj
        updated_control = initial_control.copy()
        prev_switch = 0
        self.cost = np.zeros_like(initial_control)
        tv_norm = 0
        for k in range(initial_control.shape[0]):
            for j in range(initial_control.shape[1]):
                obj = self._compute_obj_acc(k, j)
                self.cost[k, j] = obj - updated_obj
            if k == 0:
                cur_j_hat = [np.argmin(self.cost[k, :])]
            else:
                cur_cost = self.cost[k, cur_j_hat[0]]
                min_cost = min(self.cost[k, :])
                if min_cost != cur_cost and cur_cost + updated_obj > alpha * (tv_norm + 2):
                    cur_j_hat = [int(np.argmin(self.cost[k, :]))]
                    tv_norm += 2

                # cur_j_hat = np.argmin(self.cost[k, :])
            self.j_hat[k] = cur_j_hat
            if k == 0:
                self.hsequence.append(cur_j_hat)
            elif self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
            updated_control[k, :] = 0
            for singlej in cur_j_hat:
                updated_control[k, singlej] = 1
                updated_obj = self.cost[k, singlej] + updated_obj
            self._update_cur_state(cur_j_hat[0])
        # switches.binary_obj = updated_obj
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_sur(self, param):
        alpha = param['alpha']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self._initialize_matrix_exp(initial_control)
        self.continuous_obj = self._initialize_continuous_state()
        num_ctrl = initial_control.shape[1]
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = [0] * initial_control.shape[0]
        self.deviation = np.zeros_like(initial_control)
        self.deviation[0, :] = initial_control[0, :]
        cur_j_hat = np.argmax(self.deviation[0, :])
        self.j_hat[0] = cur_j_hat
        self.hsequence.append([cur_j_hat])
        self.diff_obj = 0
        prev_switch = 0
        tv_norm = 0

        binary_control = np.zeros_like(initial_control)

        for k in range(1, initial_control.shape[0]):
            self.deviation[k, :] = self.deviation[k - 1, :] + initial_control[k, :]
            self.deviation[k, int(cur_j_hat)] = self.deviation[k, int(cur_j_hat)] - 1
            self.diff_obj = max(self.diff_obj, max(abs(self.deviation[k, :] - initial_control[k, :])))

            cur_deviation = self.deviation[k, int(cur_j_hat)]
            # min_deviation = min(self.deviation[k, :])
            max_deviation = max(self.deviation[k, :])
            # if num_ctrl == 2:
            #     # if max_deviation > thre_ratio:
            #     cur_j_hat = int(np.argmax(self.deviation[k, :]))
            #     # pass
            # else:
            #     # if abs(cur_deviation - max_deviation) / abs(max_deviation - min_deviation) > thre_ratio:
            #     #     # if cur_deviation / max_deviation < 1 - thre_ratio:
            #     #     cur_j_hat = int(np.argmax(self.deviation[k, :]))
            #     if max_deviation != cur_deviation and abs(cur_deviation * self.delta_t) > tv_norm:
            #         cur_j_hat = int(np.argmax(self.deviation[k, :]))
            #         tv_norm += 0.002
            if max_deviation != cur_deviation:
                max_abs = max(abs(self.deviation[k, :]))
                test_diff = max(max_abs, abs(cur_deviation - 1), self.diff_obj) * self.delta_t
                if test_diff > alpha * (tv_norm + 2):
                    # print(test_diff, k, tv_norm)
                    cur_j_hat = int(np.argmax(self.deviation[k, :]))
                    tv_norm += 2
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append([cur_j_hat])

            binary_control[k, cur_j_hat] = 1

        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        last_diff = self.deviation[initial_control.shape[0] - 1, :].copy()
        last_diff[int(cur_j_hat)] -= 1
        self.diff_obj = max(self.diff_obj, max(abs(last_diff))) * self.delta_t
        self.tv_obj = tv_norm
        self.overall_obj = self.diff_obj + alpha * self.tv_obj
        # print("total deviation is", self.diff_obj * self.delta_t)
        # print("total objective value is", self.diff_obj * self.delta_t + tv_norm)

        # print("binary obj", self._compute_objective_value(binary_control, -1))
        return self.tau0, self.num_switches, self.hsequence

    def _initialize_matrix_exp(self, control_amps):
        self.initial_props = []
        self.diag = []
        self.decompose_matrix = []
        for hc in self.H_c:
            s, v = np.linalg.eigh(hc + self.H_d)
            self.diag.append(s)
            self.decompose_matrix.append(v)
            self.initial_props.append(
                np.dot(v.dot(np.diag(np.exp(-1j * s * self.evo_time / self.n_ts))), v.conj().T))

        self.c_matrix_exp = []
        for k in range(self.n_ts):
            self.c_matrix_exp.append(expm(-1j * (self.H_d + sum(control_amps[k, j] * self.H_c[j]
                                                                for j in range(control_amps.shape[1]))) * self.delta_t))

    def get_decomposition(self):
        return self.diag, self.decompose_matrix

    def _initialize_continuous_state(self):
        self._cur_start = self.X_0
        self._final_list = [np.identity(len(self.X_0))]
        for k in range(self.n_ts):
            _final = self._final_list[-1].dot(self.c_matrix_exp[self.n_ts - 1 - k])
            self._final_list.append(_final)
        if self.obj_type == 'fid':
            den = np.linalg.matrix_rank(self.X_targ)
            return 1 - np.abs(np.trace(
                self.X_targ.conj().T.dot(self._final_list[self.n_ts]).dot(self._cur_start))) / den
        if self.obj_type == 'energyratio':
            final_state = self._final_list[self.n_ts].dot(self._cur_start)
            return 1 - np.real(final_state.conj().T.dot(self.H_c[1]).dot(final_state)) / self.min_energy
        if self.obj_type == 'energy':
            final_state = self._final_list[self.n_ts].dot(self._cur_start)
            return np.real(final_state.conj().T.dot(self.H_c[1]).dot(final_state)) - self.min_energy

    def _update_cur_state(self, select_j):
        self._cur_start = self.initial_props[select_j].dot(self._cur_start)

    def _compute_obj_acc(self, tidx, jidx):
        if self.obj_type == 'fid':
            den = np.linalg.matrix_rank(self.X_targ)
            obj = 1 - np.abs(np.trace(self.X_targ.conj().T.dot(self._final_list[self.n_ts - 1 - tidx]).dot(
                self.initial_props[jidx]).dot(self._cur_start))) / den
        else:
            final_state = self._final_list[self.n_ts - 1 - tidx].dot(
                self.initial_props[jidx]).dot(self._cur_start)
            obj = np.real(final_state.conj().T.dot(self.H_c[1]).dot(final_state))
            if self.obj_type == 'energyratio':
                obj = 1 - obj / self.min_energy
            elif self.obj_type == 'energy':
                obj = obj - self.min_energy
        return obj

    def _update_forward(self, update_idx, select_ctrl, compute_obj=False):
        if update_idx == -1:
            self._into = [None] * (self.n_ts + 1)
            self._into[0] = self.X_0
        else:
            fwd = self.initial_props[select_ctrl].dot(self._into[update_idx])
            self._into[update_idx + 1] = fwd
        for k in range(update_idx + 1, self.n_ts):
            fwd = self.c_matrix_exp[k].dot(self._into[k])
            self._into[k + 1] = fwd
        if compute_obj:
            if self.obj_type == 'energyratio':
                return 1 - self._into[-1].conj().T.dot(self.H_c[1]).dot(self._into[-1]) / self.min_energy
            if self.obj_type == 'energy':
                return self._into[-1].conj().T.dot(self.H_c[1]).dot(self._into[-1]) - self.min_energy
            if self.obj_type == 'fid':
                den = np.linalg.matrix_rank(self.X_targ)
                return 1 - np.abs(np.trace(self.X_targ.conj().T.dot(self._into[-1]))) / den

    def _update_backward(self, update_idx, select_ctrl_list):
        if update_idx == -1:
            self._onto = []
            if self.obj_type == 'fid':
                self._onto.append(self.X_targ.conj().T)
            if self.obj_type == 'energyratio':
                self._onto.append(self._into[-1].conj().T.dot(self.H_c[1]))
            for k in range(self.n_ts):
                bwd = self._onto[k].dot(self.c_matrix_exp[self.n_ts - k - 1])
                self._onto.append(bwd)
        else:
            bwd = self._onto[self.n_ts - update_idx - 1].dot(self.initial_props[select_ctrl_list[update_idx]])
            self._onto[self.n_ts - update_idx] = bwd
            for k in range(self.n_ts - update_idx, self.n_ts):
                bwd = self._onto[k].dot(self.initial_props[select_ctrl_list[self.n_ts - k - 1]])
                self._onto[k + 1] = bwd
            # print(bwd[-1])
            # exit()

    def _compute_gradient(self, update_idx, select_ctrl_list, compute_obj=True):
        obj = self._update_forward(update_idx, select_ctrl_list[update_idx], compute_obj)
        self._update_backward(update_idx, select_ctrl_list)
        if self.obj_type == 'fid':
            temp_grad = np.zeros((self.n_ts, self.n_ctrl), dtype=complex)
            for t in range(self.n_ts):
                for j in range(self.n_ctrl):
                    grad_temp = -1j * self.delta_t * self.H_c[j]
                    g = np.trace(self._onto[self.n_ts - t - 1].dot(grad_temp).dot(self._into[t + 1]))
                    temp_grad[t, j] = g
            fid_pre = np.trace(self.X_targ.conj().T.dot(self._into[-1]))
            grad = - np.real(temp_grad * np.exp(-1j * np.angle(fid_pre)) / np.linalg.matrix_rank(self.X_targ))
        if self.obj_type == 'energyratio':
            grad = np.zeros((self.n_ts, 2))
            for k in range(self.n_ts):
                grad[k, 0] = -np.imag(self._onto[self.n_ts - k - 1].dot((-self.H_c[0]).dot(self._into[k + 1]))
                                      * self.delta_t)
                grad[k, 1] = -np.imag(self._onto[self.n_ts - k - 1].dot((-self.H_c[1]).dot(self._into[k + 1]))
                                      * self.delta_t)
            grad *= 2
            grad *= - 1 / self.min_energy
            # grad = np.expand_dims(np.array(grad), 1)
        if compute_obj:
            return grad, obj
        else:
            return grad

    def draw_extracted_control(self, fig_name):
        # if self.obj_type == "energy":
        #     n_ctrl = 2
        # else:
        #     n_ctrl = self.n_ctrl
        # control = np.zeros((self.n_ts, n_ctrl))
        #
        # for t in range(self.n_ts):
        #     control[t, int(self.j_hat[t])] = 1
        #
        # t = np.linspace(0, self.evo_time, self.n_ts + 1)
        # plt.figure(dpi=300)
        # plt.xlabel("Time")
        # plt.ylabel("Control amplitude")
        # marker_list = ['-o', '--^', '-*', '--s']
        # marker_size_list = [5, 5, 8, 5]
        # marker_step = int(self.n_ts / 5)
        # for j in range(n_ctrl):
        #     # if max(control[:, j]) > 0:
        #     plt.step(t, np.hstack((control[:, j], control[-1, j])), marker_list[j % 4],
        #              where='post', linewidth=2, label='controller ' + str(j + 1), markevery=(j, marker_step),
        #              markersize=marker_size_list[j % 4])
        # # plt.legend()
        # plt.legend(bbox_to_anchor=(1, 0.5), loc='center right', prop={'size': 6}, borderpad=0.5)
        # plt.savefig(fig_name)
        if self.obj_type in ["energy", "energyratio"]:
            n_ctrl = 2
        else:
            n_ctrl = self.n_ctrl

        draw_switchtime(self.evo_time, self.n_ts, n_ctrl, self.hsequence, self.tau0, fig_name, True)

    def draw_metric(self, fig_name, metric="gradient"):
        if self.obj_type == "energy":
            n_ctrl = 2
        else:
            n_ctrl = self.n_ctrl

        # draw the metric
        t = np.linspace(0, self.evo_time, self.n_ts + 1)
        plt.figure(dpi=300)
        plt.xlabel("Time")

        if metric == "sur":
            metric_data = self.deviation
            plt.ylabel('deviation')
        if metric == "obj":
            metric_data = self.cost
            plt.ylabel('change of objective value')

        # plt.ylabel("Metric")
        marker_list = ['-o', '--^', '-*', '--s']
        marker_size_list = [5, 5, 8, 5]
        marker_step = int(self.n_ts / 10)
        linestyle_list = ['-', '--', 'dashdot']
        for j in range(n_ctrl):
            # if max(control[:, j]) > 0:
            # plt.step(t, np.hstack((metric_data[:, j], metric_data[-1, j])), marker_list[j % 4],
            #          where='post', linewidth=2, label='controller ' + str(j + 1), markevery=(j, marker_step),
            #          markersize=marker_size_list[j % 4])
            plt.step(t, np.hstack((metric_data[:, j], metric_data[-1, j])), linestyle=linestyle_list[j % 3],
                     where='post', linewidth=2, label='controller ' + str(j + 1))
        # plt.legend(bbox_to_anchor=(1, 0.5), loc='center right', prop={'size': 6}, borderpad=0.5)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.savefig(fig_name)


if __name__ == '__main__':
    # n = 4
    # num_edges = 2
    # seed = 1
    # n_ts = 40
    # evo_time = 2
    # initial_type = "warm"
    # alpha = 0.01
    # min_up_time = 0
    # name = "EnergySTTR"
    # step = 1
    #
    # lb_threshold = 0.1
    # ub_threshold = 0.65
    #
    # if seed == 0:
    #     Jij, edges = generate_Jij_MC(n, num_edges, 100)
    #
    # else:
    #     Jij = generate_Jij(n, seed)
    #
    # C = get_ham(n, True, Jij)
    # B = get_ham(n, False, Jij)
    #
    # y0 = uniform(n)
    #
    # initial_control = "../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
    # switches = Switches(initial_control, delta_t=evo_time / n_ts)
    # warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches('naive')
    # print(ctrl_hamil_idx)
    # warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches('random', lb_threshold, ub_threshold)
    # print(ctrl_hamil_idx)
    # switches.init_gradient_computer(None, [B, C], y0[0:2 ** n], None, n_ts, evo_time, 'energy')
    # warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches('gradient')
    # print(ctrl_hamil_idx)
    d = 2
    qubit_num = 6
    molecule = "BeH2"
    target = "../result/control/Target/MoleculeVQE_BeH2_evotime20.0_n_ts200_target.csv"
    initial_type = "ave"
    Hops, H0, U0, U = generate_molecule_func(qubit_num, d, molecule)

    if target is not None:
        U = np.loadtxt(target, dtype=np.complex_, delimiter=',')
    else:
        print("Please provide the target file!")
        exit()

    n_ts = 400
    evo_time = 20

    step = 1
    alpha = 0.001
    min_up_time = 0

    name = "FakeLA"

    # The control Hamiltonians (Qobj classes)
    H_c = [Qobj(hops) for hops in Hops]
    # Drift Hamiltonian
    H_d = Qobj(H0)
    # start point for the gate evolution
    X_0 = Qobj(U0)
    # Target for the gate evolution
    X_targ = Qobj(U)

    initial_control = "../result/control/ADMM/MoleculeADMM_BeH2_evotime20.0_n_ts200_ptypeWARM_offset0.5_sum_penalty0.01_penalty0.001_ADMM_3.0_iter100.csv"
    if initial_control is None:
        print("Must provide control results of ADMM!")
        exit()

    lb_threshold = 0.1
    ub_threshold = 0.65
    threshold_ratio = 0

    switches = SwitchesAcc(initial_control, delta_t=evo_time / n_ts)
    start1 = time.time()
    switches.init_gradient_computer(H0, Hops, U0, U, n_ts, evo_time, 'fid')
    warm_start_length, num_switch, ctrl_hamil_idx = switches.obtain_switches('la', thre_ratio=threshold_ratio)
    end1 = time.time()
    print(ctrl_hamil_idx)

    fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}_extraction".format(
        name + "_" + molecule, str(evo_time), str(n_ts), str(num_switch), initial_type,
        str(min_up_time), threshold_ratio) + ".png"
    switches.draw_extracted_control(fig_name)

    fig_name = "../result/figure/SwitchTime/" + "{}_evotime_{}_n_ts{}_n_switch{}_init{}_minuptime{}_thre{}_metric".format(
        name + "_" + molecule, str(evo_time), str(n_ts), str(num_switch), initial_type,
        str(min_up_time), threshold_ratio) + ".png"
    switches.draw_metric(fig_name, "la")
