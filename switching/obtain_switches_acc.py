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


class SwitchesAcc:
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
                               min_energy=-1, max_control=None, min_control=None):
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

        if max_control:
            self.max_control = max_control
        else:
            self.max_control = len(self.H_c)
        if min_control:
            self.min_control = min_control
        else:
            self.min_control = 0

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

    def _compute_cur_cost(self, k, cur_obj):
        s = [j for j in range(len(self.H_c))]
        cur_cost = {}
        for n in range(self.min_control, self.max_control + 1):
            for i in combinations(s, n):
                cur_list = list(i)
                cur_list.sort()
                cur_cost[tuple(cur_list)] = self._compute_obj_acc(k, tuple(cur_list)) - cur_obj
        return cur_cost

    def _compute_tv(self, ctrl_1, ctrl_2):
        diff = list(set(ctrl_1) - set(ctrl_2)) + list(set(ctrl_2) - set(ctrl_1))
        return len(diff)

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
        cur_control = updated_control.copy()
        prev_switch = 0
        self.cost = np.zeros_like(initial_control)
        tv_norm = 0
        for k in range(initial_control.shape[0]):
            # cur_control = updated_control.copy()
            # for j in range(initial_control.shape[1]):
            #     cur_control[k, :] = 0
            #     cur_control[k, j] = 1
            #     obj = self._compute_obj_acc(k, j)
            #     self.cost[k, j] = obj - updated_obj
            cur_cost = self._compute_cur_cost(k, updated_obj)
            if k == 0:
                cur_j_hat = min(cur_cost, key=cur_cost.get)
            else:
                min_cost_j = min(cur_cost, key=cur_cost.get)
                min_cost = cur_cost.get(min_cost_j)
                prev_cost = cur_cost[cur_j_hat]
                tv_increase = self._compute_tv(min_cost_j, cur_j_hat)
                if min_cost != prev_cost and prev_cost + updated_obj > alpha * (tv_norm + tv_increase):
                    cur_j_hat = min_cost_j
                    tv_norm += tv_increase
            self.j_hat[k] = cur_j_hat
            # print(cur_j_hat)
            if k == 0:
                self.hsequence.append(list(cur_j_hat))
            elif self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(list(cur_j_hat))
            updated_control[k, :] = 0
            for singlej in cur_j_hat:
                updated_control[k, singlej] = 1
            updated_obj = cur_cost[cur_j_hat] + updated_obj
            self._update_cur_state(cur_j_hat)
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
        cur_j_hat = list(np.where(self.deviation[0, :] >= 0.5)[0])
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        self.diff_obj = 0
        prev_switch = 0
        tv_norm = 0

        binary_control = np.zeros_like(initial_control)
        # print(self.deviation[0, :], initial_control[0, :])
        for k in range(1, initial_control.shape[0]):
            self.deviation[k, :] = self.deviation[k - 1, :] + initial_control[k, :]
            for singlej in cur_j_hat:
                self.deviation[k, int(singlej)] = self.deviation[k, int(singlej)] - 1
            # print(self.deviation[k, :], initial_control[k, :])
            self.diff_obj = max(self.diff_obj, max(abs(self.deviation[k, :] - initial_control[k, :])))
            # print(self.diff_obj)

            update_deviation = self.deviation[k, :].copy()
            update_deviation[cur_j_hat] -= 1

            # min_deviation = min(self.deviation[k, :])
            # max_deviation = max(self.deviation[k, :])
            # if min_deviation < 0 and abs(max_deviation) > abs(min_deviation):
            test_diff = max(max(abs(update_deviation)), self.diff_obj) * self.delta_t
            select_ctrl = list(np.where(self.deviation[k, :] >= 0.5)[0])
            tv_change = self._compute_tv(cur_j_hat, select_ctrl)
            if test_diff > alpha * (tv_norm + tv_change):
                # print(test_diff, k, tv_norm)
                cur_j_hat = select_ctrl
                tv_norm += tv_change
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)

            for singlej in cur_j_hat:
                binary_control[k, singlej] = 1

        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        last_diff = self.deviation[initial_control.shape[0] - 1, :].copy()
        last_diff[cur_j_hat] -= 1
        self.diff_obj = max(self.diff_obj, max(abs(last_diff))) * self.delta_t
        self.tv_obj = tv_norm
        self.overall_obj = self.diff_obj + alpha * self.tv_obj
        # print("total deviation is", self.diff_obj * self.delta_t)
        # print("total objective value is", self.diff_obj * self.delta_t + tv_norm)

        # print("binary obj", self._compute_objective_value(binary_control, -1))
        return self.tau0, self.num_switches, self.hsequence

    def _initialize_matrix_exp(self, control_amps):
        # self.initial_props = []
        # self.diag = []
        # self.decompose_matrix = []
        # for hc in self.H_c:
        #     s, v = np.linalg.eigh(hc + self.H_d)
        #     self.diag.append(s)
        #     self.decompose_matrix.append(v)
        #     self.initial_props.append(
        #         np.dot(v.dot(np.diag(np.exp(-1j * s * self.evo_time / self.n_ts))), v.conj().T))

        self.initial_props = {}
        self.diag = {}
        self.decompose_matrix = {}
        group = [j for j in range(len(self.H_c))]
        self.ctrl_comb = []
        for n in range(len(self.H_c) + 1):
            for i in combinations(group, n):
                cur_list = list(i)
                cur_list.sort()
                s, v = np.linalg.eigh(sum(self.H_c[j] for j in cur_list) + self.H_d)
                self.diag[tuple(cur_list)] = s
                self.decompose_matrix[tuple(cur_list)] = v
                self.initial_props[tuple(cur_list)] = np.dot(
                    v.dot(np.diag(np.exp(-1j * s * self.evo_time / self.n_ts))),
                    v.conj().T)
                self.ctrl_comb.append(cur_list)
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
        # self._cur_start = self.initial_props[select_j].dot(self._cur_start)
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
