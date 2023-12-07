import sys
import time
import numpy as np
from scipy.linalg import expm
from qutip import Qobj
import qutip.control.pulseoptim as cpo
import matplotlib.pyplot as plt
import gurobipy as gb

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


class Switches:
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

    def init_gradient_computer(self, H_d, H_c, X_0, X_targ, n_ts, evo_time, obj_type='fid', phase_option="SU", min_energy=-1):
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
            self.optim = cpo.create_pulse_optimizer(Qobj(self.H_d), [Qobj(hc) for hc in self.H_c],
                                                    Qobj(self.X_0), Qobj(self.X_targ), self.n_ts, self.evo_time,
                                                    amp_lbound=0, amp_ubound=1, dyn_type='UNIT',
                                                    phase_option=self.phase_option,
                                                    init_pulse_params={"offset": 0}, gen_stats=True,
                                                    init_pulse_type="ZERO")
            self.n_ctrl = len(self.H_c)

        if self.obj_type == 'energy':
            self.n_ctrl = len(self.H_c) - 1
            self._into = None
            self._onto = None
            self.min_energy = min_energy

    def obtain_switches(self, mode='random', param=None):
        if mode == 'naive':
            return self._obtain_switches_naive()
        if mode == 'random':
            return self._obtain_switches_random(param)
        if mode == 'gradient':
            return self._obtain_switches_gradient()
        if mode == 'gradientnu':
            return self._obtain_switches_gradient_noupdate()
        if mode == 'lanu':
            return self._obtain_switches_linear_approximation_noupdate()
        if mode == 'la':
            return self._obtain_switches_linear_approximation(param)
        if mode == 'multi':
            return self._obtain_switches_without_sos1()
        if mode == 'obj':
            return self._obtain_switches_reduction(param)
        if mode == 'sur':
            return self._obtain_switches_sur(param)
        if mode == 'tv':
            # return self._obtain_switch_mip(alpha, warm_start)
            # return self._obtain_switch_absolute_diff(alpha)
            return self._obtain_switch_milp_heuristic(param)
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
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
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

    def _obtain_switches_without_sos1(self, thre_min = 0):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        num_controllers = initial_control.shape[1]
        self.tau0 = []
        self.num_switches = [0] * num_controllers
        self.hstart = [0] * num_controllers
        cur_state = [0] * num_controllers
        for j in range(num_controllers):
            # self.tau0[j] = []
            prev_switch = 0
            if initial_control[0, j] >= thre_min:
                cur_state[j] = 1
                self.hstart[j] = 1
            for k in range(1, initial_control.shape[0]):
                if initial_control[k, j] >= thre_min:
                    state = 1
                else:
                    state = 0
                if state != cur_state[j]:
                    self.tau0.append(k * self.delta_t - prev_switch)
                    prev_switch = k * self.delta_t
                    self.num_switches[j] += 1
            self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hstart

    def _obtain_switches_naive(self):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        self.j_hat[0] = np.argmax(initial_control[0, :])
        self.hsequence.append(int(self.j_hat[0]))
        prev_switch = 0
        for k in range(1, initial_control.shape[0]):
            self.j_hat[k] = np.argmax(initial_control[k, :])
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(int(self.j_hat[k]))
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_random(self, param):
        thre_min, thre_max, seed = param['thre_min'], param['thre_max'], param['seed']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        self.j_hat[0] = obtain_controller(initial_control[0, :], thre_max, thre_min, seed)
        self.hsequence.append(int(self.j_hat[0]))
        prev_switch = 0
        for k in range(1, initial_control.shape[0]):
            self.j_hat[k] = obtain_controller(initial_control[k, :], thre_max, thre_min, seed)
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(int(self.j_hat[k]))
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_gradient(self):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        updated_control = initial_control.copy()
        self.grad = np.zeros_like(initial_control)
        # compute gradient
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(updated_control)
        grad = self._compute_gradient(updated_control, -1)
        self.grad[0, :] = grad[0, :].copy()
        # cur_j_hat = int(np.argmin(grad[0, :]))
        cur_j_hat = int(np.argmax(initial_control[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        updated_control[0, :] = 0
        updated_control[0, cur_j_hat] = 1
        prev_switch = 0
        for k in range(1, initial_control.shape[0]):
            grad = self._compute_gradient(updated_control, update_idx=k - 1)
            self.grad[k, :] = grad[k, :]
            cur_j_hat = int(np.argmin(grad[k, :]))
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
            updated_control[k, :] = 0
            updated_control[k, cur_j_hat] = 1
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_gradient_noupdate(self):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        updated_control = initial_control.copy()
        # compute gradient
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(updated_control)
        grad = self._compute_gradient(updated_control, -1)
        cur_j_hat = int(np.argmin(grad[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        prev_switch = 0
        for k in range(1, initial_control.shape[0]):
            cur_j_hat = int(np.argmin(grad[k, :]))
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_linear_approximation_noupdate(self):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        updated_control = initial_control.copy()
        # compute gradient
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(updated_control)
        grad = self._compute_gradient(updated_control, -1)
        cur_j_hat = int(np.argmin(grad[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        prev_switch = 0
        for k in range(1, initial_control.shape[0]):
            cur_la = np.zeros_like(grad[k, :])
            for j in range(len(cur_la)):
                cur_la[j] = grad[k, j] - sum(grad[k, :] * updated_control[k, :])
            cur_j_hat = int(np.argmin(cur_la))
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_linear_approximation(self, param):
        alpha = param['alpha']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        num_ctrl = initial_control.shape[1]
        # num_control = initial_control.shape[1]
        # initial_control = np.zeros((400, num_control))
        # # # nonzero = [4, 14, 7, 9, 1, 3, 5, 11]
        # nonzero = [4, 10, 7, 9, 1, 3]
        # # nonzero = [2,4,5,8,10,12]
        # for idx in nonzero:
        #     initial_control[:, idx] += 1 / len(nonzero)
        # initial_control[:, 4] -= 1 / (2 * len(nonzero) - 1)
        # initial_control[0, 4] = 1
        # initial_control[1:, 3] = 1

        # print(initial_control)

        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        updated_control = initial_control.copy()
        self.la = np.zeros_like(initial_control)
        # compute gradient
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)

        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        print(self.continuous_obj)

        grad = self._compute_gradient(updated_control, -1)
        obj = self.continuous_obj
        for j in range(initial_control.shape[1]):
            self.la[0, j] = grad[0, j] - sum(grad[0, :] * updated_control[0, :])
        cur_j_hat = int(np.argmax(initial_control[0, :]))
        # cur_j_hat = int(np.argmin(self.la[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        updated_control[0, :] = 0
        updated_control[0, cur_j_hat] = 1
        prev_switch = 0
        tv_norm = 0
        for k in range(1, initial_control.shape[0]):
            grad = self._compute_gradient(updated_control, update_idx=k - 1)
            exit()
            obj = self._compute_objective_value(updated_control, update_idx=k - 1)
            # # cur_la = np.zeros_like(grad[k, :])
            # for j in range(initial_control.shape[1]):
            #     self.la[k, j] = grad[k, j] - sum(grad[k, :] * updated_control[k, :])
            #     # update tv norm
            #     if j != cur_j_hat:
            #         self.la[k, j] += alpha * 2
            # cur_la = self.la[k, int(cur_j_hat)]
            # min_la = min(self.la[k, :])
            # max_la = max(self.la[k, :])
            # if num_ctrl == 2:
            #     # if min_la < thre_ratio * obj:
            #     cur_j_hat = int(np.argmin(self.la[k, :]))
            # else:
            #     # if min_la < - thre_ratio * obj:
            #     # cur_j_hat = np.argmin(self.cost[k, :])
            #     if abs(cur_la - min_la) / abs(max_la - min_la) > thre_ratio:
            #         # cur_j_hat = int(np.argmin(self.la[k, :]))
            #         cand_cur_j_hat = np.where(self.la[k, :] == np.min(self.la[k, :]))[0]
            #         cur_j_hat = -1
            #         # for j_hat in cand_cur_j_hat:
            #         #     if j_hat not in self.hsequence:
            #         #         cur_j_hat = j_hat
            #         #         break
            #         if cur_j_hat == -1:
            #             cur_j_hat = np.argmin(self.la[k, :])
            for j in range(initial_control.shape[1]):
                self.la[k, j] = grad[k, j] - sum(grad[k, :] * updated_control[k, :])

            cur_cost = self.la[k, int(cur_j_hat)] + obj
            min_cost = min(self.la[k, :]) + obj
            if min_cost != cur_cost and cur_cost > alpha * (tv_norm + 2):
                # if min_cost != cur_cost and min_cost + updated_obj <= tv_norm + 0.002:
                cur_j_hat = int(np.argmin(self.la[k, :]))
                tv_norm += 2

            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
            updated_control[k, :] = 0
            updated_control[k, cur_j_hat] = 1
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        self.binary_obj = self._compute_objective_value(updated_control, -1)
        print(self._compute_objective_value(updated_control, -1))
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_reduction(self, param):
        alpha = param['alpha']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        num_ctrl = initial_control.shape[1]
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        updated_obj = self._compute_objective_value(initial_control, -1)
        updated_control = initial_control.copy()
        cur_control = updated_control.copy()
        prev_switch = 0
        self.cost = np.zeros_like(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        tv_norm = 0
        for k in range(initial_control.shape[0]):
            cur_control = updated_control.copy()
            for j in range(initial_control.shape[1]):
                cur_control[k, :] = 0
                cur_control[k, j] = 1
                obj = self._compute_objective_value(cur_control, update_idx=k-1)
                self.cost[k, j] = obj - updated_obj
            if k == 0:
                cur_j_hat = np.argmin(self.cost[k, :])
            else:
                cur_cost = self.cost[k, int(cur_j_hat)]
                min_cost = min(self.cost[k, :])
                # if num_ctrl == 2:
                #     # if min_cost < thre_ratio * obj:
                #     cur_j_hat = int(np.argmin(self.cost[k, :]))
                #     # pass
                # else:
                #     # if min_cost < thre_ratio * obj:
                #     # if abs(cur_cost - min_cost) / abs(max_cost - min_cost) > thre_ratio:
                #     #     cur_j_hat = int(np.argmin(self.cost[k, :]))
                #     # if cur_cost + updated_obj - c_obj > 10 / initial_control.shape[0]:
                #     #     cur_j_hat = int(np.argmin(self.cost[k, :]))
                #     if min_cost != cur_cost and cur_cost + updated_obj > tv_norm:
                #         # if min_cost != cur_cost and min_cost + updated_obj <= tv_norm + 0.002:
                #         cur_j_hat = int(np.argmin(self.cost[k, :]))
                #         tv_norm += 0.002
                if min_cost != cur_cost and cur_cost + updated_obj > alpha * (tv_norm + 2):
                    cur_j_hat = int(np.argmin(self.cost[k, :]))
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
            updated_control[k, cur_j_hat] = 1
            updated_obj = self.cost[k, cur_j_hat] + updated_obj
        # switches.binary_obj = updated_obj
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switches_sur(self, param):
        alpha = param['alpha']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        num_ctrl = initial_control.shape[1]
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(initial_control.shape[0])
        self.deviation = np.zeros_like(initial_control)
        self.deviation[0, :] = initial_control[0, :]
        cur_j_hat = np.argmax(self.deviation[0, :])
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        self.diff_obj = 0
        prev_switch = 0
        tv_norm = 0

        binary_control = np.zeros_like(initial_control)

        for k in range(1, initial_control.shape[0]):
            self.deviation[k, :] = self.deviation[k - 1, :] + initial_control[k, :]
            self.deviation[k, int(cur_j_hat)] = self.deviation[k, int(cur_j_hat)] - 1
            self.diff_obj = max(self.diff_obj, max(abs(self.deviation[k, :] - initial_control[k, :])))

            cur_deviation = self.deviation[k, int(cur_j_hat)]
            min_deviation = min(self.deviation[k, :])
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
            if max_deviation != cur_deviation and abs(cur_deviation * self.delta_t) > alpha * (tv_norm + 2):
                cur_j_hat = int(np.argmax(self.deviation[k, :]))
                tv_norm += 2
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)

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

    def _obtain_switch_mip(self, alpha=0.0, warm_start=True, warm_start_tol=0.1):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        num_ctrl = initial_control.shape[1]
        round = gb.Model()
        bin = round.addVars(self.n_ts, num_ctrl, vtype=gb.GRB.BINARY)
        # bin_val = self.b_bin.copy()
        up_diff = round.addVar()
        tv_ub = round.addVars(self.n_ts - 1, num_ctrl)

        round.addConstrs(
            gb.quicksum(initial_control[t, j] - bin[t, j] for t in range(k)) * self.delta_t + up_diff >= 0
            for j in range(num_ctrl) for k in range(1, self.n_ts + 1))
        round.addConstrs(
            gb.quicksum(initial_control[t, j] - bin[t, j] for t in range(k)) * self.delta_t - up_diff <= 0
            for j in range(num_ctrl) for k in range(1, self.n_ts + 1))

        round.addConstrs(bin[t, j] - bin[t + 1, j] + tv_ub[t, j] >= 0 for t in range(self.n_ts - 1)
                         for j in range(num_ctrl))
        round.addConstrs(bin[t, j] - bin[t + 1, j] - tv_ub[t, j] <= 0 for t in range(self.n_ts - 1)
                         for j in range(num_ctrl))

        round.addConstrs(gb.quicksum(bin[t, j] for j in range(num_ctrl)) == 1
                         for t in range(self.n_ts))
        round.setObjective(up_diff + alpha * gb.quicksum(gb.quicksum(tv_ub[t, j] for j in range(num_ctrl))
                                                         for t in range(self.n_ts - 1)))
        self.time_limit = 300
        round.Params.TimeLimit = self.time_limit
        self.mip_gap = 0.9
        self.mip_gap = 1e-3
        round.Params.MIPGap = self.mip_gap
        self.num_groups = 0

        initial_binary = np.zeros((self.n_ts, num_ctrl))
        if warm_start == 1:
            lb = 0
            ub = self.evo_time + alpha * 2 * self.n_ts
            ub_diff = self.evo_time
            lb_diff = 0
            threshold = self.evo_time
            while ub_diff - lb_diff > warm_start_tol:
                for j in range(num_ctrl):
                    cum_diff = initial_control[0, :].copy()
                    cum_diff[j] -= 1
                    if cum_diff[j] >= -threshold:
                        cum_deviation = initial_control[0, :].copy()
                        max_diff = np.max(cum_diff)
                        cur_binary = np.zeros((self.n_ts, num_ctrl))
                        cur_binary[0, j] = 1
                        cur_j = j
                        tv_norm = 0
                        for k in range(1, self.n_ts):
                            cum_diff += initial_control[k, :].copy()
                            cum_deviation += initial_control[k, :].copy()
                            if cum_diff[cur_j] - 1 < -threshold or np.max(cum_deviation) > threshold:
                                previous_j = cur_j
                                cur_j = int(np.argmax(cum_deviation))
                                if previous_j != cur_j:
                                    tv_norm += 2
                            cur_binary[k, cur_j] = 1
                            cum_diff[cur_j] -= 1
                            cum_deviation[cur_j] -= 1
                            cur_diff = np.max(np.abs(cum_diff)) * self.delta_t
                            if cur_diff > max_diff:
                                max_diff = cur_diff
                        if max_diff + alpha * tv_norm < ub:
                            print(max_diff + alpha * tv_norm, ub)
                            ub = max_diff + alpha * tv_norm
                            initial_binary = cur_binary.copy()
                        if max_diff < ub_diff:
                            ub_diff = max_diff
                            threshold = ub_diff - 0.5 * (ub_diff - lb_diff)
                    elif j == num_ctrl - 1:
                        lb_diff = threshold
                        threshold = lb_diff + 0.5 * (ub_diff - lb_diff)

            for k in range(self.n_ts):
                for j in range(num_ctrl):
                    bin[k, j].Start = initial_binary[k, j]

        round.optimize()
        self.overall_obj = round.ObjVal
        self.diff_obj = up_diff.x
        self.tv_obj = sum(sum(tv_ub[k, j].x for k in range(self.n_ts - 1)) for j in range(num_ctrl))

        bin_val = np.zeros((self.n_ts, num_ctrl))
        for j in range(num_ctrl):
            for k in range(self.n_ts):
                bin_val[k, j] = bin[k, j].x
        # bin_val = round.getAttr('X', bin)
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(self.n_ts)
        # for j in range(num_ctrl):
        #     if bin_val[0, j] == 1:
        #         cur_j_hat = j
        #         break
        cur_j_hat = int(np.argmax(bin_val[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        prev_switch = 0
        for k in range(1, self.n_ts):
            cur_j_hat = int(np.argmax(bin_val[k, :]))
            # for j in range(num_ctrl):
            #     if bin_val[k, j] == 1:
            #         cur_j_hat = j
            #         break
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        self.binary_obj = self._compute_objective_value(bin_val, -1)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switch_milp_heuristic(self, param):
        alpha = param['alpha']
        num_groups = param['ngroup']
        time_limit = param['tl']
        mip_gap = param['mipgap']
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        num_ctrl = initial_control.shape[1]
        cur_n_ts = 0
        self.overall_obj = 0
        self.diff_obj = 0
        self.tv_obj = 0
        self.mip_gap = mip_gap
        bin_val = np.zeros((self.n_ts, num_ctrl))
        self.num_groups = num_groups
        self.time_limit = time_limit / num_groups
        for group in range(num_groups):
            cur_alpha = alpha * (num_groups - group)
            sub_n_ts = int(self.n_ts / num_groups)
            if group == num_groups - 1:
                sub_n_ts = self.n_ts - cur_n_ts

            round = gb.Model()
            bin = round.addVars(sub_n_ts, num_ctrl, vtype=gb.GRB.BINARY)
            # bin_val = self.b_bin.copy()
            up_diff = round.addVar()
            tv_ub = round.addVars(sub_n_ts, num_ctrl)

            round.addConstrs((gb.quicksum(initial_control[t, j] - bin_val[t, j] for t in range(cur_n_ts))
                              + gb.quicksum(initial_control[t + cur_n_ts, j] - bin[t, j] for t in range(k)))
                             * self.delta_t + up_diff >= 0 for j in range(num_ctrl) for k in range(1, sub_n_ts + 1))
            round.addConstrs((gb.quicksum(initial_control[t, j] - bin_val[t, j] for t in range(cur_n_ts))
                              + gb.quicksum(initial_control[t + cur_n_ts, j] - bin[t, j] for t in range(k)))
                             * self.delta_t - up_diff <= 0 for j in range(num_ctrl) for k in range(1, sub_n_ts + 1))
            if group > 0:
                round.addConstrs(bin[0, j] - bin_val[cur_n_ts - 1, j] + tv_ub[0, j] >= 0 for j in range(num_ctrl))
                round.addConstrs(bin[0, j] - bin_val[cur_n_ts - 1, j] - tv_ub[0, j] <= 0 for j in range(num_ctrl))
            round.addConstrs(bin[t - 1, j] - bin[t, j] + tv_ub[t, j] >= 0 for t in range(1, sub_n_ts)
                             for j in range(num_ctrl))
            round.addConstrs(bin[t - 1, j] - bin[t, j] - tv_ub[t, j] <= 0 for t in range(1, sub_n_ts)
                             for j in range(num_ctrl))

            round.addConstrs(gb.quicksum(bin[t, j] for j in range(num_ctrl)) == 1 for t in range(sub_n_ts))
            round.setObjective(up_diff + cur_alpha * gb.quicksum(gb.quicksum(
                tv_ub[t, j] for j in range(num_ctrl)) for t in range(sub_n_ts)) + cur_alpha * self.tv_obj)
            round.Params.TimeLimit = self.time_limit
            round.Params.MIPGap = self.mip_gap
            round.optimize()

            for j in range(num_ctrl):
                for k in range(sub_n_ts):
                    bin_val[k + cur_n_ts, j] = bin[k, j].x

            cur_n_ts += sub_n_ts
            self.overall_obj = round.ObjVal
            self.diff_obj = up_diff.x
            self.tv_obj += sum(sum(tv_ub[k, j].x for k in range(sub_n_ts)) for j in range(num_ctrl))

        self.num_switches = 0

        self.tau0 = []
        self.hsequence = []
        self.j_hat = np.zeros(self.n_ts)
        # for j in range(num_ctrl):
        #     if bin_val[0, j] == 1:
        #         cur_j_hat = j
        #         break
        cur_j_hat = int(np.argmax(bin_val[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        prev_switch = 0
        for k in range(1, self.n_ts):
            cur_j_hat = int(np.argmax(bin_val[k, :]))
            # for j in range(num_ctrl):
            #     if bin_val[k, j] == 1:
            #         cur_j_hat = j
            #         break
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        # self.binary_obj = self._compute_objective_value(bin_val, -1)
        return self.tau0, self.num_switches, self.hsequence

    def _obtain_switch_absolute_diff(self, alpha=0.0):
        initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
        if self.obj_type == 'fid':
            dyn = self.optim.dynamics
            dyn.initialize_controls(initial_control)
        self.continuous_obj = self._compute_objective_value(initial_control, -1)
        num_ctrl = initial_control.shape[1]
        round = gb.Model()
        bin = round.addVars(self.n_ts, num_ctrl, vtype=gb.GRB.BINARY)
        # bin_val = self.b_bin.copy()
        up_diff = round.addVars(self.n_ts)
        tv_ub = round.addVars(self.n_ts - 1, num_ctrl)

        round.addConstrs((initial_control[k, j] - bin[k, j]) * self.delta_t + up_diff[k] >= 0
                         for j in range(num_ctrl) for k in range(self.n_ts))
        round.addConstrs((initial_control[k, j] - bin[k, j]) * self.delta_t - up_diff[k] <= 0
                         for j in range(num_ctrl) for k in range(self.n_ts))

        round.addConstrs(bin[t, j] - bin[t + 1, j] + tv_ub[t, j] >= 0 for t in range(self.n_ts - 1)
                         for j in range(num_ctrl))
        round.addConstrs(bin[t, j] - bin[t + 1, j] - tv_ub[t, j] <= 0 for t in range(self.n_ts - 1)
                         for j in range(num_ctrl))

        round.addConstrs(gb.quicksum(bin[t, j] for j in range(num_ctrl)) == 1
                         for t in range(self.n_ts))
        round.setObjective(gb.quicksum(up_diff[k] for k in range(self.n_ts)) +
                           alpha * gb.quicksum(gb.quicksum(tv_ub[t, j] for j in range(num_ctrl))
                                               for t in range(self.n_ts - 1)))
        self.mip_gap = 1e-4
        round.optimize()
        self.overall_obj = round.ObjVal
        self.diff_obj = sum(up_diff[k].x for k in range(self.n_ts))
        self.tv_obj = sum(sum(tv_ub[k, j].x for k in range(self.n_ts - 1)) for j in range(num_ctrl))

        bin_val = np.zeros((self.n_ts, num_ctrl))
        for j in range(num_ctrl):
            for k in range(self.n_ts):
                bin_val[k, j] = bin[k, j].x
        # bin_val = round.getAttr('X', bin)
        self.num_switches = 0
        self.hsequence = []
        self.tau0 = []
        self.j_hat = np.zeros(self.n_ts)
        # for j in range(num_ctrl):
        #     if bin_val[0, j] == 1:
        #         cur_j_hat = j
        #         break
        cur_j_hat = int(np.argmax(bin_val[0, :]))
        self.j_hat[0] = cur_j_hat
        self.hsequence.append(cur_j_hat)
        prev_switch = 0
        for k in range(1, self.n_ts):
            cur_j_hat = int(np.argmax(bin_val[k, :]))
            # for j in range(num_ctrl):
            #     if bin_val[k, j] == 1:
            #         cur_j_hat = j
            #         break
            self.j_hat[k] = cur_j_hat
            if self.j_hat[k] != self.j_hat[k - 1]:
                self.num_switches += 1
                self.tau0.append(k * self.delta_t - prev_switch)
                prev_switch = k * self.delta_t
                self.hsequence.append(cur_j_hat)
        self.tau0.append(initial_control.shape[0] * self.delta_t - prev_switch)
        self.binary_obj = self._compute_objective_value(bin_val, -1)
        return self.tau0, self.num_switches, self.hsequence

    def _time_evolution_energy(self, control_amps, update_idx=-1):
        if update_idx == -1:
            self._into = [self.X_0]
            for k in range(self.n_ts):
                fwd = expm(-1j * (control_amps[k, 0] * self.H_c[0] +
                                  control_amps[k, 1] * self.H_c[1]) * self.delta_t).dot(self._into[k])
                self._into.append(fwd)
        else:
            for k in range(update_idx, self.n_ts):
                fwd = expm(-1j * (control_amps[k, 0] * self.H_c[0] +
                                  control_amps[k, 1] * self.H_c[1]) * self.delta_t).dot(self._into[k])
                self._into[k + 1] = fwd

    def _back_propagation_energy(self, control_amps, update_idx=-1):
        if update_idx == -1:
            self._onto = [self._into[-1].conj().T.dot(self.H_c[1].conj().T)]
            for k in range(self.n_ts):
                bwd = self._onto[k].dot(expm(-1j * (control_amps[self.n_ts - k - 1, 0] * self.H_c[0] +
                                                    control_amps[self.n_ts - k - 1, 1] * self.H_c[1]) * self.delta_t))
                self._onto.append(bwd)
        else:
            for k in range(self.n_ts - update_idx - 1, self.n_ts):
                bwd = self._onto[k].dot(expm(-1j * (control_amps[self.n_ts - k - 1, 0] * self.H_c[0]
                                                    + control_amps[self.n_ts - k - 1, 1] * self.H_c[1]) * self.delta_t))
                self._onto[k + 1] = bwd

    def _compute_gradient(self, control_amps, update_idx=-1):
        if self.obj_type == 'fid':
            grad = self.optim.fid_err_grad_compute(control_amps.reshape(-1))
            print(self.optim.dynamics.onto_evo)
        if self.obj_type == 'energy':
            # if not (self.u == control_amps).all():
            self._time_evolution_energy(control_amps, update_idx)
            self._back_propagation_energy(control_amps, update_idx)
            # self.u = control_amps
            grad = np.zeros_like(control_amps)
            for k in range(self.n_ts):
                grad[k, 0] = -np.imag(self._onto[self.n_ts - k - 1].dot((-self.H_c[0]).dot(self._into[k + 1]))
                                      * self.delta_t)
                grad[k, 1] = -np.imag(self._onto[self.n_ts - k - 1].dot((-self.H_c[1]).dot(self._into[k + 1]))
                                      * self.delta_t)
            grad *= 2
            grad *= - 1 / self.min_energy
            # grad = np.expand_dims(np.array(grad), 1)
        return grad

    def _compute_objective_value(self, control_amps, update_idx=-1):
        if self.obj_type == 'fid':
            obj = self.optim.fid_err_func_compute(control_amps.reshape(-1))
        if self.obj_type == 'energy':
            self._time_evolution_energy(control_amps, update_idx)
            obj = 1 - np.real(self._into[-1].conj().T.dot(self.H_c[1].dot(self._into[-1]))) / self.min_energy
        return obj

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
        if self.obj_type == "energy":
            n_ctrl = 2
        else:
            n_ctrl = self.n_ctrl

        draw_switchtime(self.evo_time, n_ctrl, self.n_ts, self.hsequence, self.tau0, fig_name, True)

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
            # initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
            # deviation = np.zeros((self.n_ts, n_ctrl))
            # deviation[0, :] = initial_control[0, :]
            # for k in range(1, self.n_ts):
            #     j_hat = np.argmax(deviation[k - 1, :])
            #     deviation[k, :] = deviation[k - 1, :] + initial_control[k, :]
            #     deviation[k, j_hat] = deviation[k - 1, j_hat] + initial_control[k, j_hat] - 1
            # metric_data = deviation
            metric_data = self.deviation
            plt.ylabel('deviation')
        if metric == "gradient":
            metric_data = -self.grad
            plt.ylabel('negative gradient')
        if metric == "la":
            metric_data = -self.la
            plt.ylabel('linear approximated change of objective value')
        if metric == "obj":
            metric_data = self.cost
            plt.ylabel('change of objective value')
        if metric == "naive":
            initial_control = np.loadtxt(self.initial_control_file, delimiter=",")
            metric_data = initial_control
            plt.ylabel('Control value')

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

    switches = Switches(initial_control, delta_t=evo_time / n_ts)
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
