import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.optimize import Bounds
from scipy.optimize import minimize
from qutip import identity, sigmax, sigmaz, sigmay, tensor
from qutip.qip.operations.gates import cnot
from qutip import Qobj

sys.path.append("..")
from utils.draw_figure import draw_switchtime
from utils.auxiliary_energy import *
from utils.evolution import compute_TV_norm, compute_obj_by_switch
from utils.auxiliary_molecule import generate_molecule_func
from switching.obtain_switches_acc import SwitchesAcc


class SwitchTimeOptAcc:
    """
    class for optimization with switching time points
    """

    def __init__(self):

        self.hlist = None  # list of all control Hamiltonians
        self.ctrl_hamil_idx = None  # control sequence of the Hamiltonians

        self.tau = None  # length of time for each interval (variables to be optimized)
        self.tau0 = None  # initialize the solution
        self.switch_time = None  # switching time points

        self.x0 = None  # initial state
        self.xtarg = None  # target state

        self.evotime = None  # total evolution time
        self.num_switch = None  # number of switches

        self.forward = None
        self.backward = None
        self.hamil_expm = None
        self.final_state = None  # final state
        self.obj = None  # optimal value

        self.control = None  # descritized control results

        self.time_lb = 0
        self.time_ub = self.evotime

        self.obj_type = None
        self.phase = None

        self.lb_epsilon = 1e-04
        self.tvnorm = 0

        self.diag = []
        self.decompose_matrix = []

        self.pre_compute_time = 0
        self.evaluate_time = 0
        self.gradient_time = 0

        self.cur_x = None

    def build_optimizer(self, hlist, hd, ctrl_hamil_idx, init, x0, xtarg, evotime, num_switch, time_lb=None, time_ub=None,
                        obj_type='fid', phase="SU", decompose=None):

        self.hlist = hlist
        self.hd = hd
        if isinstance(ctrl_hamil_idx[0], int):
            ctrl_hamil_idx = [[ctrl] for ctrl in ctrl_hamil_idx]
        self.ctrl_hamil_idx = ctrl_hamil_idx
        self.hamil_list = []
        for ctrl in ctrl_hamil_idx:
            cur_hamil = sum(self.hlist[singlectrl] for singlectrl in ctrl) + self.hd
            self.hamil_list.append(cur_hamil)
        self.tau0 = init
        self.x0 = x0
        self.xtarg = xtarg
        self.evotime = evotime
        self.num_switch = num_switch
        self.phase = phase

        if time_lb:
            self.time_lb = time_lb
        if time_ub:
            self.time_ub = time_ub
        else:
            self.time_ub = self.evotime

        self.obj_type = obj_type

        if decompose is None:
            start = time.time()
            self.diag = {}
            self.decompose_matrix = {}
            decomposed_set = set()
            # for hidx in self.ctrl_hamil_idx:
            #     if hidx not in decomposed_set:
            #         s, v = np.linalg.eigh(self.hlist[hidx])
            #         self.diag[hidx], self.decompose_matrix[hidx] = s, v
            #         decomposed_set.add(hidx)
            for (ctrlidx, hamil) in zip(self.ctrl_hamil_idx, self.hamil_list):
                tuple_ctrlidx = tuple(ctrlidx)
                if tuple_ctrlidx not in decomposed_set:
                    s, v = np.linalg.eigh(hamil)
                    self.diag[tuple_ctrlidx], self.decompose_matrix[tuple_ctrlidx] = s, v
                    decomposed_set.add(tuple_ctrlidx)

            self.pre_compute_time = time.time() - start
        else:
            self.diag = dict(enumerate(decompose['diag']))
            self.decompose_matrix = dict(enumerate(decompose['matrix']))

    def obtain_state(self, x):
        """
            conduct time evolution to obtain the state
            :param x: parameters beta and gamma
            :param hamil_idx: index list of Hamiltonians
            :return: final state
            """
        # conduct time evolution
        self.forward = [self.x0]
        self.hamil_expm = []
        for k in range(self.num_switch + 1):
            # self.hamil_expm.append(expm(-1j * self.hlist[self.ctrl_hamil_idx[k]].copy() * x[k]))
            idx = self.ctrl_hamil_idx[k]
            s, v = self.diag[tuple(idx)], self.decompose_matrix[tuple(idx)]
            self.hamil_expm.append(v.dot(np.diag(np.exp(-1j * s * x[k]))).dot(v.conj().T))
            cur_state = self.hamil_expm[k].dot(self.forward[k])
            self.forward.append(cur_state)
        final_state = self.forward[-1]
        return final_state

    def objective(self, x):
        """
        compute the objective value
        :param x: parameters beta and gamma
        :param hamil_idx: index list of Hamiltonians
        :return: objective value
        """
        # conduct time evolution
        start = time.time()
        final_state = self.obtain_state(x)
        if self.obj_type == 'fid':
            den = np.linalg.matrix_rank(self.xtarg)
            if self.phase == "PSU":
                fid = np.abs(np.trace(self.xtarg.conj().T.dot(final_state))) / den
            if self.phase == "SU":
                fid = np.real(np.trace(self.xtarg.conj().T.dot(final_state))) / den
            obj = 1 - fid
        if self.obj_type == "energy":
            obj = np.real(final_state.conj().T.dot(self.hlist[1].dot(final_state)))
        self.evaluate_time += time.time() - start
        return obj

    def gradient(self, x):
        if not (self.cur_x == x).all():
            self.obtain_state(x)
            self.cur_x = x
        start = time.time()
        if self.obj_type == 'energy':
            self.backward = [self.forward[-1].conj().T.dot(self.hlist[1].conj().T)]
        if self.obj_type == 'fid':
            self.backward = [self.xtarg.conj().T]
        for k in range(self.num_switch + 1):
            bwd = self.backward[k].dot(self.hamil_expm[self.num_switch - k])
            self.backward.append(bwd)
        grad = []
        if self.obj_type == 'energy':
            for k in range(self.num_switch + 1):
                # grad += [np.imag(self.backward[self.num_switch - k].dot(
                #     self.hlist[self.ctrl_hamil_idx[k]].copy().dot(self.forward[k + 1])))]
                # grad += [np.imag(self.backward[self.num_switch - k].dot(
                #     self.hlist[self.ctrl_hamil_idx[k]].copy().dot(self.forward[k + 1]))) * 2]
                grad += [np.imag(self.backward[self.num_switch - k].dot(
                    self.hamil_list[k].copy().dot(self.forward[k + 1]))) * 2]
        if self.obj_type == 'fid':
            pre_grad = np.zeros(self.num_switch + 1, dtype=complex)
            den = np.linalg.matrix_rank(self.xtarg)
            for k in range(self.num_switch + 1):
                # grad_temp = expm_frechet(-1j * H[t] * delta_t, -1j * self.H_c[j] * delta_t, compute_expm=False)
                # grad_temp = -1j * self.hlist[self.ctrl_hamil_idx[k]].copy()
                grad_temp = -1j * self.hamil_list[k].copy()
                pre_grad[k] = np.trace(self.backward[self.num_switch - k].dot(grad_temp).dot(self.forward[k + 1]))
            if self.phase == "PSU":
                fid_pre = np.trace(self.xtarg.conj().T.dot(self.forward[-1]))
                grad = - np.real(pre_grad * np.exp(-1j * np.angle(fid_pre)) / den)
            if self.phase == "SU":
                grad = - np.real(pre_grad / den)
        self.gradient_time += time.time() - start
        return grad

    def optimize(self):
        """
        optimize the length of each control interval
        :return: result of optimization
        """
        # predefine equality constraints for parameters
        eq_cons = {'type': 'eq',
                   'fun': lambda x: sum(x) - self.evotime}

        # initialize the solution
        tau0 = self.tau0.copy()
        self.cur_x = tau0
        # set the bounds of variables
        bounds = Bounds([self.time_lb] * self.num_switch + [0], [self.time_ub] * (self.num_switch + 1))
        # minimize the objective function
        # res = minimize(self.objective, x0, method='SLSQP', constraints=eq_cons, bounds=bounds, options={'ftol': 1e-06})
        res = minimize(self.objective, tau0, method='SLSQP', constraints=eq_cons, bounds=bounds,
                       jac=self.gradient, options={'ftol': 1e-06})
        self.tau = res.x
        # retrieve the switching time points
        self.retrieve_switching_points()
        # compute the final state and optimal value
        self.final_state = self.obtain_state(self.tau)
        self.obj = res.fun

        # print(compute_obj_by_switch(self.hlist, self.hd, self.tau, self.ctrl_hamil_idx, self.x0, self.xtarg, 'fid', 'PSU'))

        return res

    def optimize_from_sur(self):
        """
        optimize the length of each control interval where the length each control interval can be zero
        :return: result of optimization
        """
        # predefine equality constraints for parameters
        eq_cons = {'type': 'eq',
                   'fun': lambda x: sum(x) - self.evotime}
        cons = [eq_cons]
        # for k in range(self.num_switch + 1):
        #     cons.append({'type': 'ineq', 'fun': lambda x: x[k] ** 2 - self.time_lb * x[k]})
        cons.append({'type': 'ineq', 'fun': lambda x: np.array([x[k] ** 2 - self.time_lb * x[k]
                                                                for k in range(self.num_switch + 1)])})

        # initialize the solution
        x0 = self.tau0.copy()
        # set the bounds of variables
        # bounds = Bounds([self.time_lb] * (self.num_switch + 1), [self.time_ub] * (self.num_switch + 1))
        bounds = Bounds([0] * (self.num_switch + 1), [self.time_ub] * (self.num_switch + 1))
        # minimize the objective function
        res = minimize(self.objective, x0, method='SLSQP', constraints=cons, bounds=bounds)
        self.tau = res.x
        # retrieve the switching time points
        self.retrieve_switching_points()
        # compute the final state and optimal value
        self.final_state = self.obtain_state(self.tau)
        self.obj = res.fun

        return res

    def retrieve_switching_points(self):
        """
        retrieve switching points from the length of control intervals
        """
        # retrieve switching time points
        if self.num_switch > 0:
            self.switch_time = np.zeros(self.num_switch)
            self.switch_time[0] = self.tau[0]
            for k in range(1, self.num_switch):
                self.switch_time[k] = self.switch_time[k - 1] + self.tau[k]

    def filter_control(self):
        ts = np.count_nonzero((self.tau > self.lb_epsilon))
        control = np.zeros((ts, len(self.hlist)))
        filtered_tau = self.tau.copy()
        filtered_tau[filtered_tau <= self.lb_epsilon] = 0
        non_zero_k = 0
        for k in range(len(self.tau)):
            if self.tau[k] > self.lb_epsilon:
                control[non_zero_k, self.ctrl_hamil_idx[k]] = 1
                non_zero_k += 1
        filtered_obj = self.objective(filtered_tau)
        self.filtered_control = control
        return control, filtered_obj

    def retrieve_control(self, num_time_step):
        """
        retrieve control results from solutions of QAOA
        :param num_time_step: number of time steps
        :return: control results
        """
        control = np.zeros((num_time_step, len(self.hlist)))
        delta_t = self.evotime / num_time_step

        cur_time_l = 0
        cur_time_r = 0
        for k in range(self.num_switch):
            cur_time_r = self.switch_time[k]
            for time_step in range(int(cur_time_l / delta_t), min(int(cur_time_r / delta_t), num_time_step)):
                for singlectrl in self.ctrl_hamil_idx[k]:
                    control[time_step, singlectrl] = 1
                # control[time_step, self.ctrl_hamil_idx[k]] = 1
            cur_time_l = cur_time_r

        for time_step in range(int(cur_time_l / delta_t), num_time_step):
            for singlectrl in self.ctrl_hamil_idx[self.num_switch]:
            # control[time_step, self.ctrl_hamil_idx[self.num_switch]] = 1
                control[time_step, singlectrl] = 1

        self.control = control

        return control

    def draw_control(self, fig_name, marker_step=250):
        """
        draw control results
        """
        # print(self.evotime, len(self.hlist), self.ctrl_hamil_idx, self.switch_time, fig_name)
        draw_switchtime(self.evotime, 1000, len(self.hlist), self.ctrl_hamil_idx, self.tau, fig_name)


    def tv_norm(self):
        """
        compute the tv norm of control results
        :return: value of tv norm
        """
        return sum(sum(abs(self.filtered_control[tstep + 1, j] - self.filtered_control[tstep, j])
                       for tstep in range(self.filtered_control.shape[0] - 1))
                   for j in range(self.filtered_control.shape[1]))
        # return 2 * np.count_nonzero((self.tau > self.lb_epsilon)) - 2


def obtain_switching_time(initial_control_file, delta_t=0.05):
    initial_control = np.loadtxt(initial_control_file, delimiter=",")
    num_switches = 0
    hsequence = []
    tau0 = []
    j_hat = np.zeros(initial_control.shape[0])
    j_hat[0] = np.argmax(initial_control[0, :])
    hsequence.append(int(j_hat[0]))
    prev_switch = 0
    for k in range(1, initial_control.shape[0]):
        j_hat[k] = np.argmax(initial_control[k, :])
        if j_hat[k] != j_hat[k - 1]:
            num_switches += 1
            tau0.append(k * delta_t - prev_switch)
            prev_switch = k * delta_t
            hsequence.append(int(j_hat[k]))
    tau0.append(initial_control.shape[0] * delta_t - prev_switch)
    return tau0, num_switches, hsequence


def obtain_switching_time_from_sur(initial_control_file, delta_t=0.05):
    initial_control = np.loadtxt(initial_control_file, delimiter=",")
    num_switches = initial_control.shape[0] - 1
    hsequence = []
    j_hat = np.zeros(initial_control.shape[0])
    for k in range(initial_control.shape[0]):
        j_hat[k] = np.argmax(initial_control[k, :])
        hsequence.append(int(j_hat[k]))
    tau0 = [delta_t] * initial_control.shape[0]
    return tau0, num_switches, hsequence


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


def obtain_switching_time_uncertain(initial_control_file, delta_t=0.05, step=1, thre_max=0.65, thre_min=0.1, seed=None):
    # compute average initial control
    initial_control_pre = np.loadtxt(initial_control_file, delimiter=",")
    n_ts = initial_control_pre.shape[0]
    initial_control = np.zeros((int(np.ceil(n_ts / step)), initial_control_pre.shape[1]))
    for k in range(0, initial_control_pre.shape[0], step):
        initial_control[int(k / step), :] = sum(initial_control_pre[k + l, :] for l in range(step)) / step

    num_switches = 0
    hsequence = []
    tau0 = []
    j_hat = np.zeros(initial_control.shape[0])
    j_hat[0] = obtain_controller(initial_control[0, :], thre_max, thre_min, seed)
    hsequence.append(int(j_hat[0]))
    prev_switch = 0
    for k in range(1, initial_control.shape[0]):
        j_hat[k] = obtain_controller(initial_control[k, :], thre_max, thre_min, seed)
        if j_hat[k] != j_hat[k - 1]:
            num_switches += 1
            tau0.append(k * delta_t - prev_switch)
            prev_switch = k * delta_t
            hsequence.append(int(j_hat[k]))
    tau0.append(initial_control.shape[0] * delta_t - prev_switch)
    print(j_hat)
    return tau0, num_switches, hsequence


def reduce_switch(max_switch, interval_len, c_sequence, epsilon=1e-5):
    num_switch = len(c_sequence) - 1
    interval_len = np.array(interval_len)
    c_sequence = np.array(c_sequence)

    # interval_len = np.array([0.15, 0.45, 0.3, 0.05, 0.5, 0.05, 0.4, 0.1])
    # c_sequence = np.array([1,0,2,0,1,0,1,0])

    while num_switch > max_switch:
        min_interval = np.where(abs(interval_len - interval_len.min()) < epsilon)[0]
        dif_switch = int(num_switch - max_switch)
        f_min_interval = min_interval
        if dif_switch < len(min_interval):
            f_min_interval = np.sort(np.random.choice(min_interval, dif_switch,
                                                      p=[1/len(min_interval)] * len(min_interval), replace=False))

        # delete the control
        for idx in range(len(f_min_interval)):
            k = f_min_interval[idx]
            if k - idx == 0:
                interval_len[k - idx + 1] += interval_len[k - idx]
            else:
                interval_len[k - idx - 1] += interval_len[k - idx]
            c_sequence = np.delete(c_sequence, k - idx)
            interval_len = np.delete(interval_len, k - idx)

        # merge the controls
        origin_num_switch = len(c_sequence)
        term_l = []
        for k in range(0, origin_num_switch):
            term_l.append(origin_num_switch)
            for l in range(k + 1, origin_num_switch):
                if c_sequence[l] == c_sequence[k]:
                    interval_len[k] += interval_len[l]
                else:
                    term_l[k] = l
                    break
        c_sequence = np.delete(c_sequence, [control for control in range(1, term_l[0])])
        interval_len = np.delete(interval_len, [control for control in range(1, term_l[0])])
        delete_ctrl = term_l[0] - 1
        for k in range(1, origin_num_switch):
            if term_l[k] != term_l[k - 1]:
                c_sequence = np.delete(c_sequence, [control for control in range(
                    k + 1 - delete_ctrl, term_l[k] - delete_ctrl)])
                interval_len = np.delete(interval_len, [control for control in range(
                    k + 1 - delete_ctrl, term_l[k] - delete_ctrl)])
                delete_ctrl += term_l[k] - k - 1

        num_switch = len(c_sequence) - 1

    return interval_len, num_switch, c_sequence




