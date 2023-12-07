import numpy as np
from itertools import combinations
from qutip import identity, sigmax, sigmaz, sigmay, tensor
from qutip.qip.operations.gates import cnot
from utils.evolution import time_evolution, compute_obj_fid


def extend_control(initial_control_name, pre_num_step, cur_num_step):
    new_control_name = initial_control_name.split(".csv")[0] + "_extend_ts_" + str(cur_num_step) + ".csv"
    initial_control = np.loadtxt(initial_control_name, delimiter=",")
    step = int(cur_num_step / pre_num_step)

    cur_control = np.zeros((cur_num_step, initial_control.shape[1]))

    for i in range(pre_num_step):
        for j in range(step):
            cur_control[i*step+j, :] = initial_control[i, :]

    np.savetxt(new_control_name, cur_control, delimiter=",")


def add_sos1(original_control):
    n_ctrl = original_control.shape[1]
    n_ts = original_control.shape[0]

    # controller 0, controller 1, controller 2, controller 1 + controller 2
    new_control = np.zeros((n_ts, 2 ** n_ctrl))

    comb_to_idx = {}
    idx = 1
    s = [j for j in range(n_ctrl)]
    for n in range(1, n_ctrl + 1):
        for i in combinations(s, n):
            cur_list = list(i)
            cur_list.sort()
            # key = ""
            # for num in cur_list:
            #     key = key + str(num)
            #     key = key +","
            comb_to_idx[str(cur_list)] = idx
            idx += 1

    # print(comb_to_idx)

    for k in range(n_ts):
        comb_list = []
        # print(np.argsort(original_control[k, :]))
        order_control_idx = np.argsort(original_control[k, :])[::-1][:n_ctrl]
        # print(order_control_idx)
        for j in range(n_ctrl - 1):
            comb_list.append(order_control_idx[j])
            comb_list.sort()
            idx = comb_to_idx[str(comb_list)]
            new_control[k, idx] = original_control[k, order_control_idx[j]] - original_control[
                k, order_control_idx[j + 1]]
        comb_list.append(order_control_idx[n_ctrl - 1])
        comb_list.sort()
        idx = comb_to_idx[str(comb_list)]
        new_control[k, idx] = original_control[k, order_control_idx[n_ctrl - 1]]
        new_control[k, 0] = 1 - original_control[k, order_control_idx[0]]
    return new_control


def add_sos1_new(original_control):
    n_ctrl = original_control.shape[1]
    n_ts = original_control.shape[0]
    new_control = np.zeros((n_ts, n_ctrl + 1))
    for k in range(n_ts):
        for idx in range(n_ctrl):
            new_control[k, idx + 1] = original_control[k, idx] / n_ctrl
        new_control[k, 0] = 1 - sum(new_control[k, 1:])
    return new_control


def remove_sos1(original_control):
    n_ctrl_sos1 = original_control.shape[1]
    n_ctrl = int(np.log2(n_ctrl_sos1))
    n_ts = original_control.shape[0]

    # controller 0, controller 1, controller 2, controller 1 + controller 2
    new_control = np.zeros((n_ts, n_ctrl))

    # print(original_control)

    idx_to_comb = {}
    idx = 1
    s = [j for j in range(n_ctrl)]
    for n in range(1, n_ctrl + 1):
        for i in combinations(s, n):
            cur_list = list(i)
            cur_list.sort()
            # key = ""
            # for num in cur_list:
            #     key = key + str(num)
            #     key = key +","
            # comb_to_idx[str(cur_list)] = idx
            idx_to_comb[idx] = cur_list
            idx += 1

    for k in range(n_ts):
        for idx in range(1, n_ctrl_sos1):
            cur_list = idx_to_comb[idx]
            for j in cur_list:
                new_control[k, j] += original_control[k, idx]

    return new_control


def remove_sos1_new(original_control):
    n_ctrl = original_control.shape[1] - 1
    n_ts = original_control.shape[0]
    new_control = np.zeros((n_ts, n_ctrl))
    for k in range(n_ts):
        for idx in range(n_ctrl):
            new_control[k, idx] = original_control[k, idx + 1] * n_ctrl
    return new_control


def add_sos1_by_file(initial_control_name):
    split1 = initial_control_name.split("/")
    split2 = split1[-1].split("_")
    new_instance_name = split2[0] + "SOS1"
    split2.pop(0)
    cur_str_1 = new_instance_name
    for s in split2:
        cur_str_1 = cur_str_1 + "_" + s
    cur_str_2 = ""
    for s in split1[:-1]:
        cur_str_2 = cur_str_2 + s + "/"
    new_control_name = cur_str_2 + cur_str_1
    # print(new_control_name)

    original_control = np.loadtxt(initial_control_name, delimiter=",")
    new_control = add_sos1(original_control)
    # print(original_control == remove_sos1(new_control))
    np.savetxt(new_control_name, new_control, delimiter=",")

    # h_d_mat = (tensor(sigmax(), sigmax()) + tensor(sigmay(), sigmay()) + tensor(sigmaz(), sigmaz())).full()
    # h_c_mat = [tensor(sigmax(), identity(2)).full(), tensor(sigmay(), identity(2)).full()]
    # X_0 = identity(4)
    # X_targ = cnot()
    #
    # q = 2
    #
    # original_result = time_evolution(
    #     h_d_mat, h_c_mat, 400, 20, original_control, X_0.full(), False, 1)
    # print(compute_obj_fid(X_targ, original_result, phase="PSU"))
    #
    # # print(X_targ)
    # print(new_control)
    #
    # new_h_c_mat = [(identity(4) * 0).full()]
    # for n in range(1, n_ctrl + 1):
    #     for i in combinations(s, n):
    #         cur_list = list(i)
    #         new_h_c_mat.append(sum(h_c_mat[j] for j in cur_list))
    # # print(h_c_mat)
    # # print(new_h_c_mat)
    #
    # print(X_0.full())
    # new_result = time_evolution(
    #     h_d_mat, new_h_c_mat, 400, 20, new_control, X_0.full(), False, 1)
    #
    # # print(X_targ)
    # print(compute_obj_fid(X_targ, new_result, phase="PSU"))


def add_sos1_by_file_new(initial_control_name):
    split1 = initial_control_name.split("/")
    split2 = split1[-1].split("_")
    new_instance_name = split2[0] + "NewSOS1"
    split2.pop(0)
    cur_str_1 = new_instance_name
    for s in split2:
        cur_str_1 = cur_str_1 + "_" + s
    cur_str_2 = ""
    for s in split1[:-1]:
        cur_str_2 = cur_str_2 + s + "/"
    new_control_name = cur_str_2 + cur_str_1
    # print(new_control_name)

    original_control = np.loadtxt(initial_control_name, delimiter=",")
    new_control = add_sos1_new(original_control)
    # print(original_control == remove_sos1(new_control))
    np.savetxt(new_control_name, new_control, delimiter=",")


if __name__ == '__main__':
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/CNOT_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv")
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/CNOT_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT.csv")
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/CNOT_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT.csv")
    #
    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/CNOT_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv")
    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/CNOT_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv")
    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/CNOT_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT_alpha0.0001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv")

    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/NOTleak_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv")
    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/NOTleak_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv")
    # add_sos1_by_file_new(
    #     "../new_result/b_control/TR+MT/NOTleak_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv")
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/NOTleak_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT.csv")
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/NOTleak_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT.csv")
    # add_sos1_by_file_new(
    #     "../new_result/c_control/Continuous/NOTleak_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv")
    extend_control("../new_result/c_control/Continuous/Energy4_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5_instance1.csv", 40, 400000)

