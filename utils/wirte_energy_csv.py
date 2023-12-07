import os
import numpy as np


def read_energy_log(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    for line in lines:
        if "ratio objective" in line:
            obj = float(line.split(" ")[-1].strip("\n"))
            # print(obj)
        if "tv norm" in line:
            tv = int(float(line.split(" ")[-1].strip("\n")))
            # print(tv)
    return obj, tv


def read_general_log(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    converted_tv = -1
    for line in lines:
        if "fun:" in line:
            obj = float(line.split(" ")[-1].strip("\n"))
            # print(obj)
        if "tv norm" in line:
            tv = int(float(line.split(" ")[-1].strip("\n")))
        if "converted tv norm" in line:
            converted_tv = int(float(line.split(" ")[-1].strip("\n")))
            # print(tv)
        if "computational time of retrieving switches" in line:
            binary_time = float(line.split(" ")[-1].strip("\n"))
        if "computational time of optimization" in line:
            compute_time = float(line.split(" ")[-1].strip("\n"))
        if "nfev:" in line:
            num_evo = int(line.split(" ")[-1].strip("\n"))
    if converted_tv != -1:
        tv = converted_tv
    return obj, tv, binary_time, compute_time, num_evo


def write_csv(alpha_list, directory):
    min_energy = [-3.96254559434774, -6.29993971669791, -4.7572977523656, -5.30539590418686, -4.10829161988264]
    obj = np.zeros((len(alpha_list), 5))
    obj_diff = np.zeros((len(alpha_list), 5))
    tv = np.zeros((len(alpha_list), 5), dtype=int)
    print([(i,alpha) for (i,alpha) in enumerate(alpha_list)])
    for filename in os.listdir(directory):
        if "EnergyAccC6obj" in filename:
            instance = int(filename.split("_")[-3][-1])
            alpha = float(filename.split("_")[-1].split(".log")[0])
            # print(instance, alpha, filename)
            if alpha in alpha_list:
                # idx = alpha_list.index(alpha)
                idx = np.where(alpha_list == alpha)[0][0]
                obj[idx, instance - 1], tv[idx, instance - 1] = read_energy_log(
                    os.path.join(directory, filename))
                obj_diff[idx, instance - 1] = (1 - obj[idx, instance - 1]) * min_energy[instance - 1] - min_energy[instance - 1]
    # print(obj)
    return obj, tv
    # np.savetxt("energy6_nts100_obj.csv", obj, delimiter=',')
    # np.savetxt("energy6_nts100_obj_diff.csv", obj_diff, delimiter=',')
    # np.savetxt("energy6_nts100_tv.csv", tv, fmt='%d', delimiter=',')


def get_all_results(instances, directory, save_file=None):
    min_energy = [[-1], [-2.511350507398168, -1.0470843748257435, -2.719582936382, -3.3330147815732003, -2.1641371229504833],
                  [-3.96254559434774, -6.29993971669791, -4.7572977523656, -5.30539590418686, -4.10829161988264]]
    obj_arr = np.zeros((len(instances), 3))
    tv_arr = np.zeros((len(instances), 3))
    binary_time_arr = np.zeros((len(instances), 3))
    opt_time_arr = np.zeros((len(instances), 3))
    num_evo_arr = np.zeros((len(instances), 3))
    for filename in os.listdir(directory):
        if ".log" not in filename:
            continue
        obj, tv, binary_time, opt_time, num_evo = read_general_log(os.path.join(directory, filename))
        print(filename, obj, tv)
        terms = filename.split("_")
        instance_group = terms[0].split("Acc")[0]
        if instance_group == "Energy":
            qubit = terms[0].split("Acc")[1][1]
            instance_name = instance_group + qubit
            if "binary" in filename:
                group_number = int(terms[-4][-1])
            else:
                group_number = int(terms[-3][-1])
            obj = 1 - obj / min_energy[int(qubit) // 2 - 1][group_number - 1]
        elif instance_group in ["NOTleak", "CNOT"]:
            evo_time = int(float(terms[2]))
            if instance_group == "NOTleak":
                instance_name = "NOT" + str(evo_time)
            else:
                instance_name = "CNOT" + str(evo_time)
        else:
            instance_name = "Circuit" + terms[1]
        instance_idx = instances.index(instance_name)
        print(instance_name)
        method_idx = 0
        if "obj" in filename:
            method_idx = 1
        elif "sur" in filename:
            method_idx = 2
        obj_arr[instance_idx, method_idx] += obj
        tv_arr[instance_idx, method_idx] += tv
        binary_time_arr[instance_idx, method_idx] += binary_time
        opt_time_arr[instance_idx, method_idx] += opt_time
        num_evo_arr[instance_idx, method_idx] += num_evo
    obj_arr[1:3, :] /= 5
    tv_arr[1:3, :] /= 5
    binary_time_arr[1:3, :] /= 5
    opt_time_arr[1:3, :] /= 5
    num_evo_arr[1:3, :] /= 5

    # print(obj_arr, tv_arr)
    results = np.concatenate((obj_arr, tv_arr, binary_time_arr, opt_time_arr, num_evo_arr), axis=1)
    if save_file:
        np.savetxt(save_file, results, delimiter=',')
    else:
        print(results)
        print("time")
        print(binary_time_arr[2, :])
        print(opt_time_arr[2, :])
        return obj_arr, tv_arr, binary_time_arr, opt_time_arr, num_evo_arr



def compute_average_acc(file_names):
    num_exp = 0
    total_time_evo = 0
    total_time = 0
    total_num_evo = 0
    for file_name in file_names:
        num_switch = int(file_name.split("switch")[1].split("_")[0])
        file = open(file_name, 'r')
        lines = file.readlines()
        for line in lines:
            if "nfev:" in line:
                num_evo = int(line.split(" ")[-1].strip("\n"))
                # print(obj)
            if "objective function evaluation time" in line:
                time_evo = float(line.split(" ")[-2].strip("\n"))
            if "computational time of optimization" in line:
                compute_time = float(line.split(" ")[-1].strip("\n"))
        total_num_evo += num_evo
        num_exp += (num_switch + 1) * num_evo
        total_time_evo += time_evo
        total_time += compute_time
    num_exp /= len(file_names)
    total_time_evo /= len(file_names)
    total_time /= len(file_names)
    total_num_evo /= len(file_names)
    print(num_exp, total_time_evo, total_time, total_num_evo)

# print(read_energy_log("result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance1_alpha_0.015.log"))
# print(read_energy_log("result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance2_alpha_0.015.log"))
# print(read_energy_log("result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch12_initwarm_instance3_alpha_0.015.log"))
# print(read_energy_log("result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance4_alpha_0.015.log"))
# print(read_energy_log("result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance5_alpha_0.015.log"))
# alpha_list = [0.001 * i for i in range(9)] + [0.009] + [0.01 * i for i in range(1, 10)] + [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.015, 0.025, 0.25, 0.075]
# alpha_list = np.array([0, 0.002, 0.006, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.14, 0.2, 0.4, 0.6, 1.2])/2
# print(sorted((alpha_list)))
# write_csv(sorted(alpha_list), "result/output/SwitchTime/alpha/")

if __name__ == '__main__':
    # get_all_results(["Energy2", "Energy4", "Energy6", "CNOT5", "CNOT10", "CNOT20", "NOT2", "NOT6", "NOT10",
    #                  "CircuitH2", "CircuitLiH", "CircuitBeH2"], "result/output/SwitchTime/fewer", None)
    # get_all_results(["Energy2", "Energy4", "Energy6", "CNOT5", "CNOT10", "CNOT20", "NOT2", "NOT6", "NOT10",
    #                  "CircuitH2", "CircuitLiH", "CircuitBeH2"], "result/output/SwitchTime",
    #                  "result/output/all_results_new.csv")
    compute_average_acc(["result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance1_alpha_0.015.log",
                         "result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance2_alpha_0.015.log",
                         "result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch12_initwarm_instance3_alpha_0.015.log",
                         "result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance4_alpha_0.015.log",
                         "result/output/SwitchTime/EnergyAccC6obj_evotime_5.0_n_ts100_n_switch10_initwarm_instance5_alpha_0.015.log"])
    compute_average_acc(["result/output/SwitchTime/EnergyAccC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance1_alpha_0.01.log",
                         "result/output/SwitchTime/EnergyAccC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance2_alpha_0.01.log",
                         "result/output/SwitchTime/EnergyAccC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance3_alpha_0.01.log",
                         "result/output/SwitchTime/EnergyAccC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance4_alpha_0.01.log",
                         "result/output/SwitchTime/EnergyAccC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance5_alpha_0.01.log"])
    # compute_average_acc(
    #     ["result/output/SwitchTime/noacc/EnergyC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance1_alpha_0.01.log",
    #      "result/output/SwitchTime/noacc/EnergyC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance2_alpha_0.01.log",
    #      "result/output/SwitchTime/noacc/EnergyC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance3_alpha_0.01.log",
    #      "result/output/SwitchTime/noacc/EnergyC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance4_alpha_0.01.log",
    #      "result/output/SwitchTime/noacc/EnergyC6sur_evotime_5.0_n_ts100_n_switch10_initwarm_instance5_alpha_0.01.log"])
    compute_average_acc(["result/output/SwitchTime/EnergyAccD6binary_evotime_5.0_n_ts100_n_switch6_initwarm_instance1_alpha_0_mintv.log",
                        "result/output/SwitchTime/EnergyAccD6binary_evotime_5.0_n_ts100_n_switch5_initwarm_instance2_alpha_0_mintv.log",
                        "result/output/SwitchTime/EnergyAccD6binary_evotime_5.0_n_ts100_n_switch6_initwarm_instance3_alpha_0_mintv.log",
                        "result/output/SwitchTime/EnergyAccD6binary_evotime_5.0_n_ts100_n_switch6_initwarm_instance4_alpha_0_mintv.log",
                        "result/output/SwitchTime/EnergyAccD6binary_evotime_5.0_n_ts100_n_switch4_initwarm_instance5_alpha_0_mintv.log"])
