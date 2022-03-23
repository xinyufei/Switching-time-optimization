import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tools.auxiliary_molecule import generate_molecule_func
from tools.auxiliary_energy import *


def draw_threshold(instance='H2', mode="sep_per"):
    threshold, obj, obj_sur, tv, tv_sur, admm_obj, max_tv = None, None, None, None, None, None, None

    if instance == 'H2':
        threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        obj = [2.5604763034259292e-06, 2.5604763034259292e-06, 4.085888842819685e-07, 4.086649161294531e-07,
               4.0658313571473315e-07, 5.159569271828701e-08, 2.8041560551361755e-07, 1.1678685185589899e-07,
               4.893899530067358e-08, 2.4012535515538502e-08]
        tv = [62, 62, 74, 74, 74, 58, 52, 48, 36, 24]
        
        obj_sur = [4.0015787350355936e-09, 4.0015787350355936e-09, 1.7960340215061876e-08, 6.788980488892093e-10,
                   4.291472838202637e-07, 7.388747391701145e-10, 7.388747391701145e-10, 9.410083823269133e-11,
                   2.5455332286483667e-09, 3.2138073535747935e-08]
        tv_sur = [72, 72, 76, 72, 72, 70, 70, 66, 62, 54]

        # admm_obj = 1.33e-05
        # max_tv = 72

    if instance == "BeH2":
        threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        obj = [7.304110769101868e-06, 8.57114753394228e-07, 3.0718463739365376e-07, 4.514792461118855e-05,
               1.8491606301740404e-07, 2.526554518933466e-07, 2.1274215824540477e-06, 4.0447738269833167e-08,
               2.2238146257791414e-07, 2.6858612756086586e-07]
        tv = [372, 312, 292, 282, 268, 274, 268, 244, 176, 128]

        obj_sur = [2.7627211252045925e-07, 3.8139089475475174e-07, 3.8139089475475174e-07, 3.8139089475475174e-07,
                   3.227216884837958e-07, 1.1590271853378908e-08, 7.613623720370555e-08, 1.4060334563303911e-08,
                   3.5548098908932957e-09, 1.366429192017904e-10]
        tv_sur = [384, 380, 380, 380, 380, 372, 374, 338, 232, 192]

        # admm_obj = 1.51e-07
        # max_tv = 384

    if instance == "LiH":
        threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999]
        obj = [0.0016300845539389819, 0.0016453065224827368, 0.0016325817157971656,
               0.001636694130179861, 0.0016406033006001186, 0.0016329232452819697, 0.001634177350660071,
               0.0016424429254702222, 0.0016512816476880188, 0.0016358850624846877, 0.0016479792765931034,
               0.0016479186490706565, 0.0016484078942400338]
               # 0.006627966615676439]
        obj_sur = [0.0015630913010576952, 0.0015608963559184952,
                   0.0015606807308062853, 0.0015560945768714474, 0.0015574275852653363, 0.001569007743440154,
                   0.0015703856822677498, 0.0015706854144919014, 0.0015701928004578924, 0.0015700497898371024]
        tv = [232, 182, 172, 170, 146, 156, 138, 112, 74, 78, 58, 30, 26]
        tv_sur = [232, 228, 226, 230, 224, 226, 218, 204, 200, 116]

    if mode == "spe_per":
        plt.figure(dpi=300)
        plt.plot(threshold[:10], [np.log10(obj_e / obj_sur[0]) for obj_e in obj[:10]], '-o', label='STR')
        plt.plot(threshold[:10], [np.log10(obj_e / obj_sur[0]) for obj_e in obj_sur], '-^', label='SUR+ST')
        # plt.plot(threshold[:10], [obj_e / obj_sur[0] for obj_e in obj[:10]], '-o', label='STR')
        # plt.plot(threshold[:10], [obj_e / obj_sur[0] for obj_e in obj_sur], '-^', label='SUR+ST')
        plt.xlabel("Threshold of selecting controller")
        plt.ylabel("Logarithm of ratio")
        # plt.ylim([-2.809, -2.776])
        plt.legend(loc='upper right')
        if instance == "H2":
            plt.savefig("../figure_paper/Molecule_H2_evotime4.0_n_ts80_obj_log10_comp_per.png")
        if instance == "LiH":
            plt.ylim([-0.004, 0.030])
            plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_obj_log10_comp_per.png")
        if instance == "BeH2":
           plt.savefig("../figure_paper/Molecule_BeH2_evotime20.0_n_ts200_obj_log10_comp_per.png")

        plt.figure(dpi=300)
        plt.plot(threshold[:10], [tv_e / tv_sur[0] for tv_e in tv[:10]], '-o', label='STR')
        plt.plot(threshold[:10], [tv_e / tv_sur[0] for tv_e in tv_sur], '-^', label='SUR+ST')
        plt.xlabel("Threshold of selecting controller")
        plt.ylabel("TV norm ratio")
        plt.legend()
        if instance == "H2":
            plt.savefig("../figure_paper/Molecule_H2_evotime4.0_n_ts80_tv_comp_per.png")
        if instance == "LiH":
            plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_tv_comp_per.png")
        if instance == "BeH2":
            plt.savefig("../figure_paper/Molecule_BeH2_evotime20.0_n_ts200_tv_comp_per.png")

    if mode == "spe_abs":
        plt.figure(dpi=300)
        plt.plot(threshold[:10], [np.log10(obj_e) for obj_e in obj], '-o', label='STR')
        plt.plot(threshold[:10], [np.log10(obj_e) for obj_e in obj_sur], '-^', label='SUR+ST')
        plt.xlabel("Threshold of selecting controller")
        plt.ylabel("Logarithm of Objective value")
        # plt.ylim([-2.809, -2.776])
        plt.legend(loc='upper right')
        if instance == "H2":
            plt.savefig("../figure_paper/Molecule_H2_evotime4.0_n_ts80_obj_log10_comp.png")
        if instance == "LiH":
            plt.ylim([-2.809, -2.776])
            plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_obj_log10_comp.png")
        if instance == "BeH2":
            plt.savefig("../figure_paper/Molecule_BeH2_evotime20.0_n_ts200_obj_log10_comp.png")

        plt.figure(dpi=300)
        plt.plot(threshold[:10], tv, '-o', label='STR')
        plt.plot(threshold[:10], tv_sur, '-^', label='SUR+ST')
        plt.xlabel("Threshold of selecting controller")
        plt.ylabel("TV norm")
        plt.legend()
        if instance == "H2":
            plt.savefig("../figure_paper/Molecule_H2_evotime4.0_n_ts80_tv_comp.png")
        if instance == "LiH":
            plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_tv_comp.png")
        if instance == "BeH2":
            plt.savefig("../figure_paper/Molecule_BeH2_evotime20.0_n_ts200_tv_comp.png")


def draw_err_bar():
    threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    num_instance = 3
    ratio_obj_all = np.zeros((num_instance, len(threshold)))
    ratio_obj_sur_all = np.zeros((num_instance, len(threshold)))
    ratio_obj_la_all = np.zeros((num_instance, len(threshold)))
    ratio_tv_all = np.zeros((num_instance, len(threshold)))
    ratio_tv_sur_all = np.zeros((num_instance, len(threshold)))
    ratio_tv_la_all = np.zeros((num_instance, len(threshold)))
    obj_h2 = [2.5604763034259292e-06, 2.5604763034259292e-06, 4.085888842819685e-07, 4.086649161294531e-07,
              4.0658313571473315e-07, 5.159569271828701e-08, 2.8041560551361755e-07, 1.1678685185589899e-07,
              4.893899530067358e-08, 2.4012535515538502e-08]
    tv_h2 = [62, 62, 74, 74, 74, 58, 52, 48, 36, 24]
    obj_sur_h2 = [4.0015787350355936e-09, 4.0015787350355936e-09, 1.7960340215061876e-08, 6.788980488892093e-10,
                  4.291472838202637e-07, 7.388747391701145e-10, 7.388747391701145e-10, 9.410083823269133e-11,
                  2.5455332286483667e-09, 3.2138073535747935e-08]
    tv_sur_h2 = [72, 72, 76, 72, 72, 70, 70, 66, 62, 54]
    obj_la_h2 = [0.0053821703228478235] * len(threshold)
    tv_la_h2 = [4] * len(threshold)
    ratio_obj_all[0, :] = np.log10(np.array([obj_e / obj_sur_h2[0] for obj_e in obj_h2]))
    ratio_obj_sur_all[0, :] = np.log10(np.array([obj_e / obj_sur_h2[0] for obj_e in obj_sur_h2]))
    ratio_obj_la_all[0, :] = np.log10(np.array([obj_e / obj_sur_h2[0] for obj_e in obj_la_h2]))
    ratio_tv_all[0, :] = np.array([tv_e / tv_sur_h2[0] for tv_e in tv_h2])
    ratio_tv_sur_all[0, :] = np.array([tv_e / tv_sur_h2[0] for tv_e in tv_sur_h2])
    ratio_tv_la_all[0, :] = np.array([tv_e / tv_sur_h2[0] for tv_e in tv_la_h2])

    obj_lih = [0.0016300845539389819, 0.0016453065224827368, 0.0016325817157971656, 0.001636694130179861,
               0.0016406033006001186, 0.0016329232452819697, 0.001634177350660071, 0.0016424429254702222,
               0.0016512816476880188, 0.0016358850624846877]
    tv_lih = [232, 182, 172, 170, 146, 156, 138, 112, 74, 78]
    obj_sur_lih = [0.0015630913010576952, 0.0015608963559184952,
                   0.0015606807308062853, 0.0015560945768714474, 0.0015574275852653363, 0.001569007743440154,
                   0.0015703856822677498, 0.0015706854144919014, 0.0015701928004578924, 0.0015700497898371024]
    tv_sur_lih = [232, 228, 226, 230, 224, 226, 218, 204, 200, 116]
    obj_la_lih = [0.0015724211989566195, 0.0015711559004031317, 0.005231124150308242, 0.0015764118661691917,
                  0.0016481498837103148, 0.0015876652356658916, 0.2197079377010338, 0.0016538190500802186,
                  0.0016538190500802186, 0.9999104954213078]
    tv_la_lih = [34, 22, 26, 22, 12, 12, 14, 4, 4, 0]
    ratio_obj_all[1, :] = np.log10(np.array([obj_e / obj_sur_lih[0] for obj_e in obj_lih]))
    ratio_obj_sur_all[1, :] = np.log10(np.array([obj_e / obj_sur_lih[0] for obj_e in obj_sur_lih]))
    ratio_obj_la_all[1, :] = np.log10(np.array([obj_e / obj_sur_lih[0] for obj_e in obj_la_lih]))
    ratio_tv_all[1, :] = np.array([tv_e / tv_sur_lih[0] for tv_e in tv_lih])
    ratio_tv_sur_all[1, :] = np.array([tv_e / tv_sur_lih[0] for tv_e in tv_sur_lih])
    ratio_tv_la_all[1, :] = np.array([tv_e / tv_sur_lih[0] for tv_e in tv_la_lih])

    obj_beh2 = [7.304110769101868e-06, 8.57114753394228e-07, 3.0718463739365376e-07, 4.514792461118855e-05,
                1.8491606301740404e-07, 2.526554518933466e-07, 2.1274215824540477e-06, 4.0447738269833167e-08,
                2.2238146257791414e-07, 2.6858612756086586e-07]
    tv_beh2 = [372, 312, 292, 282, 268, 274, 268, 244, 176, 128]
    obj_sur_beh2 = [2.7627211252045925e-07, 3.8139089475475174e-07, 3.8139089475475174e-07, 3.8139089475475174e-07,
                    3.227216884837958e-07, 1.1590271853378908e-08, 7.613623720370555e-08, 1.4060334563303911e-08,
                    3.5548098908932957e-09, 1.366429192017904e-10]
    tv_sur_beh2 = [384, 380, 380, 380, 380, 372, 374, 338, 232, 192]
    obj_la_beh2 = [1] * len(threshold)
    tv_la_beh2 = [0] * len(threshold)

    ratio_obj_all[2, :] = np.log10(np.array([obj_e / obj_sur_beh2[0] for obj_e in obj_beh2]))
    ratio_obj_sur_all[2, :] = np.log10(np.array([obj_e / obj_sur_beh2[0] for obj_e in obj_sur_beh2]))
    ratio_obj_la_all[1, :] = np.log10(np.array([obj_e / obj_sur_beh2[0] for obj_e in obj_la_beh2]))
    ratio_tv_all[2, :] = np.array([tv_e / tv_sur_beh2[0] for tv_e in tv_beh2])
    ratio_tv_sur_all[2, :] = np.array([tv_e / tv_sur_beh2[0] for tv_e in tv_sur_beh2])
    ratio_tv_la_all[1, :] = np.array([tv_e / tv_sur_beh2[0] for tv_e in tv_la_beh2])

    average_ratio_obj = np.mean(ratio_obj_all, axis=0)
    average_ratio_obj_sur = np.mean(ratio_obj_sur_all, axis=0)
    average_ratio_obj_la = np.mean(ratio_obj_la_all, axis=0)
    average_ratio_tv = np.mean(ratio_tv_all, axis=0)
    average_ratio_tv_sur = np.mean(ratio_tv_sur_all, axis=0)
    average_ratio_tv_la = np.mean(ratio_tv_la_all, axis=0)
    ratio_obj_err = np.zeros((2, len(threshold)))
    ratio_obj_err[0, :] = -np.min(ratio_obj_all, axis=0) + average_ratio_obj
    ratio_obj_err[1, :] = np.max(ratio_obj_all, axis=0) - average_ratio_obj
    ratio_obj_sur_err = np.zeros((2, len(threshold)))
    ratio_obj_sur_err[0, :] = -np.min(ratio_obj_sur_all, axis=0) + average_ratio_obj_sur
    ratio_obj_sur_err[1, :] = np.max(ratio_obj_sur_all, axis=0) - average_ratio_obj_sur
    ratio_obj_la_err = np.zeros((2, len(threshold)))
    ratio_obj_la_err[0, :] = -np.min(ratio_obj_la_all, axis=0) + average_ratio_obj_la
    ratio_obj_la_err[1, :] = np.max(ratio_obj_la_all, axis=0) - average_ratio_obj_la
    ratio_tv_err = np.zeros((2, len(threshold)))
    ratio_tv_err[0, :] = -np.min(ratio_tv_all, axis=0) + average_ratio_tv
    ratio_tv_err[1, :] = np.max(ratio_tv_all, axis=0) - average_ratio_tv
    ratio_tv_sur_err = np.zeros((2, len(threshold)))
    ratio_tv_sur_err[0, :] = -np.min(ratio_tv_sur_all, axis=0) + average_ratio_tv_sur
    ratio_tv_sur_err[1, :] = np.max(ratio_tv_sur_all, axis=0) - average_ratio_tv_sur
    ratio_tv_la_err = np.zeros((2, len(threshold)))
    ratio_tv_la_err[0, :] = -np.min(ratio_tv_la_all, axis=0) + average_ratio_tv_la
    ratio_tv_la_err[1, :] = np.max(ratio_tv_la_all, axis=0) - average_ratio_tv_la

    plt.figure(dpi=300)
    plt.errorbar(threshold, average_ratio_obj, yerr=ratio_obj_err, fmt='-o', capsize=5, label='STR')
    plt.errorbar(threshold, average_ratio_obj_la, yerr=ratio_obj_la_err, fmt='-+', capsize=5, label='STLA')
    plt.errorbar(threshold, average_ratio_obj_sur, yerr=ratio_obj_sur_err, fmt='-D', capsize=5, label='SUR+ST')
    plt.xlabel("Threshold of selecting controller")
    plt.ylabel("Logarithm of objective value ratio")
    plt.legend(loc='upper right')
    plt.savefig("../figure_paper/all_instances_obj_log10_3methods.png")

    plt.figure(dpi=300)
    plt.errorbar(threshold, average_ratio_tv, yerr=ratio_tv_err, fmt='-o', capsize=5, label='STR')
    plt.errorbar(threshold, average_ratio_tv_la, yerr=ratio_tv_la_err, fmt='-+', capsize=5, label='STLA')
    plt.errorbar(threshold, average_ratio_tv_sur, yerr=ratio_tv_sur_err, fmt='-D', capsize=5, label='SUR+ST')
    plt.xlabel("Threshold of selecting controller")
    plt.ylabel("TV norm ratio")
    plt.legend(loc='upper right')
    plt.savefig("../figure_paper/all_instances_tv_log10_3methods.png")


def draw_threshold_stla():
    threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    obj_la_lih = [0.0015724211989566195, 0.0015711559004031317, 0.005231124150308242, 0.0015764118661691917,
                  0.0016481498837103148, 0.0015876652356658916, 0.2197079377010338, 0.0016538190500802186,
                  0.0016538190500802186, 0.9999104954213078]
    tv_la_lih = [34, 22, 26, 22, 12, 12, 14, 4, 4, 0]
    plt.figure(dpi=300)
    plt.plot(threshold[:10], np.log10(np.array(obj_la_lih)), '-o')
    plt.xlabel("Threshold of selecting controller")
    plt.ylabel("Logarithm of Objective value")
    # plt.ylim([-2.809, -2.776])
    plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_obj_log10_stla.png")

    plt.figure(dpi=300)
    plt.plot(threshold[:10], tv_la_lih, '-o')
    plt.xlabel("Threshold of selecting controller")
    plt.ylabel("TV norm")
    # plt.legend()
    plt.savefig("../figure_paper/Molecule_LiH_evotime20.0_n_ts200_tv_stla.png")


def draw_obj_tv():
    threshold = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    num_instance = 3
    ratio_obj_all = np.zeros((num_instance, len(threshold)))
    ratio_obj_sur_all = np.zeros((num_instance, len(threshold)))
    ratio_tv_all = np.zeros((num_instance, len(threshold)))
    ratio_tv_sur_all = np.zeros((num_instance, len(threshold)))
    obj_h2 = [2.5604763034259292e-06, 2.5604763034259292e-06, 4.085888842819685e-07, 4.086649161294531e-07,
              4.0658313571473315e-07, 5.159569271828701e-08, 2.8041560551361755e-07, 1.1678685185589899e-07,
              4.893899530067358e-08, 2.4012535515538502e-08]
    tv_h2 = [62, 62, 74, 74, 74, 58, 52, 48, 36, 24]
    obj_sur_h2 = [4.0015787350355936e-09, 4.0015787350355936e-09, 1.7960340215061876e-08, 6.788980488892093e-10,
                  4.291472838202637e-07, 7.388747391701145e-10, 7.388747391701145e-10, 9.410083823269133e-11,
                  2.5455332286483667e-09, 3.2138073535747935e-08]
    tv_sur_h2 = [72, 72, 76, 72, 72, 70, 70, 66, 62, 54]
    ratio_obj_all[0, :] = np.log10(np.array([obj_e / obj_sur_h2[0] for obj_e in obj_h2]))
    ratio_obj_sur_all[0, :] = np.log10(np.array([obj_e / obj_sur_h2[0] for obj_e in obj_sur_h2]))
    ratio_tv_all[0, :] = np.array([tv_e / tv_sur_h2[0] for tv_e in tv_h2])
    ratio_tv_sur_all[0, :] = np.array([tv_e / tv_sur_h2[0] for tv_e in tv_sur_h2])

    obj_lih = [0.0016300845539389819, 0.0016453065224827368, 0.0016325817157971656, 0.001636694130179861,
               0.0016406033006001186, 0.0016329232452819697, 0.001634177350660071, 0.0016424429254702222,
               0.0016512816476880188, 0.0016358850624846877]
    obj_sur_lih = [0.0015630913010576952, 0.0015608963559184952,
                   0.0015606807308062853, 0.0015560945768714474, 0.0015574275852653363, 0.001569007743440154,
                   0.0015703856822677498, 0.0015706854144919014, 0.0015701928004578924, 0.0015700497898371024]
    tv_lih = [232, 182, 172, 170, 146, 156, 138, 112, 74, 78]
    tv_sur_lih = [232, 228, 226, 230, 224, 226, 218, 204, 200, 116]
    ratio_obj_all[1, :] = np.log10(np.array([obj_e / obj_sur_lih[0] for obj_e in obj_lih]))
    ratio_obj_sur_all[1, :] = np.log10(np.array([obj_e / obj_sur_lih[0] for obj_e in obj_sur_lih]))
    ratio_tv_all[1, :] = np.array([tv_e / tv_sur_lih[0] for tv_e in tv_lih])
    ratio_tv_sur_all[1, :] = np.array([tv_e / tv_sur_lih[0] for tv_e in tv_sur_lih])

    obj_beh2 = [7.304110769101868e-06, 8.57114753394228e-07, 3.0718463739365376e-07, 4.514792461118855e-05,
                1.8491606301740404e-07, 2.526554518933466e-07, 2.1274215824540477e-06, 4.0447738269833167e-08,
                2.2238146257791414e-07, 2.6858612756086586e-07]
    tv_beh2 = [372, 312, 292, 282, 268, 274, 268, 244, 176, 128]
    obj_sur_beh2 = [2.7627211252045925e-07, 3.8139089475475174e-07, 3.8139089475475174e-07, 3.8139089475475174e-07,
                    3.227216884837958e-07, 1.1590271853378908e-08, 7.613623720370555e-08, 1.4060334563303911e-08,
                    3.5548098908932957e-09, 1.366429192017904e-10]
    tv_sur_beh2 = [384, 380, 380, 380, 380, 372, 374, 338, 232, 192]
    ratio_obj_all[2, :] = np.log10(np.array([obj_e / obj_sur_beh2[0] for obj_e in obj_beh2]))
    ratio_obj_sur_all[2, :] = np.log10(np.array([obj_e / obj_sur_beh2[0] for obj_e in obj_sur_beh2]))
    ratio_tv_all[2, :] = np.array([tv_e / tv_sur_beh2[0] for tv_e in tv_beh2])
    ratio_tv_sur_all[2, :] = np.array([tv_e / tv_sur_beh2[0] for tv_e in tv_sur_beh2])

    average_ratio_obj = np.mean(ratio_obj_all, axis=0)
    average_ratio_obj_sur = np.mean(ratio_obj_sur_all, axis=0)
    average_ratio_tv = np.mean(ratio_tv_all, axis=0)
    average_ratio_tv_sur = np.mean(ratio_tv_sur_all, axis=0)

    plt.figure(dpi=300)
    plt.plot(average_ratio_tv, average_ratio_obj, 'o', label='STR')
    plt.plot(average_ratio_tv_sur, average_ratio_obj_sur, '^', label='SUR+ST')
    plt.xlabel("TV-norm ratio")
    plt.ylabel("Logarithm of objective value ratio")
    plt.legend(loc='upper right')
    plt.savefig("../figure_paper/all_instances_obj_log10_tv.png")


def draw_diff_str():
    change = [0.007315488214812316, 0.001875269953916514, 0.00019485549080711095, 0.0005134337226989638,
              0.0003158136079963736, 5.0216688023185796e-05]
    sum_change = [4.030051957309731, 2.421914577993002, 1.5029534818402843, 1.0266885773868866, 0.9614853871329476,
                  0.7503757548819144]
    max_change = [0.14860567218377674, 0.04071606059538557, 0.018703382188763662, 0.010397513838499406,
                  0.006543987744289881, 0.004391594411919053]
    max_chosen_change = [0.007416306107414172, 0.0020338362101174345, 0.0010098776378355545, 0.00048749522195989936,
                         0.0004005764389356514, 0.0002813090392926876]
    time_steps = [80, 160, 240, 320, 400, 480]

    # plt.figure(dpi=300)
    # plt.plot(time_steps, np.log10(np.array(max_chosen_change)), '-o', label='Maximum chosen change')
    # plt.plot(time_steps, -2 * np.log10(np.array(time_steps)) + 2, '-^', label='-2logT+C')
    # plt.plot(time_steps, np.array(max_chosen_change), '-o', label='Maximum chosen change')
    # plt.plot(time_steps, 2 / np.array(time_steps), '-^', label='C/T^2')
    # plt.xlabel('Time steps')
    # plt.legend()
    # plt.ylabel('Common Logarithm')
    # plt.savefig("../figure_paper/MoleculeVQEADMMSTCost_H2_evotime4.0_maxchosen_ts_nolog.png")

    # exit()

    plt.figure(dpi=300)
    # plt.plot(time_steps, np.log10(np.array(change)), '-o', label='Difference')
    plt.plot(time_steps, np.array(change), '-o', label='Difference')
    plt.xlabel('Time steps')
    plt.legend()
    plt.ylabel('Common Logarithm')
    plt.savefig("../figure_paper/MoleculeVQEADMMSTCost_H2_evotime4.0_diff_ts_nolog.png")
    exit()

    plt.figure(dpi=300)
    plt.plot(time_steps, np.log10(np.array(sum_change)), '-o', label='Summation of changes')
    plt.plot(time_steps, -np.log10(np.array(time_steps)) + 2.7, '-^', label='-logT+C')
    plt.xlabel('Time steps')
    plt.legend()
    # plt.ylabel('Objective value')
    plt.savefig("../figure_paper/MoleculeVQEADMMSTCost_H2_evotime4.0_sumchange_ts.png")

    plt.figure(dpi=300)
    plt.plot(time_steps, np.log10(np.array(max_change)), '-o', label='Max change')
    plt.plot(time_steps, -2 * np.log10(np.array(time_steps)) + 3.1, '-^', label='-2logT+C')
    plt.xlabel('Time steps')
    plt.legend()
    # plt.ylabel('Objective value')
    plt.savefig("../figure_paper/MoleculeVQEADMMSTCost_H2_evotime4.0_maxchange_ts.png")


def draw_str_diff_ub_molecule():
    change = [0.007315488214812316, 0.001875269953916514, 0.00019485549080711095, 0.0005134337226989638,
              0.0003158136079963736, 5.0216688023185796e-05]
    ts_list = [80, 160, 240, 320, 400, 480]
    d = 2
    qubit_num = 2
    molecule = "H2"
    Hops, H0, U0, U = generate_molecule_func(qubit_num, d, molecule)
    evo_time = 4
    max_sigma = 0
    for H in Hops:
        u, s, vh = np.linalg.svd(H, full_matrices=True)
        if max_sigma < s[0]:
            max_sigma = s[0]

    # print(max_sigma)

    ub_list = []
    ub_1_list = []
    ub_2_list = []
    for n_ts in ts_list:
        c_control_name = "../example/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts" + str(n_ts) + \
                         "_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter100.csv"
        c_control = np.loadtxt(c_control_name, delimiter=',')

        # print(abs(np.sum(c_control, axis=1)))
        epsilon = np.max(abs(np.sum(c_control, axis=1) - 1))
        epsilon_sum = np.sum(abs(np.sum(c_control, axis=1) - 1))
        delta_t = evo_time / n_ts

        print(epsilon, epsilon_sum * delta_t)

        C1 = (1 + epsilon) ** 2 * max_sigma ** 2
        C2 = (1 + epsilon) * max_sigma
        C0 = 2 ** qubit_num * (1 + epsilon) * max_sigma * epsilon_sum / (1 - epsilon) * delta_t

        ub_1 = 2 * C1 * np.exp(C2 * delta_t) * delta_t
        ub_2 = C0
        ub = ub_1 + ub_2
        ub_1_list.append(ub_1)
        ub_2_list.append(ub_2)
        ub_list.append(ub)

    # draw the figure
    print(ub_list)
    print(ub_1_list)
    print(ub_2_list)
    plt.figure(dpi=300)
    plt.xlabel("Time steps")
    plt.ylabel("Common logarithm")
    plt.plot(ts_list, np.log10(change), '-o', label='difference')
    plt.plot(ts_list, np.log10(ub_list), linestyle="-", marker='s',
             label="upper bound")
    plt.plot(ts_list, np.log10(ub_1_list), linestyle="--", marker='^', label="first term of the upper bound")
    plt.plot(ts_list, np.log10(ub_2_list), linestyle="--", marker='+', markersize='8',
             label="second term of the upper bound")
    # plt.plot(delta_t_list, integral_err, label='Maximum integral error')
    plt.legend(prop={'size': 6})
    # plt.legend()
    plt.savefig("../figure_paper/MoleculeNew_H2_evotime4.0_str_error_delta_t_log10.png")
    

def draw_str_diff_ub_energy():
    change = [0.0022654349904387416, 0.00037150083126546996, -3.3161086668842543e-06, -1.5042721027147543e-05,
              5.654499404983415e-05, 5.400374985864431e-06]
    ts_list = [40, 80, 120, 160, 200, 240]
    n = 2
    num_edges = 1
    Jij, edges = generate_Jij_MC(n, num_edges, 100)

    # n = 4
    # Jij = generate_Jij(n, 1)

    C = get_ham(n, True, Jij)
    B = get_ham(n, False, Jij)

    Hops = [B, C]
    evo_time = 2
    max_sigma = 0
    for H in Hops:
        u, s, vh = np.linalg.svd(H, full_matrices=True)
        if max_sigma < s[0]:
            max_sigma = s[0]

    u, s, vh = np.linalg.svd(C, full_matrices=True)
    sum_sigma = sum(s)

    # print(max_sigma)

    Emin = 1
    # Emin = -2.511

    ub_list = []
    ub_1_list = []
    ub_2_list = []
    for n_ts in ts_list:
        delta_t = evo_time / n_ts

        C1 = 2 ** n * sum_sigma * (2 ** n * max_sigma ** 2 + 2 + 2 ** n * max_sigma) / abs(Emin)
        C2 = max_sigma

        ub_1 = 2 * C1 * np.exp(C2 * delta_t) * delta_t
        ub_2 = 0
        ub = ub_1 + ub_2
        ub_1_list.append(ub_1)
        ub_2_list.append(ub_2)
        ub_list.append(ub)

    # draw the figure
    print(ub_list)
    print(ub_1_list)
    print(ub_2_list)
    plt.figure(dpi=300)
    plt.xlabel("Time steps")
    plt.ylabel("Objective value")
    plt.plot(ts_list, change, '-o', label='difference')
    plt.plot(ts_list, ub_list, linestyle="-", marker='s', label="upper bound")
    # plt.plot(ts_list, np.log10(ub_1_list), linestyle="--", marker='^', label="first term of the upper bound")
    # plt.plot(ts_list, np.log10(ub_2_list), linestyle="--", marker='+', markersize='8',
    #          label="second term of the upper bound")
    # plt.plot(delta_t_list, integral_err, label='Maximum integral error')
    plt.legend(prop={'size': 6})
    # plt.legend()
    # plt.savefig("../figure_paper/EnergyADMM2_evotime2.0_str_error_delta_t_log10.png")
    plt.savefig("../figure_paper/EnergyADMM2_evotime2.0_str_error_delta_t.png")


def draw_str_error_ub_molecule_log2():
    change = [0.041733463100893764, 0.010212036771837019, 0.0027534119707577354, 0.0007731417198063584,
              0.00013517161104248387, 4.905932617771391e-05, 1.5732045017924357e-05]
    ts_list = [32, 64, 128, 256, 512, 1024, 2048]
    d = 2
    qubit_num = 2
    molecule = "H2"
    Hops, H0, U0, U = generate_molecule_func(qubit_num, d, molecule)
    evo_time = 4
    max_sigma = 0
    for H in Hops:
        u, s, vh = np.linalg.svd(H, full_matrices=True)
        if max_sigma < s[0]:
            max_sigma = s[0]

    # print(max_sigma)

    ub_list = []
    ub_1_list = []
    ub_2_list = []
    for n_ts in ts_list:
        c_control_name = "../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts" + str(n_ts) + \
                         "_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter50.csv"
        c_control = np.loadtxt(c_control_name, delimiter=',')

        # print(abs(np.sum(c_control, axis=1)))
        epsilon = np.max(abs(np.sum(c_control, axis=1) - 1))
        epsilon_sum = np.sum(abs(np.sum(c_control, axis=1) - 1))
        delta_t = evo_time / n_ts

        print(epsilon, epsilon_sum * delta_t)

        C1 = (1 + epsilon) ** 2 * max_sigma ** 2
        C2 = (1 + epsilon) * max_sigma
        C0 = (1 + epsilon) * max_sigma * epsilon_sum / np.abs(1 - epsilon) * delta_t

        ub_1 = 2 * C1 * np.exp(C2 * delta_t) * delta_t
        ub_2 = C0
        ub = ub_1 + ub_2
        ub_1_list.append(ub_1)
        ub_2_list.append(ub_2)
        ub_list.append(ub)

    # draw the figure
    print(ub_list)
    print(ub_1_list)
    print(ub_2_list)
    plt.figure(dpi=300)
    plt.xlabel("Binary logarithm of time steps")
    plt.ylabel("Binary logarithm")
    plt.plot(np.log2(ts_list), np.log2(change), '-o', label='difference')
    plt.plot(np.log2(ts_list), np.log2(ub_list), linestyle="-", marker='s', label="upper bound")
    # plt.plot(np.log2(ts_list), np.log2(ub_1_list), linestyle="--", marker='^', label="first term of the upper bound")
    # plt.plot(np.log2(ts_list), np.log2(ub_2_list), linestyle="--", marker='+', markersize='8',
    #          label="second term of the upper bound")
    # plt.plot(delta_t_list, integral_err, label='Maximum integral error')
    plt.legend(prop={'size': 6})
    # plt.legend()
    plt.savefig("../figure_paper/MoleculeNew_H2_evotime4.0_str_error_delta_t_log2_one_ub.png")


def draw_str_error_ub_energy_log2():
    change = [0.015805083194930347, 0.0010995902816599568, -0.0005207824336386224, -0.00025172322191790997,
              -0.00022925154556485694, -0.0023870906889555954, -0.003598010584822564]
    ts_list = [32, 64, 128, 256, 512, 1024, 2048]
    n = 4
    num_edges = 2
    # Jij, edges = generate_Jij_MC(n, num_edges, 100)

    # n = 4
    Jij = generate_Jij(n, 1)

    C = get_ham(n, True, Jij)
    B = get_ham(n, False, Jij)

    Hops = [B, C]
    evo_time = 2
    max_sigma = 0
    for H in Hops:
        u, s, vh = np.linalg.svd(H, full_matrices=True)
        if max_sigma < s[0]:
            max_sigma = s[0]

    u, s, vh = np.linalg.svd(C, full_matrices=True)
    max_sigma_energy = s[0]

    # print(max_sigma)

    # Emin = 1
    Emin = -2.511

    ub_list = []
    ub_1_list = []
    ub_2_list = []
    for n_ts in ts_list:
        delta_t = evo_time / n_ts

        C1 = max_sigma_energy * (max_sigma ** 2 + 2 + evo_time * max_sigma) / abs(Emin)
        C2 = max_sigma

        ub_1 = 2 * C1 * np.exp(C2 * delta_t) * delta_t
        ub_2 = 0
        ub = ub_1 + ub_2
        ub_1_list.append(ub_1)
        ub_2_list.append(ub_2)
        ub_list.append(ub)

    # draw the figure
    print(ub_list)
    print(ub_1_list)
    print(ub_2_list)
    plt.figure(figsize=(8, 6), dpi=300)
    plt.xlabel("Binary logarithm of time steps")
    plt.ylabel("Difference")
    plt.plot(np.log2(ts_list), change, '-o', label='difference')
    # plt.plot(np.log2(ts_list), ub_list, linestyle="-", marker='s', label="upper bound")
    # plt.plot(np.log2(ts_list), [0] * len(ts_list), linestyle='-', marker='s', label="")
    plt.legend(prop={'size': 6})
    ax = plt.gca()
    ax.spines['right'].set_color("None")
    ax.spines['top'].set_color("None")
    ax.spines['bottom'].set_position(("data", 0))
    # plt.legend()
    # plt.savefig("../figure_paper/EnergyADMM2_evotime2.0_str_error_delta_t_log10.png")
    plt.savefig("../figure_paper/EnergyADMM4_evotime2.0_str_error_delta_t_zero.png")


if __name__ == '__main__':
    # draw_threshold("H2", "spe_per")
    # draw_threshold("LiH", "spe_per")
    # draw_threshold("BeH2", per=True)
    # draw_err_bar()
    draw_threshold_stla()
    # draw_obj_tv()
    # draw_diff_str()
    # draw_str_diff_ub_molecule()
    # draw_str_diff_ub_energy()
    # draw_str_error_ub_molecule_log2()
    # draw_str_error_ub_energy_log2()
