#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
#python Molecule_acc.py --name=MoleculeAccD --molecule=LiH --qubit_num=4 --evo_time=20 --n_ts=200 \
#  --initial_type='ave' --extract="binary" \
#  --target="../new_result/c_control/Target/LiH_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1.csv" \
#  --b_control="../new_result/b_control/TR+MT/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"

# 2nd choice: 0.035+ave
#for i in 0.03;
#do
#python Molecule_acc.py --name=MoleculeAccC --molecule=LiH --qubit_num=4 --evo_time=20 --n_ts=200 \
#  --initial_type='warm' --extract="obj" --alpha=$i \
#  --target="../new_result/c_control/Target/LiH_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1.csv";
##python Molecule_acc.py --name=MoleculeAccC --molecule=LiH --qubit_num=4 --evo_time=20 --n_ts=200 \
##  --initial_type='ave' --extract="obj" --alpha=$i \
##  --target="../new_result/c_control/Target/LiH_target.csv" \
##  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1.csv";
#done

for i in 0.06;
do
python Molecule_acc.py --name=MoleculeAccC --molecule=LiH --qubit_num=4 --evo_time=20 --n_ts=200 \
  --initial_type='warm' --extract="sur" --alpha=$i \
  --target="../new_result/c_control/Target/LiH_target.csv" \
  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1.csv";
#python Molecule_acc.py --name=MoleculeAccC --molecule=LiH --qubit_num=4 --evo_time=20 --n_ts=200 \
#  --initial_type='ave' --extract="sur" --alpha=$i \
#  --target="../new_result/c_control/Target/LiH_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_LiH_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.1.csv";
done
cd ../script/