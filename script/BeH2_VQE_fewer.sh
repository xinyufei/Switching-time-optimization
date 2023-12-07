#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=40 \
  --initial_type='warm' --extract="obj" --alpha=0.035 \
  --target="../new_result/c_control/Target/BeH2_target.csv" \
  --c_control="../new_result/c_control/Continuous/MoleculeVQE_BeH2_evotime20.0_n_ts40_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv"

#python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=40 \
#  --initial_type='ave' --extract="sur" --alpha=0.2 \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts40_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv"

cd ../new_script/