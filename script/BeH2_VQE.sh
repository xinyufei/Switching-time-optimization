#!/user/bin/env bash
#conda activate qcopt
cd ../acc_example/
#for i in 0.2;
#do
##python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
##  --initial_type='warm' --extract="sur" --alpha=$i \
##  --target="../new_result/c_control/Target/BeH2_target.csv" \
##  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv";
##python Molecule.py --name=MoleculeC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
##  --initial_type='warm' --extract="sur" --alpha=$i \
##  --target="../new_result/c_control/Target/BeH2_target.csv" \
##  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv";
#done

#python Molecule.py --name=MoleculeC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
#  --initial_type='warm' --extract="la" --alpha=0.03 \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv"
#python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
#  --initial_type='warm' --extract="la" --alpha=0.03 \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv"

for i in 0.03;
do
#python Molecule.py --name=MoleculeC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
#  --initial_type='ave' --extract="obj" --alpha=0.03 \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv";
python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
  --initial_type='ave' --extract="obj" --alpha=$i \
  --target="../new_result/c_control/Target/BeH2_target.csv" \
  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv";
done

#python Molecule.py --name=MoleculeD --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
#  --initial_type='ave' --extract="binary" \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv" \
#  --b_control="../new_result/b_control/TR+MT/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"
#python Molecule_acc.py --name=MoleculeAccD --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
#  --initial_type='ave' --extract="binary" \
#  --target="../new_result/c_control/Target/BeH2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv" \
#  --b_control="../new_result/b_control/TR+MT/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"

cd ../new_script/