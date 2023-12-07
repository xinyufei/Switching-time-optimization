#conda activate qcopt
cd ../acc_example/
#python Molecule_acc.py --name=MoleculeAccD --molecule=H2 --qubit_num=2 --evo_time=10 --n_ts=100 \
#  --initial_type='ave' --extract="binary" \
#  --target="../new_result/c_control/Target/H2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_H2_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty1.0.csv" \
#  --b_control="../new_result/b_control/TR+MT/MoleculeVQENew_H2_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty1.0_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"

#for i in 0.045;
#do
#python Molecule_acc.py --name=MoleculeAccC --molecule=H2 --qubit_num=2 --evo_time=10 --n_ts=100 \
#  --initial_type='warm' --extract="obj" --alpha=$i \
#  --target="../new_result/c_control/Target/H2_target.csv" \
#  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_H2_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty1.0.csv";
#done

for i in 0.01;
do
python Molecule_acc.py --name=MoleculeAccC --molecule=H2 --qubit_num=2 --evo_time=10 --n_ts=100 \
  --initial_type='warm' --extract="sur" --alpha=$i \
  --target="../new_result/c_control/Target/H2_target.csv" \
  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_H2_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty1.0.csv";
done
cd ../new_script/