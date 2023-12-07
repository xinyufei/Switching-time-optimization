#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
#python CNOT_acc.py --name=CNOTAccC --evo_time=5 --n_ts=100 --initial_type='warm' --extract="obj" --alpha=0.01 \
#  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv"
#python CNOT_acc.py --name=CNOTAccC --evo_time=10 --n_ts=200 --initial_type='ave' --extract="obj" --alpha=0.003 \
#  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT.csv"
#python CNOT_acc.py --name=CNOTAccC --evo_time=20 --n_ts=400 --initial_type='ave' --extract="obj" --alpha=0.01 \
#  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT.csv"

#python CNOT_acc.py --name=CNOTAccC --evo_time=5 --n_ts=100 --initial_type='ave' --extract="sur" --alpha=0.02 \
#  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv"
#python CNOT_acc.py --name=CNOTAccC --evo_time=10 --n_ts=200 --initial_type='warm' --extract="sur" --alpha=0.007 \
#  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT.csv"
for i in 0.015;
do
python CNOT_acc.py --name=CNOTAccC --evo_time=20 --n_ts=400 --initial_type='warm' --extract="sur" --alpha=$i \
  --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT.csv";
done

#python CNOT_acc.py --name=CNOTAccD --evo_time=5 --n_ts=100 --initial_type='ave' --extract='binary' \
#    --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv" \
#    --b_control="../new_result/b_control/TR+MT/CNOTSOS1_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python CNOT_acc.py --name=CNOTAccD --evo_time=10 --n_ts=200 --initial_type='ave' --extract='binary' \
#    --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT.csv" \
#    --b_control="../new_result/b_control/TR+MT/CNOTSOS1_evotime10.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python CNOT_acc.py --name=CNOTAccD --evo_time=20 --n_ts=400 --initial_type='ave' --extract='binary' \
#    --c_control="../new_result/c_control/Continuous/CNOTSOS1_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT.csv" \
#    --b_control="../new_result/b_control/TR+MT/CNOTSOS1_evotime20.0_n_ts400_ptypeCONSTANT_offset0.5_objUNIT_alpha0.0001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"

cd ../new_script/