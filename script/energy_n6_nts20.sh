#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
#python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=1 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance1.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance1_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=2 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance2.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance2_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=3 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance3.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance3_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=4 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance4.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance4_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=5 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance5.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance5_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"

#for i in 0.015;
#do
#  python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=1 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="obj" --alpha=$i\
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance1.csv";
#  python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=2 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="obj" --alpha=$i\
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance2.csv";
#  python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=3 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="obj" --alpha=$i\
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance3.csv";
#  python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=4 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="obj" --alpha=$i\
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance4.csv";
#  python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=5 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="obj" --alpha=$i\
#  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance5.csv";
#done

for i in 0.01;
do
python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=1 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i\
  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance1.csv";
python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=2 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i\
  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance2.csv";
python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=3 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i\
  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance3.csv";
python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=4 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i\
  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance4.csv";
python energy_acc.py --name=EnergyAccC --n=6 --num_edge=3 --rgraph=1 --seed=5 --evo_time=5 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i\
  --c_control="../new_result/c_control/Continuous/Energy6_evotime5.0_n_ts20_ptypeCONSTANT_offset0.5_instance5.csv";
done


cd ../new_script/