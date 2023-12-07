#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
#python NOTleak_acc.py --name=NOTleakAccD --evo_time=2 --n_ts=20 --initial_type='warm' --extract='binary' \
# --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT.csv" \
# --b_control="../new_result/b_control/TR+MT/NOTleakSOS1_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"
#python NOTleak_acc.py --name=NOTleakAccD --evo_time=6 --n_ts=60 --initial_type='warm' --extract='binary' \
# --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT.csv" \
# --b_control="../new_result/b_control/TR+MT/NOTleakSOS1_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"
#python NOTleak_acc.py --name=NOTleakAccD --evo_time=10 --n_ts=100 --initial_type='warm' --extract='binary' \
# --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv" \
# --b_control="../new_result/b_control/TR+MT/NOTleakSOS1_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT_alpha0.001_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup5_0_sigma0.25_eta0.001_threshold30_iter100_typeminup_time5.csv"

#for i in 0.01;
#do
#python NOTleak_acc.py --name=NOTleakAccC --evo_time=2 --n_ts=20 --initial_type='ave' --extract="obj" --alpha=$i \
#  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT.csv";
#done
#for i in 0.0006;
#do
#python NOTleak_acc.py --name=NOTleakAccC --evo_time=6 --n_ts=60 --initial_type='ave' --extract="obj" --alpha=$i \
#  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT.csv";
#done
#for i in 0.004;
#do
#python NOTleak_acc.py --name=NOTleakAccC --evo_time=10 --n_ts=100 --initial_type='warm' --extract="obj" --alpha=$i \
#  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv";
#done

#for i in 0.05;
#do
#python NOTleak_acc.py --name=NOTleakAccC --evo_time=2 --n_ts=20 --initial_type='warm' --extract="sur" --alpha=$i \
#  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime2.0_n_ts20_ptypeCONSTANT_offset0.5_objUNIT.csv";
#done
for i in 0.006;
do
python NOTleak_acc.py --name=NOTleakAccC --evo_time=6 --n_ts=60 --initial_type='warm' --extract="sur" --alpha=$i \
  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime6.0_n_ts60_ptypeCONSTANT_offset0.5_objUNIT.csv";
done
#for i in 0.03;
#do
#python NOTleak_acc.py --name=NOTleakAccC --evo_time=10 --n_ts=100 --initial_type='warm' --extract="sur" --alpha=$i \
#  --c_control="../new_result/c_control/Continuous/NOTleakSOS1_evotime10.0_n_ts100_ptypeCONSTANT_offset0.5_objUNIT.csv";
#done
cd ../acc_script/