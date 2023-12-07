#!/user/bin/env bash

#conda activate qcopt
cd ../acc_example/
#python energy.py --name=Energytrmt --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/Continuous/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5.csv" \
#  --b_control="../new_result/b_control/TR+MT/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#python energy.py --name=Energymintv --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv" \
#  --b_control="../new_result/b_control/mintv/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5_minup10_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
#
#python energy.py --name=Energybest --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv" \
#  --b_control="../new_result/b_control/bestsol/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_maxswitch5_sigma0.25_eta0.001_threshold30_iter100_typemaxswitch_switch5.csv"
#
#python energy.py --name=Energysurtv --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="binary" \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv" \
#  --b_control="../new_result/b_control/tvmodel/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_SUR_sigma0.25_eta0.001_threshold30_iter100_typetv.csv"

#python energy.py --name=Energy --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='ave' --extract="obj" --alph=0.01 \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv"
#python energy.py --name=Energy --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='ave' --extract="la" --alph=0.01 \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv"
#python energy.py --name=Energy --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='ave' --extract="sur" --alph=0.01 \
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv"
#python energy.py --name=Energy --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="tv" --alph=0.01 --ngroup=16\
#  --c_control="../new_result/c_control/ADMM/EnergyADMM2_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100.csv"

# for obj and la, 0.1 is the best
#python energy.py --name=EnergyC --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="obj" --alpha=0.1 \
#  --c_control="../new_result/c_control/Continuous/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5.csv"
#python energy.py --name=EnergyC --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="la" --alpha=0.1 \
#  --c_control="../new_result/c_control/Continuous/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5.csv"
# for sur and mw, 0.03 is the best
for i in 0.075;
do
python energy_acc.py --name=EnergyAccC --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='ave' --extract="sur" --alpha=$i \
  --c_control="../new_result/c_control/Continuous/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5.csv";
done
#python energy.py --name=EnergyC --n=2 --num_edge=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="tv" --alpha=0.03 --ngroup=1 \
#  --c_control="../new_result/c_control/Continuous/Energy2_evotime2.0_n_ts40_ptypeCONSTANT_offset0.5.csv"
cd ../new_script/