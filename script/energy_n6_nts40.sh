#!/user/bin/env bash

#conda activate qcopt
cd ../example/
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="sur" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM6_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="la" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM6_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM6_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
cd ../script/