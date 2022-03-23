#!/user/bin/env bash

#conda activate qcopt
cd ../example/
#python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="sur" --thre_ratio=0 \
#  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
#python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="la" --thre_ratio=0 \
#  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
#python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=32 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts32_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=64 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts64_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=128 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts128_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=256 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts256_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=512 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts512_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=1024 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts1024_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
#python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=2048 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts2048_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter30_instance1.csv"
cd ../script/