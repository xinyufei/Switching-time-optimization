#!/user/bin/env bash

#conda activate qcopt
cd ../example/
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=80 --initial_type='warm' --extract="sur" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts80_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter100.csv"
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=80 --initial_type='warm' --extract="la" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts80_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter100.csv"
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=80 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts80_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter100.csv"
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=512 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts512_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter50.csv"
python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
  --evo_time=4 --n_ts=2048 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts2048_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter50.csv"
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=128 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts128_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter50.csv"
#python Molecule.py --name=MoleculeNew --molecule=H2 --qubit_num=2 \
#  --evo_time=4 --n_ts=256 --initial_type='warm' --extract="obj" --thre_ratio=0 \
#  --target="../result/control/Target/MoleculeNEW_H2_evotime4.0_n_ts80_target.csv" \
#  --c_control="../result/control/ADMM/MoleculeADMMNew_H2_evotime4.0_n_ts256_ptypeWARM_offset0.5_sum_penalty1.0_penalty0.001_ADMM_0.5_iter50.csv"
cd ../script/