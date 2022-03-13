#!/user/bin/env bash

#conda activate qcopt
cd ../example/
python Molecule.py --name=MoleculeVQE --molecule=LiH --qubit_num=4 \
  --evo_time=20 --n_ts=200 --initial_type='warm' --extract="sur" --thre_ratio=0 \
  --target="../result/control/Target/MoleculeVQE_LiH_evotime20.0_n_ts200_target.csv" \
  --c_control="../result/control/ADMM/MoleculeVQEADMM_LiH_evotime20.0_n_ts200_ptypeWARM_offset0.5_sum_penalty0.1_penalty0.001_ADMM_3.0_iter100.csv"
python Molecule.py --name=MoleculeVQE --molecule=LiH --qubit_num=4 \
  --evo_time=20 --n_ts=200 --initial_type='warm' --extract="la" --thre_ratio=0 \
  --target="../result/control/Target/MoleculeVQE_LiH_evotime20.0_n_ts200_target.csv" \
  --c_control="../result/control/ADMM/MoleculeVQEADMM_LiH_evotime20.0_n_ts200_ptypeWARM_offset0.5_sum_penalty0.1_penalty0.001_ADMM_3.0_iter100.csv"
python Molecule.py --name=MoleculeVQE --molecule=LiH --qubit_num=4 \
  --evo_time=20 --n_ts=200 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --target="../result/control/Target/MoleculeVQE_LiH_evotime20.0_n_ts200_target.csv" \
  --c_control="../result/control/ADMM/MoleculeVQEADMM_LiH_evotime20.0_n_ts200_ptypeWARM_offset0.5_sum_penalty0.1_penalty0.001_ADMM_3.0_iter100.csv"
cd ../script/