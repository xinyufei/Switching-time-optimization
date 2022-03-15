# Switching Time Optimization for Binary Control Pulse in Quantum Systems
This repository contains the source code used in the computational experiments of the paper: 
**Switching Time Optimization for Binary Control Pulse in Quantum Systems**. 

In our paper, we first propose a switching time optimization model for the binary quantum 
control problem with given controller sequences and a gradient method to solve it. 
Then we introduce three algorithms to derive controller sequences from discretized continuous controls. 
The algorithms include:
* Extraction based on linear approximated objective value (ELA)
* Extraction based on objective values (EOV)
* Extraction based on SUR (ESUR)

## Test instances
There are two test instances in the paper:
* Energy minimization problem
* Circuit compilation problem

For each instance, we apply the above three extraction algorithms to obtain controller
sequences and then solve the switching time optimization model. 

## Installation
### Requirements
* Python >= 3.8
* qiskit >= 0.29.0, scipy >= 1.6.2, qutip >= 4.6

## Usage
### Stored results
All the control results are stored in the folder ```result/control/SwitchTime/```. All the output control figures are stored in 
```result/figure/SwitchTime```. The output files are stored in ```result/output/Switchtime/```. One can change the 
paths in files to change the positions. 

**Before starting your own experiments, we suggest deleting the above three folders to clear all the existing results.** 

### Required files
* Discretized continuous control files (input as parameter ```c_control```).
* Target unitary operator for circuit compilation problem (input as parameter ```target```).

Both required files should be **.csv** files. 

### Example code
There are three options for the extraction algorithms. 
* ELA: ```--extract='la'```
* EOV: ```--extract='obj'```
* SUR: ```--extract='sur```

To run an energy minimization problem with 4 qubits, randomly generated graph for Hamiltonian controllers, 
evolution time as 2, time steps as 40, extracting controller sequences by EOV with selecting threshold 0:
```shell
python energy.py --name=Energy --n=4 --num_edge=2 --rgraph=1 --seed=1 --evo_time=2 --n_ts=40 --initial_type='warm' \
  --extract="obj" --thre_ratio=0 \
  --c_control="../result/control/ADMM/EnergyADMM4_evotime2.0_n_ts40_ptypeWARM_offset0.5_penalty0.01_ADMM_10.0_iter100_instance1.csv"
```

To run circuit compilation problem on molecule LiH with evolution time 20, time steps 200, extracting controller
sequences by EOV with selecting threshold 0:
```shell
python Molecule.py --name=MoleculeVQE --molecule=LiH --qubit_num=4 \
  --evo_time=20 --n_ts=200 --initial_type='warm' --extract="obj" --thre_ratio=0 \
  --target="../result/control/Target/MoleculeVQE_LiH_evotime20.0_n_ts200_target.csv" \
  --c_control="../result/control/ADMM/MoleculeVQEADMM_LiH_evotime20.0_n_ts200_ptypeWARM_offset0.5_sum_penalty0.1_penalty0.001_ADMM_3.0_iter100.csv"
```

## Acknowledgement
We thank Dr.Lucas Brady for providing the code used in the paper **Optimal Protocols in Quantum Annealing and 
QAOA Problems** (https://arxiv.org/pdf/2003.08952.pdf).

We refer to the paper **Partial Compilation of Variational Algorithms for 
Noisy Intermediate-Scale Quantum Machines** (https://arxiv.org/pdf/1909.07522.pdf) for generating the 
circuit compilation problem. Their code is presented at https://github.com/epiqc/PartialCompilation.

We refer to the paper **Binary Control Pulse Optimization for Quantum Systems** and solve the model with 
TV regularizer to obtain discretized continuous control for extracting controller sequences. Their code is 
presented at https://github.com/xinyufei/Quantum-Control-qutip.