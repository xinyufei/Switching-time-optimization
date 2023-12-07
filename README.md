# Switching Time Optimization for Binary Control Pulse in Quantum Systems
This repository contains the source code used in the computational experiments of the paper: 
[**Switching Time Optimization for Binary Control Pulse in Quantum Systems**](https://arxiv.org/pdf/2308.03132.pdf). 

In our paper, we first introduce two algorithms to derive controller sequences from discretized continuous controls. 
The algorithms include:
* Extraction based on objective values (Obj)
* Extraction based on cumulative difference (SUR)

Then we propose a switching time optimization model for the binary quantum 
control problem with given controller sequences and a gradient method to solve it. In addition, we design a 
  pre-computing acceleration technique to accelerate the time evolution simulation.

## Citation
If you use our code in your research, please cite our paper:

> [**Switching Time Optimization for Binary Control Pulse in Quantum Systems**](https://arxiv.org/pdf/2308.03132.pdf) <br />
> Xinyu Fei, Lucas T. Brady, Jeffrey Larson, Sven Leyffer, Siqian Shen <br />
> ```
> @article{fei2023switching,
>   title={Switching Time Optimization for Binary Quantum Optimal Control},
>   author={Fei, Xinyu and Brady, Lucas T and Larson, Jeffrey and Leyffer, Sven and Shen, Siqian},
>   journal={arXiv preprint arXiv:2308.03132},
>   year={2023}
> }
> ```

## Installation
### Prerequisites
* Python >= 3.8
* qiskit >= 0.29.0, scipy >= 1.6.2, qutip >= 4.6

### Installation Steps
1. Clone the repository
2. Set up a virtual environment (optional)
3. Install the above dependencies by 
```shell
pip install qiskit scipy qutip
```


## Test instances
In our paper, we test our algorithms on four quantum control instances:
* Energy minimization problem
* CNOT gate estimation problem
* NOT gate estimation problem with energy leakage
* Circuit compilation problem

For each instance, we solve the switching time optimization model with controller 
sequences obtained from the following three methods respectively:
* Discretized control in the paper [**Binary Control Pulse Optimization for Quantum Systems**](https://quantum-journal.org/papers/q-2023-01-04-892/pdf/) <br />
* Extraction based on objective values (Obj)
* Extraction based on cumulative difference (SUR)

## Usage
### Parameter lists
**Required parameters**
* ```--name```: name of the instance
* ```--n```: number of qubits
* ```--evo_time```: evolution time
* ```--n_ts```: number of time steps
* ```--initial_type```: initial type of the control time intervals, including 'warm', 'average', and 'random'
  * 'warm': the initial control time intervals are set as the value of obtained discretized binary control intervals
  * 'average': the initial control time intervals are set as evolution time / number of intervals
  * 'random': the initial control time intervals are set as random values between 0 and evolution time
* ```--extract```: extraction algorithm, including 'binary', 'obj', and 'sur'
    * 'binary': discretized control in the paper [**Binary Control Pulse Optimization for Quantum Systems**](https://quantum-journal.org/papers/q-2023-01-04-892/pdf/)
    * 'obj': extraction based on objective values
    * 'sur': extraction based on cumulative difference
* ```--alpha```: parameter for the extraction algorithm Obj/SUR
* ```--c_control```: file path of the discretized continuous control

**Energy minimization problem only parameters**
* ```--num_edge```: number of edges in the randomly generated graph for Hamiltonian controllers
* ```--rgraph```: whether to generate a random graph for Hamiltonian controllers 
* ```--seed```: random seed for generating the random graph

**Other instances only parameters**
* ```--target```: file path of the target unitary operator
* ```--phase```: phase parameter to compute infidelity objective function

**Molecule compilation problem only parameters**
* ```--molecule```: molecule name, including 'BeH2', 'LiH', 'H2'

**Discretized control parameters**
* ```--b_control```: file path of the discretized binary control for controller sequences

### Stored results
All the control results are stored in the folder ```result/control/SwitchTime/```. All the output control figures are stored in 
```result/figure/SwitchTime```. The output files are stored in ```result/output/Switchtime/```. One can change the 
paths in files to change the positions. 

**Before starting your own experiments, we suggest deleting the above three folders to clear all the existing results.** 

### Example data files
* Discretized continuous control files (input as parameter ```result/control/c_control```).
* Target unitary operator for circuit compilation problem (input as parameter ```result/control/Target```).

Both required files should be **.csv** files. 

### Usage examples
To run an energy minimization problem with 6 qubits, randomly generated graph for Hamiltonian controllers, 
evolution time as 5, time steps as 100, extracting controller sequences by previous discretized binary solution:
```shell
python energy_acc.py --name=EnergyAccD --n=6 --num_edge=3 --rgraph=1 --seed=1 --evo_time=5 --n_ts=100 --initial_type='warm' --extract="binary" \
  --c_control="../result/c_control/Continuous/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance1.csv" \
  --b_control="..result/b_control/TR+MT/Energy6_evotime5.0_n_ts100_ptypeCONSTANT_offset0.5_instance1_alpha0.01_sigma0.25_eta0.001_threshold30_iter100_typetvc_minup10_1_sigma0.25_eta0.001_threshold30_iter100_typeminup_time10.csv"
```

To run circuit compilation problem on molecule BeH2 with evolution time 20, time steps 200, extracting controller
sequences by Method 2 (Obj) with switching penalty 0:
```shell
python Molecule_acc.py --name=MoleculeAccC --molecule=BeH2 --qubit_num=6 --evo_time=20 --n_ts=200 \
  --initial_type='warm' --extract="obj" --alpha=0 \
  --target="../new_result/c_control/Target/BeH2_target.csv" \
  --c_control="../new_result/c_control/Continuous/MoleculeVQENew_BeH2_evotime20.0_n_ts200_ptypeCONSTANT_offset0.5_objUNIT_sum_penalty0.01.csv"
```

## Important modules
* ```example```: different quantum control examples
* ```script```: testing script examples
* ```switchint```: algorithm to obtain controller sequences and solve
switching time optimization model
* ```utils```: utility functions
## Acknowledgement
We thank Dr. Lucas Brady for providing the code used in the paper [**Optimal Protocols in Quantum Annealing and 
QAOA Problems**](https://arxiv.org/pdf/2003.08952.pdf).

We refer to the paper [**Partial Compilation of Variational Algorithms for 
Noisy Intermediate-Scale Quantum Machines**](https://arxiv.org/pdf/1909.07522.pdf) for generating the 
circuit compilation problem. Their code is presented at https://github.com/epiqc/PartialCompilation.

We refer to the paper [**Binary Control Pulse Optimization for Quantum Systems**](https://quantum-journal.org/papers/q-2023-01-04-892/pdf/) and solve the model with 
TV regularizer to obtain discretized continuous control for extracting controller sequences. Their code is 
presented at https://github.com/xinyufei/Quantum-Control-qutip.

## References
[1] Brady, Lucas T., et al. "Optimal protocols in quantum annealing and quantum approximate optimization algorithm problems." _Physical Review Letters_ 126.7 (2021): 070505.

[2] Gokhale, Pranav, et al. "Partial compilation of variational algorithms for noisy intermediate-scale quantum machines." _Proceedings of the 52nd Annual IEEE/ACM International Symposium on Microarchitecture_. 2019.

[3] Fei, Xinyu, et al. "Binary control pulse optimization for quantum systems." _Quantum_ 7 (2023): 892.

## Developers
Xinyu Fei (xinyuf@umich.edu)

## Contact
Xinyu Fei (xinyuf@umich.edu)

Siqian Shen (siqian@umich.edu)