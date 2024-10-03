
******************************
progressive pnj-BKZ
******************************

Please use gcc-8.5 to compile the files, and firstly run

.. code-block:: bash

   git clone https://github.com/cr-marcstevens/parallel-hashmap
   
   PYTHON=python3 ./bootstrap.sh #threads for compiling
   ./rebuild.sh --noyr -j 30 #threads for compiling

Before implement our code, please follow the compile guidance in the topic **G6K - GPU Tensor** We add some files in `G6K - GPU Tensor`(https://github.com/WvanWoerden/G6K-GPU-Tensor) to run a two-step mode for solving u-SVP problem in G6K-GPU with a blocksize selection method. One can generate the reduction strategy and solve the LWE instance through the reduction strategy by running the following command:


.. code-block:: bash

    source ./activate
    python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "enumbs" --load_lwe "lwe_challenge" --max_RAM 43
 
It means that use blocksize strategy generation algorithm EnumBS(three choices: enumbs, bssav1, bssav2), 32 threads and 2 gpus to solve LWE challenge with `(n,alpha) = (40,0.025)` in float_type "dd" with maximal memory 43 log2(bit). Since the gaps of $n$ and $\alpha$ between two LWE challenges are  $5$ and $0.005$, to test more LWE instance, we could randomly generate an LWE instance in specific $n$ and $\alpha$ by setting `--load_lwe ` as `lwe_instance`. 



Besides, we give the guidance of our experiments in the article https://eprint.iacr.org/2022/1343.pdf. 

**Experiment Environment. ** Intel Xeon 5128 16c 32@2.3GHz, 1.48T RAM and NVIDIA Geforce RTX 3090 * 2, Ubuntu 20.04, gcc8, cuda 11.4, NVIDIA-SMI 470.199.02.

By the way, for the practical test, the blocksize and jump strategy is only optimized in our machine with Intel Xeon 5128 16c 32@2.3GHz, 1.48T RAM and NVIDIA Geforce RTX 3090 * 2. If you want to generate your optimized blocksize and jump strategy (but also has a speedup effect) in your machine, you should test the practical cost model for your machine and modify the data in the function `refined_practical_pump_cost_dd` and `refined_PnJBKZ_cost_model_dd` in `cost.cpp` file. For convenience, one can also generate the strategy using a theoretical cost model, just set the `cost_model = 1`.

Now we explain how to generate the experimental data we provided in the supplement materials. One should download the open source code https://github.com/Summwer/pro-pnj-bkz, first compile it following the guidance, then run the corresponding codes in the main directory.



## Figure 5 & Table 2

We run the experiment in `implement_lwechal_forall.sh`(dd float_type) to solve LWE challenge with strategy in default g6k, bssa, or enumbs and obtain the cost information in Figure 2(a) and Figure 2(b). It stores the test result in the folder  `lwechal-test`. We also run an experiment in `implement_lwe_instance_forall.sh` to test the cost of LWE instances among the above three solvers while with growth of $n$ in "dd" float_type, the experiment result is stored in the folder `lwe-instance-test`. It is the source data of Figure 2(c) and Figure 2(d). All the above result shows the cost of T(ProPnJBKZ(EnumBS))<T(ProPnJBKZ(BSSA))<T(default G6K).

Beside in each log of test, we also print the cost for each strategy generation and the detailed strategy. We also list some of them in `Table 2`.

All the data is stored in the folder `Fig.5 cost-comparison`.




## **Figure4&9~22.**

This is the code description document for all verification experiments of the accuracy of PnJBKZ Simulator in Section 4.1.2 (Section: Performance of PnJBKZ simulator). 

One should download the open source code https://github.com/Summwer/pro-pnj-bkz

First, according to the beta, jump, and tours parameters you want to verify, modify the corresponding parameters in the `lwe_challenge_gen_rr.py` file in https://github.com/Summwer/pro-pnj-bkz, and run it 20 times or any number of times you want to experiment by the following command like:

```bash
python lwe_challenge_gen_rr.py 75 --lwe/alpha 0.005 --bkz/jump 9 --pump/down_sieve True --bkz/blocksizes "[95,95,95,95,95,95,95,95,95,95,95,95]" --gpus 2 --threads 32
```

It will record the actual rr value obtained by the lattice basis reduction of PnJBKZ-$(\beta,J)$ at different numbers of tours in the folder `simulator-test`, and then we can draw a figure like Figure 6. 

After obtaining the $\mathsf{rr}$ value, it stores the values in the folder `simulator-test` . Then  enter the folder `simulator-test` and run

```bash
python Fig7_Fig8_Verification_Experiments_of_PnJBKZ_Simulator_2024.py
```

(We've pre-stored the generated rr in `simulator-test`, one can implement the command above directly, we also give the test data and implemented code in the folder `Figure4&9~22`) It will print and output the result of calculating the error between the PnJBKZ simulator simulation value and the actual reduced rr value, which is shown in Fig9. At the same time, the program will draw Fig10 to verify the accuracy of the PnJBKZ simulator under the corresponding reduction parameters. Fig11~Fig22 give more tests about PnJBKZ simulator. For the limit of supplementary materials, we only put partial test results in the folder. For the entire test results, please see  https://github.com/Summwer/pro-pnj-bkz/simulator-test.



## Figure8

To test the difference of failure probability of the Pump Dimension Estimation used in default G6K and our work, we generate 100 randomly LWE instances for each $(n,\alpha)$ by running the command:

```bash
python PumpDimEst_comparison.py
```
and obtain `Figure8`.



## Table3,9-11

In `Table3,9-11`. We compare the actual running process of LWE challenge and our simulated results. One can obtain the simulated results by running

```python
python strategy_simulation.py
```

and the running log of each LWE challenge is in the folder `lwechal-test`.







## Table5. NIST-est

All the logs mentioned have been stored in the folder `Table5. NIST-est`. One should download the open source code https://github.com/Summwer/lwe-estimator-with-PnJBKZ, enter the folder `cpp` first compile it by `rebuild.sh`, then run the following command to obatain a nist-est generated by enumbs:

```bash
./implement_all_NIST_schemes.sh
```

Through the above command, we'll get the file `enumbs(cumprob+prob)+list_decoding[MATZOV22].log`  whose list-decoding complexity is from [[MATZOV22]](MATZOV, “Report on the Security of LWE: Improved Dual Lattice Attack,” Apr.
2022.) in the folder `nist-round3-est-result`.  Besides, we also give the log in the folder `Table5. NIST-est`. It gives the estimation process in detail and give the blocksize and jump strategy generated by EnumBS for each NIST scheme we estimate.



For the column of "Previous", we use adapt https://github.com/lducas/leaky-LWE-Estimator developed by Leo Ducas to a cumulated Gate algorithm as the paper `A Refined Hardness Estimation of LWE in Two-Step Mode` said. Since we need the gate count in matzov22, so one can obtain the result by running the command in folder `sage`:

```bash
sage NIST-pro-bkz-avgbeta-matzov.sage
```

in sage environment.



For the column of "Two-step", one should download the open source code https://github.com/Summwer/lwe-estimator-with-PnJBKZ, and enter the folder `sage/NIST-round3`. Initiate the sage environment by running `sage`, then run

```
load("NIST-two-step-matzov.sage")    
```

and obatain the estimated result in two-step mode with trival reduction strategy and [[MATZOV22]](MATZOV, “Report on the Security of LWE: Improved Dual Lattice Attack,” Apr.
2022.) version.







## Fig. 24~28. Practical Cost Model

In the folder `Fig. 21~25. Practical Cost Model`, we give the cost test result in it. We also give the code for constructing the practical cost model of Pump and PnJBKZ here. Especially, one can test their own practical cost data for cost model by running the command:

```bash
python practical_pump_cost_test.py | tee practical_pump_cost.log
```

 for pump cost model and 

```bash
python practical_PnJBKZ_cost_test.py
```

for PnJBKZ cost model, which will stores in the folder `each_pump_cost_in_pnjbkz` and one can construct the PnJBKZ cost model using the data and running the python file first. Then, add the generated data into  `pnjbkz_cost_model.py` and run it.

Besides, one can construct the practical Pump cost model by running

```bash
python practical_pump_cost.py
```





## Table7. optimize m

We also try to optimize number of LWE samples in our work. One can run the code 

```bash
./optimize_m.sh
```

to obtain the optimized $m$ as Table 9 shown. The log result is stored in the folder in the folder `optimize-m`.

## Table 8 & Table 4

To generate the strategy for our updated LWE records to run them in Table 4, one can run the command

```bash
./implement_unsolved_lwechal.sh
```

and the results are stored in the folder  `strategy-for-unsolved-lwechal`.