
******************************
progressive pnj-BKZ
******************************

Please use gcc-8.5 to compile the files, and firstly run

.. code-block:: bash

   git clone https://github.com/cr-marcstevens/parallel-hashmap
   python=PYTHON3 ./bootstrap.sh -j 30 #threads for compiling
   ./rebuild.sh --noyr -j 30 #threads for compiling

Before implement our code, please follow the compile guidance in the topic **G6K - GPU Tensor** We add some files in `G6K - GPU Tensor`(https://github.com/WvanWoerden/G6K-GPU-Tensor) to run a two-step mode for solving u-SVP problem in G6K-GPU with a blocksize selection method. One can generate the reduction strategy and solve the LWE instance through the reduction strategy by running the following command:


.. code-block:: bash

    source ./activate
    python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "enumbs" --load_lwe "lwe_challenge"

It means that use blocksize strategy generation algorithm EnumBS(three choices: enumbs, bssav1, bssav2), 32 threads and 2 gpus to solve LWE challenge with `(n,alpha) = (40,0.025)` in float_type "dd". 

By the way, for the practical test, the blocksize and jump strategy is only optimized in our machine with Intel Xeon 5128 16c 32@2.3GHz, 1.48T RAM and NVIDIA Geforce RTX 3090 * 2. If you want to generate your optimized blocksize and jump strategy (but also has a speedup effect) in your machine, you should test the practical cost model for you machine and modify the data in the function `get_k1_k2_pump` and `get_k1_k2_pnj` in `cost.cpp` file. For convenience, one can also generate the strategy using a theoretical cost model, just set the `cost_model = 1`.





We run the experiment in `implement_low_dim_qd.sh`(qd float_type) and `implement_low_dim_dd.sh`(dd float_type) to solve LWE challenge with strategy in default g6k, bssa, or enumbs and obtain the cost information in Table 1. It stores the test result in the folder `lwechal-test-qd` and `lwechal-test-dd`. We also run an experiment in `implement_lwe_instance.sh` to test the cost of LWE instances among the above three solvers while with growth of $n$ in "dd" float_type, the experiment result is stored in the folder `lwechal-instance-test-dd`. All the above result shows the cost of T(ProPnjBKZ(EnumBS))<=T(ProPnjBKZ(BSSA))<T(default G6K).




The LWE instance in detail, including the actual basis quality(we regard "slope" as a basis quality representation) and actual cost(in practical cost model with threads = 32 and gpus = 2 * (RTX3090 GPU) (sec)) after each pnj-BKZ(beta,J,tours) reduction, also obatained from the above code. We select the process information of LWE challenge ($40,0.030$)  as 'Practical' column in Table 5. 

.. code-block:: bash

    ./implement_low_dim_qd.sh


Call
.. code-block:: bash

    python module_comparison.py | tee module_comparison.log

could output the norm of $\lVert {\bf{b}}_0 \rVert$ and PSC estimation after pnj-BKZ/pump reduction in the same time cost, we've stored them in the module_comparison.log, it shows the Table 3 and Table 4 in https://eprint.iacr.org/2022/1343.



One can set beta, jump, and tours parameters you want to verify, modify the corresponding parameters in the `lwe_challenge_gen_rr.py` file in https://github.com/Summwer/pro-pnj-bkz, and run it 20 times or any number of times you want to experiment by the following command like:

.. code-block:: bash

python lwe_challenge_gen_rr.py 75 --lwe/alpha 0.005 --bkz/jump 9 --pump/down_sieve True --bkz/blocksizes "[95,95,95,95,95,95,95,95,95,95,95,95]" --gpus 2 --threads 32


It will record the actual rr value obtained by the lattice basis reduction of PnjBKZ-$(\beta,J)$ at different numbers of tours in the folder `simulator-test/gs-lengths-simulator`, and then we can draw Figure4. 

We also implement the experiment about the fitness of pnj-BKZ simulator in `simulator-test` folder, the file `Fig5_Fig6_Verification_Experiments_of_PnjBKZ_Simulator_2024.py` is used to draw the figure about the comparison between actual reduced gs-lengths and our simulated gs-lengths. One can run 

.. code-block:: bash

python Fig5_Fig6_Verification_Experiments_of_PnjBKZ_Simulator_2024.py


in the folder `simulator-test`. (We've pre-stored the generated rr in `simulator-test`, one can implement the command above directly, we also give the test data and implemented code in the folder `Figure4&5&6&10~18`) It will print and output the result of calculating the error between the PnjBKZ simulator simulation value and the actual reduced rr value, which is shown in Fig5. At the same time, the program will draw Fig6 to verify the accuracy of the PnjBKZ simulator under the corresponding reduction parameters. Fig10~Fig18 give more tests about PnjBKZ simulator.


To test the difference of failure probability of the Pump Dimension Estimate used in default G6K and our work, we generate 100 randomly LWE instances for each $(n,\alpha)$ by running the command:

```bash
python PumpDimEst_comparison.py
```

and obtain Figure7.


One can test the cost of each Pump/PnjBKZ by implementing the file `practical_cost_test.py` by running 

.. code-block:: bash

python practical_cost_test.py

The default setting is gpus=2, threads=32. 









******************************
G6K - GPU Tensor
******************************

G6K is an open-source C++ and Python (2) library that implements several Sieve algorithms to be used in more advanced lattice reduction tasks. It follows the stateful machine framework from: 

Martin R. Albrecht and Léo Ducas and Gottfried Herold and Elena Kirshanova and Eamonn W. Postlethwaite and Marc Stevens, 
The General Sieve Kernel and New Records in Lattice Reduction.

The main source is available in `fplll/g6k <https://github.com/fplll/g6k>`__

This fork expands the G6K implementation with GPU, and in particular Tensor Core, accelerated sieves, and is accompanied by the work:

Léo Ducas, Marc Stevens, Wessel van Woerden,
Advanced Lattice Sieving on GPUs, with Tensor Cores, 
Eurocrypt 2021 (`eprint <https://eprint.iacr.org/2021/141.pdf>`__).

Note the this fork has been expanded from a `pretty old commit <https://github.com/fplll/g6k/commit/11e202967bf16ce5fe40258597fed54849e10a69>`__.

The CPU-only version of the BDGL-like sieve has been integrated into the `main g6k repository <https://github.com/fplll/g6k>`__, with further improvements, and we aim for long term maintenance. 
The GPU implementation has been made public in this repository, but with a lower commitment to quality, documentation and maintenance. Nevertheless feel free to create issues in this repository.

Building the library
====================

The code has only been tested on the NVIDIA Turing generation, and might not work on more recent GPUs.

You will need the current master of FPyLLL and a recent version of the CUDA Toolkit. See ``bootstrap.sh`` for creating all dependencies from scratch except for the CUDA Toolkit:

.. code-block:: bash

    ./bootstrap.sh                # once only: creates local python env, builds fplll, fpylll and G6K
    source g6k-env/bin/activate   # for every new shell: activates local python env
    ./rebuild.sh -f -y            # whenever you want to rebuild G6K

Otherwise, you will need fplll and fpylll already installed and build the G6K Cython extension **in place** like so:

.. code-block:: bash

    pip install Cython
    pip install -r requirements.txt
    ./rebuild.sh -f -y

Remove ``-f`` option to compile faster (fewer optimisations). 
The ``-y`` option significantly reduces the memory footprint, but disables the standard cpu-only sieves. See ``rebuild.sh`` for more options.


Code examples
=============

You can run a single svp-challenge instance on a multiple cores and multiple GPUs, for example:

.. code-block:: bash

    ./svp_challenge.py 100 --threads 4 --gpus 1 --verbose

Will run a svp-challenge using 4 CPU threads and a single GPU.

For more details on the parameters used for the `SVP records <https://www.latticechallenge.org/svp-challenge/halloffame.php>`__ see Section 7.2 of the `paper <https://eprint.iacr.org/2021/141.pdf>`__ or ``runchal2.sh``.

BDGL-sieve
----------

The BDGL-like GPU sieve can be enabled by running

.. code-block:: bash

    ./svp_challenge.py 100 --threads 4 --gpus 1 --gpu_bucketer bdgl --verbose

Acknowledgements
================

This project was supported through the European Union PROMETHEUS project (Horizon 2020 Research and Innovation Program, grant 780701), ERC-StGARTICULATE project (no. 947821), and the RCADG-ALGSTRONGCRYPTO project (no. 740972).
