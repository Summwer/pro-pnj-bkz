
******************************
progressive pnj-BKZ
******************************

Please use gcc-8.5 to compile the files, and firstly run

.. code-block:: bash

   git clone https://github.com/cr-marcstevens/parallel-hashmap
   ./rebuild.sh

Before implement our code, please follow the compile guidance in the topic **G6K - GPU Tensor** We add some files in `G6K - GPU Tensor`(https://github.com/WvanWoerden/G6K-GPU-Tensor) to run a two-step mode for solving u-SVP problem in G6K-GPU with a blocksize selection method. One should first generate the blocksize strategy in https://github.com/Summwer/lwe-estimator-with-pnjbkz/tree/main/cpp. Then run the following command:

.. code-block:: bash

    source ./activate
    python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,8,1),(104,8,1)]"

It means that use blocksize strategy `(beta,J,tours) in [(91,8,1),(104,8,1)]`, 32 threads and 2 gpus to solve LWE challenge with `(n,alpha) = (40,0.025)`. 


We've given some instances in `implement_low_dim.sh` to solve LWE challenge with strategy in default g6k, bssa, or enumbs and obtain the cost information in Table 4 in https://eprint.iacr.org/2022/1343.

.. code-block:: bash

    ./implement_low_dim.sh



We test an LWE instance in detail and return the actual basis quality(we regard "slope" as a basis quality representation) and actual cost(in practical cost model with threads = 32 and gpus = 2 * (RTX3090 GPU) (sec)) after each pnj-BKZ(beta,J,tours) reduction, one can run it by the following command

.. code-block:: bash

    ./quality-predict-test.sh

and get the 'Simulation' column for LWE challenge ($40,0.035$) in Table 5 in https://eprint.iacr.org/2022/1343.


Call
.. code-block:: bash

    python module_comparison.py

could output the norm of ||b0|| and PSC estimation after pnj-BKZ/pump reduction in the same time cost and get the Table 1 and Table 2 in https://eprint.iacr.org/2022/1343.


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
