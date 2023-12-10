mkdir "lwechal-test"
cd "lwechal-test"
mkdir "default_g6k_main(32+2gpus)"
mkdir "d4f-default-g6k"
cd "d4f-default-g6k"
mkdir "bssa(32+2gpus)"
mkdir "enumbs(32+2gpus)"
cd ..
cd ..
###########################################
# LWE challenge solved by default G6K test
#
###########################################

python lwe_challenge.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/40-025.log

python lwe_challenge.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True  | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/45-020.log

python lwe_challenge.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/50-015.log

python lwe_challenge.py 55 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/55-010.log

python lwe_challenge.py 60 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/60-010.log

python lwe_challenge.py 70 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/70-005.log

python lwe_challenge.py 75 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/75-005.log

python lwe_challenge.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/80-005.log

python lwe_challenge.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/40-030.log


###########################################
# EnumBS strategy test:  d4f-default-g6k
#
###########################################

<<<<<<< HEAD

python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(96,8,2)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/40-025.log

python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(89,8,2),(107,8,1)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/45-020.log

python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(89,8,2),(105,8,1)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/50-015.log


python lwe_challenge_last_pump.py 55 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(92,8,2)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/55-010.log


python lwe_challenge_last_pump.py 60 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(91,8,3),(117,8,1)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/60-010.log

python lwe_challenge_last_pump.py 70 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(92,8,3)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/70-005.log

python lwe_challenge_last_pump.py 75 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(109,8,3)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/75-005.log

python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(79,8,1),( 90,8,2),(117,8,2),(117,4,1)]" | tee "lwechal-test"/"bssa(32+2gpus)-d4f-default-g6k"/80-005.log
=======
python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(89,9,1),(114,10,1)]" --pump/dsvp 138 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025.log


python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(89,9,1),(117,10,1),(119,10,1)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030.log


python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(90,9,1),(117,10,1)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020.log

>>>>>>> 96f7bfed91ccf56534cd646ca2de6d56a3d619d8

python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(90,9,1),(114,10,1)]" --pump/dsvp 140 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015.log


###########################################
# BSSA strategy test: d4f-default-g6k
#
###########################################

<<<<<<< HEAD
# ###########################################
# # EnumBS strategy test
# #
# ###########################################


# #J=8,max_loop=5
python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,8,1),(104,8,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-default-g6k"/"J=8,max_loop=5"/40-025.log

python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(79,8,1),(91,8,1),(108,8,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-default-g6k"/"J=8,max_loop=5"/45-020.log

python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(79,8,1),(91,8,1),(106,8,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-default-g6k"/"J=8,max_loop=5"/50-015.log


python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,9,1),(111,11,2),(117,11,2),(120,11,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-default-g6k"/"J=100,max_loop=10"/80-005.log


python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(89,9,1),(90,9,1),(117,10,3),(119,10,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-default-g6k"/"J=100,max_loop=10"/80-005.log
=======
>>>>>>> 96f7bfed91ccf56534cd646ca2de6d56a3d619d8

python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,3)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025.log

python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,7,3),(79,4,1),(104,4,2)]" --pump/dsvp 153 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030.log

<<<<<<< HEAD
python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(106,11,1),(117,12,1),(117,4,1)]" | tee "lwechal-test"/"enumbs(32+2gpus)-d4f-theo2"/"J=100,max_loop=10"/40-030.log
=======
python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 145 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020.log
>>>>>>> 96f7bfed91ccf56534cd646ca2de6d56a3d619d8

python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 144 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015.log

