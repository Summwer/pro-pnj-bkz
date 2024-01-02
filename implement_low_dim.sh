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

python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(89,9,1),(114,10,1)]" --pump/dsvp 138 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025.log


python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(89,9,1),(117,10,1),(119,10,1)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030.log



python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(90,9,1),(117,10,1)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020.log


python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(73,8,1),(90,9,1),(114,10,1)]" --pump/dsvp 143 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015.log


#test 40-035
python lwe_challenge_last_pump.py 40 --lwe/alpha 0.035 --threads 55 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(76,8,1),(91,9,1),(117,10,1),(117,4,1),(132,4,1),(141,4,1)]" | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-035.log


###########################################
# BSSA strategy test: d4f-default-g6k
#
###########################################


python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,3)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025.log

python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,7,3),(79,4,1),(104,4,2)]" --pump/dsvp 153 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030.log

python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 145 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020.log

python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 144 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015.log

