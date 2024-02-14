mkdir "lwechal-test-dd"
cd "lwechal-test-dd"
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

python lwe_challenge.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd"| tee "lwechal-test-dd"/"default_g6k_main(32+2gpus)"/40-025.log

python lwe_challenge.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd"| tee "lwechal-test-dd"/"default_g6k_main(32+2gpus)"/45-020.log

python lwe_challenge.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd"| tee "lwechal-test-dd"/"default_g6k_main(32+2gpus)"/50-015.log


python lwe_challenge.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd"| tee "lwechal-test-dd"/"default_g6k_main(32+2gpus)"/40-030.log



# ###########################################
# # EnumBS strategy test:  d4f-default-g6k
# #
# ###########################################

python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.025 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --load_lwe "lwe_challenge" --float_type "dd" --strategy_method "enumbs" | tee "lwechal-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025.log

python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --strategy_method "enumbs"  --load_lwe "lwe_challenge" --float_type "dd" | tee "lwechal-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030.log


python ProPnjBKZ_for_lwe.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --load_lwe "lwe_challenge" --float_type "dd" --strategy_method "enumbs" | tee "lwechal-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020.log


python ProPnjBKZ_for_lwe.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --load_lwe "lwe_challenge" --float_type "dd" --strategy_method "enumbs" | tee "lwechal-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015.log




# ###########################################
# # BSSA strategy test: d4f-default-g6k
# #
# ###########################################


python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "bssav1" --load_lwe "lwe_challenge" | tee "lwechal-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025.log

python ProPnjBKZ_for_lwe.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "bssav1" --load_lwe "lwe_challenge" | tee "lwechal-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030.log

python ProPnjBKZ_for_lwe.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "bssav1" --load_lwe "lwe_challenge" | tee "lwechal-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020.log

python ProPnjBKZ_for_lwe.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --float_type "dd" --strategy_method "bssav1" --load_lwe "lwe_challenge" | tee "lwechal-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015.log



