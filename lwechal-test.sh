

# mkdir "enumbs(32+2gpus)-d4f-theo2"
# cd "enumbs(32+2gpus)-d4f-theo2"
# mkdir "J=100,max_loop=10"
# cd "J=100,max_loop=10"
# # mkdir "40-030"
# cd ..
# cd ..



# python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,9,1),(111,11,2),(117,11,2),(120,11,1)]" | tee "enumbs(32+2gpus)-d4f-default-g6k"/"enumbs(32+2gpus)-d4f-theo2"/80-005.log


# for i in $(seq 1 10) 
#     do 
#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(106,11,1),(117,12,1),(117,4,1)]" | tee "enumbs(32+2gpus)-d4f-theo2"/"J=100,max_loop=10"/40-030/${i+1}.log
#     done

# python lwe_challenge_last_pump.py 45 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(106,11,1),(116,12,1),(117,4,1),(120,4,1)]" | tee "enumbs(32+2gpus)-d4f-theo2"/"J=100,max_loop=10"/45-025.log


# mkdir "bssa(32+2gpus)-d4f-theo2"
# cd "bssa(32+2gpus)-d4f-theo2"
# mkdir "J=100,max_loop=10"
# cd ..

# python lwe_challenge_last_pump.py 45 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(113,12,2),(117,4,1),(121,4,1)]" | tee "bssa(32+2gpus)-d4f-theo2"/"J=100,max_loop=10"/45-025.log




# cd "enumbs(32+2gpus)-d4f-default-g6k"
# cd "J=100,max_loop=10"
# mkdir "40-030"
# cd "40-030"
# mkdir  "last-pump-theo2"
# cd ..
# cd ..

# for i in $(seq 10);
#     do 
#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(116,12,2),(120,12,1)]" | tee "enumbs(32+2gpus)-d4f-default-g6k"/"J=100,max_loop=10"/40-030/last-pump-theo2/$i.log;
#     done


python lwe_challenge.py 45 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True | tee "default_g6k_main(32+2gpus)"/45-025-to-d.log