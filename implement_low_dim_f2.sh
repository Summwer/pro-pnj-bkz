mkdir "lwechal-test"
cd "lwechal-test"
mkdir "default_g6k_main(32+2gpus)"
mkdir "d4f-default-g6k"
cd "d4f-default-g6k"
mkdir "bssa(32+2gpus)"
mkdir "enumbs(32+2gpus)"
cd ..
# mkdir "d4f-theo2"
# cd "d4f-theo2"
# mkdir "bssa(32+2gpus)"
# mkdir "enumbs(32+2gpus)"
# cd ..
cd ..
###########################################
# LWE challenge solved by default G6K test
#
###########################################

# python lwe_challenge.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/40-025.log

# python lwe_challenge.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True  | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/45-020.log

# python lwe_challenge.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/50-015.log

# python lwe_challenge.py 55 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/55-010.log

# python lwe_challenge.py 60 --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/60-010.log

# python lwe_challenge.py 70 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/70-005.log

# python lwe_challenge.py 75 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/75-005.log

# python lwe_challenge.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/80-005.log

# python lwe_challenge.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True | tee "lwechal-test"/"default_g6k_main(32+2gpus)"/40-030.log



# ###########################################
# # BSSA strategy test: d4f-default-g6k
# #
# ###########################################

# for i in $(seq 5);
#     do
#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[True]" --pump/dsvp 143 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025-$i.log
    

#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(113,4,3)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030-$i.log

#     python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(77,4,4)]" --pump/dsvp 144 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020-$i.log

#     python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(77,4,4)]" --pump/dsvp 143 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015-$i.log

#     python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,7,5),(77,7,2),(113,4,3)]" --pump/dsvp 150 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/80-005-$i.log


#     ###########################################
#     # EnumBS strategy test:  d4f-default-g6k
#     #
#     ###########################################

#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,9,1),(113,11,1)]" --pump/dsvp 138 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025-$i.log


#     python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(104,10,1),(117,11,1),(119,11,1)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030-$i.log



#     python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(112,11,2)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020-$i.log


#     python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(77,7,1),(101,10,1),(111,11,1)]" --pump/dsvp 139 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015-$i.log


#     python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,7,4),(77,7,2),(111,11,1),(117,11,2),(119,11,1)]" --pump/dsvp 152  | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/80-005-$i.log

#     done





###########################################
# BSSA strategy test: d4f-default-g6k
#
###########################################

for i in $(seq 5);
    do
    python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,3)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025-$i.log

    python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,7,3),(79,4,1),(104,4,2)]" --pump/dsvp 153 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030-$i.log

    python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 145 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020-$i.log

    python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(79,8,5)]" --pump/dsvp 144 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015-$i.log

    python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(79,8,4),(89,9,5),(116,4,2)]" --pump/dsvp 151 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/80-005-$i.log


    ###########################################
    # EnumBS strategy test:  d4f-default-g6k
    #
    ###########################################

    python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,9,1),(114,10,1)]" --pump/dsvp 138 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025-$i.log


    python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,8,1),(89,9,1),(117,10,1),(119,10,1)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030-$i.log



    python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,8,1),(90,9,1),(117,10,1)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020-$i.log


    python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,8,1),(90,9,1),(114,10,1)]" --pump/dsvp 140 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015-$i.log


    python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(89,9,1),(90,9,1),(117,10,3),(119,10,1)]" --pump/dsvp 152  | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/80-005-$i.log

    done





# ###########################################
# # BSSA strategy test: theo2
# #
# ###########################################


# # python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[True]" --pump/dsvp 143 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-025.log

# # python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(113,4,3)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/40-030.log

# # python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(77,4,4)]" --pump/dsvp 144 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/45-020.log

# # python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "skip" --bkz/blocksizes "[(77,4,4)]" --pump/dsvp 143 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/50-015.log

# # python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,7,5),(77,7,2),(113,4,3)]" --pump/dsvp 150 | tee "lwechal-test"/"d4f-default-g6k"/"bssa(32+2gpus)"/80-005.log


# ###########################################
# # EnumBS strategy test:  d4f-default-g6k
# #
# ###########################################

# # python lwe_challenge_last_pump.py 40 --lwe/alpha 0.025 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(91,9,1),(113,11,1)]" --pump/dsvp 138 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-025.log


# python lwe_challenge_last_pump.py 40 --lwe/alpha 0.030 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(104,10,1),(117,11,1),(119,11,1)]" --pump/dsvp 152 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/40-030.log



# # python lwe_challenge_last_pump.py 45 --lwe/alpha 0.020 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(112,11,2)]" --pump/dsvp 141 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/45-020.log


# # python lwe_challenge_last_pump.py 50 --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(77,7,1),(101,10,1),(111,11,1)]" --pump/dsvp 139 | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/50-015.log


# # python lwe_challenge_last_pump.py 80 --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(73,7,4),(77,7,2),(111,11,1),(117,11,2),(119,11,1)]" --pump/dsvp 152  | tee "lwechal-test"/"d4f-default-g6k"/"enumbs(32+2gpus)"/80-005.log



