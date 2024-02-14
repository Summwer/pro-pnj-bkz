mkdir "lwechal-instance-test-dd"
cd "lwechal-instance-test-dd"
mkdir "default_g6k_main(32+2gpus)"
mkdir "d4f-default-g6k"
cd "d4f-default-g6k"
mkdir "bssa(32+2gpus)"
mkdir "enumbs(32+2gpus)"
cd ..
cd ..


for j in $(seq 5 6) #do
do

    # for i in $(seq 49 49) 
    # do

    # # # # ###########################################
    # # # # # BSSA strategy test: d4f-default-g6k
    # # # # #
    # # # # ###########################################

    # python ProPnjBKZ_for_lwe.py ${i} --lwe/alpha 0.015 --gpus 2 --threads 32 --verbose True --pump/down_sieve True  --strategy_method "bssav1" --load_lwe "lwe_instance" --float_type "dd" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/${i}-015-${j}.log


    # # # ###########################################
    # # # # EnumBS strategy test:  d4f-default-g6k
    # # # #
    # # # ###########################################

    # python ProPnjBKZ_for_lwe.py ${i} --lwe/alpha 0.015 --gpus 2 --threads 32 --verbose True --pump/down_sieve True  --strategy_method "enumbs" --load_lwe "lwe_instance" --float_type "dd"| tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-015-${j}.log

    # # ###########################################
    # # # LWE instance solved by default G6K test
    # # #
    # # ###########################################

    # python lwe_instance.py ${i} --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True  | tee "lwechal-instance-test-dd"/"default_g6k_main(32+2gpus)"/${i}-015-${j}.log

    # done

    ###########################################
    # LWE instance solved by default G6K test
    #
    ###########################################


    for i in $(seq 56 56) #do
    do
    
    # python lwe_instance.py $i --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True --pump/down_sieve True  | tee "lwechal-instance-test-dd"/"default_g6k_main(32+2gpus)"/$i-010-${j}.log

    # # ###########################################
    # # # EnumBS strategy test:  d4f-default-g6k
    # # #
    # # ###########################################

    python ProPnjBKZ_for_lwe.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True  --strategy_method "enumbs" --load_lwe "lwe_instance" --float_type "dd"| tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-010-${j}.log
   
    # # # ###########################################
    # # # # BSSA strategy test: d4f-default-g6k
    # # # #
    # # # ###########################################

    # python ProPnjBKZ_for_lwe.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True  --strategy_method "bssav1" --load_lwe "lwe_instance" --float_type "dd"| tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"bssa(32+2gpus)"/${i}-010-${j}.log

    done

done