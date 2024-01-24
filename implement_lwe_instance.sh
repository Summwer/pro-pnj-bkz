mkdir "lwechal-instance-test-dd"
cd "lwechal-instance-test-dd"
mkdir "default_g6k_main(32+2gpus)"
mkdir "d4f-default-g6k"
cd "d4f-default-g6k"
mkdir "bssa(32+2gpus)"
mkdir "enumbs(32+2gpus)"
cd ..
cd ..

for j in $(seq 1 10) #do
do
    ###########################################
    # LWE instance solved by default G6K test
    #
    ###########################################
    for i in $(seq 41 50) #do
    do
        python lwe_instance.py ${i} --lwe/alpha 0.015 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" | tee "lwechal-instance-test-dd"/"default_g6k_main(32+2gpus)"/${i}-015-${j}.log
    done


    for i in $(seq 51 60) #do
    do
        python lwe_instance.py $i --lwe/alpha 0.010 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" | tee "lwechal-instance-test-dd"/"default_g6k_main(32+2gpus)"/$i-010-${j}.log
    done


    for i in $(seq 71 80) 
    do
        python lwe_instance.py "$i" --lwe/alpha 0.005 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" | tee "lwechal-instance-test-dd"/"default_g6k_main(32+2gpus)"/"$i"-005-${j}.log
    done


    # ###########################################
    # # EnumBS strategy test:  d4f-default-g6k
    # #
    # ###########################################


    for i in $(seq 41 50) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.015 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "enumbs" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-015-${j}.log
    done


    for i in $(seq 51 60) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "enumbs" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-010-${j}.log
    done

    for i in $(seq 71 80) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "enumbs" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-005-${j}.log
    done



    # # ###########################################
    # # # BSSA strategy test: d4f-default-g6k
    # # #
    # # ###########################################


    for i in $(seq 41 50) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.015 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "bssa" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-015-${j}.log
    done


    for i in $(seq 51 60) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "bssa" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-010-${j}.log
    done

    for i in $(seq 71 80) 
    do
    python ProPnjBKZ_for_lwe_instance.py ${i} --lwe/alpha 0.010 --gpus 2 --threads 32 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --strategy_method "bssa" | tee "lwechal-instance-test-dd"/"d4f-default-g6k"/"enumbs(32+2gpus)"/${i}-005-${j}.log
    done

done