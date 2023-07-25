
cd "lwechal-test"
mkdir "quality-predict-test(32+2gpus)"
cd ..



python lwe_challenge_last_pump.py 40 --lwe/alpha 0.035 --threads 32 --gpus 2 --verbose True --pump/down_sieve True --pump/saturation_error "ignore" --bkz/blocksizes "[(83,8,1),(93,8,1),(108,8,1),(117, 8,1),(119,4,1),(133,4,1)]" | tee "lwechal-test"/"quality-predict-test(32+2gpus)"/40-035.log
