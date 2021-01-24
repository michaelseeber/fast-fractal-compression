# Meetings

previous on Telegram


## 16.04.2020
C straight forward implementation - DONE

Simplifiying assumptions:
* Greyscale
* Image size divisible by 8
* ~~Squared images~~ (seems to work with rectangular images after all)

Other input restrictions:
* ~~1024 maximal size~~ (we will most likely not use the old implementation of find_optimal_mappings() which was responsible for allocating too much memory for large images)

Benchmarking Infrastucture:
* ~~blockwise_sum_of_xmuly - TODO~~ (will not be used since it was merged back into the nested loop for better optization flexibility in commit SHA: fe017cea8723689d60a3beb169bf71f2b50fda94)
* ~~find_optimal_mappings - TODO~~ (done in this commit)

Cost measure:
* flops
* sqrt
* div