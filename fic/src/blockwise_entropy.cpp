#include <iostream>
#include <map>
#include <math.h>
#include <algorithm>

float round_prec(float v, int precision){
    int r = pow(10.0, precision);
    float val = (int)(v * r );
    return (float)val / r;
}

float calculate_entropy(const int pixels_per_block, float * block){
    float entropy = 0.0;
    std::map<float, int> counts;

    for(int i = 0; i < pixels_per_block; ++i){
        counts[round_prec(block[i], 1)]++;
    }
    std::map<float, int>::iterator it = counts.begin();

    while(it !=counts.end()){
        float p_x = (float)it->second/pixels_per_block;
        if(p_x > 0) entropy-=p_x * log(p_x)/log(2);
        else std::cout << "Impossible negative probability" << std::endl;
        it++;
    }
    return entropy;

}

// Calculate the entropy of each block and store it in entropies
// len(entropies) should be num_blocks
float blockwise_entropy(const int pixels_per_block, const int num_blocks, float *img, float * entropies, float percentage_discard){
    float tot_entropy = 0.0;
    for(int i = 0; i < num_blocks; i++){
        int block_offset = i * pixels_per_block;
        float entropy = calculate_entropy(pixels_per_block, &img[block_offset]);
        entropies[i] = entropy;
        tot_entropy += entropy;
    }
    std::sort(entropies, entropies + num_blocks);
    float threshold = entropies[(int((1 -percentage_discard) * num_blocks))];
    // Return the threshold that limits where a certain percentage will be discarded
    return threshold;
}

float blockwise_entropy_standard_deviation(float mean, float * entropies, const int num_blocks){
    float sum = 0.0;
    for (int i = 0; i < num_blocks; i++){
        float diff = mean - entropies[i];
        sum += diff * diff;
    }
    return sum / num_blocks;
}