
#ifndef FIC_BLOCKWISE_ENTROPY_H
#define FIC_BLOCKWISE_ENTROPY_H

float calculate_entropy(const int block_size, float * block);
float blockwise_entropy(const int pixels_per_block, const int num_blocks, float *img, float * entropies, float percentage_discard);
float blockwise_entropy_standard_deviation(float mean, float * entropies, const int num_blocks);
#endif //FIC_BLOCKWISE_ENTROPY_H