//namespace reference {
    void blockwise_sum(const int num_blocks, const int block_num_pixels, float *img_blocks, float *img_blockwise_sum) {
        for (int i = 0; i < num_blocks; i++) {
            img_blockwise_sum[i] = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                img_blockwise_sum[i] += img_blocks[i * block_num_pixels + j];
            }
        }
    }
//}
