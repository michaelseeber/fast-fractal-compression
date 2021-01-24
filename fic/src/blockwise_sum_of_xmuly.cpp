namespace reference {
    void blockwise_sum_of_xmuly(const int x_num_blocks, const int y_num_blocks, const int block_num_pixels, double *x_blocks, double **y_blocks_rotations, double **blockwise_sum_xy) {
        for (int i = 0; i < x_num_blocks; i++) {
            for (int j = 0; j < y_num_blocks; j++) {
                blockwise_sum_xy[0][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[1][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[2][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[3][i * y_num_blocks + j] = 0.0;
                for (int p = 0; p < block_num_pixels; p++) {
                    blockwise_sum_xy[0][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[0][j * block_num_pixels + p];
                    blockwise_sum_xy[1][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[1][j * block_num_pixels + p];
                    blockwise_sum_xy[2][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[2][j * block_num_pixels + p];
                    blockwise_sum_xy[3][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[3][j * block_num_pixels + p];
                }
            }
        }
    }
}

namespace fast {
    void blockwise_sum_of_xmuly(const int x_num_blocks, const int y_num_blocks, const int block_num_pixels, double *x_blocks, double **y_blocks_rotations, double **blockwise_sum_xy) {
        for (int i = 0; i < x_num_blocks; i++) {
            for (int j = 0; j < y_num_blocks; j++) {
                blockwise_sum_xy[0][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[1][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[2][i * y_num_blocks + j] = 0.0;
                blockwise_sum_xy[3][i * y_num_blocks + j] = 0.0;
                for (int p = 0; p < block_num_pixels; p++) {
                    blockwise_sum_xy[0][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[0][j * block_num_pixels + p];
                    blockwise_sum_xy[1][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[1][j * block_num_pixels + p];
                    blockwise_sum_xy[2][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[2][j * block_num_pixels + p];
                    blockwise_sum_xy[3][i * y_num_blocks + j] += x_blocks[i * block_num_pixels + p] * y_blocks_rotations[3][j * block_num_pixels + p];
                }
            }
        }
    }
}
