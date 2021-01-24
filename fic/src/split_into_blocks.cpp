//namespace reference {
//correct
    void split_into_blocks(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;

        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i++) {
                    int img_j = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i * width + img_j];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot90(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                int img_i_rot = block_idx_col * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = block_idx_row * block_num_rows;
                    int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot180(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot270(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = block_idx_col * block_num_cols;
                int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    int img_j_rot = block_idx_row * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;

        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i++) {
                    int img_j = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i * width + img_j];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot90(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                int img_i_rot = block_idx_col * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = block_idx_row * block_num_rows;
                    int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot180(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot270(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = block_idx_col * block_num_cols;
                int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    int img_j_rot = block_idx_row * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot90(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot180(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot270(const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = block_idx_col * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = block_idx_row * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    //--------------implementations of splitting image into overlapping blocks = using stride argument that was previously fixed to 16------------
    void split_into_blocks(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;

        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i = block_idx_row * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i++) {
                    int img_j = block_idx_col * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i * width + img_j];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot90(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                int img_i_rot = block_idx_col * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = block_idx_row * block_num_rows;
                    int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot180(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot270(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = block_idx_col * block_num_cols;
                int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    int img_j_rot = block_idx_row * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = block_idx_row * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot90(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot180(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * stride;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    int img_j_rot = block_idx_col * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    
    void split_into_blocks_hflip_rot270(const int stride, const int block_num_rows, const int block_num_cols, const int height, const int width, float *img, float *img_blocks) {
        //const int img_num_blocks_vertical = height / stride - (block_num_cols / stride) + 1;
        //const int img_num_blocks_horizontal = width / stride - (block_num_rows / stride) + 1;
        const int img_num_blocks_vertical =   stride < block_num_cols ? height / stride - (block_num_cols / stride) + 1 : height / stride;
        const int img_num_blocks_horizontal = stride < block_num_rows ? width / stride - (block_num_rows / stride) + 1  : width / stride;
        const int block_num_pixels = block_num_rows * block_num_cols;

        //std::cout << "stride = " << stride << std::endl;
        //std::cout << "img_num_blocks_vertical = " << img_num_blocks_vertical << std::endl;
        //std::cout << "img_num_blocks_horizontal = " << img_num_blocks_horizontal << std::endl;
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            //std::cout << "\tblock_idx_row = " << block_idx_row << std::endl;
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i_rot = block_idx_col * stride;
                //std::cout << "\tblock_offset = " << block_offset << std::endl;
                //std::cout << "\t\tblock_idx_col = " << block_idx_col << std::endl;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //std::cout << "\t\timg_i_rot = " << img_i_rot << std::endl;
                    int img_j_rot = block_idx_row * stride;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        //std::cout << "\t\t\tblock_i = " << block_i << std::endl;
                        //std::cout << "\t\t\tblock_j = " << block_j << std::endl;
                        //std::cout << "\t\t\timg_j_rot = " << img_j_rot << std::endl;
                        //std::cout << "\t\t\timg_blocks[" << block_offset + block_i * block_num_rows + block_j << "] = img[" << img_i_rot * width + img_j_rot << "];" << std::endl;
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
/*
//--------------------WRONG implements same wrong compression as in mannyray repository----------------------------------------------------------
    void split_into_blocks(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;

        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i++) {
                    int img_j = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i * width + img_j];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot90(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                int img_i_rot = block_idx_col * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = block_idx_row * block_num_rows;
                    //int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    int img_j_rot = block_idx_row * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot180(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                int img_i_rot = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                    int img_j_rot = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
    void split_into_blocks_rot270(const int block_num_rows, const int block_num_cols, const int height, const int width, double *img, double *img_blocks) {
        const int img_num_blocks_vertical = height / block_num_cols;
        const int img_num_blocks_horizontal = width / block_num_rows;
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        for(int block_idx_row = 0; block_idx_row < img_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                //int img_i_rot = block_idx_col * block_num_cols;
                //int img_i_rot = (img_num_blocks_horizontal-1 - block_idx_col) * block_num_cols;
                int img_i_rot = block_idx_col * block_num_cols;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i_rot++) {
                    //int img_j_rot = (img_num_blocks_vertical-1 - block_idx_row) * block_num_rows;
                    int img_j_rot = block_idx_row * block_num_rows;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j_rot++) {
                        img_blocks[block_offset + block_i * block_num_rows + block_j] = img[img_i_rot * width + img_j_rot];
                    }
                }
            }
        }
    }
*/
//}
