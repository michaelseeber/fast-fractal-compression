#include <fic/common.h>
#include <fic/blockwise_sum.h>
#include <fic/blockwise_sum_of_squares.h>
#include <fic/blockwise_sum_of_xmuly.h>
#include <fic/blockwise_entropy.h>
#include <fic/difference_norm.h>
#include <immintrin.h>

#define MAPSTORE 5
#define DISCARD_PERCENTAGE 0.9

namespace reference {
    void find_optimal_mappings(const int block_num_rows,
            const int block_num_cols,
            const int img_gray_num_blocks,
            const int img_small_num_blocks,
            float *img_gray_blocks,
            float **img_small_blocks_rotations,
            float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        //double *gray_blockwise_sumofsquares  = new double[img_gray_num_blocks];
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        //blockwise_sum_of_squares(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sumofsquares);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);
        
        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            for(int i = 0; i < img_small_num_blocks; i++) {    
                float denom = block_num_pixels * small_blockwise_sumofsquares[i] - std::pow(small_blockwise_sums[i], 2);
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums[i] * gray_blockwise_sums[j]) / denom;
                    float brightness = (gray_blockwise_sums[j] * small_blockwise_sumofsquares[i]  -  small_blockwise_sums[i] * graymulsmall_blockwise_sum) / denom;
                    
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }
                    
                    float error;
                    reference::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    //fast::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    
                    if(mappings[j*MAPSTORE + 0] > error) {
                        mappings[j*MAPSTORE + 0] = error;
                        mappings[j*MAPSTORE + 1] = contrast;
                        mappings[j*MAPSTORE + 2] = brightness;
                        mappings[j*MAPSTORE + 3] = rot;
                        mappings[j*MAPSTORE + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;
        
        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_flip(const int block_num_rows,
            const int block_num_cols,
            const int img_gray_num_blocks,
            const int img_small_num_blocks,
            float *img_gray_blocks,
            float **img_small_blocks_rotations,
            float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        //double *gray_blockwise_sumofsquares  = new double[img_gray_num_blocks];
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        //blockwise_sum_of_squares(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sumofsquares);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);
        
        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            for(int i = 0; i < img_small_num_blocks; i++) {    
                float denom = block_num_pixels * small_blockwise_sumofsquares[i] - std::pow(small_blockwise_sums[i], 2);
                for(int rot_flip = 0; rot_flip < 8; rot_flip++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot_flip][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    
                    float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums[i] * gray_blockwise_sums[j]) / denom;
                    float brightness = (gray_blockwise_sums[j] * small_blockwise_sumofsquares[i]  -  small_blockwise_sums[i] * graymulsmall_blockwise_sum) / denom;
                    
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot_flip][img_small_blocks_offset + p] + brightness;
                    }
                    
                    float error;
                    reference::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    
                    if(mappings[j*MAPSTORE + 0] > error) {
                        mappings[j*MAPSTORE + 0] = error;
                        mappings[j*MAPSTORE + 1] = contrast;
                        mappings[j*MAPSTORE + 2] = brightness;
                        mappings[j*MAPSTORE + 3] = rot_flip;
                        mappings[j*MAPSTORE + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;
        
        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }



    void find_optimal_mappings_with_entropy(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        // If the difference in entropies is higher than *threshold*, don't calculate and compare error differences
        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);
        std::cout << "Threshold " << threshold << std::endl;
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            //float grey_entropy = gray_entropy[j];
            for(int i = 0; i < img_small_num_blocks; i++) {

                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                        //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                       continue;

                }
                float denom = block_num_pixels * small_blockwise_sumofsquares[i] - std::pow(small_blockwise_sums[i], 2);
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }

                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums[i] * gray_blockwise_sums[j]) / denom;
                    float brightness = (gray_blockwise_sums[j] * small_blockwise_sumofsquares[i]  -  small_blockwise_sums[i] * graymulsmall_blockwise_sum) / denom;

                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }

                    float error;
                    reference::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    //fast::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);

                    if(mappings[j*MAPSTORE + 0] > error) {
                        mappings[j*MAPSTORE + 0] = error;
                        mappings[j*MAPSTORE + 1] = contrast;
                        mappings[j*MAPSTORE + 2] = brightness;
                        mappings[j*MAPSTORE + 3] = rot;
                        mappings[j*MAPSTORE + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
}

namespace fast {
    void find_optimal_mappings_scalar_replacement(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        //double *gray_blockwise_sumofsquares  = new double[img_gray_num_blocks];
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        //blockwise_sum_of_squares(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sumofsquares);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {

            int img_small_blocks_offset = 0;
            float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i++) {

                float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                float small_blockwise_sums_i = small_blockwise_sums[i];

                float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }

                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;

                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }

                    float error;
                    fast::difference_norm_scalar_replacement(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    //fast::difference_norm(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);

                    if(mappings[j*MAPSTORE + 0] > error) {
                        mappings[j*MAPSTORE + 0] = error;
                        mappings[j*MAPSTORE + 1] = contrast;
                        mappings[j*MAPSTORE + 2] = brightness;
                        mappings[j*MAPSTORE + 3] = rot;
                        mappings[j*MAPSTORE + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_scalar_replacement_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }
                    float error;
                    fast::difference_norm_scalar_replacement_modified(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_flip_scalar_replacement_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_transformations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_transformations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_transformations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int transformation = 0; transformation < 8; transformation++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_transformations[transformation][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_transformations[transformation][img_small_blocks_offset + p] + brightness;
                    }
                    float error;
                    fast::difference_norm_scalar_replacement_modified(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = transformation;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_inline_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //inline difference_norm
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    float error = 0.f;
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_flip_inline_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //inline difference_norm
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int rot_flip = 0; rot_flip < 8; rot_flip++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot_flip][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    float error = 0.f;
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot_flip][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot_flip;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_strength_reduction_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //inline difference_norm
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                //const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                const float inv_denom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) * inv_denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) * inv_denom;
                    float error = 0.f;
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_code_motion_modified(const int block_num_rows,
            const int block_num_cols,
            const int img_gray_num_blocks,
            const int img_small_num_blocks,
            float *img_gray_blocks,
            float **img_small_blocks_rotations,
            float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;
        
        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);
        
        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float gray_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {    
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                const float temp0 = block_num_pixels / denom;
                const float temp1 = gray_blockwise_sums_j / denom;
                const float temp2 = small_blockwise_sums_i / denom;
                const float temp3 = temp1 * small_blockwise_sums_i;
                const float temp4 = temp1 * small_blockwise_sumofsquares_i;
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = temp0 * graymulsmall_blockwise_sum - temp3;
                    const float brightness = temp4  - temp2 * graymulsmall_blockwise_sum;
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }
                    float error = 0.f;
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*MAPSTORE + 0] > error) {
                        mappings[j*MAPSTORE + 0] = error;
                        mappings[j*MAPSTORE + 1] = contrast;
                        mappings[j*MAPSTORE + 2] = brightness;
                        mappings[j*MAPSTORE + 3] = rot;
                        mappings[j*MAPSTORE + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;
        
        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_ilp_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_ilp_modified_reordered(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);
        
        int num_rotations = 4;
        float *img_small_blocks_rotations_reordered = new float[img_small_num_blocks * block_num_pixels * num_rotations];
        for (int i = 0; i < img_small_num_blocks; ++i) {
            for (int p = 0; p < block_num_pixels; ++p) {
                img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + p * num_rotations + 0] = img_small_blocks_rotations[0][i * block_num_pixels + p];
                img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + p * num_rotations + 1] = img_small_blocks_rotations[1][i * block_num_pixels + p];
                img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + p * num_rotations + 2] = img_small_blocks_rotations[2][i * block_num_pixels + p];
                img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + p * num_rotations + 3] = img_small_blocks_rotations[3][i * block_num_pixels + p];
            }
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                int p2 = 0;
                for(int p = 0; p < block_num_pixels; p++, p2 += 4) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 0] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 1] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 2] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 3] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                p2 = 0;
                for(int p = 0; p < block_num_pixels; p++, p2 += 4){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 0] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 1] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 2] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations_reordered[img_small_blocks_offset + p2 + 3] + brightness_rotation3;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels * num_rotations;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_flip_ilp_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations[4][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations[5][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations[6][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations[7][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations[4][img_small_blocks_offset + p] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations[5][img_small_blocks_offset + p] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations[6][img_small_blocks_offset + p] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations[7][img_small_blocks_offset + p] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_flip_ilp_modified_splitloop(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }
                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations[4][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations[5][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations[6][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations[7][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations[4][img_small_blocks_offset + p] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations[5][img_small_blocks_offset + p] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations[6][img_small_blocks_offset + p] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations[7][img_small_blocks_offset + p] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation4) {
                    mappings[j*MAPSTORE + 0] = error_rotation4;
                    mappings[j*MAPSTORE + 1] = contrast_rotation4;
                    mappings[j*MAPSTORE + 2] = brightness_rotation4;
                    mappings[j*MAPSTORE + 3] = 4;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation5) {
                    mappings[j*MAPSTORE + 0] = error_rotation5;
                    mappings[j*MAPSTORE + 1] = contrast_rotation5;
                    mappings[j*MAPSTORE + 2] = brightness_rotation5;
                    mappings[j*MAPSTORE + 3] = 5;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation6) {
                    mappings[j*MAPSTORE + 0] = error_rotation6;
                    mappings[j*MAPSTORE + 1] = contrast_rotation6;
                    mappings[j*MAPSTORE + 2] = brightness_rotation6;
                    mappings[j*MAPSTORE + 3] = 6;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation7) {
                    mappings[j*MAPSTORE + 0] = error_rotation7;
                    mappings[j*MAPSTORE + 1] = contrast_rotation7;
                    mappings[j*MAPSTORE + 2] = brightness_rotation7;
                    mappings[j*MAPSTORE + 3] = 7;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_flip_ilp_modified_splitloop2(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }
                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations[4][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations[5][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations[6][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations[7][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations[4][img_small_blocks_offset + p] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations[5][img_small_blocks_offset + p] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations[6][img_small_blocks_offset + p] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations[7][img_small_blocks_offset + p] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_flip_ilp_modified_reordered(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);
        
        int num_rotations = 8;
        float *img_small_blocks_rotations_reordered = new float[img_small_num_blocks * block_num_pixels * num_rotations];
        for (int i = 0; i < img_small_num_blocks; ++i) {
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    for (int rot = 0; rot < num_rotations; ++rot) {
                        img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row * 8 + col) * num_rotations + rot] =
                            img_small_blocks_rotations[rot][i * block_num_pixels + row * 8 + col];
                    }
                }
            }
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 0] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 1] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 2] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 3] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 4] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 5] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 6] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 7] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 0] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 1] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 2] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 3] + brightness_rotation3;
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 4] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 5] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 6] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations) + 7] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation4) {
                    mappings[j*MAPSTORE + 0] = error_rotation4;
                    mappings[j*MAPSTORE + 1] = contrast_rotation4;
                    mappings[j*MAPSTORE + 2] = brightness_rotation4;
                    mappings[j*MAPSTORE + 3] = 4;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation5) {
                    mappings[j*MAPSTORE + 0] = error_rotation5;
                    mappings[j*MAPSTORE + 1] = contrast_rotation5;
                    mappings[j*MAPSTORE + 2] = brightness_rotation5;
                    mappings[j*MAPSTORE + 3] = 5;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation6) {
                    mappings[j*MAPSTORE + 0] = error_rotation6;
                    mappings[j*MAPSTORE + 1] = contrast_rotation6;
                    mappings[j*MAPSTORE + 2] = brightness_rotation6;
                    mappings[j*MAPSTORE + 3] = 6;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation7) {
                    mappings[j*MAPSTORE + 0] = error_rotation7;
                    mappings[j*MAPSTORE + 1] = contrast_rotation7;
                    mappings[j*MAPSTORE + 2] = brightness_rotation7;
                    mappings[j*MAPSTORE + 3] = 7;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels * num_rotations;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }
    
    void find_optimal_mappings_GCC_vectorization_modified(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        //inline difference_norm
        const int block_num_pixels = block_num_rows * block_num_cols;
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                #pragma GCC unroll 4
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    #pragma GCC ivdep
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    float error = 0.f;
                    #pragma GCC ivdep
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    //horizontal_sum implementations
    //"best" horizontal sum implementation
    float horizontal_sum(__m256 v)
    {
        __m256 x = _mm256_permute2f128_ps(v, v, 1);
        __m256 y = _mm256_add_ps(v, x);
        x = _mm256_shuffle_ps(y, y, _MM_SHUFFLE(2, 3, 0, 1));
        x = _mm256_add_ps(x, y);
        y = _mm256_shuffle_ps(x, x, _MM_SHUFFLE(1, 0, 3, 2));
        return _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(x, y)));
    }

    // //similar to above regarding speed but more instructionsr
    float horizontal_sum1(__m256 x) {
        const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
        const __m128 loQuad = _mm256_castps256_ps128(x);
        const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
        const __m128 loDual = sumQuad;
        const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
        const __m128 sumDual = _mm_add_ps(loDual, hiDual);
        const __m128 lo = sumDual;
        const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
        const __m128 sum = _mm_add_ss(lo, hi);
        return _mm_cvtss_f32(sum);
    }

    // //straight forward initial implementation but slowest because of hadd
    float horizontal_sum2(__m256 a) {
        __m256 t1 = _mm256_hadd_ps(a,a);
        __m256 t2 = _mm256_hadd_ps(t1,t1);
        __m128 t3 = _mm256_extractf128_ps(t2,1);
        __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2),t3);
        return _mm_cvtss_f32(t4);        
    }

#if __AVX2__
    void find_optimal_mappings_avx(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                
                }

                // //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);

                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_avx_splitloop(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                
                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                
                }

                // //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
            
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
            


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
            
                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);


                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);


                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);


                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);

                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                
                }

                // //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
    
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
        
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


              
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

               
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

            
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

            
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
            
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);

            
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_avx_hsuminlined(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                
                }

                //inlined horizontal sum
                __m256 hsum0_x = _mm256_permute2f128_ps(graymulsmall_blockwise_sum_rotation0_row, graymulsmall_blockwise_sum_rotation0_row, 1);
                __m256 hsum0_y = _mm256_add_ps(graymulsmall_blockwise_sum_rotation0_row, hsum0_x);
                hsum0_x = _mm256_shuffle_ps(hsum0_y, hsum0_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum0_x = _mm256_add_ps(hsum0_x, hsum0_y);
                hsum0_y = _mm256_shuffle_ps(hsum0_x, hsum0_x, _MM_SHUFFLE(1, 0, 3, 2));
                float graymulsmall_blockwise_sum_rotation0 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum0_x, hsum0_y)));

                __m256 hsum1_x = _mm256_permute2f128_ps(graymulsmall_blockwise_sum_rotation1_row, graymulsmall_blockwise_sum_rotation1_row, 1);
                __m256 hsum1_y = _mm256_add_ps(graymulsmall_blockwise_sum_rotation1_row, hsum1_x);
                hsum1_x = _mm256_shuffle_ps(hsum1_y, hsum1_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum1_x = _mm256_add_ps(hsum1_x, hsum1_y);
                hsum1_y = _mm256_shuffle_ps(hsum1_x, hsum1_x, _MM_SHUFFLE(1, 0, 3, 2));
                float graymulsmall_blockwise_sum_rotation1 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum1_x, hsum1_y)));

                __m256 hsum2_x = _mm256_permute2f128_ps(graymulsmall_blockwise_sum_rotation2_row, graymulsmall_blockwise_sum_rotation2_row, 1);
                __m256 hsum2_y = _mm256_add_ps(graymulsmall_blockwise_sum_rotation2_row, hsum2_x);
                hsum2_x = _mm256_shuffle_ps(hsum2_y, hsum2_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum2_x = _mm256_add_ps(hsum2_x, hsum2_y);
                hsum2_y = _mm256_shuffle_ps(hsum2_x, hsum2_x, _MM_SHUFFLE(1, 0, 3, 2));
                float graymulsmall_blockwise_sum_rotation2 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum2_x, hsum2_y)));

                __m256 hsum3_x = _mm256_permute2f128_ps(graymulsmall_blockwise_sum_rotation3_row, graymulsmall_blockwise_sum_rotation3_row, 1);
                __m256 hsum3_y = _mm256_add_ps(graymulsmall_blockwise_sum_rotation3_row, hsum3_x);
                hsum3_x = _mm256_shuffle_ps(hsum3_y, hsum3_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum3_x = _mm256_add_ps(hsum3_x, hsum3_y);
                hsum3_y = _mm256_shuffle_ps(hsum3_x, hsum3_x, _MM_SHUFFLE(1, 0, 3, 2));
                float graymulsmall_blockwise_sum_rotation3 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum3_x, hsum3_y)));
                        
                        

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

    
                //inlined horizontal sum
                hsum0_x = _mm256_permute2f128_ps(error_rotation0_vec, error_rotation0_vec, 1);
                hsum0_y = _mm256_add_ps(error_rotation0_vec, hsum0_x);
                hsum0_x = _mm256_shuffle_ps(hsum0_y, hsum0_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum0_x = _mm256_add_ps(hsum0_x, hsum0_y);
                hsum0_y = _mm256_shuffle_ps(hsum0_x, hsum0_x, _MM_SHUFFLE(1, 0, 3, 2));
                float error_rotation0 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum0_x, hsum0_y)));

                hsum1_x = _mm256_permute2f128_ps(error_rotation1_vec, error_rotation1_vec, 1);
                hsum1_y = _mm256_add_ps(error_rotation1_vec, hsum1_x);
                hsum1_x = _mm256_shuffle_ps(hsum1_y, hsum1_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum1_x = _mm256_add_ps(hsum1_x, hsum1_y);
                hsum1_y = _mm256_shuffle_ps(hsum1_x, hsum1_x, _MM_SHUFFLE(1, 0, 3, 2));
                float error_rotation1 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum1_x, hsum1_y)));

                hsum2_x = _mm256_permute2f128_ps(error_rotation2_vec, error_rotation2_vec, 1);
                hsum2_y = _mm256_add_ps(error_rotation2_vec, hsum2_x);
                hsum2_x = _mm256_shuffle_ps(hsum2_y, hsum2_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum2_x = _mm256_add_ps(hsum2_x, hsum2_y);
                hsum2_y = _mm256_shuffle_ps(hsum2_x, hsum2_x, _MM_SHUFFLE(1, 0, 3, 2));
                float error_rotation2 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum2_x, hsum2_y)));

                hsum3_x = _mm256_permute2f128_ps(error_rotation3_vec, error_rotation3_vec, 1);
                hsum3_y = _mm256_add_ps(error_rotation3_vec, hsum3_x);
                hsum3_x = _mm256_shuffle_ps(hsum3_y, hsum3_y, _MM_SHUFFLE(2, 3, 0, 1));
                hsum3_x = _mm256_add_ps(hsum3_x, hsum3_y);
                hsum3_y = _mm256_shuffle_ps(hsum3_x, hsum3_x, _MM_SHUFFLE(1, 0, 3, 2));
                float error_rotation3 = _mm_cvtss_f32(_mm256_castps256_ps128(_mm256_add_ps(hsum3_x, hsum3_y)));


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_avx_althsum(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum1(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum1(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum1(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum1(graymulsmall_blockwise_sum_rotation3_row);

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum1(error_rotation0_vec);
                float error_rotation1 = horizontal_sum1(error_rotation1_vec);
                float error_rotation2 = horizontal_sum1(error_rotation2_vec);
                float error_rotation3 = horizontal_sum1(error_rotation3_vec);


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_flip_avx(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation4_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation5_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation6_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation7_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                    graymulsmall_blockwise_sum_rotation4_row =  _mm256_fmadd_ps(img_small_blocks_rotation4_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation4_row);
                    graymulsmall_blockwise_sum_rotation5_row =  _mm256_fmadd_ps(img_small_blocks_rotation5_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation5_row);
                    graymulsmall_blockwise_sum_rotation6_row =  _mm256_fmadd_ps(img_small_blocks_rotation6_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation6_row);
                    graymulsmall_blockwise_sum_rotation7_row =  _mm256_fmadd_ps(img_small_blocks_rotation7_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation7_row);                
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);
                float graymulsmall_blockwise_sum_rotation4 = horizontal_sum(graymulsmall_blockwise_sum_rotation4_row);
                float graymulsmall_blockwise_sum_rotation5 = horizontal_sum(graymulsmall_blockwise_sum_rotation5_row);
                float graymulsmall_blockwise_sum_rotation6 = horizontal_sum(graymulsmall_blockwise_sum_rotation6_row);
                float graymulsmall_blockwise_sum_rotation7 = horizontal_sum(graymulsmall_blockwise_sum_rotation7_row);
                

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);
                __m256 contrast_rotation4_vec = _mm256_set1_ps(contrast_rotation4);
                __m256 contrast_rotation5_vec = _mm256_set1_ps(contrast_rotation5);
                __m256 contrast_rotation6_vec = _mm256_set1_ps(contrast_rotation6);
                __m256 contrast_rotation7_vec = _mm256_set1_ps(contrast_rotation7);
                

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);
                __m256 brightness_rotation4_vec = _mm256_set1_ps(brightness_rotation4);
                __m256 brightness_rotation5_vec = _mm256_set1_ps(brightness_rotation5);
                __m256 brightness_rotation6_vec = _mm256_set1_ps(brightness_rotation6);
                __m256 brightness_rotation7_vec = _mm256_set1_ps(brightness_rotation7);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation4_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation5_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation6_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation7_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);
                    __m256 x1_4 = _mm256_fmadd_ps(img_small_blocks_rotation4_row, contrast_rotation4_vec, brightness_rotation4_vec);
                    __m256 x1_5 = _mm256_fmadd_ps(img_small_blocks_rotation5_row, contrast_rotation5_vec, brightness_rotation5_vec);
                    __m256 x1_6 = _mm256_fmadd_ps(img_small_blocks_rotation6_row, contrast_rotation6_vec, brightness_rotation6_vec);
                    __m256 x1_7 = _mm256_fmadd_ps(img_small_blocks_rotation7_row, contrast_rotation7_vec, brightness_rotation7_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);
                    __m256 difference_4_vec = _mm256_sub_ps(x1_4, x2);
                    __m256 difference_5_vec = _mm256_sub_ps(x1_5, x2);
                    __m256 difference_6_vec = _mm256_sub_ps(x1_6, x2);
                    __m256 difference_7_vec = _mm256_sub_ps(x1_7, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                    error_rotation4_vec = _mm256_fmadd_ps(difference_4_vec, difference_4_vec, error_rotation4_vec);
                    error_rotation5_vec = _mm256_fmadd_ps(difference_5_vec, difference_5_vec, error_rotation5_vec);
                    error_rotation6_vec = _mm256_fmadd_ps(difference_6_vec, difference_6_vec, error_rotation6_vec);
                    error_rotation7_vec = _mm256_fmadd_ps(difference_7_vec, difference_7_vec, error_rotation7_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);
                float error_rotation4 = horizontal_sum(error_rotation4_vec);
                float error_rotation5 = horizontal_sum(error_rotation5_vec);
                float error_rotation6 = horizontal_sum(error_rotation6_vec);
                float error_rotation7 = horizontal_sum(error_rotation7_vec);


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_flip_avx_splitloop(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);
                   
                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);
                

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);
                

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);


                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);;

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation4_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation5_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation6_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation7_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation4_row =  _mm256_fmadd_ps(img_small_blocks_rotation4_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation4_row);
                    graymulsmall_blockwise_sum_rotation5_row =  _mm256_fmadd_ps(img_small_blocks_rotation5_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation5_row);
                    graymulsmall_blockwise_sum_rotation6_row =  _mm256_fmadd_ps(img_small_blocks_rotation6_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation6_row);
                    graymulsmall_blockwise_sum_rotation7_row =  _mm256_fmadd_ps(img_small_blocks_rotation7_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation7_row);                
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation4 = horizontal_sum(graymulsmall_blockwise_sum_rotation4_row);
                float graymulsmall_blockwise_sum_rotation5 = horizontal_sum(graymulsmall_blockwise_sum_rotation5_row);
                float graymulsmall_blockwise_sum_rotation6 = horizontal_sum(graymulsmall_blockwise_sum_rotation6_row);
                float graymulsmall_blockwise_sum_rotation7 = horizontal_sum(graymulsmall_blockwise_sum_rotation7_row);
                

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;

                __m256 contrast_rotation4_vec = _mm256_set1_ps(contrast_rotation4);
                __m256 contrast_rotation5_vec = _mm256_set1_ps(contrast_rotation5);
                __m256 contrast_rotation6_vec = _mm256_set1_ps(contrast_rotation6);
                __m256 contrast_rotation7_vec = _mm256_set1_ps(contrast_rotation7);
                

                __m256 brightness_rotation4_vec = _mm256_set1_ps(brightness_rotation4);
                __m256 brightness_rotation5_vec = _mm256_set1_ps(brightness_rotation5);
                __m256 brightness_rotation6_vec = _mm256_set1_ps(brightness_rotation6);
                __m256 brightness_rotation7_vec = _mm256_set1_ps(brightness_rotation7);

                __m256 error_rotation4_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation5_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation6_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation7_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    __m256 x1_4 = _mm256_fmadd_ps(img_small_blocks_rotation4_row, contrast_rotation4_vec, brightness_rotation4_vec);
                    __m256 x1_5 = _mm256_fmadd_ps(img_small_blocks_rotation5_row, contrast_rotation5_vec, brightness_rotation5_vec);
                    __m256 x1_6 = _mm256_fmadd_ps(img_small_blocks_rotation6_row, contrast_rotation6_vec, brightness_rotation6_vec);
                    __m256 x1_7 = _mm256_fmadd_ps(img_small_blocks_rotation7_row, contrast_rotation7_vec, brightness_rotation7_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_4_vec = _mm256_sub_ps(x1_4, x2);
                    __m256 difference_5_vec = _mm256_sub_ps(x1_5, x2);
                    __m256 difference_6_vec = _mm256_sub_ps(x1_6, x2);
                    __m256 difference_7_vec = _mm256_sub_ps(x1_7, x2);

                    error_rotation4_vec = _mm256_fmadd_ps(difference_4_vec, difference_4_vec, error_rotation4_vec);
                    error_rotation5_vec = _mm256_fmadd_ps(difference_5_vec, difference_5_vec, error_rotation5_vec);
                    error_rotation6_vec = _mm256_fmadd_ps(difference_6_vec, difference_6_vec, error_rotation6_vec);
                    error_rotation7_vec = _mm256_fmadd_ps(difference_7_vec, difference_7_vec, error_rotation7_vec);
                
                }

                //ugly horizontal sum now again:)
                float error_rotation4 = horizontal_sum(error_rotation4_vec);
                float error_rotation5 = horizontal_sum(error_rotation5_vec);
                float error_rotation6 = horizontal_sum(error_rotation6_vec);
                float error_rotation7 = horizontal_sum(error_rotation7_vec);

                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_flip_avx_splitloop_with_entropy(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                        //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                       continue;

                }

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);
                   
                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);
                

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);
                

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);


                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);;

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i+=1) {

                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                        //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                       continue;

                }

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation4_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation5_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation6_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation7_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation4_row =  _mm256_fmadd_ps(img_small_blocks_rotation4_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation4_row);
                    graymulsmall_blockwise_sum_rotation5_row =  _mm256_fmadd_ps(img_small_blocks_rotation5_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation5_row);
                    graymulsmall_blockwise_sum_rotation6_row =  _mm256_fmadd_ps(img_small_blocks_rotation6_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation6_row);
                    graymulsmall_blockwise_sum_rotation7_row =  _mm256_fmadd_ps(img_small_blocks_rotation7_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation7_row);                
                }

                //ugly horizontal sum now :) 
                float graymulsmall_blockwise_sum_rotation4 = horizontal_sum(graymulsmall_blockwise_sum_rotation4_row);
                float graymulsmall_blockwise_sum_rotation5 = horizontal_sum(graymulsmall_blockwise_sum_rotation5_row);
                float graymulsmall_blockwise_sum_rotation6 = horizontal_sum(graymulsmall_blockwise_sum_rotation6_row);
                float graymulsmall_blockwise_sum_rotation7 = horizontal_sum(graymulsmall_blockwise_sum_rotation7_row);
                

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;

                __m256 contrast_rotation4_vec = _mm256_set1_ps(contrast_rotation4);
                __m256 contrast_rotation5_vec = _mm256_set1_ps(contrast_rotation5);
                __m256 contrast_rotation6_vec = _mm256_set1_ps(contrast_rotation6);
                __m256 contrast_rotation7_vec = _mm256_set1_ps(contrast_rotation7);
                

                __m256 brightness_rotation4_vec = _mm256_set1_ps(brightness_rotation4);
                __m256 brightness_rotation5_vec = _mm256_set1_ps(brightness_rotation5);
                __m256 brightness_rotation6_vec = _mm256_set1_ps(brightness_rotation6);
                __m256 brightness_rotation7_vec = _mm256_set1_ps(brightness_rotation7);

                __m256 error_rotation4_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation5_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation6_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation7_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation4_row = _mm256_load_ps(&img_small_blocks_rotations[4][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation5_row = _mm256_load_ps(&img_small_blocks_rotations[5][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation6_row = _mm256_load_ps(&img_small_blocks_rotations[6][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation7_row = _mm256_load_ps(&img_small_blocks_rotations[7][img_small_blocks_offset+row_offset]);

                    __m256 x1_4 = _mm256_fmadd_ps(img_small_blocks_rotation4_row, contrast_rotation4_vec, brightness_rotation4_vec);
                    __m256 x1_5 = _mm256_fmadd_ps(img_small_blocks_rotation5_row, contrast_rotation5_vec, brightness_rotation5_vec);
                    __m256 x1_6 = _mm256_fmadd_ps(img_small_blocks_rotation6_row, contrast_rotation6_vec, brightness_rotation6_vec);
                    __m256 x1_7 = _mm256_fmadd_ps(img_small_blocks_rotation7_row, contrast_rotation7_vec, brightness_rotation7_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_4_vec = _mm256_sub_ps(x1_4, x2);
                    __m256 difference_5_vec = _mm256_sub_ps(x1_5, x2);
                    __m256 difference_6_vec = _mm256_sub_ps(x1_6, x2);
                    __m256 difference_7_vec = _mm256_sub_ps(x1_7, x2);

                    error_rotation4_vec = _mm256_fmadd_ps(difference_4_vec, difference_4_vec, error_rotation4_vec);
                    error_rotation5_vec = _mm256_fmadd_ps(difference_5_vec, difference_5_vec, error_rotation5_vec);
                    error_rotation6_vec = _mm256_fmadd_ps(difference_6_vec, difference_6_vec, error_rotation6_vec);
                    error_rotation7_vec = _mm256_fmadd_ps(difference_7_vec, difference_7_vec, error_rotation7_vec);
                
                }

                //ugly horizontal sum now again:)
                float error_rotation4 = horizontal_sum(error_rotation4_vec);
                float error_rotation5 = horizontal_sum(error_rotation5_vec);
                float error_rotation6 = horizontal_sum(error_rotation6_vec);
                float error_rotation7 = horizontal_sum(error_rotation7_vec);

                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_avx_reordered(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        //blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        //blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        //blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        int num_rotations = 4;
        float *img_small_blocks_rotations_reordered = new float[img_small_num_blocks * block_num_pixels * num_rotations];
        for (int i = 0; i < img_small_num_blocks; ++i) {
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    for (int rot = 0; rot < num_rotations; ++rot) {
                        img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row * 8 + col) * num_rotations + rot] =
                            img_small_blocks_rotations[rot][i * block_num_pixels + row * 8 + col];
                    }
                }
            }
        }

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            //NOT USED (unroll twice to utilize 8 floats)
            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

         
                __m256 graymulsmall_blockwise_sum_rotation_vec0 = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
                __m256 graymulsmall_blockwise_sum_rotation_vec1 = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
                __m256 graymulsmall_blockwise_sum_rotation_vec2 = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
                __m256 graymulsmall_blockwise_sum_rotation_vec3 = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

                for(int row = 0; row < 8; row++)
                {
                    int row_offset = row*8;

                    __m256 img_gray_blocks_vec0 = _mm256_set_ps(img_gray_blocks[img_gray_blocks_offset + row_offset + 1], img_gray_blocks[img_gray_blocks_offset + row_offset + 1], img_gray_blocks[img_gray_blocks_offset + row_offset + 1], img_gray_blocks[img_gray_blocks_offset + row_offset + 1], img_gray_blocks[img_gray_blocks_offset + row_offset + 0], img_gray_blocks[img_gray_blocks_offset + row_offset + 0], img_gray_blocks[img_gray_blocks_offset + row_offset + 0], img_gray_blocks[img_gray_blocks_offset + row_offset + 0]);
                    __m256 img_gray_blocks_vec1 = _mm256_set_ps(img_gray_blocks[img_gray_blocks_offset + row_offset + 3], img_gray_blocks[img_gray_blocks_offset + row_offset + 3], img_gray_blocks[img_gray_blocks_offset + row_offset + 3], img_gray_blocks[img_gray_blocks_offset + row_offset + 3], img_gray_blocks[img_gray_blocks_offset + row_offset + 2], img_gray_blocks[img_gray_blocks_offset + row_offset + 2], img_gray_blocks[img_gray_blocks_offset + row_offset + 2], img_gray_blocks[img_gray_blocks_offset + row_offset + 2]);
                    __m256 img_gray_blocks_vec2 = _mm256_set_ps(img_gray_blocks[img_gray_blocks_offset + row_offset + 5], img_gray_blocks[img_gray_blocks_offset + row_offset + 5], img_gray_blocks[img_gray_blocks_offset + row_offset + 5], img_gray_blocks[img_gray_blocks_offset + row_offset + 5], img_gray_blocks[img_gray_blocks_offset + row_offset + 4], img_gray_blocks[img_gray_blocks_offset + row_offset + 4], img_gray_blocks[img_gray_blocks_offset + row_offset + 4], img_gray_blocks[img_gray_blocks_offset + row_offset + 4]);
                    __m256 img_gray_blocks_vec3 = _mm256_set_ps(img_gray_blocks[img_gray_blocks_offset + row_offset + 7], img_gray_blocks[img_gray_blocks_offset + row_offset + 7], img_gray_blocks[img_gray_blocks_offset + row_offset + 7], img_gray_blocks[img_gray_blocks_offset + row_offset + 7], img_gray_blocks[img_gray_blocks_offset + row_offset + 6], img_gray_blocks[img_gray_blocks_offset + row_offset + 6], img_gray_blocks[img_gray_blocks_offset + row_offset + 6], img_gray_blocks[img_gray_blocks_offset + row_offset + 6]);


                    __m256 img_small_blocks_vec0 = _mm256_load_ps(&img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row_offset + 0) * num_rotations]);
                    __m256 img_small_blocks_vec1 = _mm256_load_ps(&img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row_offset + 2) * num_rotations]);
                    __m256 img_small_blocks_vec2 = _mm256_load_ps(&img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row_offset + 4) * num_rotations]);
                    __m256 img_small_blocks_vec3 = _mm256_load_ps(&img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row_offset + 6) * num_rotations]);

                    graymulsmall_blockwise_sum_rotation_vec0 =  _mm256_fmadd_ps(img_small_blocks_vec0, img_gray_blocks_vec0, graymulsmall_blockwise_sum_rotation_vec0);
                    graymulsmall_blockwise_sum_rotation_vec1 =  _mm256_fmadd_ps(img_small_blocks_vec1, img_gray_blocks_vec1, graymulsmall_blockwise_sum_rotation_vec1);
                    graymulsmall_blockwise_sum_rotation_vec2 =  _mm256_fmadd_ps(img_small_blocks_vec2, img_gray_blocks_vec2, graymulsmall_blockwise_sum_rotation_vec2);
                    graymulsmall_blockwise_sum_rotation_vec3 =  _mm256_fmadd_ps(img_small_blocks_vec3, img_gray_blocks_vec3, graymulsmall_blockwise_sum_rotation_vec3);
                }
                graymulsmall_blockwise_sum_rotation_vec0 = _mm256_add_ps(graymulsmall_blockwise_sum_rotation_vec0, graymulsmall_blockwise_sum_rotation_vec1);
                graymulsmall_blockwise_sum_rotation_vec2 = _mm256_add_ps(graymulsmall_blockwise_sum_rotation_vec2, graymulsmall_blockwise_sum_rotation_vec3);

                graymulsmall_blockwise_sum_rotation_vec0 = _mm256_add_ps(graymulsmall_blockwise_sum_rotation_vec0, graymulsmall_blockwise_sum_rotation_vec2);

                __m128 graymulsmall_blockwise_sum_rotation_vec01 = _mm256_castps256_ps128(graymulsmall_blockwise_sum_rotation_vec0);
                __m128 graymulsmall_blockwise_sum_rotation_vec23 = _mm256_extractf128_ps(graymulsmall_blockwise_sum_rotation_vec0, 1);
                __m128 graymulsmall_blockwise_sum = _mm_add_ps(graymulsmall_blockwise_sum_rotation_vec01, graymulsmall_blockwise_sum_rotation_vec23);
                float *graysums_tmp  = new float[4]; //y
                _mm_store_ps(graysums_tmp, graymulsmall_blockwise_sum);
                //std::cout << "graysums_tmp: " << graysums_tmp[0] << ", " << graysums_tmp[1] << ", " << graysums_tmp[2] << ", " << graysums_tmp[3] << std::endl;
                //return;
                float graymulsmall_blockwise_sum_rotation0 = graysums_tmp[0];
                float graymulsmall_blockwise_sum_rotation1 = graysums_tmp[1];
                float graymulsmall_blockwise_sum_rotation2 = graysums_tmp[2];
                float graymulsmall_blockwise_sum_rotation3 = graysums_tmp[3];

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {   
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);
                
                }

                //ugly horizontal sum now again:) 
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);


                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_flip_avx_all_transformations_in_vector(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        //blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        //blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        //blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        int num_rotations = 8;
        float *img_small_blocks_rotations_reordered = new float[img_small_num_blocks * block_num_pixels * num_rotations];
        for (int i = 0; i < img_small_num_blocks; ++i) {
            for (int row = 0; row < 8; ++row) {
                for (int col = 0; col < 8; ++col) {
                    for (int rot = 0; rot < num_rotations; ++rot) {
                        img_small_blocks_rotations_reordered[i * num_rotations * block_num_pixels + (row * 8 + col) * num_rotations + rot] =
                            img_small_blocks_rotations[rot][i * block_num_pixels + row * 8 + col];
                    }
                }
            }
        }


        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - (small_blockwise_sums_i*small_blockwise_sums_i);

                //Holds all 8 transformations
                __m256 graymulsmall_blockwise_sum_rotation_vec = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

                for(int p = 0; p < block_num_pixels; p++)
                {
                    __m256 img_gray_blocks_vec =  _mm256_set1_ps(img_gray_blocks[img_gray_blocks_offset + p]);

                    //load all 8 transforms
                    __m256 img_small_blocks_vec = _mm256_load_ps(&img_small_blocks_rotations_reordered[img_small_blocks_offset + (p * num_rotations)]);

                    graymulsmall_blockwise_sum_rotation_vec =  _mm256_fmadd_ps(img_small_blocks_vec, img_gray_blocks_vec, graymulsmall_blockwise_sum_rotation_vec);
                }
                
                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                __m256 denom_vec = _mm256_set1_ps(1/denom);

                __m256 block_num_pixels_vec = _mm256_set1_ps(block_num_pixels);
                __m256 second_term = _mm256_set1_ps(-1*small_blockwise_sums_i*grey_blockwise_sums_j);
                __m256 temp1 =  _mm256_fmadd_ps(block_num_pixels_vec,graymulsmall_blockwise_sum_rotation_vec, second_term);
                __m256 contrast_rotation_vec = _mm256_mul_ps(temp1, denom_vec);

                __m256 small_blockwise_sums_vec = _mm256_set1_ps(small_blockwise_sums_i);
                __m256 first_term = _mm256_set1_ps(-1*grey_blockwise_sums_j*small_blockwise_sumofsquares_i);
                __m256 temp2 =  _mm256_fmadd_ps(small_blockwise_sums_vec, graymulsmall_blockwise_sum_rotation_vec, first_term) ;
                __m256 brightness_rotation_vec = _mm256_mul_ps(temp2, denom_vec);


                __m256 error_rotation_vec = _mm256_set1_ps(0.f);

                for(int p = 0; p < block_num_pixels; p++)
                {   
                    __m256 img_small_blocks_rotations_vec = _mm256_load_ps(&img_small_blocks_rotations_reordered[img_small_blocks_offset+(p * num_rotations)]);
        
                    __m256 x1 = _mm256_fmadd_ps(img_small_blocks_rotations_vec, contrast_rotation_vec, brightness_rotation_vec);
                    __m256 x2 = _mm256_set1_ps(img_gray_blocks[img_gray_blocks_offset + p]);
                    __m256 difference_vec = _mm256_sub_ps(x1, x2);

                    error_rotation_vec = _mm256_fmadd_ps(difference_vec, difference_vec, error_rotation_vec);
                
                }
                //Can we avxerize this?
                float *error  = static_cast<float*>(_mm_malloc(8*sizeof(float), 32));
                float *contrast = static_cast<float*>(_mm_malloc(8*sizeof(float), 32));
                float *brightness = static_cast<float*>(_mm_malloc(8*sizeof(float), 32));
                
                _mm256_store_ps(error, error_rotation_vec);
                _mm256_store_ps(contrast, contrast_rotation_vec);
                _mm256_store_ps(brightness, brightness_rotation_vec);

                // __m256 mappings_vec = _mm256_set1_ps(mappings[j*MAPSTORE + 0]);
                // __m256 comparison_results = _mm256_cmp_ps(mappings_vec,error_rotation_vec, _CMP_GT_OQ )
                // __m256 set_vec = comparison_results

                if(mappings[j*MAPSTORE + 0] > error[0]) {
                    mappings[j*MAPSTORE + 0] = error[0];
                    mappings[j*MAPSTORE + 1] = contrast[0];
                    mappings[j*MAPSTORE + 2] = brightness[0];
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[1]) {
                    mappings[j*MAPSTORE + 0] = error[1];
                    mappings[j*MAPSTORE + 1] = contrast[1];
                    mappings[j*MAPSTORE + 2] = brightness[1];
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[2]) {
                    mappings[j*MAPSTORE + 0] = error[2];
                    mappings[j*MAPSTORE + 1] = contrast[2];
                    mappings[j*MAPSTORE + 2] = brightness[2];
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[3]) {
                    mappings[j*MAPSTORE + 0] = error[3];
                    mappings[j*MAPSTORE + 1] = contrast[3];
                    mappings[j*MAPSTORE + 2] = brightness[3];
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[4]) {
                    mappings[j*MAPSTORE + 0] = error[4];
                    mappings[j*MAPSTORE + 1] = contrast[4];
                    mappings[j*MAPSTORE + 2] = brightness[4];
                    mappings[j*MAPSTORE + 3] = 4;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[5]) {
                    mappings[j*MAPSTORE + 0] = error[5];
                    mappings[j*MAPSTORE + 1] = contrast[5];
                    mappings[j*MAPSTORE + 2] = brightness[5];
                    mappings[j*MAPSTORE + 3] = 5;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[6]) {
                    mappings[j*MAPSTORE + 0] = error[6];
                    mappings[j*MAPSTORE + 1] = contrast[6];
                    mappings[j*MAPSTORE + 2] = brightness[6];
                    mappings[j*MAPSTORE + 3] = 6;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error[7]) {
                    mappings[j*MAPSTORE + 0] = error[7];
                    mappings[j*MAPSTORE + 1] = contrast[7];
                    mappings[j*MAPSTORE + 2] = brightness[7];
                    mappings[j*MAPSTORE + 3] = 7;
                    mappings[j*MAPSTORE + 4] = i;
                }

                
                
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    void find_optimal_mappings_avx_with_entropy(const int block_num_rows,
                               const int block_num_cols,
                               const int img_gray_num_blocks,
                               const int img_small_num_blocks,
                               float *img_gray_blocks,
                               float **img_small_blocks_rotations,
                               float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        // blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        // blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        // blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        for (int i = 0; i < img_gray_num_blocks; i++) {
            float s = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                s += img_gray_blocks[i * block_num_pixels + j];
            }
            gray_blockwise_sums[i] = s;
        }

        for (int i = 0; i < img_small_num_blocks; i++) {
            float s = 0.0;
            float s_sq = 0.0;
            float val = 0.0;
            for (int j = 0; j < block_num_pixels; j++) {
                val = img_small_blocks_rotations[0][i * block_num_pixels + j];
                s += val;
                s_sq += val*val;
            }
            small_blockwise_sums[i] = s;
            small_blockwise_sumofsquares[i] = s_sq;
        }

        std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];

            for(int i = 0; i < img_small_num_blocks; i+=1) {
                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                        //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                       continue;

                }

                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - std::pow(small_blockwise_sums_i, 2);

                __m256 graymulsmall_blockwise_sum_rotation0_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation1_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation2_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                __m256 graymulsmall_blockwise_sum_rotation3_row = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                //__m256 graymulsmall_blockwise_sum_rotation0_row_second = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                //__m256 graymulsmall_blockwise_sum_rotation1_row_second = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                //__m256 graymulsmall_blockwise_sum_rotation2_row_second = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);
                //__m256 graymulsmall_blockwise_sum_rotation3_row_second = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,0.f);

                for(int row = 0; row < 8; row++)
                {
                    int row_offset = row*8;

                    __m256 img_gray_blocks_row = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);


                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    graymulsmall_blockwise_sum_rotation0_row =  _mm256_fmadd_ps(img_small_blocks_rotation0_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation0_row);
                    graymulsmall_blockwise_sum_rotation1_row =  _mm256_fmadd_ps(img_small_blocks_rotation1_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation1_row);
                    graymulsmall_blockwise_sum_rotation2_row =  _mm256_fmadd_ps(img_small_blocks_rotation2_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation2_row);
                    graymulsmall_blockwise_sum_rotation3_row =  _mm256_fmadd_ps(img_small_blocks_rotation3_row, img_gray_blocks_row, graymulsmall_blockwise_sum_rotation3_row);

                }

                // //ugly horizontal sum now :)
                float graymulsmall_blockwise_sum_rotation0 = horizontal_sum(graymulsmall_blockwise_sum_rotation0_row);
                float graymulsmall_blockwise_sum_rotation1 = horizontal_sum(graymulsmall_blockwise_sum_rotation1_row);
                float graymulsmall_blockwise_sum_rotation2 = horizontal_sum(graymulsmall_blockwise_sum_rotation2_row);
                float graymulsmall_blockwise_sum_rotation3 = horizontal_sum(graymulsmall_blockwise_sum_rotation3_row);

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;


                __m256 contrast_rotation0_vec = _mm256_set1_ps(contrast_rotation0);
                __m256 contrast_rotation1_vec = _mm256_set1_ps(contrast_rotation1);
                __m256 contrast_rotation2_vec = _mm256_set1_ps(contrast_rotation2);
                __m256 contrast_rotation3_vec = _mm256_set1_ps(contrast_rotation3);

                __m256 brightness_rotation0_vec = _mm256_set1_ps(brightness_rotation0);
                __m256 brightness_rotation1_vec = _mm256_set1_ps(brightness_rotation1);
                __m256 brightness_rotation2_vec = _mm256_set1_ps(brightness_rotation2);
                __m256 brightness_rotation3_vec = _mm256_set1_ps(brightness_rotation3);

                __m256 error_rotation0_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation1_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation2_vec = _mm256_set1_ps(0.f);
                __m256 error_rotation3_vec = _mm256_set1_ps(0.f);


                for(int row = 0; row < 8; row++)
                {
                    int row_offset = row*8;

                    __m256 img_small_blocks_rotation0_row = _mm256_load_ps(&img_small_blocks_rotations[0][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation1_row = _mm256_load_ps(&img_small_blocks_rotations[1][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation2_row = _mm256_load_ps(&img_small_blocks_rotations[2][img_small_blocks_offset+row_offset]);
                    __m256 img_small_blocks_rotation3_row = _mm256_load_ps(&img_small_blocks_rotations[3][img_small_blocks_offset+row_offset]);

                    __m256 x1_0 = _mm256_fmadd_ps(img_small_blocks_rotation0_row, contrast_rotation0_vec, brightness_rotation0_vec);
                    __m256 x1_1 = _mm256_fmadd_ps(img_small_blocks_rotation1_row, contrast_rotation1_vec, brightness_rotation1_vec);
                    __m256 x1_2 = _mm256_fmadd_ps(img_small_blocks_rotation2_row, contrast_rotation2_vec, brightness_rotation2_vec);
                    __m256 x1_3 = _mm256_fmadd_ps(img_small_blocks_rotation3_row, contrast_rotation3_vec, brightness_rotation3_vec);

                    __m256 x2 = _mm256_load_ps(&img_gray_blocks[img_gray_blocks_offset + row_offset]);

                    __m256 difference_0_vec = _mm256_sub_ps(x1_0, x2);
                    __m256 difference_1_vec = _mm256_sub_ps(x1_1, x2);
                    __m256 difference_2_vec = _mm256_sub_ps(x1_2, x2);
                    __m256 difference_3_vec = _mm256_sub_ps(x1_3, x2);

                    error_rotation0_vec = _mm256_fmadd_ps(difference_0_vec, difference_0_vec, error_rotation0_vec);
                    error_rotation1_vec = _mm256_fmadd_ps(difference_1_vec, difference_1_vec, error_rotation1_vec);
                    error_rotation2_vec = _mm256_fmadd_ps(difference_2_vec, difference_2_vec, error_rotation2_vec);
                    error_rotation3_vec = _mm256_fmadd_ps(difference_3_vec, difference_3_vec, error_rotation3_vec);

                }

                //ugly horizontal sum now again:)
                float error_rotation0 = horizontal_sum(error_rotation0_vec);
                float error_rotation1 = horizontal_sum(error_rotation1_vec);
                float error_rotation2 = horizontal_sum(error_rotation2_vec);
                float error_rotation3 = horizontal_sum(error_rotation3_vec);

                if(mappings[j*MAPSTORE + 0] > error_rotation0) {
                    mappings[j*MAPSTORE + 0] = error_rotation0;
                    mappings[j*MAPSTORE + 1] = contrast_rotation0;
                    mappings[j*MAPSTORE + 2] = brightness_rotation0;
                    mappings[j*MAPSTORE + 3] = 0;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation1) {
                    mappings[j*MAPSTORE + 0] = error_rotation1;
                    mappings[j*MAPSTORE + 1] = contrast_rotation1;
                    mappings[j*MAPSTORE + 2] = brightness_rotation1;
                    mappings[j*MAPSTORE + 3] = 1;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation2) {
                    mappings[j*MAPSTORE + 0] = error_rotation2;
                    mappings[j*MAPSTORE + 1] = contrast_rotation2;
                    mappings[j*MAPSTORE + 2] = brightness_rotation2;
                    mappings[j*MAPSTORE + 3] = 2;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation3) {
                    mappings[j*MAPSTORE + 0] = error_rotation3;
                    mappings[j*MAPSTORE + 1] = contrast_rotation3;
                    mappings[j*MAPSTORE + 2] = brightness_rotation3;
                    mappings[j*MAPSTORE + 3] = 3;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;

    }


#endif

    void find_optimal_mappings_scalar_replacement_modified_entropy(const int block_num_rows,
                                                                   const int block_num_cols,
                                                                   const int img_gray_num_blocks,
                                                                   const int img_small_num_blocks,
                                                                   float *img_gray_blocks,
                                                                   float **img_small_blocks_rotations,
                                                                   float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *vec1 = new float[block_num_pixels];
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);
        std::cout << "Threshold " << threshold << std::endl;

        //std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                    //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                    continue;

                }
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    for (int p = 0; p < block_num_pixels; p++) {
                        vec1[p] = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                    }
                    float error;
                    fast::difference_norm_scalar_replacement_modified(block_num_cols, img_gray_blocks + img_gray_blocks_offset, vec1, &error);
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        //std::cout << "\t\tDone." << std::endl;

        delete[] vec1;
        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    static inline void find_optimal_mappings_inline_modified_with_entropy(const int block_num_rows,
                                                            const int block_num_cols,
                                                            const int img_gray_num_blocks,
                                                            const int img_small_num_blocks,
                                                            float *img_gray_blocks,
                                                            float **img_small_blocks_rotations,
                                                            float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);

        //inline difference_norm
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                    //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                    continue;

                }
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                for(int rot = 0; rot < 4; rot++) {
                    float graymulsmall_blockwise_sum = 0.0;
                    for(int p = 0; p < block_num_pixels; p++) {
                        graymulsmall_blockwise_sum += img_small_blocks_rotations[rot][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    }
                    //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                    const float contrast = (block_num_pixels * graymulsmall_blockwise_sum - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                    const float brightness = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum) / denom;
                    float error = 0.f;
                    for(int p = 0; p < block_num_pixels; p++){
                        const float x1 = contrast * img_small_blocks_rotations[rot][img_small_blocks_offset + p] + brightness;
                        const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                        const float difference = x1 - x2;
                        const float difference_squared = difference * difference;
                        error += difference_squared;
                    }
                    if(mappings[j*5 + 0] > error) {
                        mappings[j*5 + 0] = error;
                        mappings[j*5 + 1] = contrast;
                        mappings[j*5 + 2] = brightness;
                        mappings[j*5 + 3] = rot;
                        mappings[j*5 + 4] = i;
                    }
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_flip_ilp_modified_with_entropy(const int block_num_rows,
                                                 const int block_num_cols,
                                                 const int img_gray_num_blocks,
                                                 const int img_small_num_blocks,
                                                 float *img_gray_blocks,
                                                 float **img_small_blocks_rotations,
                                                 float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                    //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                    continue;

                }
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations[4][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations[5][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations[6][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations[7][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations[4][img_small_blocks_offset + p] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations[5][img_small_blocks_offset + p] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations[6][img_small_blocks_offset + p] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations[7][img_small_blocks_offset + p] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation4) {
                    mappings[j*5 + 0] = error_rotation4;
                    mappings[j*5 + 1] = contrast_rotation4;
                    mappings[j*5 + 2] = brightness_rotation4;
                    mappings[j*5 + 3] = 4;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation5) {
                    mappings[j*5 + 0] = error_rotation5;
                    mappings[j*5 + 1] = contrast_rotation5;
                    mappings[j*5 + 2] = brightness_rotation5;
                    mappings[j*5 + 3] = 5;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation6) {
                    mappings[j*5 + 0] = error_rotation6;
                    mappings[j*5 + 1] = contrast_rotation6;
                    mappings[j*5 + 2] = brightness_rotation6;
                    mappings[j*5 + 3] = 6;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation7) {
                    mappings[j*5 + 0] = error_rotation7;
                    mappings[j*5 + 1] = contrast_rotation7;
                    mappings[j*5 + 2] = brightness_rotation7;
                    mappings[j*5 + 3] = 7;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }


    void find_optimal_mappings_flip_ilp_modified_splitloop_with_entropy(const int block_num_rows,
                                                           const int block_num_cols,
                                                           const int img_gray_num_blocks,
                                                           const int img_small_num_blocks,
                                                           float *img_gray_blocks,
                                                           float **img_small_blocks_rotations,
                                                           float *mappings) {
        const int block_num_pixels = block_num_rows * block_num_cols;

        // Calculate entropies
        float *small_entropy = new float[img_small_num_blocks];
        float threshold = blockwise_entropy(block_num_pixels, img_small_num_blocks, img_small_blocks_rotations[0], small_entropy, DISCARD_PERCENTAGE);

        //reconsturct loop to improve ILP, locality
        float *gray_blockwise_sums  = new float[img_gray_num_blocks]; //y
        float *small_blockwise_sums = new float[img_small_num_blocks]; //x
        float *small_blockwise_sumofsquares = new float[img_small_num_blocks];
        blockwise_sum(img_gray_num_blocks,  block_num_pixels, img_gray_blocks,               gray_blockwise_sums);
        blockwise_sum(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sums);
        blockwise_sum_of_squares(img_small_num_blocks, block_num_pixels, img_small_blocks_rotations[0], small_blockwise_sumofsquares);

        std::cout << "\t\tStart: nested loops that solve least squares problems and blockwise_sum_of_xmuly" << std::endl;
        int img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                // Don't compare rotation errors if the difference in entropies is higher than threshold
                //if(small_entropy[i]  > threshold + 3 * standard_deviation){
                if(small_entropy[i]  > threshold ){
                    //std::cout << "The entropy is too different " << small_entropy[i] - grey_entropy << std::endl;
                    continue;

                }
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation0 = 0.f;
                float graymulsmall_blockwise_sum_rotation1 = 0.f;
                float graymulsmall_blockwise_sum_rotation2 = 0.f;
                float graymulsmall_blockwise_sum_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation0 += img_small_blocks_rotations[0][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation1 += img_small_blocks_rotations[1][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation2 += img_small_blocks_rotations[2][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation3 += img_small_blocks_rotations[3][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }
                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation0 = (block_num_pixels * graymulsmall_blockwise_sum_rotation0 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation1 = (block_num_pixels * graymulsmall_blockwise_sum_rotation1 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation2 = (block_num_pixels * graymulsmall_blockwise_sum_rotation2 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation3 = (block_num_pixels * graymulsmall_blockwise_sum_rotation3 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation0 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation0) / denom;
                const float brightness_rotation1 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation1) / denom;
                const float brightness_rotation2 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation2) / denom;
                const float brightness_rotation3 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i - small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation3) / denom;
                float error_rotation0 = 0.f;
                float error_rotation1 = 0.f;
                float error_rotation2 = 0.f;
                float error_rotation3 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_0 = contrast_rotation0 * img_small_blocks_rotations[0][img_small_blocks_offset + p] + brightness_rotation0;
                    const float x1_1 = contrast_rotation1 * img_small_blocks_rotations[1][img_small_blocks_offset + p] + brightness_rotation1;
                    const float x1_2 = contrast_rotation2 * img_small_blocks_rotations[2][img_small_blocks_offset + p] + brightness_rotation2;
                    const float x1_3 = contrast_rotation3 * img_small_blocks_rotations[3][img_small_blocks_offset + p] + brightness_rotation3;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_0 = x1_0 - x2;
                    const float difference_1 = x1_1 - x2;
                    const float difference_2 = x1_2 - x2;
                    const float difference_3 = x1_3 - x2;
                    const float difference_0_squared = difference_0 * difference_0;
                    const float difference_1_squared = difference_1 * difference_1;
                    const float difference_2_squared = difference_2 * difference_2;
                    const float difference_3_squared = difference_3 * difference_3;
                    error_rotation0 += difference_0_squared;
                    error_rotation1 += difference_1_squared;
                    error_rotation2 += difference_2_squared;
                    error_rotation3 += difference_3_squared;
                }
                if(mappings[j*5 + 0] > error_rotation0) {
                    mappings[j*5 + 0] = error_rotation0;
                    mappings[j*5 + 1] = contrast_rotation0;
                    mappings[j*5 + 2] = brightness_rotation0;
                    mappings[j*5 + 3] = 0;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation1) {
                    mappings[j*5 + 0] = error_rotation1;
                    mappings[j*5 + 1] = contrast_rotation1;
                    mappings[j*5 + 2] = brightness_rotation1;
                    mappings[j*5 + 3] = 1;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation2) {
                    mappings[j*5 + 0] = error_rotation2;
                    mappings[j*5 + 1] = contrast_rotation2;
                    mappings[j*5 + 2] = brightness_rotation2;
                    mappings[j*5 + 3] = 2;
                    mappings[j*5 + 4] = i;
                }
                if(mappings[j*5 + 0] > error_rotation3) {
                    mappings[j*5 + 0] = error_rotation3;
                    mappings[j*5 + 1] = contrast_rotation3;
                    mappings[j*5 + 2] = brightness_rotation3;
                    mappings[j*5 + 3] = 3;
                    mappings[j*5 + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }

        img_gray_blocks_offset = 0;
        for(int j = 0; j < img_gray_num_blocks; j++) {
            int img_small_blocks_offset = 0;
            const float grey_blockwise_sums_j = gray_blockwise_sums[j];
            for(int i = 0; i < img_small_num_blocks; i++) {
                const float small_blockwise_sumofsquares_i = small_blockwise_sumofsquares[i];
                const float small_blockwise_sums_i = small_blockwise_sums[i];
                const float denom = block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i;
                //const float INVdenom = 1.f / (block_num_pixels * small_blockwise_sumofsquares_i - small_blockwise_sums_i * small_blockwise_sums_i);
                float graymulsmall_blockwise_sum_rotation4 = 0.f;
                float graymulsmall_blockwise_sum_rotation5 = 0.f;
                float graymulsmall_blockwise_sum_rotation6 = 0.f;
                float graymulsmall_blockwise_sum_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++) {
                    graymulsmall_blockwise_sum_rotation4 += img_small_blocks_rotations[4][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation5 += img_small_blocks_rotations[5][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation6 += img_small_blocks_rotations[6][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                    graymulsmall_blockwise_sum_rotation7 += img_small_blocks_rotations[7][img_small_blocks_offset + p] * img_gray_blocks[img_gray_blocks_offset + p];
                }

                //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
                const float contrast_rotation4 = (block_num_pixels * graymulsmall_blockwise_sum_rotation4 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation5 = (block_num_pixels * graymulsmall_blockwise_sum_rotation5 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation6 = (block_num_pixels * graymulsmall_blockwise_sum_rotation6 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float contrast_rotation7 = (block_num_pixels * graymulsmall_blockwise_sum_rotation7 - small_blockwise_sums_i * grey_blockwise_sums_j) / denom;
                const float brightness_rotation4 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation4) / denom;
                const float brightness_rotation5 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation5) / denom;
                const float brightness_rotation6 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation6) / denom;
                const float brightness_rotation7 = (grey_blockwise_sums_j * small_blockwise_sumofsquares_i  -  small_blockwise_sums_i * graymulsmall_blockwise_sum_rotation7) / denom;
                float error_rotation4 = 0.f;
                float error_rotation5 = 0.f;
                float error_rotation6 = 0.f;
                float error_rotation7 = 0.f;
                for(int p = 0; p < block_num_pixels; p++){
                    const float x1_4 = contrast_rotation4 * img_small_blocks_rotations[4][img_small_blocks_offset + p] + brightness_rotation4;
                    const float x1_5 = contrast_rotation5 * img_small_blocks_rotations[5][img_small_blocks_offset + p] + brightness_rotation5;
                    const float x1_6 = contrast_rotation6 * img_small_blocks_rotations[6][img_small_blocks_offset + p] + brightness_rotation6;
                    const float x1_7 = contrast_rotation7 * img_small_blocks_rotations[7][img_small_blocks_offset + p] + brightness_rotation7;
                    const float x2 = img_gray_blocks[img_gray_blocks_offset + p];
                    const float difference_4 = x1_4 - x2;
                    const float difference_5 = x1_5 - x2;
                    const float difference_6 = x1_6 - x2;
                    const float difference_7 = x1_7 - x2;
                    const float difference_4_squared = difference_4 * difference_4;
                    const float difference_5_squared = difference_5 * difference_5;
                    const float difference_6_squared = difference_6 * difference_6;
                    const float difference_7_squared = difference_7 * difference_7;
                    error_rotation4 += difference_4_squared;
                    error_rotation5 += difference_5_squared;
                    error_rotation6 += difference_6_squared;
                    error_rotation7 += difference_7_squared;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation4) {
                    mappings[j*MAPSTORE + 0] = error_rotation4;
                    mappings[j*MAPSTORE + 1] = contrast_rotation4;
                    mappings[j*MAPSTORE + 2] = brightness_rotation4;
                    mappings[j*MAPSTORE + 3] = 4;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation5) {
                    mappings[j*MAPSTORE + 0] = error_rotation5;
                    mappings[j*MAPSTORE + 1] = contrast_rotation5;
                    mappings[j*MAPSTORE + 2] = brightness_rotation5;
                    mappings[j*MAPSTORE + 3] = 5;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation6) {
                    mappings[j*MAPSTORE + 0] = error_rotation6;
                    mappings[j*MAPSTORE + 1] = contrast_rotation6;
                    mappings[j*MAPSTORE + 2] = brightness_rotation6;
                    mappings[j*MAPSTORE + 3] = 6;
                    mappings[j*MAPSTORE + 4] = i;
                }
                if(mappings[j*MAPSTORE + 0] > error_rotation7) {
                    mappings[j*MAPSTORE + 0] = error_rotation7;
                    mappings[j*MAPSTORE + 1] = contrast_rotation7;
                    mappings[j*MAPSTORE + 2] = brightness_rotation7;
                    mappings[j*MAPSTORE + 3] = 7;
                    mappings[j*MAPSTORE + 4] = i;
                }
                img_small_blocks_offset += block_num_pixels;
            }
            img_gray_blocks_offset += block_num_pixels;
        }
        std::cout << "\t\tDone." << std::endl;

        delete[] gray_blockwise_sums;
        delete[] small_blockwise_sums;
        delete[] small_blockwise_sumofsquares;
    }

    
    void find_optimal_mappings(const int block_num_rows,
            const int block_num_cols,
            const int img_gray_num_blocks,
            const int img_small_num_blocks,
            float *img_gray_blocks,
            float **img_small_blocks_rotations,
            float *mappings) {
        
        //here we can select the function (e.g. the fastest find_optimal_mappings function) we want to use in src/compress_fast.cpp
        //reference::find_optimal_mappings(
#ifdef _SCALAR_REPLACEMENT_
        //find_optimal_mappings_scalar_replacement(
        find_optimal_mappings_scalar_replacement_modified(
#elif defined _INLINE_
        find_optimal_mappings_inline_modified(
#elif defined _ILP_
        //find_optimal_mappings_ilp_modified_reordered(
        find_optimal_mappings_ilp_modified(
#elif defined _MANUAL_AVX_
        find_optimal_mappings_avx(
#else
        //here select function that you want to use in executable: fic_compress_fast
        //find_optimal_mappings_vectorized_modified(
        //find_optimal_mappings_avx(
        //find_optimal_mappings_avx_with_entropy(
        //reference::find_optimal_mappings(
        find_optimal_mappings_ilp_modified(
#endif
                block_num_rows,
                block_num_cols,
                img_gray_num_blocks,
                img_small_num_blocks,
                img_gray_blocks,
                img_small_blocks_rotations,
                mappings);
    }
    
    void find_optimal_mappings_flip(const int block_num_rows,
            const int block_num_cols,
            const int img_gray_num_blocks,
            const int img_small_num_blocks,
            float *img_gray_blocks,
            float **img_small_blocks_transformations,
            float *mappings) {
        
        //-----------------------------------only for fic compression with rotation AND flip-------------------------------------
        
        //here we can select the function (e.g. the fastest find_optimal_mappings function) we want to use in src/compress_fast.cpp
//#ifdef _SCALAR_REPLACEMENT_
//        find_optimal_mappings_flip_scalar_replacement_modified(
//#elif defined _INLINE_
//        find_optimal_mappings_flip_inline_modified(
//#elif defined _ILP_
//        find_optimal_mappings_flip_ilp_modified(
//#elif defined _AVX_
//        find_optimal_mappings_flip_avx_splitloop(
//#else
//        //here select function that you want to use in executable: fic_compress_flip_fast
//        //find_optimal_mappings_flip_avx_splitloop(
//        //reference::find_optimal_mappings_with_entropy(
//        //reference::find_optimal_mappings_flip(
//        //find_optimal_mappings_flip_ilp_modified_with_entropy(
//        find_optimal_mappings_flip_avx_splitloop_with_entropy(
//#endif
#ifdef _FAST_REFERENCE_
        reference::find_optimal_mappings_flip(
#elif defined _MODIFIED_
        find_optimal_mappings_flip_modified(
#elif defined _SCALAR_REPLACEMENT_
        find_optimal_mappings_flip_scalar_replacement_modified(
#elif defined _INLINE_
        find_optimal_mappings_flip_inline_modified(
#elif defined _ILP_
        //find_optimal_mappings_flip_ilp_modified(
        find_optimal_mappings_flip_ilp_modified_splitloop(
#elif defined _MANUAL_AVX_
        //find_optimal_mappings_flip_avx(
        find_optimal_mappings_flip_avx_splitloop(
#elif defined _ENTROPY_
        //reference::find_optimal_mappings_with_entropy(
        //find_optimal_mappings_flip_ilp_modified_with_entropy(
        //find_optimal_mappings_flip_ilp_modified_splitloop_with_entropy(
        find_optimal_mappings_flip_avx_splitloop_with_entropy(
#else
        //here select function that you want to use in executable: fic_compress_flip_fast
        //find_optimal_mappings_flip_avx_splitloop(
        
        //fast unoptimized
        //reference::find_optimal_mappings_flip(
        
        //fast optimized
        find_optimal_mappings_flip_avx_splitloop(
#endif
                block_num_rows,
                block_num_cols,
                img_gray_num_blocks,
                img_small_num_blocks,
                img_gray_blocks,
                img_small_blocks_transformations,
                mappings);
    }



}
