#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <fic/common.h>
#include <fic/resize.h>
#include <fic/rotate.h>
#include <fic/split_into_blocks.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include <string>
#include <fstream>

constexpr int num_iterations = 20;
//constexpr int blockDimension = 8;
constexpr int block_num_rows = 8;
constexpr int block_num_cols = 8;
constexpr int block_num_pixels = block_num_rows * block_num_cols;
//constexpr int mapStore = 5;

int main(int argc, char *argv[]){
    if ( argc != 3 ){
        std::cout<<"usage: "<< argv[0] <<" <compression_file> <image_output_name_extension>\n";
        return -1;
    }

    std::cout << "Start: decompress image" << std::endl;
    //read the compression file
    std::string compression_file = argv[1];
    std::string output_file_name = argv[2];
    std::ifstream infile1(compression_file.c_str());
    int stride;
    int blockDimension;
    int img_gray_num_blocks_horizontal;
    int img_gray_num_blocks_vertical;
    infile1 >> stride >> blockDimension >> img_gray_num_blocks_horizontal >> img_gray_num_blocks_vertical;
    std::cout << "\tstride = " << stride << std::endl;
    const int stride_small = stride / 2;
    const int img_gray_num_blocks = img_gray_num_blocks_vertical * img_gray_num_blocks_horizontal;
    std::cout << "\timg_gray_num_blocks = " << img_gray_num_blocks << std::endl;
    float *mappings = (float *)malloc(img_gray_num_blocks*4*sizeof(float));
    for(int i = 0; i < img_gray_num_blocks; i++)
        for(int j = 0; j < 4; j++)
            infile1 >> mappings[i * 4 + j];
    const int height = img_gray_num_blocks_vertical * blockDimension;
    const int width = img_gray_num_blocks_horizontal * blockDimension;
    std::cout << "\theight = " << height << "\twidth = " << width << std::endl;
    const int height_small = height / 2;
    const int width_small = width / 2;
    //const int img_small_num_blocks_vertical = height_small / stride_small;
    //const int img_small_num_blocks_horizontal = width_small / stride_small;
    const int img_small_num_blocks_vertical =   stride_small < block_num_cols ? height_small / stride_small - (block_num_cols / stride_small) + 1 : height_small / stride_small;
    const int img_small_num_blocks_horizontal = stride_small < block_num_rows ? width_small / stride_small - (block_num_rows / stride_small) + 1  : width_small / stride_small;
    const int img_small_num_blocks = img_small_num_blocks_vertical * img_small_num_blocks_horizontal;
    
    float *img_gray               = (float *)calloc(height*width, sizeof(float));
    float *img_small              = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot90        = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot180       = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot270       = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_hflip        = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_hflip_rot90  = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_hflip_rot180 = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_hflip_rot270 = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_gray_blocks             = new float[img_gray_num_blocks * block_num_pixels];
    float **img_small_blocks_transformations = new float *[8];
    img_small_blocks_transformations[0]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[1]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[2]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[3]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[4]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[5]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[6]      = new float[img_small_num_blocks * block_num_pixels];
    img_small_blocks_transformations[7]      = new float[img_small_num_blocks * block_num_pixels];
    
    std::cout << "\tStart: use optimal mappings to iteratively reconstruct the image" << std::endl;
    for (int i = 0; i < num_iterations; i++) {
        //------------uncomment in order to visualize the individual decompression iterations--------------------
        std::string output_file_name_current_iter = output_file_name + "_iter" + std::to_string(i);
        std::cout << "Start: save image as fic/images/" << output_file_name_current_iter << ".png" << std::endl;
        //can be used to resolve relative paths
        filesystem::resolver resolver;
        resolver.append("");
        resolver.append("../");
        resolver.append("fic/");
        unsigned char *img_gray_int = (unsigned char *)malloc(height*width*sizeof(unsigned char));
        for (int p = 0; p < height*width; p++) {
            if (img_gray[p] < 0.f)
                img_gray[p] = 0.f;
            else if (img_gray[p] > 255.f)
                img_gray[p] = 255.f;
            img_gray_int[p] = (uint8_t)img_gray[p];
        }
        stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path(output_file_name_current_iter + ".png")).str().c_str(), width, height, 1, img_gray_int, width);
        delete[] img_gray_int;
        //-----------------------------------------------------------------------------------------------------
        
        resize(height, width, img_gray, img_small);
        //rotate(height_small, width_small, img_small, img_small_rot90, img_small_rot180, img_small_rot270);
        rotate_and_hflip(height_small, width_small, img_small,
                img_small_rot90,
                img_small_rot180,
                img_small_rot270,
                img_small_hflip,
                img_small_hflip_rot90,
                img_small_hflip_rot180,
                img_small_hflip_rot270);
        split_into_blocks(                           block_num_rows, block_num_cols, height,       width,       img_gray,               img_gray_blocks);
        split_into_blocks(             stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small,              img_small_blocks_transformations[0]);
        //split_into_blocks_rot90(       stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot90,        img_small_blocks_transformations[1]);
        split_into_blocks_rot90(       stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_rot90,        img_small_blocks_transformations[1]);
        split_into_blocks_rot180(      stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot180,       img_small_blocks_transformations[2]);
        //split_into_blocks_rot270(      stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot270,       img_small_blocks_transformations[3]);
        split_into_blocks_rot270(      stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_rot270,       img_small_blocks_transformations[3]);
        split_into_blocks_hflip(       stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_hflip,        img_small_blocks_transformations[4]);
        //split_into_blocks_hflip_rot90( stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_hflip_rot90,  img_small_blocks_transformations[5]);
        split_into_blocks_hflip_rot90( stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_hflip_rot90,  img_small_blocks_transformations[5]);
        split_into_blocks_hflip_rot180(stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_hflip_rot180, img_small_blocks_transformations[6]);
        //split_into_blocks_hflip_rot270(stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_hflip_rot270, img_small_blocks_transformations[7]);
        split_into_blocks_hflip_rot270(stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_hflip_rot270, img_small_blocks_transformations[7]);

        int img_gray_blocks_offset = 0;
        for (int b = 0; b < img_gray_num_blocks; b++) {
            const float contrast =       mappings[b * 4 + 0];
            const float brightness =     mappings[b * 4 + 1];
            const int   transformation = mappings[b * 4 + 2];
            const int   block_idx =      mappings[b * 4 + 3];
            for (int p = 0; p < block_num_pixels; p++)
                img_gray_blocks[img_gray_blocks_offset + p] = contrast * img_small_blocks_transformations[transformation][block_idx * block_num_pixels + p] + brightness;
            img_gray_blocks_offset += block_num_pixels;
        }
        
        // undo split
        for(int block_idx_row = 0; block_idx_row < img_gray_num_blocks_vertical; block_idx_row++) {
            for(int block_idx_col = 0; block_idx_col < img_gray_num_blocks_horizontal; block_idx_col++) {
                const int block_idx = block_idx_row * img_gray_num_blocks_horizontal + block_idx_col;
                const int block_offset = block_idx * block_num_pixels;
                int img_i = block_idx_row * block_num_rows;
                for (int block_i = 0; block_i < block_num_rows; block_i++, img_i++) {
                    int img_j = block_idx_col * block_num_cols;
                    for (int block_j = 0; block_j < block_num_cols; block_j++, img_j++) {
                        img_gray[img_i * width + img_j] = img_gray_blocks[block_offset + block_i * block_num_rows + block_j];
                        if (img_gray[img_i * width + img_j] < 0.f)
                            img_gray[img_i * width + img_j] = 0.f;
                        else if (img_gray[img_i * width + img_j] > 255.f)
                            img_gray[img_i * width + img_j] = 255.f;
                    }
                }
            }
        }

    }
    std::cout << "\tDone." << std::endl;
    
    delete[] mappings;
    delete[] img_small;
    delete[] img_small_blocks_transformations[0];
    delete[] img_small_blocks_transformations[1];
    delete[] img_small_blocks_transformations[2];
    delete[] img_small_blocks_transformations[3];
    delete[] img_small_blocks_transformations[4];
    delete[] img_small_blocks_transformations[5];
    delete[] img_small_blocks_transformations[6];
    delete[] img_small_blocks_transformations[7];
    delete[] img_small_blocks_transformations;
    delete[] img_gray_blocks;
    std::cout << "Done." << std::endl;

    std::cout << "Start: save image as fic/images/" << output_file_name << ".png" << std::endl;
    //can be used to resolve relative paths
    filesystem::resolver resolver;
    resolver.append("");
    resolver.append("../");
    resolver.append("fic/");
    unsigned char *img_gray_int = (unsigned char *)malloc(height*width*sizeof(unsigned char));
    for (int p = 0; p < height*width; p++) {
        if (img_gray[p] < 0.f)
            img_gray[p] = 0.f;
        else if (img_gray[p] > 255.f)
            img_gray[p] = 255.f;
        img_gray_int[p] = (uint8_t)img_gray[p];
    }
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path(output_file_name + ".png")).str().c_str(), width, height, 1, img_gray_int, width);
    delete[] img_gray;
    delete[] img_gray_int;
    std::cout << "Done" << std::endl;
}
