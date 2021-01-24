#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <fic/common.h>
#include <fic/resize.h>
#include <fic/rotate.h>
#include <fic/split_into_blocks.h>
#include <fic/find_optimal_mappings.h>
#include <cstdlib>
#include <fstream>
#include <chrono>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb/stb_image_write.h"
#include "fic/tsc_x86.h"

//using namespace reference;
using namespace fast;

int main(int argc, char **argv) {
    //can be used to resolve relative paths
    filesystem::resolver resolver;
    resolver.append("");
    resolver.append("../");
    resolver.append("fic/");
    
    //load image
    int height, width, channels;
    unsigned char *img;
    int stride = 16;
    if (argc == 2) {
        // image file is passed as argument, handle it
        std::string input_filename = argv[1];
        filesystem::path input_path(input_filename);
        if (input_path.extension() == "png" ||
            input_path.extension() == "jpg" ||
            input_path.extension() == "jpeg" ||
            input_path.extension() == "gif" ||
            input_path.extension() == "bmp") {
            // read the image file
            img = stbi_load(input_path.make_absolute().str().c_str(), &width, &height, &channels, 0);
            if (img == NULL) {
                std::cerr << "Error: failed to load the image" << std::endl;
                return -1;
            }
        } else {
            std::cerr << "Error: unknown file \"" << input_filename
            << "\", expected an extension of type .png, .gif, .jpg, .jpeg or .bmp" << std::endl;
            return -1;
        }
    } else if (argc == 3) {
        stride = atoi(argv[2]);
        // image file is passed as argument, handle it
        std::string input_filename = argv[1];
        filesystem::path input_path(input_filename);
        if (input_path.extension() == "png" ||
            input_path.extension() == "jpg" ||
            input_path.extension() == "jpeg" ||
            input_path.extension() == "gif" ||
            input_path.extension() == "bmp") {
            // read the image file
            img = stbi_load(input_path.make_absolute().str().c_str(), &width, &height, &channels, 0);
            if (img == NULL) {
                std::cerr << "Error: failed to load the image" << std::endl;
                return -1;
            }
        } else {
            std::cerr << "Error: unknown file \"" << input_filename
            << "\", expected an extension of type .png, .gif, .jpg, .jpeg or .bmp" << std::endl;
            return -1;
        }
    } else {
        std::cout << "Loading lena_color.gif. If you want to use other images:\n" <<
            "usage: relative/path/to/fic <relative/path/to/custom_image>\n" <<
            "And note: lena_gray.gif looks grayscale but has 4 channels. Check for yourself by using relative/path/to/fic <relative/path/to/lena_gray.gif>\n" << std::endl;
        // read the GIF image file
        img = stbi_load((resolver.resolve("images").make_absolute()/filesystem::path("lena_color.gif")).str().c_str(), &width, &height, &channels, 0);
        if (img == NULL) {
            std::cerr << "Error: failed to load the GIF image" << std::endl;
            return -1;
        }
    }
    std::cout << "Image loaded successfully: height = " << height << "\twidth = " << width << "\tchannels = " << channels << "\n" << std::endl;
    
    //convert to float array, ensure grayscale
    float *img_gray = (float *)malloc(height*width*sizeof(float));
    if (channels >= 3) {
        //convert color image to grayscale
        std::cout << "Start: convert color image to grayscale" << std::endl;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                img_gray[i * width + j] = \
                    0.299 * (uint8_t)img[i * channels * width + j * channels + 0] +
                    0.587 * (uint8_t)img[i * channels * width + j * channels + 1] +
                    0.114 * (uint8_t)img[i * channels * width + j * channels + 2];
            }
        }
        std::cout << "Done.\n" << std::endl;
    } else {
        for (int idx = 0; idx < height*width; idx++) {
            img_gray[idx] = (float)img[idx];
        }
    }
    stbi_image_free(img);


    //calculate flops
    const int block_num_rows_c = 8;
    const int block_num_cols_c = 8;

    const int block_num_pixels_c = block_num_cols_c * block_num_rows_c; 
    const int img_gray_num_blocks_vertical_c = height / block_num_cols_c;
    const int img_gray_num_blocks_horizontal_c = width / block_num_rows_c;
    const int img_gray_num_blocks_c = img_gray_num_blocks_vertical_c * img_gray_num_blocks_horizontal_c;
    const int height_small_c = height / 2;
    const int width_small_c = width / 2;
    const int img_small_num_blocks_vertical_c = height_small_c / block_num_cols_c;
    const int img_small_num_blocks_horizontal_c = width_small_c / block_num_rows_c;
    const int img_small_num_blocks_c = img_small_num_blocks_vertical_c * img_small_num_blocks_horizontal_c;
    const int num_rot = 4; // because no flip

    // nested loop
    long num_flops = 3; // compute denom
    num_flops = num_flops + 8 * 8 * num_rot * 2; // compute graymulsmall...
    //num_flops = num_flops + 3 * 8 * num_rot; // compute horizontal_sum
#ifdef _AVX_
    num_flops = num_flops + 7 * num_rot; // compute horizontal_sum
#endif
    num_flops = num_flops + 3 * num_rot; // compute contrast
    num_flops = num_flops + 3 * num_rot; // compute brightness
    num_flops = num_flops + 8 * (8 * 2 * num_rot + 8 * num_rot + 8 * 2 * num_rot); // compute error
    //num_flops = num_flops + 3 * 8 * num_rot; // compute horizontal_sum again
    num_flops = num_flops + 7 * num_rot; // compute horizontal_sum again
    num_flops = num_flops * img_small_num_blocks_c; // 2nd layer loop
    num_flops = num_flops * img_gray_num_blocks_c; // 1st layer loop

    // inlined blockwise sum
    num_flops = num_flops + img_gray_num_blocks_c * block_num_pixels_c; // compute gray_blockwise_sums
    num_flops = num_flops + 3 * img_small_num_blocks_c * block_num_pixels_c; // compute small_blockwise_sums and small_blockwise_sumofsquares

    // stuff outside find_optimal_mappings
    num_flops = num_flops + height * width * 3 / 4; // resize

    // number of divisions
    long num_divs = num_rot; // contrast
    num_divs = num_divs + num_rot; // brightness
    num_divs = num_divs * img_small_num_blocks_c * img_gray_num_blocks_c; // inner 2 loops
    num_divs = num_divs + height * width / 4; // resize

    // number of square roots
    long num_sqrts = 0;

    // compute total number of flops
    long total_num_flops = num_flops + num_divs + num_sqrts;
    
    //compress
    myInt64 start, end;

    std::cout << "Start: compress image" << std::endl;

    auto t0 = std::chrono::steady_clock::now();
    start = start_tsc();

    constexpr int mapStore = 5;
    constexpr int block_num_cols = 8; 
    constexpr int block_num_rows = 8; 
    constexpr int block_num_pixels = block_num_cols * block_num_rows; 
    const int height_small = height / 2;
    const int width_small = width / 2;
    const int stride_small = stride / 2;
    const int img_gray_num_blocks_vertical = height / block_num_cols;
    const int img_gray_num_blocks_horizontal = width / block_num_rows;
    const int img_gray_num_blocks = img_gray_num_blocks_vertical * img_gray_num_blocks_horizontal;
    //const int img_small_num_blocks_vertical = height_small / stride_small;
    //const int img_small_num_blocks_horizontal = width_small / stride_small;
    const int img_small_num_blocks_vertical =   stride_small < block_num_cols ? height_small / stride_small - (block_num_cols / stride_small) + 1 : height_small / stride_small;
    const int img_small_num_blocks_horizontal = stride_small < block_num_rows ? width_small / stride_small - (block_num_rows / stride_small) + 1  : width_small / stride_small;
    const int img_small_num_blocks = img_small_num_blocks_vertical * img_small_num_blocks_horizontal;
    
    float *img_small        = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot90  = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot180 = (float *)malloc(height_small*width_small*sizeof(float));
    float *img_small_rot270 = (float *)malloc(height_small*width_small*sizeof(float));

    float *img_gray_blocks         =  static_cast<float*>(_mm_malloc(img_gray_num_blocks * block_num_pixels * sizeof(float), 32));
    // float *img_gray_blocks         = new float[img_gray_num_blocks * block_num_pixels];
    float **img_small_blocks_rotations = new float *[4];

    img_small_blocks_rotations[0]    = static_cast<float*>(_mm_malloc(img_small_num_blocks * block_num_pixels * sizeof(float), 32));
    img_small_blocks_rotations[1]    = static_cast<float*>(_mm_malloc(img_small_num_blocks * block_num_pixels * sizeof(float), 32));
    img_small_blocks_rotations[2]    = static_cast<float*>(_mm_malloc(img_small_num_blocks * block_num_pixels * sizeof(float), 32));
    img_small_blocks_rotations[3]    = static_cast<float*>(_mm_malloc(img_small_num_blocks * block_num_pixels * sizeof(float), 32));
    //img_small_blocks_rotations[0]    = new float[img_small_num_blocks * block_num_pixels];
    //img_small_blocks_rotations[1]    = new float[img_small_num_blocks * block_num_pixels];
    //img_small_blocks_rotations[2]    = new float[img_small_num_blocks * block_num_pixels];
    //img_small_blocks_rotations[3]    = new float[img_small_num_blocks * block_num_pixels];
    float *mappings = new float[img_gray_num_blocks * mapStore];
    
    for(int i = 0; i < img_gray_num_blocks * mapStore; i++){
        mappings[i] = std::numeric_limits<float>::max();
    }
    
    std::cout << "\tStart: resize" << std::endl;
    resize(height, width, img_gray, img_small);
    std::cout << "\tDone." << std::endl;
    std::cout << "\tStart: rotate" << std::endl;
    rotate(height_small, width_small, img_small, img_small_rot90, img_small_rot180, img_small_rot270);
    std::cout << "\tDone." << std::endl;
    std::cout << "\tStart: split image into blocks with stride = " << stride << " resulting in " << img_small_num_blocks << " blocks" << std::endl;
    split_into_blocks(       block_num_rows, block_num_cols, height,       width,       img_gray,         img_gray_blocks);
    split_into_blocks(       stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small,        img_small_blocks_rotations[0]);
    //split_into_blocks_rot90( stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot90,  img_small_blocks_rotations[1]);
    split_into_blocks_rot90( stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_rot90,  img_small_blocks_rotations[1]);
    split_into_blocks_rot180(stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot180, img_small_blocks_rotations[2]);
    //split_into_blocks_rot270(stride_small, block_num_rows, block_num_cols, height_small, width_small, img_small_rot270, img_small_blocks_rotations[3]);
    split_into_blocks_rot270(stride_small, block_num_rows, block_num_cols, width_small, height_small, img_small_rot270, img_small_blocks_rotations[3]);
    free(img_gray);
    free(img_small);
    free(img_small_rot90);
    free(img_small_rot180);
    free(img_small_rot270);
    std::cout << "\tDone." << std::endl;
    
    std::cout << "\tStart: find optimal mapping" << std::endl;
    find_optimal_mappings(block_num_rows,
            block_num_cols,
            img_gray_num_blocks,
            img_small_num_blocks,
            img_gray_blocks,
            img_small_blocks_rotations,
            mappings);       
    end = stop_tsc(start);
    auto t1 = std::chrono::steady_clock::now();
    double ns =  (double)std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    std::cout << "\tDone." << std::endl;
    std::cout << "Done.\t" << 1e-9 * ns << " s\n" << std::endl;


    //output flops/cycle
    double perf = round((100.0 * total_num_flops) / (double)end) / 100.0;
    std::cout << "Running: compress_fast " << perf << " flops/ cycles" << std::endl;

    std::cout << "Start: write fic_compress_fast_output" << std::endl;
    std::ofstream ofile((resolver.resolve("images").make_absolute()/filesystem::path("fic_compress_fast_output")).str());
    ofile << stride << " " << block_num_cols << " " << img_gray_num_blocks_horizontal << " " << img_gray_num_blocks_vertical << std::endl;
    for(int i = 0; i < img_gray_num_blocks; i++){
        for(int j = 1; j < mapStore; j++){
            ofile << mappings[i * mapStore + j] << " ";
        }
        ofile << std::endl;
    }
    ofile.close();
    std::cout << "Done.\n" << std::endl;

    //-------------------------------------DEBUGGING-----------------------------------------------------------------------------
    //let us find the absolute path of the projects image folder using the resolver class
    /*
    std::cout << "Start: write images for debug purposes" << std::endl;
    unsigned char *img_small_int               = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_rot90_int         = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_rot180_int        = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_rot270_int        = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_gray_blocks_int         = (unsigned char *)malloc(height*width*sizeof(unsigned char));
    unsigned char *img_small_blocks_int        = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_blocks_rot90_int  = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_blocks_rot180_int = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    unsigned char *img_small_blocks_rot270_int = (unsigned char *)malloc(height_small*width_small*sizeof(unsigned char));
    for (int idx = 0; idx < height_small*width_small; idx++) {
        img_small_int[idx]               = (uint8_t)img_small[idx];
        img_small_rot90_int[idx]         = (uint8_t)img_small_rot90[idx];
        img_small_rot180_int[idx]        = (uint8_t)img_small_rot180[idx];
        img_small_rot270_int[idx]        = (uint8_t)img_small_rot270[idx];
        img_small_blocks_int[idx]        = (uint8_t)img_small_blocks_rotations[0][idx];
        img_small_blocks_rot90_int[idx]  = (uint8_t)img_small_blocks_rotations[1][idx];
        img_small_blocks_rot180_int[idx] = (uint8_t)img_small_blocks_rotations[2][idx];
        img_small_blocks_rot270_int[idx] = (uint8_t)img_small_blocks_rotations[3][idx];
    }
    for (int idx = 0; idx < height*width; idx++) {
        img_gray_blocks_int[idx]        = (uint8_t)img_gray_blocks[idx];
    }
    //stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_rot90.png")).str().c_str(), width, height, 1, img_gray_int, width);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_rot0.png")).str().c_str(),   width_small,  height_small, 1, img_small_int,        width_small);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_rot90.png")).str().c_str(),  height_small, width_small,  1, img_small_rot90_int,  height_small);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_rot180.png")).str().c_str(), width_small,  height_small, 1, img_small_rot180_int, width_small);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_rot270.png")).str().c_str(), height_small, width_small,  1, img_small_rot270_int, height_small);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_gray_blocks.png")).str().c_str(), block_num_cols, img_gray_num_blocks * block_num_rows, 1, img_gray_blocks_int, block_num_cols);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_small_blocks.png")).str().c_str(), block_num_cols, img_small_num_blocks * block_num_rows, 1, img_small_blocks_int, block_num_cols);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_small_blocks_rot90.png")).str().c_str(), block_num_cols, img_small_num_blocks * block_num_rows, 1, img_small_blocks_rot90_int, block_num_cols);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_small_blocks_rot180.png")).str().c_str(), block_num_cols, img_small_num_blocks * block_num_rows, 1, img_small_blocks_rot180_int, block_num_cols);
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("test_small_blocks_rot270.png")).str().c_str(), block_num_cols, img_small_num_blocks * block_num_rows, 1, img_small_blocks_rot270_int, block_num_cols);
    free(img_small_int);
    free(img_small_rot90_int);
    free(img_small_rot180_int);
    free(img_small_rot270_int);
    free(img_gray_blocks_int);
    free(img_small_blocks_int);
    free(img_small_blocks_rot90_int);
    free(img_small_blocks_rot180_int);
    free(img_small_blocks_rot270_int);
    std::cout << "Done.\n" << std::endl;
    */
    //-----------------------------------DEBUGGING--------------------------------------------------------------------------------
    
    //free allocated memory
    //free(img_gray);
    delete[] mappings;
    delete[] img_gray_blocks;
    delete[] img_small_blocks_rotations[0];
    delete[] img_small_blocks_rotations[1];
    delete[] img_small_blocks_rotations[2];
    delete[] img_small_blocks_rotations[3];
    delete[] img_small_blocks_rotations;
    
    return 0;
}
