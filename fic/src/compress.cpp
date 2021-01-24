#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <fic/common.h>
#include <fic/block_store.h>
#include <fic/blur.h>
#include <fic/difference_norm.h>
#include <fic/map_one_matrix_to_another.h>
#include <fic/mapping_between_two_matrices.h>
#include <fstream>
#include <chrono>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb/stb_image_write.h"
#include "fic/tsc_x86.h"

using namespace reference;
//using namespace fast;

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
            << "\", expected an extension of type .png, .gif, .jpg, or .bmp" << std::endl;
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
    
    //convert to double array, ensure grayscale
    double *img_gray = (double *)malloc(height*width*sizeof(double));
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
            img_gray[idx] = (double)img[idx];
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

    //long num_flops = (long) height * height * width * width / (4 * block_num_pixels_c * block_num_pixels_c) // outer 2 loops of nested loop
    //    * (long) (3 + num_rot * (7 * block_num_pixels_c + 6)) // 3rd nested loop
    //    + (long) 7 * height * width / 4; // outside the nested loop
    //long num_divs = (long) height * height * width * width / (num_rot * block_num_pixels_c * block_num_pixels_c) * 16;
    //long num_sqrts = (long) height * height * width * width / (num_rot * block_num_pixels_c * block_num_pixels_c) * 4;
    //num_flops = num_flops + num_divs + num_sqrts;

    long long int num_flops = 0;
    // mappings_between_two_matrices
    num_flops = num_flops + block_num_pixels_c * 8; // loop n
    num_flops = num_flops + 3; // compute denom
    num_flops = num_flops + 3 + 3; // compute map[0] and map[1]
    num_flops = num_flops + block_num_pixels_c * 2; // map_one_matrix_to_another
    num_flops = num_flops + block_num_pixels_c * 3; // difference_norm
    num_flops = num_flops * num_rot; // because it's in a loop

    // other stuff inside nested loop
    num_flops = num_flops + block_num_pixels_c * 2; // map_one_matrix_to_another
    num_flops = num_flops + block_num_pixels_c * 3; // difference_norm

    // multiplications
    num_flops = num_flops * img_small_num_blocks_c; // 2nd outermost loop
    num_flops = num_flops * img_gray_num_blocks_c; // 3rd outermost loop

    long long int num_divs = 0;
    // mappings_between_two_matrices
    num_divs = num_divs + 2; // compute denom
    num_divs = num_divs + 2; // div-rence norm
    num_divs = num_divs * num_rot; // because it's in a loop

    // other stuff inside nested loop
    num_divs = num_divs + 2; // div-rence norm

    // multiplications
    num_divs = num_divs * img_small_num_blocks_c; // 2nd outermost loop
    num_divs = num_divs * img_gray_num_blocks_c; // 3rd outermost loop

    long long int num_sqrts = 0;
    // mappings_between_two_matrices
    num_sqrts = num_sqrts + 1; // div-rence norm
    num_sqrts = num_sqrts * num_rot; // because it's in a loop

    // other stuff inside nested loop
    num_sqrts = num_sqrts + 1; // div-rence norm
    
    // multiplications
    num_sqrts = num_sqrts * img_small_num_blocks_c; // 2nd outermost loop
    num_sqrts = num_sqrts * img_gray_num_blocks_c; // 3rd outermost loop


    // compute total number of flops
    long long int total_num_flops = num_flops + num_divs + num_sqrts;

    //compress
    myInt64 start, end;

    std::cout << "Start: compress image" << std::endl;

    auto t0 = std::chrono::steady_clock::now();
    start = start_tsc();

    int blockDimension8 = 8;
    int blockDimension16 = 16;
    int blockRowCount8, blockColCount8, blockRowCount16, blockColCount16;
    double **blockStore8, **blockStore16;

    read_matrix(height, width, blockDimension8, img_gray, &blockRowCount8, &blockColCount8, &blockStore8);
    read_matrix(height, width, blockDimension16, img_gray, &blockRowCount16, &blockColCount16, &blockStore16);
    int blockCount8 = blockRowCount8 * blockColCount8;
    int blockCount16 = blockRowCount16 * blockColCount16;
    blur(blockCount16, &blockDimension16, blockStore16);
    
    int blockMaps = blockCount8;
    int mapStore = 5;
    int mapSize = blockMaps;
    double **mappingData = new double *[blockMaps];
    for(int i = 0; i < blockMaps; i++){
        mappingData[i] = new double[mapStore];
        for(int j = 0; j < mapStore; j++){
            //for comparison purposes later, if you are trying to come up with some sort of minimal mapping
            //then it makes sense to store some variables as worst case possible maxs values
            mappingData[i][j] = std::numeric_limits<double>::max();
        }
    }
    
    //find optimal mapping 
    std::cout << "blockCount = " << blockCount8 << std::endl;
    std::cout << "\tStart: find optimal mapping" << std::endl;
    double *dataArray = new double[blockDimension8 * blockDimension8];
    double *mapArray = new double[mapStore - 2];
    for(int i = 0; i < blockCount16; i++) {    
        for(int j = 0; j < blockCount8; j++) {
            mapping_between_two_matrices(blockDimension8, blockStore16[i], blockStore8[j], mapArray);
            map_one_matrix_to_another(blockDimension8, blockStore16[i], mapArray, dataArray);
            double error;
            difference_norm(blockDimension8, blockStore8[j], dataArray, &error);
            if(mappingData[j][0] > error){
                mappingData[j][0] = error;
                for(int k = 1; k < 1 + mapStore - 2; k++){
                    mappingData[j][k] = mapArray[k-1];
                }
                mappingData[j][mapStore-1] = i;
            }
        }
    }

    end = stop_tsc(start);
    auto t1 = std::chrono::steady_clock::now();
    double ns =  (double)std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    delete [] mapArray;
    delete [] dataArray;
    std::cout << "\tDone." << std::endl;
    //std::cout << "Done.\n" << std::endl;
    std::cout << "Done.\t" << 1e-9 * ns << " s\n" << std::endl;

    //output flops/cycle
    double perf = round((100.0 * total_num_flops) / (double)end) / 100.0;
    std::cout << "Running: compress_reference " << perf << " flops/ cycles" << std::endl;

    std::cout << "Start: write fic_compress_output" << std::endl;
    std::ofstream ofile((resolver.resolve("images").make_absolute()/filesystem::path("fic_compress_output")).str());
    ofile << stride << " " << blockDimension8 << " " << blockColCount8 << " " << blockRowCount8 << std::endl;
    for(int i = 0; i < mapSize; i++){
        for(int j = 1; j < mapStore; j++){
            ofile << mappingData[i][j] << " ";
        }
        ofile << std::endl;
    }
    ofile.close();
    std::cout << "Done.\n" << std::endl;

    /*
    //write PNGs
    //let us find the absolute path of the projects image folder using the resolver class
    std::cout << "Start: write fic_output.png" << std::endl;
    unsigned char *img_gray_int = (unsigned char *)malloc(height*width*sizeof(unsigned char));
    for (int idx = 0; idx < height*width; idx++) {
        img_gray_int[idx] = (uint8_t)img_gray[idx];
    }
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("fic_output.png")).str().c_str(), width, height, 1, img_gray_int, width);
    free(img_gray_int);
    std::cout << "Done.\n" << std::endl;
    */
    
    //free allocated memory
    free(img_gray);
    
    return 0;
}
