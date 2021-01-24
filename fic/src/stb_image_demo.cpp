#include <filesystem/path.h>
#include <filesystem/resolver.h>
#include <fic/common.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

int main(int argc, char **argv) {
    //can be used to resolve relative paths
    filesystem::resolver resolver;
    resolver.append("");
    resolver.append("../");
    resolver.append("fic/");
    
    //load image
    int height, width, channels;
    unsigned char *img;
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
    
    //modify image (example)
    std::cout << "Start: modify color image" << std::endl;
    //color 10th-100th row red
    for (int i = 10; i < 100; i++) {
        for (int j = 0; j < width; j++) {
            img[i * channels * width + j * channels + 0] = 255;//r
            img[i * channels * width + j * channels + 1] = 0;//g
            img[i * channels * width + j * channels + 2] = 0;//b
        }
    }
    //color 50th-100th column blue
    for (int i = 0; i < height; i++) {
        for (int j = 50; j < 100; j++) {
            img[i * channels * width + j * channels + 0] = 0;//r
            img[i * channels * width + j * channels + 1] = 0;//g
            img[i * channels * width + j * channels + 2] = 255;//b
        }
    }
    std::cout << "Done.\n" << std::endl;
    
    //convert color image to grayscale
    std::cout << "Start: convert color image to grayscale" << std::endl;
    unsigned char *img_gray = (unsigned char *)malloc(height*width*sizeof(unsigned char));
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            img_gray[i * width + j] = \
                0.299 * (uint8_t)img[i * channels * width + j * channels + 0] +
                0.587 * (uint8_t)img[i * channels * width + j * channels + 1] +
                0.114 * (uint8_t)img[i * channels * width + j * channels + 2];
        }
    }
    std::cout << "Done.\n" << std::endl;

    //write PNGs
    //let us find the absolute path of the projects image folder using the resolver class
    std::cout << "Start: write demo_gray_output.png" << std::endl;
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("demo_gray_output.png")).str().c_str(), width, height, 1, img_gray, width);
    std::cout << "Done.\n" << std::endl;
    std::cout << "Start: write demo_color_output.png" << std::endl;
    stbi_write_png((resolver.resolve("images").make_absolute()/filesystem::path("demo_color_output.png")).str().c_str(), width, height, channels, img, width * channels);
    std::cout << "Done.\n" << std::endl;
    
    //free allocated memory
    stbi_image_free(img);
    free(img_gray);
    
    return 0;
}
