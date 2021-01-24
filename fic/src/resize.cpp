//namespace reference {
    void resize(const int height, const int width, double *image, double *image_small) {
        const int height_small = height / 2;
        const int width_small = width / 2;
        for(int i = 0; i < height_small; i++){
            for(int j = 0; j < width_small; j++){
                image_small[i * width_small + j]  = image[2*i   * width + 2*j    ];
                image_small[i * width_small + j] += image[2*i   * width + 2*j + 1];
                image_small[i * width_small + j] += image[(2*i) * width + 2*j    ];
                image_small[i * width_small + j] += image[(2*i) * width + 2*j + 1];
                image_small[i * width_small + j] /= 4;
            }
        }
    }
    void resize(const int height, const int width, float *image, float *image_small) {
        const int height_small = height / 2;
        const int width_small = width / 2;
        for(int i = 0; i < height_small; i++){
            for(int j = 0; j < width_small; j++){
                image_small[i * width_small + j]  = image[2*i   * width + 2*j    ];
                image_small[i * width_small + j] += image[2*i   * width + 2*j + 1];
                image_small[i * width_small + j] += image[(2*i) * width + 2*j    ];
                image_small[i * width_small + j] += image[(2*i) * width + 2*j + 1];
                image_small[i * width_small + j] /= 4;
            }
        }
    }
//}
