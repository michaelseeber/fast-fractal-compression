//namespace reference {
//correct
    void rotate(const int height, const int width, double *img, double *img_rot90, double *img_rot180, double *img_rot270) {
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                //img_rot90[ j * height + i] = img[(height-1- i)  * width + j            ];
                img_rot90[ j * height + (height-1 - i)]            = img[i * width + j];
                //img_rot180[i * width + j] = img[(height-1 - i) * width  + (width-1 - j)];
                img_rot180[(height-1 - i) * width + (width-1 - j)] = img[i * width + j];
                //img_rot270[j * height + i] = img[i             * width  + (width-1 - j)];
                img_rot270[(width-1 - j) * height + i]             = img[i * width + j];
            }
        }
    }
    void rotate(const int height, const int width, float *img, float *img_rot90, float *img_rot180, float *img_rot270) {
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                //img_rot90[ j * height + i] = img[(height-1- i)  * width + j            ];
                img_rot90[ j * height + (height-1 - i)]            = img[i * width + j];
                //img_rot180[i * width + j] = img[(height-1 - i) * width  + (width-1 - j)];
                img_rot180[(height-1 - i) * width + (width-1 - j)] = img[i * width + j];
                //img_rot270[j * height + i] = img[i             * width  + (width-1 - j)];
                img_rot270[(width-1 - j) * height + i]             = img[i * width + j];
            }
        }
    }
/*
//--------------------WRONG implements same wrong compression as in mannyray repository----------------------------------------------------------
    void rotate(const int height, const int width, double *img, double *img_rot90, double *img_rot180, double *img_rot270) {
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                //img_rot90[ j * height + i] = img[(height-1- i)  * width + j            ];
                img_rot90[ j * height + i] = img[i * width + j];
                //img_rot180[i * width + j] = img[(height-1 - i) * width  + (width-1 - j)];
                img_rot180[i * width + j] = img[i * width + j];
                //img_rot270[j * height + i] = img[i             * width  + (width-1 - j)];
                img_rot270[j * height + i] = img[i * width + j];
            }
        }
    }
*/
    void rotate_and_hflip(const int height, const int width, float *img,
            float *img_rot90,
            float *img_rot180,
            float *img_rot270,
            float *img_hflip,
            float *img_hflip_rot90,
            float *img_hflip_rot180,
            float *img_hflip_rot270) {
        // TODO: horizontal fip = flip around vertical axis?
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                img_rot90[ j                    * height + (height-1 - i)] = img[i * width + j];
                img_rot180[(height-1 - i)       * width  + (width-1 - j)]  = img[i * width + j];
                img_rot270[(width-1 - j)        * height + i]              = img[i * width + j];
                img_hflip[i                     * width  + (width-1 - j)]  = img[i * width + j];
                img_hflip_rot90[(width-1 - j)   * height + (height-1 - i)] = img[i * width + j];
                img_hflip_rot180[(height-1 - i) * width  + j]              = img[i * width + j];
                img_hflip_rot270[j              * height + i]              = img[i * width + j];
            }
        }
    }
//}
