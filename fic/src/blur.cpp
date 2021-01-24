#include <iostream>

void blur_matrix(const int blockDimensionOriginal, double **target) {
    int blockDimension = blockDimensionOriginal/2;
    double *a = new double[blockDimension*blockDimension];
    for(int i = 0; i < blockDimension; i++){
        for(int j = 0; j < blockDimension; j++){
            a[i*blockDimension + j]  = (*target)[i*4*    blockDimension + j*2];
            a[i*blockDimension + j] += (*target)[i*4*    blockDimension + j*2 + 1] ;
            a[i*blockDimension + j] += (*target)[(i*4+2)*blockDimension + j*2];
            a[i*blockDimension + j] += (*target)[(i*4+2)*blockDimension + j*2 + 1];
            a[i*blockDimension + j] /= 4;
        }
    }
    delete [] *target;
    *target = a;
}

void blur(const int blockCount, int *blockDimension, double **blockStore) {
    for(int i = 0; i < blockCount; i++){
        //double * ap = blur_matrix(blockDimension, blockStore[i]);
        //delete [] blockStore[i];
        //blockStore[i] = ap;
        blur_matrix(*blockDimension, &(blockStore[i]));
        //blur_matrix(*blockDimension, blockStore+i);
    }
    *blockDimension /= 2;
    //column count and row count should still stay the same
}
