//void read_matrix(const int height, const int width, const int blockDimension, const double *image, int *blockCount, double ***blockStore) {
void read_matrix(const int height, const int width, const int blockDimension, const double *image, int *blockRowCount, int *blockColCount, double ***blockStore) {
    //int blockRowCount = height / blockDimension;
    //int blockColCount = width / blockDimension;
    //*blockCount = blockRowCount * blockColCount;
    //*blockStore = new double *[*blockCount];
    //for(int i = 0; i < *blockCount; i++){
    *blockRowCount = height / blockDimension;
    *blockColCount = width / blockDimension;
    int blockCount = *blockRowCount * *blockColCount;
    *blockStore = new double *[blockCount];
    for(int i = 0; i < blockCount; i++){
        (*blockStore)[i] = new double[blockDimension*blockDimension];
    }
    int idx = 0;
    //for each 'block row'
    for(int i = 0; i < *blockColCount; i++){
        //each 'block row' consists of blockDimension number of pixel rows
        for(int j = 0; j< blockDimension; j++){
            //each 'block collumn' section within a row
            for(int k = 0; k < *blockRowCount;k++){
                //each 'block collumn is blockDimension pixels wide
                for(int m = 0; m < blockDimension; m++){
                    //populate the array correctly in order that it is read
                    (*blockStore)[i*(*blockRowCount) + k][j*blockDimension + m] = image[idx];
                    idx++;
                }
            }
        }
    }
    
}
