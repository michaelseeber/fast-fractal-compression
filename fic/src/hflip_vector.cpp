namespace reference {
    void hflip_vector(const int blockDimension, double *vector) {
        //in memory rotation of matrix not implemented
        double **rows = new double*[blockDimension];
        int index = 0;
        for(int i = 0; i < blockDimension; i++){
            rows[i] = new double[blockDimension];
            for(int j = 0; j < blockDimension; j++){
                rows[i][j] = vector[index];
                index++;
            }
        }
    
        //time to do hflip
        for(int i = 0; i < blockDimension; i++){//column
            for(int j = 0; j < blockDimension; j++){//row
                vector[i*blockDimension + (blockDimension-1 - j)] = rows[i][j];
            }
        }
    
        for(int i = 0; i < blockDimension; i++){
            delete [] rows[i];
        }
        delete [] rows;
    }
}

