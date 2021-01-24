namespace reference {
    void rotate_vector(const int blockDimension, double *vector) {
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
    
        //time to do rotate
        for(int i = 0; i < blockDimension; i++){//column
            for(int j = 0; j < blockDimension; j++){//row
                vector[j*blockDimension + (blockDimension-1 - i)] = rows[i][j];
            }
        }
    
        for(int i = 0; i < blockDimension; i++){
            delete [] rows[i];
        }
        delete [] rows;
    }
}

namespace fast {
    void rotate_vector(const int blockDimension, double *vector) {
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
    
        //time to do rotate
        for(int i = 0; i < blockDimension; i++){//column
            for(int j = 0; j < blockDimension; j++){//row
                //vector[i*blockDimension + j] = rows[j][i];
                vector[j*blockDimension + (blockDimension-1 - i)] = rows[i][j];
            }
        }
    
        for(int i = 0; i < blockDimension; i++){
            delete [] rows[i];
        }
        delete [] rows;
    }
}
