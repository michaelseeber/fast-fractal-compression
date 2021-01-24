//#include <fic/common.h>
#include <cmath>
#include <limits>
#include <fic/map_one_matrix_to_another.h>
#include <fic/difference_norm.h>
#include <fic/rotate_vector.h>

namespace reference {
    void mapping_between_two_matrices(const int blockDimension, double *x, double *y, double *mapArray) {
        double min_error =  std::numeric_limits<double>::max();
        double *vec1 = new double[blockDimension*blockDimension];
    
        double *map = new double[3];
        map[2] = 0;//the rotation bit
    
        for(int i = 0; i < 4; i++){
            //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
            int n = blockDimension * blockDimension;
            double sumx = 0.0;
            double sumx2 = 0.0;
            double sumxy = 0.0;
            double sumy = 0.0;
            double sumy2 = 0.0;
    
            for(int i=0;i<n;i++){
                sumx  += x[i];
                sumy  += y[i];
                sumx2 += std::pow((x[i]),2);
                sumxy += x[i] * y[i];
                sumy2 += std::pow((y[i]),2);
            }
    
            double denom = (n * sumx2 - std::pow(sumx,2));
    
            map[0] = (n * sumxy  -  sumx * sumy) / denom;
            map[1] = (sumy * sumx2  -  sumx * sumxy) / denom;
    
            map_one_matrix_to_another(blockDimension, x, map, vec1);
    
            double error;
            difference_norm(blockDimension, y, vec1, &error);
    
            if(min_error > error){
                min_error = error;
                mapArray[0] = map[0];
                mapArray[1] = map[1];
                mapArray[2] = i;
            }
    
            rotate_vector(blockDimension, x);
        }

        delete []vec1;
        delete []map;
    }
}

namespace fast {
    void mapping_between_two_matrices(const int blockDimension, double *x, double *y, double *mapArray) {
        double min_error =  std::numeric_limits<double>::max();
        double *vec1 = new double[blockDimension*blockDimension];
    
        double *map = new double[3];
        map[2] = 0;//the rotation bit
    
        for(int i = 0; i < 4; i++){
            //http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
            int n = blockDimension * blockDimension;
            double sumx = 0.0;
            double sumx2 = 0.0;
            double sumxy = 0.0;
            double sumy = 0.0;
            double sumy2 = 0.0;
    
            for(int i=0;i<n;i++){
                sumx  += x[i];
                sumy  += y[i];
                sumx2 += std::pow((x[i]),2);
                sumxy += x[i] * y[i];
                sumy2 += std::pow((y[i]),2);
            }
    
            double denom = (n * sumx2 - std::pow(sumx,2));
    
            map[0] = (n * sumxy  -  sumx * sumy) / denom;
            map[1] = (sumy * sumx2  -  sumx * sumxy) / denom;
    
            map_one_matrix_to_another(blockDimension, x, map, vec1);
    
            double error;
            difference_norm(blockDimension, y, vec1, &error);
    
            if(min_error > error){
                min_error = error;
                mapArray[0] = map[0];
                mapArray[1] = map[1];
                mapArray[2] = i;
            }
    
            rotate_vector(blockDimension, x);
        }
    }
}
