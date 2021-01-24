#include <fic/rotate_vector.h>
#include <fic/hflip_vector.h>

namespace reference {
    void map_one_matrix_to_another_with_flip(const int blockDimension, double *one, double *result_array, double *two) {//TODO one & two in right order?
        if (result_array[2] < 4) {
            //rotate the correct amount of times
            for(int i = 0; i < result_array[2]; i++) {
                reference::rotate_vector(blockDimension, one);
            }
    
            for(int i = 0; i < blockDimension*blockDimension; i++) {
                two[i] = one[i] * result_array[0] + result_array[1];
            }
            
            //always rotate back to original
            if(result_array[2] == 0) {
                return;
            }
            
            for(int i = result_array[2]; i < 4; i++) {
                reference::rotate_vector(blockDimension, one);
            } 
            
        } else {
            reference::hflip_vector(blockDimension, one);
            
            //rotate the correct amount of times
            for(int i = 4; i < result_array[2]; i++) {
                reference::rotate_vector(blockDimension, one);
            }
    
            for(int i = 0; i < blockDimension*blockDimension; i++) {
                two[i] = one[i] * result_array[0] + result_array[1];
            }
            
            //always rotate back to original
            if(result_array[2] == 4) {
                reference::hflip_vector(blockDimension, one);
                return;
            }
            
            for(int i = result_array[2]; i < 8; i++) {
                reference::rotate_vector(blockDimension, one);
            } 
            reference::hflip_vector(blockDimension, one);
        }
    }
}

//namespace fast {
//}
