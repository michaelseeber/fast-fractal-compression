#if !defined(__FIC_SPLIT_INTO_BLOCKS_H)
#define __FIC_SPLIT_INTO_BLOCKS_H

//split into non-overlapping blocks
void split_into_blocks(             const int, const int, const int, const int, double *, double *);
void split_into_blocks_rot90(       const int, const int, const int, const int, double *, double *);
void split_into_blocks_rot180(      const int, const int, const int, const int, double *, double *);
void split_into_blocks_rot270(      const int, const int, const int, const int, double *, double *);
void split_into_blocks(             const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot90(       const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot180(      const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot270(      const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip(       const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot90( const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot180(const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot270(const int, const int, const int, const int, float *,  float *);
//split into overlapping blocks using the stride argument
void split_into_blocks(             const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot90(       const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot180(      const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_rot270(      const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip(       const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot90( const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot180(const int, const int, const int, const int, const int, float *,  float *);
void split_into_blocks_hflip_rot270(const int, const int, const int, const int, const int, float *,  float *);


#endif // __FIC_SPLIT_INTO_BLOCKS_H
