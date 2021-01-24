#if !defined(__FIC_BLOCKWISE_SUM_OF_XMULY_H)
#define __FIC_BLOCKWISE_SUM_OF_XMULY_H

namespace reference {
    void blockwise_sum_of_xmuly(const int, const int, const int, double *, double **, double **);
}

namespace fast {
    void blockwise_sum_of_xmuly(const int, const int, const int, double *, double **, double **);
}


#endif // __FIC_BLOCKWISE_SUM_OF_XMULY_H
