#ifndef BASIC_H
#define BASIC_H

/* Basic data types */
#define Index matlib_index


/* Counting time */
#define START_TIMMING(tb)                            \
    do {                                             \
        clock_gettime(CLOCK_REALTIME, &tb);          \
                                                     \
    } while (0)                            

#define GET_DURATION(tb, te, dt)                    \
    do {                                            \
        clock_gettime(CLOCK_REALTIME, &te);         \
        dt = (double)(te.tv_sec-tb.tv_sec)*1.0e3 +  \
             (double)(te.tv_nsec-tb.tv_nsec)/1.0e6; \
                                                    \
    } while (0)                                     


#endif
