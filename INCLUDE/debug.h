#ifndef DEBUG_H
#define DEBUG_H

/* NDEBUG is recognized by assert.h
 * */

#ifdef NDEBUG
    #define DEBUG_TEST 0
#else
    #define DEBUG_TEST 1
#endif

#define debug_enter(fmt, ...)                   \
        do { if (DEBUG_TEST)                    \
                fprintf( stderr,                \
                         "%s:%s: >(" fmt ")\n", \
                         __FILE__,              \
                         __func__,              \
                         __VA_ARGS__);          \
           } while (0)


#define debug_exit(fmt, ...)                    \
        do { if (DEBUG_TEST)                    \
                fprintf( stderr,                \
                         "%s:%s: <(" fmt ")\n", \
                         __FILE__,              \
                         __func__,              \
                         __VA_ARGS__);          \
           } while (0)


#define debug_body(fmt, ...)                       \
        do { if (DEBUG_TEST)                       \
                fprintf( stderr,                   \
                         "%s:%d:%s: -(" fmt ")\n", \
                         __FILE__,                 \
                         __LINE__,                 \
                         __func__,                 \
                         __VA_ARGS__);             \
           } while (0)


#define debug_print(fmt, ...)                  \
        do {                                   \
            fprintf( stderr,                   \
                     "%s:%d:%s: -(" fmt ")\n", \
                     __FILE__,                 \
                     __LINE__,                 \
                     __func__,                 \
                     __VA_ARGS__);             \
           } while (0)


/* Writing code blocks */
#define BEGIN_DEBUG if(DEBUG_TEST){
#define END_DEBUG }


#endif
