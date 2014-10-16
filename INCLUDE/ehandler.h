#ifndef EHANDLER_H
#define EHANDLER_H
/* Note:
 * Too much error checking slows down the prgram; 
 * hence, invoke error checking only in higher-level functions
 * or wrapper functions.
 *
 * */

/*
 * Use stringification (#) to log a warning.
 * */ 
#define warn_if(EXPR, fmt,...)                   \
        do { if(EXPR)                            \
            fprintf( stderr,                     \
                     "%s:%s: WARNING: (%s)=TRUE:"\
                     fmt "\n",                   \
                     __FILE__,                   \
                     __func__,                   \
                    #EXPR,                       \
                    __VA_ARGS__);                \
           } while (0)


/* Prints formatted error msg to stderr 
 * */
#define eprint(fmt, ...)                        \
        do {                                    \
            fprintf( stderr,                    \
                     "%s:%s: ERROR: " fmt "\n", \
                     __FILE__,                  \
                     __func__,                  \
                     __VA_ARGS__);              \
           } while (0)

#define eprintb(fmt, ...)                          \
        do {                                       \
            fprintf( stderr,                       \
                     "%s:%d:%s: ERROR: " fmt "\n", \
                     __FILE__,                     \
                     __LINE__,                     \
                     __func__,                     \
                     __VA_ARGS__);                 \
           } while (0)

/* Terminate execution of the program and print error msg to stderr 
 * */
#define term_exec(fmt, ...)                     \
        do {                                    \
            fprintf( stderr,                    \
                     "%s:%s: ERROR: " fmt "\n", \
                     __FILE__,                  \
                     __func__,                  \
                     __VA_ARGS__);              \
            exit(EXIT_FAILURE);                 \
           } while (0)

/* Provides line number in the error message 
 * */
#define term_execb(fmt, ...)                       \
        do {                                       \
            fprintf( stderr,                       \
                     "%s:%d:%s: ERROR: " fmt "\n", \
                     __FILE__,                     \
                     __LINE__,                     \
                     __func__,                     \
                     __VA_ARGS__);                 \
            exit(EXIT_FAILURE);                    \
           } while (0)


#define err_check(EXPR, LABEL, fmt,...)          \
        do { if(EXPR){                           \
            fprintf( stderr,                     \
                     "%s:%s: Error: (%s)=TRUE:"  \
                     fmt "\n",                   \
                     __FILE__,                   \
                     __func__,                   \
                    #EXPR,                       \
                    __VA_ARGS__);                \
            goto LABEL;}                         \
           } while (0)

#endif

