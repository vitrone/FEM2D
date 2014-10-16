/*============================================================================+/
 | pthpool.c
 | Defines library functions for creating thread pool. 
/+============================================================================*/
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#define NDEBUG
#include "assert.h"
#include "pthpool.h"
/*============================================================================*/

static void* pthpool_schedule_task(void *mp);
/*============================================================================*/

static void* pthpool_schedule_task(void *mp)
{
    pthpool_data_t* pth = (pthpool_data_t*) mp;
    debug_enter("Executing thread index: %d", pth->thread_index);
    pthpool_task_t task; 

    clock_t start, duration;
    
    while(1) /* Keep the thread waiting in an infinite loop */ 
    {
        debug_body("Entering infinite loop (thread: %d)", pth->thread_index);
        pthread_mutex_lock(&(pth->lock));
        while((pth->action) == PTHPOOL_WAIT)
        {
            debug_body("waiting (thread: %d)", pth->thread_index);
            /* Keep waiting on the condition variable until the thread is awaken */ 
            /* Here the mutex gets unlocked */ 
            pthread_cond_wait(&(pth->notify), &(pth->lock)); 
            /* on returning we regain the lock */ 
            debug_body("notice recieved (thread: %d)", pth->thread_index);
        }
        if((pth->action)==PTHPOOL_EXIT)
        {
            debug_body("exit request made for thread: %d", pth->thread_index);
            break;
        }

        /* do the task */
        debug_body("performing task with thread: %d", pth->thread_index);
        task = *(pth->task);
        pthread_mutex_unlock(&(pth->lock));
        task.function(task.argument);

        /* Check if any threads are waiting for finishing the task */ 
        debug_body("task completed with thread: %d", pth->thread_index);
        pthread_mutex_lock(&(pth->lock));
        if((pth->action) == PTHPOOL_NOTIFY)
        {
            debug_body("notify thread: %d", pth->thread_index);
            (pth->action) = PTHPOOL_WAIT;
            debug_body("thread index: %d, action : WAIT", pth->thread_index);
            pthread_cond_signal(&(pth->notify)); 
            pthread_mutex_unlock(&(pth->lock));
        }
        else
        {
            (pth->action) = PTHPOOL_WAIT;
            debug_body("thread index: %d, action : WAIT", pth->thread_index);
            pthread_mutex_unlock(&(pth->lock));
        }
    }
    pthread_mutex_unlock(&(pth->lock));
    debug_exit("exiting thread: %d", pth->thread_index);
    pthread_exit(NULL);
    return(NULL);
}
/*============================================================================*/

void pthpool_create_threads
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp
)
{
    debug_enter("Number of threads: %d", num_threads);
    matlib_index i, j;
    pthread_attr_t attr;

    /* initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int pthread_r;

    for(i=0; i<num_threads; i++)
    {
        mp[i].thread_index = i;
        /* set the CPU affinity */ 
        CPU_ZERO(&mp[i].cpu);
        CPU_SET( i%MAX_NUM_CPU, &mp[i].cpu);
        pthread_r = pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &mp[i].cpu);
        
        pthread_mutex_init(&(mp[i].lock), NULL);
        pthread_cond_init(&(mp[i].notify), NULL);
        mp[i].action = PTHPOOL_WAIT;
        debug_body("thread index: %d, action : WAIT", mp[i].thread_index);
        pthread_r = pthread_create( &(mp[i].thread), 
                                    &attr, 
                                    pthpool_schedule_task, 
                                    (void *) &mp[i]); 

        //pthread_setaffinity_np(mp[i].thread, sizeof(cpu_set_t), &mp[i].cpu);
        if(pthread_r)
        {
            term_exec("failed to create the pthread (return value: %d)", pthread_r);
        }

        BEGIN_DEBUG
            cpu_set_t cpuset_tmp;

            pthread_getaffinity_np(mp[i].thread, sizeof(cpu_set_t), &cpuset_tmp);
            if (CPU_ISSET(i, &cpuset_tmp))
            {
                debug_print("thread: %d, CPU affinity: %d", i, i);
            }
        END_DEBUG
    }
    debug_exit("%s", "");
}
/*============================================================================*/

void pthpool_exec_task
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp, 
    pthpool_task_t* task
)
{
    debug_enter("Number of threads: %d", num_threads);
    matlib_index i;
    for(i=0; i<num_threads; i++)
    {
        debug_enter("pthpool_data_t*: %p", &mp[i]);
        pthread_mutex_lock(&(mp[i].lock));
        debug_body("thread id: %d", i);

        mp[i].task = &task[i];
        (mp[i].action) = PTHPOOL_PERFORM;
        debug_body("thread index: %d, action : PERFORM", mp[i].thread_index);

        pthread_cond_signal(&(mp[i].notify));
        pthread_mutex_unlock(&(mp[i].lock));
        
    }
    //sleep(1);
    for(i=0; i<num_threads; i++)
    {
        pthread_mutex_lock(&(mp[i].lock));
        if((mp[i].action) == PTHPOOL_PERFORM)
        {
            debug_body("thread id: %d", i);
            (mp[i].action) = PTHPOOL_NOTIFY;
            debug_body("thread index: %d, action : NOTIFY", mp[i].thread_index);
            debug_body("Going to wait for thread: %d", mp[i].thread_index);
            pthread_cond_wait(&(mp[i].notify), &(mp[i].lock));
            debug_body("completed thread: %d", mp[i].thread_index);
        }
        debug_body("Thread id: %d, unlocking", i); 
        pthread_mutex_unlock(&(mp[i].lock));
        
    }
    debug_exit("%s", "");
}
/*============================================================================*/

void pthpool_exec_task_nosync
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp, 
    pthpool_task_t* task
)
{
    debug_enter("Number of threads: %d", num_threads);
    matlib_index i;
    for(i=0; i<num_threads; i++)
    {
        debug_enter("pthpool_data_t*: %p", &mp[i]);
        pthread_mutex_lock(&(mp[i].lock));
        debug_body("thread id: %d", i);

        mp[i].task = &task[i];
        (mp[i].action) = PTHPOOL_PERFORM;
        debug_body("thread index: %d, action : PERFORM", mp[i].thread_index);

        pthread_cond_signal(&(mp[i].notify));
        pthread_mutex_unlock(&(mp[i].lock));
        
    }
    debug_exit("%s", "");
}

void pthpool_sync_threads
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp 
)
{
    debug_enter("Number of threads: %d", num_threads);
    matlib_index i;
    for(i=0; i<num_threads; i++)
    {
        pthread_mutex_lock(&(mp[i].lock));
        debug_body("thread id: %d", i);

        if((mp[i].action) == PTHPOOL_PERFORM)
        {
            (mp[i].action) = PTHPOOL_NOTIFY;
            debug_body("thread index: %d, action : NOTIFY", mp[i].thread_index);
            debug_body("Going to wait for thread: %d", mp[i].thread_index);
            pthread_cond_wait(&(mp[i].notify), &(mp[i].lock));
        }
        debug_body("completed thread: %d", mp[i].thread_index);
        pthread_mutex_unlock(&(mp[i].lock));
    }
    debug_exit("%s", "");
}

/*============================================================================*/

void pthpool_destroy_threads
( 
    matlib_index    num_threads, 
    pthpool_data_t* mp 
)
{
    debug_enter("number of threads: %d", num_threads);
    matlib_index pthread_r, i;
    void* status;
    for( i=0; i<num_threads; i++)
    {
        debug_body("thread id: %d", i);
        pthread_mutex_lock(&(mp[i].lock));
        if((mp[i].action) == PTHPOOL_PERFORM)
        {
            (mp[i].action) = PTHPOOL_NOTIFY;
            debug_body("thread index: %d, action : NOTIFY", mp[i].thread_index);
            debug_body("Going to wait for task to finish (thread index: %d)", i);
            pthread_cond_wait(&(mp[i].notify), &(mp[i].lock)); 
            (mp[i].action) = PTHPOOL_EXIT;
            debug_body("thread index: %d, action : EXIT", mp[i].thread_index);
        }
        else if((mp[i].action) == PTHPOOL_WAIT)
        {
            (mp[i].action) = PTHPOOL_EXIT;
            debug_body("thread index: %d, action : EXIT", mp[i].thread_index);
            pthread_cond_signal(&(mp[i].notify));
        }
        pthread_mutex_unlock(&(mp[i].lock));
    }
    for(i=0; i<num_threads; i++)
    {
        pthread_r = pthread_join((mp[i].thread), &status); 
    }
    for(i=0; i<num_threads; i++)
    {
        pthread_r = pthread_mutex_destroy(&(mp[i].lock));
        pthread_r = pthread_cond_destroy(&(mp[i].notify));
    }
    debug_exit("%s", "");
}

/*============================================================================*/

void pthpool_func
(
    matlib_index*   Np,
    void**          shared_data,
    void*           thfunc,
    matlib_index    num_threads,
    pthpool_data_t* mp

)
{
    debug_enter("Np: [%d, %d]", Np[0], Np[1]);

    matlib_index i; 
    /* define the shared data */ 

    matlib_index nsdata[num_threads][2];

    pthpool_arg_t  arg[num_threads];
    pthpool_task_t task[num_threads];


    for(i=0; i<num_threads-1; i++)
    {
        nsdata[i][0] = i*Np[0];
        nsdata[i][1] = (i+1)*Np[0];
        arg[i].shared_data    = shared_data; 
        arg[i].nonshared_data = (void**)&nsdata[i];
        arg[i].thread_index   = i;
        /* Define the task */ 
        task[i].function  = (void*)thfunc;
        task[i].argument  = &arg[i];
    }
    nsdata[i][0] = i*Np[0];
    nsdata[i][1] = Np[1];
    arg[i].shared_data    = shared_data; 
    arg[i].nonshared_data = (void**)&nsdata[i];
    arg[i].thread_index   = i;
    /* Define the task */ 
    task[i].function  = (void*)thfunc;
    task[i].argument  = &arg[i];

    debug_body("%s", "created task");

    pthpool_exec_task(num_threads, mp, task);
    
    debug_exit("%s", "");
}

/*============================================================================*/



