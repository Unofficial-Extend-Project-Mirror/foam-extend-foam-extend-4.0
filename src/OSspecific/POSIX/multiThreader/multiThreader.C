/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    multiThreader

Description
    Implementation of the multiThreader class

Author
    Sandeep Menon
    University of Massachusetts Amherst

\*----------------------------------------------------------------------------*/

#include "multiThreader.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IOmultiThreader, 0);

bool multiThreader::debug = false;
bool Mutex::debug = false;
bool rwMutex::debug = false;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiThreader::multiThreader(int numThreads)
:
    maxQueueSize_(10),
    poolInfo_(NULL)
{
    if (numThreads > 0)
    {
        numThreads_ = numThreads;

        if (debug)
        {
            Info << "Initializing threading environment with "
                 << numThreads_ << " threads." << endl;
        }
    }
    else
    {
        // Default number of threads at one (single-threaded)
        numThreads_ = 1;

        if (debug)
        {
            Info << "Defaulting threading environment to one thread." << endl;
        }
    }

    // Initialize the thread pool
    initializeThreadPool();
}


Mutex::Mutex()
{
    // Set attributes based on debug flag
    pthread_mutexattr_t attribute;
    pthread_mutexattr_init(&attribute);

    if (debug)
    {
        pthread_mutexattr_settype(&attribute, PTHREAD_MUTEX_ERRORCHECK);
    }
    else
    {
        pthread_mutexattr_settype(&attribute, PTHREAD_MUTEX_NORMAL);
    }

    if (pthread_mutex_init(&lock_, &attribute))
    {
        FatalErrorIn("multiThreader::Mutex::Mutex()")
            << "Unable to initialize mutex"
            << abort(FatalError);
    }

    // Destroy the attribute
    pthread_mutexattr_destroy(&attribute);
}


rwMutex::rwMutex()
{
    // Set attributes for the mutex
    pthread_rwlockattr_t attribute;
    pthread_rwlockattr_init(&attribute);

    // Set the attribute type
    // Writer-preferred appears to be the default,
    // but set it explicitly anyway.
#ifndef darwin
    pthread_rwlockattr_setkind_np(&attribute, PTHREAD_RWLOCK_PREFER_WRITER_NP);
#endif

    if (pthread_rwlock_init(&lock_, &attribute))
    {
        FatalErrorIn("multiThreader::rwMutex::rwMutex()")
            << "Unable to initialize read-write mutex"
            << abort(FatalError);
    }

    // Destroy the attribute
    pthread_rwlockattr_destroy(&attribute);
}


Conditional::Conditional()
{
    if (pthread_cond_init(&condition_, NULL))
    {
        FatalErrorIn("multiThreader::Conditional::Conditional()")
            << "Unable to initialize condition"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

multiThreader::~multiThreader()
{
    destroyThreadPool();
}


Mutex::~Mutex()
{
    if (pthread_mutex_destroy(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::~Mutex()")
            << "Unable to destroy mutex"
            << abort(FatalError);
    }
}


rwMutex::~rwMutex()
{
    if (pthread_rwlock_destroy(&lock_))
    {
        FatalErrorIn("multiThreader::rwMutex::~rwMutex()")
            << "Unable to destroy read-write mutex"
            << abort(FatalError);
    }
}


Conditional::~Conditional()
{
    if (pthread_cond_destroy(&condition_))
    {
        FatalErrorIn("multiThreader::Conditional::~Conditional()")
            << "Unable to destroy condition"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void multiThreader::initializeThreadPool()
{
    // Initialize threads only if multi-threaded
    if (multiThreaded())
    {
        // Allocate the threadPool structure
        poolInfo_ = new threadPool;

        // Initialize fields
        poolInfo_->threader = this;
        poolInfo_->numThreads = numThreads_;
        poolInfo_->queueSize = 0;
        poolInfo_->threads = new pthread_t[numThreads_];
        poolInfo_->head = NULL;
        poolInfo_->tail = NULL;

        // Initialize flags
        poolInfo_->queueClosed = false;
        poolInfo_->shutDown = false;

        // Initialize thread attributes
        pthread_attr_init(&(poolInfo_->attr));
        pthread_attr_setdetachstate
        (
            &(poolInfo_->attr),
            PTHREAD_CREATE_JOINABLE
        );

        // Create worker threads and have them wait for jobs
        for (int tIndex = 0; tIndex < numThreads_; tIndex++)
        {
            int status = pthread_create
                         (
                             &(poolInfo_->threads[tIndex]),
                             &(poolInfo_->attr),
                             reinterpret_cast<externThreadFunctionType>
                             (
                                 poolThread
                             ),
                             reinterpret_cast<void *>
                             (
                                 poolInfo_
                             )
                         );

            if (status != 0)
            {
                FatalErrorIn("multiThreader::initializeThreadPool()")
                    << "pthread_create could not initialize thread: "
                    << tIndex
                    << abort(FatalError);
            }
        }
    }
}


threadReturnType multiThreader::poolThread(void *arg)
{
    // Typecast the argument into the required structure
    threadPool *poolInfo = reinterpret_cast<threadPool *>(arg);

    // Work queue loop
    while (true)
    {
        // Lock the work queue
        poolInfo->queueLock.lock();

        // Wait for work to arrive in the queue
        while ((poolInfo->queueSize == 0) && (!poolInfo->shutDown))
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info << "poolThread::Wait on queueNotEmpty." << endl;
            }
#           endif
            poolInfo->threader->waitForCondition
                                (
                                    poolInfo->queueNotEmpty,
                                    poolInfo->queueLock
                                );
        }

        // Check for shutdown
        if (poolInfo->shutDown)
        {
            poolInfo->queueLock.unlock();
            pthread_exit(NULL);
        }

        // Pick an item off the queue, and get to work
        workQueueItem *myWorkItem = poolInfo->head;
        poolInfo->queueSize--;
        if (poolInfo->queueSize == 0)
        {
            poolInfo->head = poolInfo->tail = NULL;
        }
        else
        {
            poolInfo->head = myWorkItem->next;
        }

        // Handle a waiting destructor
        if (poolInfo->queueSize == 0)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info << "poolThread::Signaling: Empty queue." << endl;
            }
#           endif
            poolInfo->threader->signal(poolInfo->queueEmpty);
        }

        // Unlock the work queue
        poolInfo->queueLock.unlock();

        // Perform the work
        myWorkItem->function(myWorkItem->arg);

        // Free up the work item
        delete myWorkItem;
    }

    return threadReturnValue;
}


void multiThreader::addToWorkQueue
(
    void (*tFunction)(void*),
    void *arg
) const
{
    if (singleThreaded())
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "addToWorkQueue:: Not in multiThreaded mode." << endl;
        }
#       endif

        return;
    }

    // Lock the work queue
    poolInfo_->queueLock.lock();

    // If occupied, wait for the queue to free-up
    while
    (
         (poolInfo_->queueSize == maxQueueSize_)
      && (!(poolInfo_->shutDown || poolInfo_->queueClosed))
    )
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info << "addToWorkQueue:: Wait on queueNotFull." << endl;
        }
#       endif
        waitForCondition(poolInfo_->queueNotFull, poolInfo_->queueLock);
    }

    // Is the pool in the process of being destroyed?
    // Unlock the mutex and return to caller.
    if (poolInfo_->shutDown || poolInfo_->queueClosed)
    {
        poolInfo_->queueLock.unlock();
        return;
    }

    // Allocate a new work structure
    workQueueItem *newWorkItem = new workQueueItem;
    newWorkItem->function = tFunction;
    newWorkItem->arg = arg;
    newWorkItem->next = NULL;

    // Add new work structure to the queue
    if (poolInfo_->queueSize == 0)
    {
        poolInfo_->tail = poolInfo_->head = newWorkItem;
        broadCast(poolInfo_->queueNotEmpty);
    }
    else
    {
        poolInfo_->tail->next = newWorkItem;
        poolInfo_->tail = newWorkItem;
    }

    poolInfo_->queueSize++;

    // Unlock the work queue
    poolInfo_->queueLock.unlock();
}


void multiThreader::destroyThreadPool()
{
    // Destroy threads only if multi-threaded
    if (multiThreaded())
    {
        // Lock the work queue
        poolInfo_->queueLock.lock();

        // Is a shutdown already in progress?
        if (poolInfo_->queueClosed || poolInfo_->shutDown)
        {
            // Unlock the mutex and return
            poolInfo_->queueLock.unlock();
            return;
        }

        poolInfo_->queueClosed = true;

        // Wait for workers to drain the queue
        while (poolInfo_->queueSize != 0)
        {
            waitForCondition(poolInfo_->queueEmpty, poolInfo_->queueLock);
        }

        poolInfo_->shutDown = true;

        // Unlock the work queue
        poolInfo_->queueLock.unlock();

        // Wake up workers so that they check the shutdown flag
        broadCast(poolInfo_->queueNotEmpty);
        broadCast(poolInfo_->queueNotFull);

        // Wait for all workers to exit
        for(int i=0; i < numThreads_; i++)
        {
            if (pthread_join(poolInfo_->threads[i],NULL))
            {
                FatalErrorIn("multiThreader::destroyThreadPool()")
                    << "pthread_join failed."
                    << abort(FatalError);
            }
        }

        // Destroy the attribute
        pthread_attr_destroy(&(poolInfo_->attr));

        // Deallocate the work-queue and pool structure
        delete [] poolInfo_->threads;

        workQueueItem *currentNode;
        while(poolInfo_->head != NULL)
        {
            currentNode = poolInfo_->head->next;
            poolInfo_->head = poolInfo_->head->next;
            delete currentNode;
        }

        delete poolInfo_;
    }
}


void multiThreader::waitForCondition
(
    Conditional& condition,
    Mutex& mutex
) const
{
    if (pthread_cond_wait(condition(),mutex()))
    {
        FatalErrorIn("multiThreader::waitForCondition()")
            << "Conditional wait failed."
            << abort(FatalError);
    }
}


void multiThreader::broadCast(Conditional& condition) const
{
    if (pthread_cond_broadcast(condition()))
    {
        FatalErrorIn("multiThreader::broadCast()")
            << "Unable to broadcast."
            << abort(FatalError);
    }
}


void multiThreader::signal(Conditional& condition) const
{
    if (pthread_cond_signal(condition()))
    {
        FatalErrorIn("multiThreader::signal()")
            << "Unable to signal."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the number of threads
int multiThreader::getNumThreads() const
{
    return numThreads_;
}


//- Obtain the thread ID for a given index
pthread_t multiThreader::getID(int index) const
{
    if (multiThreaded())
    {
        if (poolInfo_ && index > -1 && index < numThreads_)
        {
            return poolInfo_->threads[index];
        }
        else
        {
            FatalErrorIn("multiThreader::getID(int index)")
                << "Invalid request for ID."
                << abort(FatalError);
        }
    }

    return pthread_self();
}


//- Return true if the number of threads is equal to one.
bool multiThreader::singleThreaded() const
{
    return (numThreads_ == 1);
}


//- Return true if the number of threads is more than one.
bool multiThreader::multiThreaded() const
{
    return (numThreads_ > 1);
}


//- Return the maxQueueSize
int multiThreader::getMaxQueueSize() const
{
    return maxQueueSize_;
}


//- Set the maxQueueSize
void multiThreader::setMaxQueueSize(int size) const
{
    if (size > 0)
    {
        maxQueueSize_ = size;
    }
    else
    {
        FatalErrorIn("multiThreader::setMaxQueueSize(int size)")
            << "Improper value for MaxQueueSize."
            << abort(FatalError);
    }
}


void Mutex::lock() const
{
    if (pthread_mutex_lock(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::lock()")
            << "Unable to lock mutex."
            << abort(FatalError);
    }
}


bool Mutex::tryLock() const
{
    label retVal;

    if ((retVal = pthread_mutex_trylock(&lock_)) != 0)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            if (retVal == EINVAL)
            {
                FatalErrorIn("multiThreader::Mutex::trylock()")
                    << "Mutex returned EINVAL."
                    << abort(FatalError);
            }
            if (retVal == EFAULT)
            {
                FatalErrorIn("multiThreader::Mutex::trylock()")
                    << "Mutex returned EFAULT."
                    << abort(FatalError);
            }
        }
#       endif
    }

    return retVal;
}


void Mutex::unlock() const
{
    if (pthread_mutex_unlock(&lock_))
    {
        FatalErrorIn("multiThreader::Mutex::unlock()")
            << "Unable to unlock the mutex."
            << abort(FatalError);
    }
}


void rwMutex::lock(const lockType lType) const
{
    if (lType == READ_LOCK)
    {
        if (pthread_rwlock_rdlock(&lock_))
        {
            FatalErrorIn("multiThreader::rwMutex::lock()")
                << "Unable to read lock the mutex."
                << abort(FatalError);
        }
    }
    else
    if (lType == WRITE_LOCK)
    {
        if (pthread_rwlock_wrlock(&lock_))
        {
            FatalErrorIn("multiThreader::rwMutex::lock()")
                << "Unable to write lock the mutex."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("multiThreader::rwMutex::lock()")
            << "Undefined mutex type."
            << abort(FatalError);
    }
}


bool rwMutex::tryLock(const lockType lType) const
{
    label retVal = -1;

    if (lType == READ_LOCK)
    {
        if ((retVal = pthread_rwlock_tryrdlock(&lock_)) != 0)
        {
            if (retVal == EINVAL)
            {
                FatalErrorIn("multiThreader::rwMutex::tryLock()")
                    << "Read mutex returned EINVAL."
                    << abort(FatalError);
            }

            if (retVal == EFAULT)
            {
                FatalErrorIn("multiThreader::rwMutex::tryLock()")
                    << "Read mutex returned EFAULT."
                    << abort(FatalError);
            }
        }
    }
    else
    if (lType == WRITE_LOCK)
    {
        if ((retVal = pthread_rwlock_trywrlock(&lock_)) != 0)
        {
            if (retVal == EINVAL)
            {
                FatalErrorIn("multiThreader::rwMutex::tryLock()")
                    << "Write mutex returned EINVAL."
                    << abort(FatalError);
            }

            if (retVal == EFAULT)
            {
                FatalErrorIn("multiThreader::rwMutex::tryLock()")
                    << "Write mutex returned EFAULT."
                    << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn("multiThreader::rwMutex::tryLock()")
            << "Undefined mutex type."
            << abort(FatalError);
    }

    return retVal;
}


void rwMutex::unlock() const
{
    if (pthread_rwlock_unlock(&lock_))
    {
        FatalErrorIn("multiThreader::rwMutex::unlock()")
            << "Unable to unlock the read-write mutex."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void multiThreader::operator=(const multiThreader& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("multiThreader::operator=(const multiThreader&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
