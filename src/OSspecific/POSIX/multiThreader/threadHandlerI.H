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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor for threadHandler
template<class T>
threadHandler<T>::threadHandler
(
    T& tRef,
    const multiThreader& threader
)
:
    tRef_(tRef),
    threader_(threader),
    argList_(0),
    nThreads_(threader.getNumThreads()),
    pthreadID_(pthread_self()),
    master_(false),
    predicate_(false)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

template<class T>
threadHandler<T>::~threadHandler()
{
    // Null the argument list pointer copies
    // to avoid multiple deallocations.
    forAll(argList_, indexI)
    {
        argList_[indexI] = NULL;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a reference to the template class
template<class T>
inline T& threadHandler<T>::reference()
{
    return tRef_;
}


// Set a size for the argument list
template<class T>
inline void threadHandler<T>::setSize
(
    const label size
)
{
    argList_.setSize(size, NULL);
}


// Set a argument pointer for a particular index
template<class T>
inline void threadHandler<T>::set
(
    const label index,
    void* argPtr
)
{
    if (!argList_.size())
    {
        FatalErrorIn("threadHandler<T>::set")
            << "Attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else
    if (index < 0 || index >= argList_.size())
    {
        FatalErrorIn("threadHandler<T>::set")
            << "index " << index << " out of range 0 ... "
            << argList_.size()-1
            << abort(FatalError);
    }

    argList_[index] = argPtr;
}


// Return a reference to the multiThreader
template<class T>
inline const multiThreader& threadHandler<T>::threader() const
{
    return threader_;
}


// Return the number of threads
template<class T>
inline label threadHandler<T>::nThreads() const
{
    return nThreads_;
}


// Designate as master thread
template<class T>
inline void threadHandler<T>::setMaster() const
{
    master_ = true;
}


// Designate as slave thread
template<class T>
inline void threadHandler<T>::setSlave() const
{
    master_ = false;
}


// Is this a master thread?
template<class T>
inline bool threadHandler<T>::master() const
{
    return (master_ == true);
}


// Is this a slave thread?
template<class T>
inline bool threadHandler<T>::slave() const
{
    return !master();
}


// Lock this thread
template<class T>
inline void threadHandler<T>::lock
(
    const signalType sType
) const
{
    if (sType == START)
    {
        startMutex_.lock();
    }
    else
    if (sType == STOP)
    {
        stopMutex_.lock();
    }
}


// Unlock this thread
template<class T>
inline void threadHandler<T>::unlock
(
    const signalType sType
) const
{
    if (sType == START)
    {
        startMutex_.unlock();
    }
    else
    if (sType == STOP)
    {
        stopMutex_.unlock();
    }
}


// Send signal to a waiting conditional
template<class T>
inline void threadHandler<T>::sendSignal
(
    const signalType sType
) const
{
    lock(sType);

    if (predicate(sType))
    {
        InfoIn("threadHandler::sendSignal()")
            << "Predicate is already set."
            << endl;
    }
    else
    {
        // Set predicate before signalling
        setPredicate(sType);
    }

    if (sType == START)
    {
        threader().signal
        (
            startConditional_
        );
    }
    else
    if (sType == STOP)
    {
        threader().signal
        (
            stopConditional_
        );
    }
    else
    {
        FatalErrorIn("threadHandler::sendSignal()")
            << "Undefined enumerant."
            << abort(FatalError);
    }

    unlock(sType);
}


// Wait for signal
template<class T>
inline void threadHandler<T>::waitForSignal
(
    const signalType sType
) const
{
    if (sType == START)
    {
        threader().waitForCondition
        (
            startConditional_,
            startMutex_
        );
    }
    else
    if (sType == STOP)
    {
        threader().waitForCondition
        (
            stopConditional_,
            stopMutex_
        );
    }
    else
    {
        FatalErrorIn("threadHandler::waitForSignal()")
            << "Undefined enumerant."
            << abort(FatalError);
    }

    if (!predicate(sType))
    {
        FatalErrorIn("threadHandler::waitForSignal()")
            << "Spurious wake-up."
            << abort(FatalError);
    }

    unsetPredicate(sType);

    // Unlock the acquired mutex
    unlock(sType);
}


// Return state of the predicate variable
template<class T>
inline bool threadHandler<T>::predicate
(
    const signalType sType
) const
{
    return predicate_[sType];
}


// Set the predicate variable
template<class T>
inline void threadHandler<T>::setPredicate
(
    const signalType sType
) const
{
    predicate_[sType] = true;
}


// Unset the predicate variable
template<class T>
inline void threadHandler<T>::unsetPredicate
(
    const signalType sType
) const
{
    predicate_[sType] = false;
}


// Return the ID
template<class T>
inline void threadHandler<T>::setID(const pthread_t& pt)
{
    pthreadID_ = pt;
}


// Return the ID
template<class T>
inline pthread_t threadHandler<T>::ID() const
{
    return pthreadID_;
}


// Does the calling thread correspond to this handler?
template<class T>
inline bool threadHandler<T>::self() const
{
    return pthread_equal(ID(), pthread_self());
}


// Return an argument pointer at a particular index
template<class T>
inline void * threadHandler<T>::operator()
(
    const label index
)
{
    if (!argList_.size())
    {
        FatalErrorIn("threadHandler<T>::operator()")
            << "Attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else
    if (index < 0 || index >= argList_.size())
    {
        FatalErrorIn("threadHandler<T>::operator()")
            << "index " << index << " out of range 0 ... "
            << argList_.size()-1
            << abort(FatalError);
    }
    else
    if (!argList_[index])
    {
        FatalErrorIn("threadHandler<T>::operator()")
            << "Hanging pointer. This is not allowed."
            << abort(FatalError);
    }

    return argList_[index];
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Lock all threads provided by sequence
template <class T>
void lockThreads
(
    const List<label>& sequence,
    const PtrList<threadHandler<T> >& handler
)
{
    forAll(sequence, i)
    {
        handler[sequence[i]].lock(threadHandler<T>::START);
        handler[sequence[i]].lock(threadHandler<T>::STOP);

        handler[sequence[i]].unsetPredicate(threadHandler<T>::START);
        handler[sequence[i]].unsetPredicate(threadHandler<T>::STOP);
    }
}


// Synchronize all threads provided by sequence
template <class T>
void synchronizeThreads
(
    const List<label>& sequence,
    const PtrList<threadHandler<T> >& handler
)
{
    forAll(sequence, i)
    {
        // Wait for a signal from this thread before moving on.
        handler[sequence[i]].waitForSignal(threadHandler<T>::STOP);
    }
}


// Execute threads for the submitted static function by sequence
template <class T>
void executeThreads
(
    const List<label>& sequence,
    PtrList<threadHandler<T> >& handler,
    void (*tFunction)(void*)
)
{
    if (!handler.size())
    {
        FatalErrorIn
        (
            "\n\n"
            "template <class T>\n"
            "void executeThreads\n"
            "(\n"
            "    const List<label>& sequence,\n"
            "    PtrList<threadHandler<T> >& handler,\n"
            "    void (*tFunction)(void*)\n"
            ")\n"
        )
            << "Empty handler list."
            << abort(FatalError);
    }

    // Fetch threader reference
    const multiThreader& threader = handler[0].threader();

    // Lock slave threads by sequence
    lockThreads(sequence, handler);

    forAll(sequence, i)
    {
        // Submit jobs to the work queue
        threader.addToWorkQueue(tFunction, &(handler[sequence[i]]));

        // Wait for a signal from this thread before moving on.
        handler[sequence[i]].waitForSignal(threadHandler<T>::START);
    }

    // Synchronize threads
    synchronizeThreads(sequence, handler);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
