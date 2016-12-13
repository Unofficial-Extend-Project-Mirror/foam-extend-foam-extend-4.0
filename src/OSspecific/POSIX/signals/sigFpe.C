/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include <stdint.h>
#include "error.H"
#include "sigFpe.H"

#include "JobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#ifdef LINUX_GNUC

#   ifndef __USE_GNU
#       define __USE_GNU
#   endif

#   include <fenv.h>
#   include <malloc.h>

#elif defined(sgiN32) || defined(sgiN32Gcc)

#   include <sigfpe.h>

#elif defined(__APPLE__)

// #   include <fenv.h>
#include <xmmintrin.h>
#include <mach/mach.h>

#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigFpe::oldAction_;


#if defined(LINUX)

void *(*Foam::sigFpe::old_malloc_hook)(size_t, const void *) = NULL;

void* Foam::sigFpe::my_malloc_hook(size_t size, const void *caller)
{
    void *result;

    // Restore all old hooks
    __malloc_hook = old_malloc_hook;

    // Call recursively
    result = malloc (size);

    // initialize to signalling nan
#   ifdef WM_SP

    const uint32_t sNAN = 0x7ff7fffflu;

    int nScalars = size / sizeof(scalar);

    uint32_t* dPtr = reinterpret_cast<uint32_t*>(result);

    for (int i = 0; i < nScalars; i++)
    {
        *dPtr++ = sNAN;
    }

#   else

    const uint64_t sNAN = 0x7ff7ffffffffffffllu;

    int nScalars = size/sizeof(scalar);

    uint64_t* dPtr = reinterpret_cast<uint64_t*>(result);

    for (int i = 0; i < nScalars; i++)
    {
        *dPtr++ = sNAN;
    }

#   endif

    // Restore our own hooks
    __malloc_hook = my_malloc_hook;

    return result;
}

#elif defined(__APPLE__)

void *(*Foam::sigFpe::system_malloc_)(malloc_zone_t *zone, size_t size)=NULL;

void* Foam::sigFpe::nan_malloc_(malloc_zone_t *zone, size_t size)
{
    void *result=system_malloc_(zone,size);

    // initialize to signalling NaN
#   ifdef WM_SP

    const uint32_t sNAN = 0x7ff7fffflu;
    uint32_t* dPtr = reinterpret_cast<uint32_t*>(result);

#   else

    const uint64_t sNAN = 0x7ff7ffffffffffffllu;
    uint64_t* dPtr = reinterpret_cast<uint64_t*>(result);

#   endif

    const size_t nScalars = size/sizeof(scalar);
    for (size_t i = 0; i < nScalars; ++i)
    {
        *dPtr++ = sNAN;
    }

    return result;
}

#endif


#if defined(LINUX_GNUC) || defined(__APPLE__)

void Foam::sigFpe::sigFpeHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, NULL) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigSegv::sigFpeHandler()"
        )   << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}

#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    oldAction_.sa_handler = NULL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    if (env("FOAM_SIGFPE"))
    {
#       ifdef LINUX_GNUC

        // Reset signal
        if (oldAction_.sa_handler && sigaction(SIGFPE, &oldAction_, NULL) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigFpe::~sigFpe()"
            )   << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }

#       endif
    }

    if (env("FOAM_SETNAN"))
    {
#       ifdef LINUX_GNUC

        // Reset to standard malloc
        if (oldAction_.sa_handler)
        {
            __malloc_hook = old_malloc_hook;
        }

            #       elif defined(__APPLE__)

        if(system_malloc_!=NULL) {
            malloc_zone_t *zone = malloc_default_zone();
            if(zone==NULL) {
                FatalErrorIn("Foam__sigFpe::set")
                    << "Could not get malloc_default_zone()." << endl
                        << "Seems like this version of Mac OS X doesn't support FOAM_SETNAN"
                        << endl
                        << exit(FatalError);

            }

            if(zone->version>=8)
            {
                vm_protect(
                    mach_task_self(),
                    (uintptr_t)zone,
                    sizeof(malloc_zone_t),
                    0,
                    VM_PROT_READ | VM_PROT_WRITE
                );//remove the write protection
            }
            zone->malloc=system_malloc_;
            system_malloc_=NULL;
            if(zone->version==8)
            {
                vm_protect(
                    mach_task_self(),
                    (uintptr_t)zone,
                    sizeof(malloc_zone_t),
                    0,
                    VM_PROT_READ
                );//put the write protection back
            }

        }

#       endif
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::set(const bool verbose)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorIn
        (
            "Foam::sigFpe::set()"
        )   << "Cannot call sigFpe::set() more than once"
            << abort(FatalError);
    }

    if (env("FOAM_SIGFPE"))
    {
        if (verbose)
        {
            Info<< "SigFpe   : Enabling floating point exception trapping"
                << " (FOAM_SIGFPE)." << endl;
        }

#       ifdef LINUX_GNUC

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        struct sigaction newAction;
        newAction.sa_handler = sigFpeHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigFpe::set()"
            )   << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }


#       elif defined(sgiN32) || defined(sgiN32Gcc)

        sigfpe_[_DIVZERO].abort=1;
        sigfpe_[_OVERFL].abort=1;
        sigfpe_[_INVALID].abort=1;

        sigfpe_[_DIVZERO].trace=1;
        sigfpe_[_OVERFL].trace=1;
        sigfpe_[_INVALID].trace=1;

        handle_sigfpes
        (
            _ON,
            _EN_DIVZERO
          | _EN_INVALID
          | _EN_OVERFL,
            0,
            _ABORT_ON_ERROR,
            NULL
        );

#       elif defined(__APPLE__)

        struct sigaction newAction;
        newAction.sa_handler = sigFpeHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigFpe::set()"
            )   << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_DIV_ZERO);

        _mm_setcsr( _MM_MASK_MASK &~
        (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );

#       endif
    }


    if (env("FOAM_SETNAN"))
    {
        if (verbose)
        {
            Info<< "SetNaN   : Initialising allocated memory to NaN"
                << " (FOAM_SETNAN)." << endl;
        }

#       ifdef LINUX_GNUC

        // Set our malloc
        __malloc_hook = Foam::sigFpe::my_malloc_hook;

#elif defined(__APPLE__)

        if(system_malloc_!=NULL) {
            FatalErrorIn("Foam__sigFpe::set")
                << "system_malloc_ already reset." << endl
                    << "This should never happen"
                    << endl
                    << exit(FatalError);
        }

        malloc_zone_t *zone = malloc_default_zone();
        if(zone==NULL) {
            FatalErrorIn("Foam__sigFpe::set")
                << "Could not get malloc_default_zone()." << endl
                    << "Seems like this version of Mac OS X doesn't support FOAM_SETNAN"
                    << endl
                    << exit(FatalError);
        }
        // According to http://bkdc.ubiquity.ro/2011/07/how-to-set-malloc-hooks-in-osx-lion-107.html
        if(zone->version>=8)
        {
            vm_protect(
                mach_task_self(),
                (uintptr_t)zone,
                sizeof(malloc_zone_t),
                0,
                VM_PROT_READ | VM_PROT_WRITE
            );//remove the write protection
        }
        system_malloc_=zone->malloc;
        zone->malloc=Foam::sigFpe::nan_malloc_;
        if(zone->version==8)
        {
            vm_protect(
                mach_task_self(),
                (uintptr_t)zone,
                sizeof(malloc_zone_t),
                0,
                VM_PROT_READ
            );//put the write protection back
        }

#       endif
    }
}


// ************************************************************************* //
