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

Class
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Controlled by two env vars:
    @param FOAM_SIGFPE \n
        exception trapping
    @param FOAM_SETNAN \n
        initialization of all malloced memory to NaN. If FOAM_SIGFPE
        also set, this will cause usage of uninitialized scalars to trigger
        an abort.

SourceFiles
    sigFpe.C

\*---------------------------------------------------------------------------*/

#ifndef sigFpe_H
#define sigFpe_H

#include "OSspecific.H"
#include <signal.h>

#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64) || defined(linuxARM7)
#    define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
#    define LINUX_GNUC
#endif

#ifdef __APPLE__
#include <malloc/malloc.h>

// these are defined by the mach-headers and would break compilation of Switch.H
#undef TRUE
#undef FALSE

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private data

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

#       ifdef LINUX

            //- Saved old malloc
            static void *(*old_malloc_hook)(size_t, const void *);

            //- nan malloc function. From malloc_hook manpage.
            static void* my_malloc_hook(size_t size, const void *caller);

#elif defined (__APPLE__)

    //- pointer to the original malloc that is overrided
    static void *(*system_malloc_)(malloc_zone_t *zone, size_t size);

    //- the overriding handler
    static void* nan_malloc_(malloc_zone_t *zone, size_t size);

#       endif


    // Static data members

#   if defined(LINUX_GNUC) || defined(__APPLE__)

        //- Handler for caught signals
        static void sigFpeHandler(int);

#   endif


public:


    // Constructors

        sigFpe();


    // Destructor

        ~sigFpe();


    // Member functions

        void set(const bool verbose);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
