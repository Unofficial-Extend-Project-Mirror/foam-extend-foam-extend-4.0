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
    sigFpe

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
    Copyright            : (C) 2011 Symscape
    Website              : www.symscape.com
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "sigFpe.H"

#include "JobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#include <float.h>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

__p_sig_fn_t Foam::sigFpe::oldAction_ = SIG_DFL;

static unsigned int fpOld_ = 0;


static void clearFpe()
{
    //_clearfp();
    //_controlfp(fpOld_, 0xFFFFFFFF);
}


void Foam::sigFpe::sigFpeHandler(int)
{
    const __p_sig_fn_t success = ::signal(SIGFPE, oldAction_);

    // Reset old handling
    if( SIG_ERR == success )
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

    clearFpe();

    // Throw signal (to old handler)
    ::raise(SIGFPE);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    oldAction_ = SIG_DFL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    if( env("FOAM_SIGFPE") )
    {
        clearFpe();

        // Reset signal
        const __p_sig_fn_t success = ::signal(SIGFPE, oldAction_);
        oldAction_ = SIG_DFL;

        if( SIG_ERR == success )
        {
            FatalErrorIn
            (
                "Foam::sigFpe::~sigFpe()"
            )   << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }
    }

    if( env("FOAM_SETNAN") )
    {
        WarningIn("Foam::sigFpe::~sigFpe()")
            << "FOAM_SETNAN not supported under MSwindows "
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::set(const bool verbose)
{
    if( SIG_DFL != oldAction_ )
    {
        FatalErrorIn
        (
            "Foam::sigFpe::set()"
        )   << "Cannot call sigFpe::set() more than once"
            << abort(FatalError);
    }

    if( env("FOAM_SIGFPE") )
    {
        if( verbose )
        {
            Info<< "SigFpe : Enabling floating point exception trapping"
                << " (FOAM_SIGFPE)." << endl;
        }
/*
        fpOld_ = _controlfp(0, 0);
        const unsigned int fpNew =
          fpOld_ & ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW);
        _controlfp(fpNew, _MCW_EM);

        oldAction_ = ::signal(SIGFPE, &Foam::sigFpe::sigFpeHandler);

        if( SIG_ERR == oldAction_ )
        {
            oldAction_ = SIG_DFL;

            FatalErrorIn
            (
                "Foam::sigFpe::set()"
            )   << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }
*/
    }


    if( env("FOAM_SETNAN") )
    {
        if( verbose )
        {
            WarningIn("Foam::sigFpe::set()")
              << "FOAM_SETNAN not supported under MSwindows "
              << endl;
        }
    }
}


// ************************************************************************* //
