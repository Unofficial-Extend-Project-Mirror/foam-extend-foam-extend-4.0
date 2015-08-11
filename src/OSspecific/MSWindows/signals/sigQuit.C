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

#include "error.H"
#include "sigQuit.H"
#include "JobInfo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

__p_sig_fn_t Foam::sigQuit::oldAction_ = SIG_DFL;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigQuit::sigQuitHandler(int)
{
    // Reset old handling
    const __p_sig_fn_t success = ::signal(SIGBREAK, oldAction_);

    if( SIG_ERR == success )
    {
        FatalErrorIn
        (
            "Foam::sigQuit::sigQuitHandler()"
        )   << "Cannot reset SIGQUIT trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    ::raise(SIGBREAK);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigQuit::sigQuit()
{
    oldAction_ = SIG_DFL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigQuit::~sigQuit()
{
    // Reset old handling
    if( SIG_DFL != oldAction_ )
    {
        const __p_sig_fn_t success = ::signal(SIGBREAK, oldAction_);
        oldAction_ = SIG_DFL;

        if( SIG_ERR == success )
        {
            FatalErrorIn("Foam::sigQuit::~sigQuit()")
                << "Cannot reset SIGBREAK trapping"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigQuit::set(const bool verbose)
{
    if( SIG_DFL != oldAction_ )
    {
        FatalErrorIn
        (
            "Foam::sigQuit::set()"
        )   << "Cannot call sigQuit::set() more than once"
            << abort(FatalError);
    }

    oldAction_ = ::signal(SIGBREAK, &Foam::sigQuit::sigQuitHandler);

    if( SIG_ERR == oldAction_ )
    {
        oldAction_ = SIG_DFL;

        FatalErrorIn
        (
            "Foam::sigQuit::set()"
        )   << "Cannot set SIGQUIT trapping"
            << abort(FatalError);
    }
}


// ************************************************************************* //
