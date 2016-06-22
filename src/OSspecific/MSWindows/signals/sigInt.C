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

#include "error.H"
#include "sigInt.H"
#include "JobInfo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

__p_sig_fn_t Foam::sigInt::oldAction_ = SIG_DFL;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigInt::sigIntHandler(int)
{
    // Reset old handling
    const __p_sig_fn_t success = ::signal(SIGINT, oldAction_);

    if( SIG_ERR == success )
    {
        FatalErrorIn
        (
            "Foam::sigInt::sigIntHandler()"
        )   << "Cannot reset SIGINT trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    // Throw signal (to old handler)
    raise(SIGINT);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigInt::sigInt()
{
    oldAction_ = SIG_DFL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigInt::~sigInt()
{
    // Reset old handling
    const __p_sig_fn_t success = ::signal(SIGINT, oldAction_);
    oldAction_ = SIG_DFL;

    if( SIG_ERR == success )
    {
        FatalErrorIn
        (
            "Foam::sigInt::~sigInt()"
        )   << "Cannot reset SIGINT trapping"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigInt::set(const bool verbose)
{
    if( SIG_DFL != oldAction_ )
    {
        FatalErrorIn
        (
            "Foam::sigInt::set()"
        )   << "Cannot call sigInt::set() more than once"
            << abort(FatalError);
    }

    oldAction_ = ::signal(SIGINT, &Foam::sigInt::sigIntHandler);

    if( SIG_ERR == oldAction_ )
    {
        oldAction_ = SIG_DFL;

        FatalErrorIn
        (
            "Foam::sigInt::set()"
        )   << "Cannot set SIGINT trapping"
            << abort(FatalError);
    }
}


// ************************************************************************* //
