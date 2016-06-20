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

Description

\*---------------------------------------------------------------------------*/

#include <unistd.h>

#include "error.H"
#include "timer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::timer, 0);

jmp_buf Foam::timer::envAlarm;

struct sigaction Foam::timer::oldAction_;

unsigned int Foam::timer::oldTimeOut_ = 0;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::timer::signalHandler(int)
{
    if (debug)
    {
        Info<< "Foam::timer::signalHandler(int sig) : "
            << " timed out. Jumping."
            << endl;
    }
    longjmp(envAlarm, 1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
Foam::timer::timer(const unsigned int newTimeOut)
:
    newTimeOut_(newTimeOut)
{

    if (newTimeOut > 0)
    {
        // Is singleton since handler is static function
        if (oldTimeOut_ != 0)
        {
            FatalErrorIn
            (
                "Foam::timer::timer(const unsigned int)"
            )   << "timer already used."
                << abort(FatalError);
        }

        // Install alarm signal handler:
        // - do not block any signals while in it
        // - clear list of signals to mask
        struct sigaction newAction;
        newAction.sa_handler = timer::signalHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);

        if (sigaction(SIGALRM, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::timer::timer(const unsigned int)"
            )   << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }

        oldTimeOut_ = ::alarm(newTimeOut);

        if (debug)
        {
            Info<< "Foam::timer::timer(const unsigned int) : "
                << " installing timeout " << int(newTimeOut_)
                << " seconds"
                << " (overriding old timeout " << int(oldTimeOut_)
                << ")." << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timer::~timer()
{
    if (newTimeOut_ > 0)
    {
        if (debug)
        {
            Info<< "Foam::timer::~timer(const unsigned int) : timeOut="
                << int(newTimeOut_)
                << " : resetting timeOut to " << int(oldTimeOut_) << endl;
        }

        // Reset timer
        ::alarm(oldTimeOut_);
        oldTimeOut_ = 0;

        // Restore signal handler
        if (sigaction(SIGALRM, &oldAction_, NULL) < 0)
        {
            FatalErrorIn
            (
                "Foam::timer::~timer(const struct sigaction&"
                "const struct sigaction&)"
            )   << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }
    }
}

// ************************************************************************* //
