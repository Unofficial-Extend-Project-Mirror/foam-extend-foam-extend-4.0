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

#include "scalar.H"
#include "MSwindows.H"

#define WINVER 0x0500 // To access CreateTimerQueueTimer
#include <windows.h>

#define SIGALRM 14

#include "error.H"
#include "timer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::timer, 0);

jmp_buf Foam::timer::envAlarm;

__p_sig_fn_t Foam::timer::oldAction_ = SIG_DFL;

static HANDLE hTimer_ = NULL;

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


static VOID CALLBACK timerExpired(PVOID lpParam, BOOLEAN TimerOrWaitFired)
{
    ::raise(SIGALRM);
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
        if( NULL != hTimer_ )
        {
            FatalErrorIn
            (
                "Foam::timer::timer(const unsigned int)"
            )   << "timer already used."
                << abort(FatalError);
        }

        // Install alarm signal handler:
        oldAction_ = ::signal(SIGALRM, &Foam::timer::signalHandler);

        if( SIG_ERR == oldAction_ )
        {
            oldAction_ = SIG_DFL;

            FatalErrorIn
            (
                "Foam::timer::timer(const unsigned int)"
            )   << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }

        if (debug)
        {
            Info<< "Foam::timer::timer(const unsigned int) : "
                << " installing timeout " << int(newTimeOut_)
                << " seconds." << endl;
        }

        const bool success =
          ::CreateTimerQueueTimer(&hTimer_,
                                  NULL,
                                  (WAITORTIMERCALLBACK)timerExpired,
                                  NULL ,
                                  newTimeOut * 1000,
                                  0, 0);

        if (!success)
        {
            hTimer_ = NULL;
            FatalErrorIn
            (
                "Foam::timer::timer(const unsigned int)"
            )   << "CreateTimerQueueTimer, "
                << MSwindows::getLastError()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timer::~timer()
{
    if (newTimeOut_ > 0)
    {
        // Reset timer
        const bool timerSuccess =
          ::DeleteTimerQueueTimer(NULL, hTimer_, NULL);
        hTimer_ = NULL;

        if (!timerSuccess)
        {
            FatalErrorIn
            (
                "Foam::timer::~timer() "
            )   << "DeleteTimerQueueTimer, "
                << MSwindows::getLastError()
                << abort(FatalError);
        }

        if (debug)
        {
            Info<< "Foam::timer::~timer() timeOut="
                << int(newTimeOut_) << endl;
        }

        const __p_sig_fn_t signalSuccess = signal(SIGALRM, oldAction_);
        oldAction_ = SIG_DFL;

        // Restore signal handler
        if (SIG_ERR == signalSuccess)
        {
            FatalErrorIn
            (
                "Foam::timer::~timer()"
            )   << "sigaction(SIGALRM) error"
                << abort(FatalError);
        }
    }
}

// ************************************************************************* //
