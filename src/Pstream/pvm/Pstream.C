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

#include "Pstream.H"
#include "PstreamReduceOps.H"

#include <cstring>
#include <cstdlib>
#include <csignal>

#include <pvm3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Pstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
}


bool Pstream::init(int& argc, char**& argv)
{
    // Set the comunications options
    pvm_setopt(PvmRoute, PvmRouteDirect);

    // Get the ID of this processor
    int mytid = pvm_mytid();

#ifdef USECRAYSHMEM

    // Get the size of the NULL group
    procIDs_.setSize(pvm_gsize(NULL));

    // For each processor of the NULL group get its ID
    for (int proci=0; proci<ProcIDs.size(); proci++)
    {
        procIDs_[proci] = pvm_gettid(NULL, proci);
    }

#else

    // Initialisation message type
    int initMsgType = 0;

    // If this is not a slave then it must be the master.
    // Master spawns the rest of the child processes in the same manner as MPI
    if (string(argv[argc-1]) != "-slave")
    {
        // Last argument is number of processors in parallel run
        int nProcs = atoi(argv[argc-1]);

        // If it is less than 2 this is not a parallel run!
        if (nProcs < 2)
        {
            FatalErrorIn("Pstream::init(int& argc, char**& argv)")
                << "Attempt to run parallel on < 2 processors ... stopping."
                << abort(FatalError);
        }


        Info<< "Starting parallel run on " << nProcs << " processors ... "
            << nl << endl;


        // set size of ID list
        procIDs_.setSize(nProcs);
        procIDs_ = 0;

        // I am the master
        myProcNo_ = 1;

        // Put my ID in the list
        procIDs_[0] = mytid;

        // Setup arguments of children
        typedef char* charPtr;
        char** Argv = new charPtr[argc + 1];

        for (int i=0; i<argc-1; i++)
        {
            Argv[i] = new char[strlen(argv[i+1] + 1)];
            strcpy(Argv[i], argv[i+1]);
        }

        Argv[argc-1] = new char[7];
        strcpy(Argv[argc-1], "-slave");

        Argv[argc] = NULL;

        // Spawn children as copies of me
        if
        (
            pvm_spawn
            (
                argv[0],
                Argv,
                PvmTaskDefault,
                "",
                nProcs-1,
                &(procIDs_[1])
            ) != nProcs-1
        )
        {
            FatalErrorIn("Pstream::init(int& argc, char**& argv)")
                << "Unable to spawn processes ... stopping."
                << abort(FatalError);
        }


        // Broadcast task IDs to all children
        pvm_setopt(PvmRoute, PvmRouteDirect);
        pvm_initsend(PvmDataDefault);
        pvm_pkint((int*)(&nProcs), 1, 1);
        pvm_pkint(procIDs_.begin(), nProcs, 1);
        pvm_mcast(procIDs_.begin(), nProcs, initMsgType);


        Info<< "nProcs : " << nProcs << endl;
        Info<< "TIDS   : ";
        for (int proci=0; proci<procIDs_.size(); proci++)
        {
            cout<< hex << procIDs_[proci] << ' ';
        }
        cout<< dec << nl << std::endl;
    }
    else
    {
        // Receive processor data from master
        pvm_recv(-1, initMsgType);

        // Should have received the number of processors in the run
        int nProcs;
        pvm_upkint(&nProcs, 1, 1);

        // ... set size of ID list
        procIDs_.setSize(nProcs);

        // ... and unpack the processor IDs
        pvm_upkint(procIDs_.begin(), nProcs, 1);
    }

#endif

    // Find which processor number this is
    for (int proci=0; proci<procIDs_.size(); proci++)
    {
        if (procIDs_[proci] == mytid)
        {
            break;
        }
    }

    // Set the processor numbers to start from 1
    myProcNo_ = proci + 1;

    /*
    if (pvm_joingroup("foam") < 0)
    {
        FatalErrorIn("Pstream::init(int& argc, char**& argv)")
            << "Pstream::init(int*, char **[]) : "
            << "could not join group ... stopping."
            << abort(FatalError);
    }

    pvm_barrier("foam", nProcs());
    */

    // Setup signal handler to catch an interupt (^C) and abort the run
    // This doesn't work, it causes
    // libpvm [t40003]: pvm_sendsig(): Not implemented
    // libpvm [t40003]: pvm_kill(): Not implemented
    // messages
    //signal(SIGINT, stop);

    if (master())
    {
        Sout<< "Master started successfully." << nl << endl;
    }
    else
    {
        Sout<< "Child " << myProcNo_ << " started successfully." << nl << endl;
    }

    setParRun();

    // Everything is OK
    return true;
}


void Pstream::exit(int errnum)
{
    //pvm_lvgroup("foam");

    if (errnum != 0)
    {
        for (int proci=1; proci<=procIDs_.size(); proci++)
        {
            if (proci != myProcNo())
            {
                pvm_kill(procID(proci));
            }
        }
    }

    pvm_exit();
    ::exit(errnum);
}


void Pstream::abort()
{
    for (int proci=1; proci<=procIDs_.size(); proci++)
    {
        if (proci != myProcNo())
        {
            pvm_kill(procID(proci));
        }
    }

    pvm_exit();
    //::abort();
}


void reduce(scalar& Value, const sumOp<scalar>& bop)
{
    if (Pstream::parRun())
    {
#       ifdef PVM_REDUCE
        if
        (
            pvm_reduce
            (
                PvmSum,
                &Value,
                1,
                PVM_DOUBLE,
                Pstream::msgType(),
                "foam",
                0
            ) != PvmOk
        )
        {
            FatalErrorIn
            (
                "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
            )   << "pvm_reduce failed"
                << abort(FatalError);
        }
#       endif

        if (Pstream::master())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                scalar value;
                int atid, atag, alen;

                if
                (
                    pvm_precv
                    (
                        Pstream::procID(slave),
                        Pstream::msgType(),
                        &value,
                        1,
                        PVM_DOUBLE,
                        &atid, &atag, &alen
                    ) != PvmOk
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "pvm_precv failed"
                        << abort(FatalError);
                }

                Value = bop(Value, value);
            }
        }
        else
        {
            if
            (
                pvm_psend
                (
                    Pstream::procID(Pstream::masterNo()),
                    Pstream::msgType(),
                    &Value,
                    1,
                    PVM_DOUBLE
                ) != PvmOk
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "pvm_psend failed"
                    << abort(FatalError);
            }
        }


        if (Pstream::master())
        {
            pvm_initsend(PvmDataDefault);
            pvm_pkdouble(&Value, 1, 1);

            if
            (
                pvm_mcast
                (
                    (int*)Pstream::procIDs().begin(),
                    Pstream::nProcs(),
                    Pstream::msgType()
                ) != PvmOk
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "pvm_mcast failed"
                    << abort(FatalError);
            }
        }
        else
        {
            if
            (
                pvm_recv
                (
                    Pstream::procID(Pstream::masterNo()),
                    Pstream::msgType()
                ) <= 0
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "pvm_psend failed"
                    << abort(FatalError);
            }

            pvm_upkdouble(&Value, 1, 1);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
