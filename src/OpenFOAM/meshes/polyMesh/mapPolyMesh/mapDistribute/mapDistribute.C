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

#include "mapDistribute.H"
#include "commSchedule.H"
#include "HashSet.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapDistribute::calcSchedule() const
{
    // Communications: send and receive processor
    List<labelPair> allComms;

    {
        HashSet<labelPair, labelPair::Hash<> > commsSet(Pstream::nProcs());

        // Find what communication is required
        forAll(subMap_, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                if (subMap_[procI].size() > 0)
                {
                    // I need to send to procI
                    commsSet.insert(labelPair(Pstream::myProcNo(), procI));
                }
                if (constructMap_[procI].size() > 0)
                {
                    // I need to receive from procI
                    commsSet.insert(labelPair(procI, Pstream::myProcNo()));
                }
            }
        }
        allComms = commsSet.toc();
    }


    // Reduce
    if (Pstream::master())
    {
        // Receive and merge
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(Pstream::blocking, slave);
            List<labelPair> nbrData(fromSlave);

            forAll(nbrData, i)
            {
                if (findIndex(allComms, nbrData[i]) == -1)
                {
                    label sz = allComms.size();
                    allComms.setSize(sz+1);
                    allComms[sz] = nbrData[i];
                }
            }
        }
        // Send back
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(Pstream::blocking, slave);
            toSlave << allComms;
        }
    }
    else
    {
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());
            toMaster << allComms;
        }
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
            fromMaster >> allComms;
        }
    }


    // Determine my schedule.
    labelList mySchedule
    (
        commSchedule
        (
            Pstream::nProcs(),
            allComms
        ).procSchedule()[Pstream::myProcNo()]
    );

    // Processors involved in my schedule
    schedulePtr_.reset
    (
        new List<labelPair>
        (
            IndirectList<labelPair>(allComms, mySchedule)
        )
    );


    //if (debug)
    //{
    //    Pout<< "I need to:" << endl;
    //    const List<labelPair>& comms = schedule();
    //    forAll(comms, i)
    //    {
    //        const labelPair& twoProcs = comms[i];
    //        label sendProc = twoProcs[0];
    //        label recvProc = twoProcs[1];
    //
    //        if (recvProc == Pstream::myProcNo())
    //        {
    //            Pout<< "    receive from " << sendProc << endl;
    //        }
    //        else
    //        {
    //            Pout<< "    send to " << recvProc << endl;
    //        }
    //    }
    //}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap
)
:
    constructSize_(constructSize),
    subMap_(subMap),
    constructMap_(constructMap),
    schedulePtr_()
{}


//- (optionally destructively) construct from components
Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    labelListList& subMap,
    labelListList& constructMap,
    const bool reUse                // clone or reuse
)
:
    constructSize_(constructSize),
    subMap_(subMap, reUse),
    constructMap_(constructMap, reUse),
    schedulePtr_()
{}


Foam::mapDistribute::mapDistribute
(
    const labelList& sendProcs,
    const labelList& recvProcs
)
:
    constructSize_(sendProcs.size()),
    schedulePtr_()
{
    if (sendProcs.size() != recvProcs.size())
    {
        FatalErrorIn
        (
            "mapDistribute::mapDistribute(const labelList&, const labelList&)"
        )   << "The send and receive data is not the same length. sendProcs:"
            << sendProcs.size() << " recvProcs:" << recvProcs.size()
            << abort(FatalError);
    }

    // Per processor the number of samples we have to send/receive.
    labelList nSend(Pstream::nProcs(), 0);
    labelList nRecv(Pstream::nProcs(), 0);

    forAll(sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        // Note that also need to include local communication (both
        // RecvProc and sendProc on local processor)

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Count destination processor.
            nSend[recvProc]++;
        }
        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            nRecv[sendProc]++;
        }
    }

    subMap_.setSize(Pstream::nProcs());
    constructMap_.setSize(Pstream::nProcs());
    forAll(nSend, procI)
    {
        subMap_[procI].setSize(nSend[procI]);
        constructMap_[procI].setSize(nRecv[procI]);
    }
    nSend = 0;
    nRecv = 0;

    forAll(sendProcs, sampleI)
    {
        label sendProc = sendProcs[sampleI];
        label recvProc = recvProcs[sampleI];

        if (Pstream::myProcNo() == sendProc)
        {
            // I am the sender. Store index I need to send.
            subMap_[recvProc][nSend[recvProc]++] = sampleI;
        }
        if (Pstream::myProcNo() == recvProc)
        {
            // I am the receiver.
            constructMap_[sendProc][nRecv[sendProc]++] = sampleI;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
