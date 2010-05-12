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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Distribute list.
template<class T>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field
)
{
    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                OPstream toNbr(Pstream::blocking, domain);
                toNbr << IndirectList<T>(field, subMap[domain])();
            }
        }

        // Subset myself
        List<T> subField(IndirectList<T>(field, subMap[Pstream::myProcNo()]));

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        forAll(map, i)
        {
            field[map[i]] = subField[i];
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                IPstream fromNbr(Pstream::blocking, domain);
                List<T> subField(fromNbr);

                const labelList& map = constructMap[domain];

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize);

        // Subset myself
        List<T> subField(IndirectList<T>(field, subMap[Pstream::myProcNo()]));

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            newField[map[i]] = subField[i];
        }

        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am sender. Send to recvProc.
                OPstream toNbr(Pstream::scheduled, recvProc);
                toNbr << IndirectList<T>(field, subMap[recvProc])();
            }
            else
            {
                // I am receiver. Receive from sendProc.
                IPstream fromNbr(Pstream::scheduled, sendProc);
                List<T> subField(fromNbr);

                const labelList& map = constructMap[sendProc];

                forAll(map, i)
                {
                    newField[map[i]] = subField[i];
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        List<T> newField(constructSize);

        // Subset myself
        List<T> subField(IndirectList<T>(field, subMap[Pstream::myProcNo()]));

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            newField[map[i]] = subField[i];
        }

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                OPstream toNbr(Pstream::nonBlocking, domain);
                toNbr << IndirectList<T>(field, subMap[domain])();
            }
        }


        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            if (domain != Pstream::myProcNo())
            {
                IPstream fromNbr(Pstream::nonBlocking, domain);
                List<T> subField(fromNbr);

                const labelList& map = constructMap[domain];

                forAll(map, i)
                {
                    newField[map[i]] = subField[i];
                }
            }
        }
        OPstream::waitRequests();
        IPstream::waitRequests();

        field.transfer(newField);
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


// ************************************************************************* //
