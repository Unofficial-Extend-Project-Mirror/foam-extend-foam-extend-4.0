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
    List<T>& field,
    const int tag
)
{
    if (!Pstream::parRun())
    {
        // Do only me to me.

        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        forAll(map, i)
        {
            field[map[i]] = subField[i];
        }

        return;
    }

    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain, 0, tag);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

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
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

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
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            newField[map[i]] = subField[i];
        }

        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr(Pstream::scheduled, recvProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[recvProc]);
                }
                {
                    IPstream fromNbr(Pstream::scheduled, recvProc, 0, tag);
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        newField[map[i]] = subField[i];
                    }
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr(Pstream::scheduled, sendProc, 0, tag);
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        newField[map[i]] = subField[i];
                    }
                }
                {
                    OPstream toNbr(Pstream::scheduled, sendProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[sendProc]);
                }
            }
        }

        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    OPstream toDomain(Pstream::nonBlocking, domain, 0, tag);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving. Do not block.

            {
                // Set up 'send' to myself
                const labelList& mySubMap = subMap[Pstream::myProcNo()];
                List<T> mySubField(mySubMap.size());
                forAll(mySubMap, i)
                {
                    mySubField[i] = field[mySubMap[i]];
                }
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        field[map[i]] = mySubField[i];
                    }
                }
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    IPstream str(Pstream::nonBlocking, domain, 0, tag);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    forAll(map, i)
                    {
                        field[map[i]] = recvField[i];
                    }
                }
            }
        }
        else // contiguous data
        {
            // Set up sends to neighbours

            List<List<T > > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.byteSize(),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    IPstream::read
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].byteSize(),
                        tag
                    );
                }
            }


            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());

                forAll(map, i)
                {
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);

            // Receive sub field from myself (sendFields[Pstream::myProcNo()])
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
            }


            // Wait for all to finish
            Pstream::waitRequests(nOutstanding);

            // Collect neighbour fields
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    forAll(map, i)
                    {
                        field[map[i]] = subField[i];
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule "
            << Pstream::commsTypeNames[commsType]
            << abort(FatalError);
    }
}


// Distribute list.
template<class T, class CombineOp>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field,
    const CombineOp& cop,
    const T& nullValue,
    const int tag
)
{
    if (!Pstream::parRun())
    {
        // Do only me to me.

        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        forAll(map, i)
        {
            cop(field[map[i]], subField[i]);
        }
        return;
    }

    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain, 0, tag);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        forAll(map, i)
        {
            cop(field[map[i]], subField[i]);
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain, 0, tag);
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize, nullValue);

        // Subset myself
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            cop(newField[map[i]], subField[i]);
        }

        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr(Pstream::scheduled, recvProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[recvProc]);
                }
                {
                    IPstream fromNbr(Pstream::scheduled, recvProc, 0, tag);
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        cop(newField[map[i]], subField[i]);
                    }
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr(Pstream::scheduled, sendProc, 0, tag);
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    forAll(map, i)
                    {
                        cop(newField[map[i]], subField[i]);
                    }
                }
                {
                    OPstream toNbr(Pstream::scheduled, sendProc, 0, tag);
                    toNbr << UIndirectList<T>(field, subMap[sendProc]);
                }
            }
        }

        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        label nOutstanding = Pstream::nRequests();

        if (!contiguous<T>())
        {
            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    OPstream toDomain(Pstream::nonBlocking, domain, 0, tag);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving. Do not block.

            {
                // Set up 'send' to myself
                List<T> mySubField(field, subMap[Pstream::myProcNo()]);
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                field = nullValue;
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        cop(field[map[i]], mySubField[i]);
                    }
                }
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    IPstream str(Pstream::nonBlocking, domain, 0, tag);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    forAll(map, i)
                    {
                        cop(field[map[i]], recvField[i]);
                    }
                }
            }
        }
        else // Contiguous data
        {
            // Set up sends to neighbours

            List<List<T> > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.size()*sizeof(T),
                        tag
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    IPstream::read
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].size()*sizeof(T),
                        tag
                    );
                }
            }

            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);
            field = nullValue;

            // Receive sub field from myself (subField)
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
            }

            // Wait for all to finish
            Pstream::waitRequests(nOutstanding);

            // Collect neighbour fields
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    forAll(map, i)
                    {
                        cop(field[map[i]], subField[i]);
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


// ************************************************************************* //
