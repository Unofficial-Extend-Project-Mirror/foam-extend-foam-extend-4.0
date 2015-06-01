/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "regionCouplePolyPatch.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::regionCouplePolyPatch::fastExpand
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    // with communication

    // Algorithm:
    // 1) Master processor holds maps of all zone addressing (data provided)
    // and all remote zone addressing (data required)
    // 2) Each processor will send the locally active data to the master
    // 3) Master assembles all the data
    // 4) Master sends to all processors the data they need to receive
    //
    // Notes:
    // A) If the size of zone addressing is zero, data is not sent
    // B) Communicated data on each processor has the size of live faces
    // C) Expanded data will be equal to actual data from other processors
    //    only for the faces marked in remote; for other faces, it will be
    //    equal to zero
    // D) On processor zero, complete data is available
    // HJ, 4/Jun/2011
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > regionCouplePolyPatch::fastExpand\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (localParallel())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > regionCouplePolyPatch::fastExpand"
            "("
            "    const Field<Type>& ff"
            ") const"
        )   << "Requested expand on local parallel.  This is not allowed"
            << abort(FatalError);
    }

    // This function requires send-receive-addressing, but usage is not
    // symmetric across processors. Hence trigger re-calculate at this point
    // HR, 10/Jul/2013
    if (Pstream::parRun() && !localParallel())
    {
        receiveAddr();
        shadow().receiveAddr();
    }

    // Expand the field to zone size
    tmp<Field<Type> > texpandField
    (
        new Field<Type>(zone().size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField();

    if (Pstream::master())
    {
        // Insert master processor
        const labelList& za = zoneAddressing();

        forAll (za, i)
        {
            expandField[za[i]] = ff[i];
        }

        // Master receives and inserts data from all processors for which
        // receiveAddr contains entries
        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            const labelList& curRAddr = receiveAddr()[procI];

            if (!curRAddr.empty())
            {
                Field<Type> receiveBuf(curRAddr.size());

                // Opt: reconsider mode of communication
                IPstream::read
                (
                    Pstream::blocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize()
                );

                // Insert received information
                forAll (curRAddr, i)
                {
                    expandField[curRAddr[i]] = receiveBuf[i];
                }
            }
        }

        // Expanded field complete, send required data to other processors
        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            const labelList& curSAddr = shadow().sendAddr()[procI];

            if (!curSAddr.empty())
            {
                Field<Type> sendBuf(curSAddr.size());

                forAll (curSAddr, i)
                {
                    sendBuf[i] = expandField[curSAddr[i]];
                }

                // Opt: reconsider mode of communication
                OPstream::write
                (
                    Pstream::blocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf.begin()),
                    sendBuf.byteSize()
                );
            }
        }
    }
    else
    {
        // Send local data to master and receive remote data
        // If patch is empty, communication is avoided
        // HJ, 4/Jun/2011
        if (size())
        {
            // Opt: reconsider mode of communication
            OPstream::write
            (
                Pstream::blocking,
                Pstream::masterNo(),
                reinterpret_cast<const char*>(ff.begin()),
                ff.byteSize()
            );
        }

        // Prepare to receive remote data
        const labelList& rza = shadow().remoteZoneAddressing();

        if (!rza.empty())
        {
            Field<Type> receiveBuf(rza.size());

            // Opt: reconsider mode of communication
            IPstream::read
            (
                Pstream::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(receiveBuf.begin()),
                receiveBuf.byteSize()
            );

            // Insert the data into expanded field
            forAll (rza, i)
            {
                expandField[rza[i]] = receiveBuf[i];
            }
        }
    }

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::regionCouplePolyPatch::interpolate
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != shadow().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > regionCouplePolyPatch::interpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect slave patch field size.  Field size: "
            << ff.size() << " patch size: " << shadow().size()
            << abort(FatalError);
    }

    if (localParallel())
    {
        // No expansion or filtering needed.  HJ, 4/Jun/2011

        if (empty())
        {
            // Patch empty, no interpolation
            return tmp<Field<Type> >(new Field<Type>());
        }

        // Interpolate field
        if (master())
        {
            return patchToPatch().slaveToMaster(ff);
        }
        else
        {
            return patchToPatch().masterToSlave(ff);
        }
    }
    else
    {
        // Expand shadow
        Field<Type> expandField = shadow().fastExpand(ff);

        tmp<Field<Type> > tresult(new Field<Type>(size()));
        Field<Type>& result = tresult();

        if (master())
        {
            patchToPatch().maskedSlaveToMaster
            (
                expandField,
                result,
                zoneAddressing()
            );
        }
        else
        {
            patchToPatch().maskedMasterToSlave
            (
                expandField,
                result,
                zoneAddressing()
            );
        }

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::regionCouplePolyPatch::interpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = interpolate(tff());
    tff.clear();
    return tint;
}


template<class Type>
void Foam::regionCouplePolyPatch::bridge
(
    const Field<Type>& bridgeField,
    Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > regionCouplePolyPatch::bridge\n"
            "(\n"
            "    const Field<Type>& ff,\n"
            "    Field<Type>& ff\n"
            ") const"
        )   << "Incorrect patch field size for bridge.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (bridgeOverlap())
    {
        if (empty())
        {
            // Patch empty, no bridging
            return;
        }

        if (localParallel())
        {
            if (master())
            {
                patchToPatch().bridgeMaster(bridgeField, ff);
            }
            else
            {
                patchToPatch().bridgeSlave(bridgeField, ff);
            }
        }
        else
        {
            // Note: since bridging is only a local operation
            if (master())
            {
                patchToPatch().maskedBridgeMaster
                (
                    bridgeField,
                    ff,
                    zoneAddressing()
                );
            }
            else
            {
                patchToPatch().maskedBridgeSlave
                (
                    bridgeField,
                    ff,
                    zoneAddressing()
                );
            }
        }
    }
}


// ************************************************************************* //
