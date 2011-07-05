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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "ggiPolyPatch.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// New.  HJ, 12/Jun/2011
template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::fastExpand
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
    // C) Expanded data will be equal to actual data fronm other processors
    //    only for the faces marked in remote; for other faces, it will be
    //    equal to zero
    // D) On processor zero, complete data is available
    // HJ, 4/Jun/2011
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::fastExpand\n"
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
            "tmp<Field<Type> > ggiPolyPatch::fastExpand"
            "("
            "    const Field<Type>& ff"
            ") const"
        )   << "Requested expand on local parallel.  This is not allowed"
            << abort(FatalError);
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


// Obsolete.  HJ, 12/Jun/2011
template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::expand
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::expand"
            "("
            "    const Field<Type>& ff"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    // Expand the field to zone size
    tmp<Field<Type> > texpandField
    (
        new Field<Type>(zone().size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField();

    const labelList& addr = zoneAddressing();

    forAll (addr, i)
    {
        expandField[addr[i]] = ff[i];
    }

    // Parallel data exchange: update surface field on all processors
    // Operation is performed if the patch is not localised
    // HJ, 8/Dec/2009

    if (!localParallel())
    {
        reduce(expandField, sumOp<Field<Type> >());
    }

    return texpandField;
}


// Obsolete.  HJ, 12/Jun/2011
template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::filter
(
    const Field<Type>& ef
) const
{
    // Check and expand the field from patch size to zone size
    if (ef.size() != zone().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::filter"
            "("
            "    const Field<Type>& ef"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << ef.size() << " patch size: " << zone().size()
            << abort(FatalError);
    }

    // Filter the field to zone size
    tmp<Field<Type> > tfilterField
    (
        new Field<Type>(size(), pTraits<Type>::zero)
    );

    Field<Type>& filterField = tfilterField();

    const labelList& addr = zoneAddressing();

    forAll (addr, i)
    {
        filterField[i] = ef[addr[i]];
    }

    return tfilterField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::interpolate
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != shadow().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::interpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect slave patch field size.  Field size: "
            << ff.size() << " patch size: " << shadow().size()
            << abort(FatalError);
    }

#   if 1

    // New.  HJ, 12/Jun/2011
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

#   else

    // Old treatment
    // Obsolete.  HJ, 12/Jun/2011

    // Expand the field to zone size
    Field<Type> expandField = shadow().expand(ff);

    Field<Type> zoneField;

    // Interpolate expanded field
    if (master())
    {
        zoneField = patchToPatch().slaveToMaster(expandField);
    }
    else
    {
        zoneField = patchToPatch().masterToSlave(expandField);
    }

    return this->filter(zoneField);

#   endif
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::interpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = interpolate(tff());
    tff.clear();
    return tint;
}


template<class Type>
void Foam::ggiPolyPatch::bridge
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
            "tmp<Field<Type> > ggiPolyPatch::bridge\n"
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
#       if 1

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

#       else

        // Expand the field to zone size
        Field<Type> expandBridge = this->expand(bridgeField);
        Field<Type> expandField = this->expand(ff);

        if (master())
        {
            patchToPatch().bridgeMaster(expandBridge, expandField);
        }
        else
        {
            patchToPatch().bridgeSlave(expandBridge, expandField);
        }

        // Contract the field from zone size to patch size
        const labelList& addr = zoneAddressing();

        forAll (addr, i)
        {
            ff[i] = expandField[addr[i]];
        }

#       endif
    }
}


// ************************************************************************* //
