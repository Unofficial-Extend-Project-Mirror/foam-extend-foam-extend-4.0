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

#include "ggiGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > ggiGAMGInterface::fastReduce(const UList<Type>& ff) const
{
    // Algorithm
    // Local processor contains faceCells part of the zone and requires
    // zoneAddressing part.
    // For fast communications, each processor will send the faceCells and
    // zoneAddressing to the master.  Master will assemble global zone
    // and send off messages to all processors containing only
    // the required data
    // HJ, 24/Jun/2011

    if (ff.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiGAMGInterface::fastReduce"
            "("
            "    const UList<Type>& ff"
            ") const"
        )   << "Wrong field size.  ff: " << ff.size()
            << " interface: " << this->size()
            << abort(FatalError);
    }

    if (localParallel() || !Pstream::parRun())
    {
        // Field remains identical: no parallel communications required
        tmp<Field<Type> > tresult(new Field<Type>(ff));

        return tresult;
    }

    // Execute reduce if not already done
    if (!initReduce_)
    {
        initFastReduce();
    }

    if (Pstream::master())
    {
        // Master collects information and distributes data.
        Field<Type> expandField(zoneSize(), pTraits<Type>::zero);

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
            const labelList& curRAddr = receiveAddr_[procI];

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
            const labelList& curSAddr = sendAddr_[procI];

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

        // Note: different from ggi patch: field reduction happens within
        // fastReduce.  HJ, 26/Jun/2011
        const labelList& sza = shadowInterface().zoneAddressing();

        tmp<Field<Type> > tredField
        (
            new Field<Type>(sza.size(), pTraits<Type>::zero)
        );
        Field<Type>& redField = tredField();

        // Select elements from shadow zone addressing
        forAll (sza, i)
        {
            redField[i] = expandField[sza[i]];
        }

        return tredField;
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
        const labelList& sza = shadowInterface().zoneAddressing();

        tmp<Field<Type> > treceiveBuf
        (
            new Field<Type>(sza.size(), pTraits<Type>::zero)
        );
        Field<Type>& receiveBuf = treceiveBuf();

        if (!sza.empty())
        {
            // Opt: reconsider mode of communication
            IPstream::read
            (
                Pstream::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(receiveBuf.begin()),
                receiveBuf.byteSize()
            );

            // Note: different from ggi patch: field reduction happens within
            // fastReduce.  HJ, 26/Jun/2011
        }

        return treceiveBuf;
    }
}


template<class Type>
tmp<Field<Type> > ggiGAMGInterface::expand(const UList<Type>& pd) const
{
    // Expand interface data to complete zone
    tmp<Field<Type> > ted
    (
        new Field<Type>(zoneSize(), pTraits<Type>::zero)
    );
    Field<Type>& ed = ted();

    const labelList& za = zoneAddressing();

    forAll (za, i)
    {
        ed[za[i]] = pd[i];
    }

    // Reduce zone data if the patch is distributed
    if (!localParallel())
    {
        reduce(ed, sumOp<Field<Type> >());
    }

    return ted;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
