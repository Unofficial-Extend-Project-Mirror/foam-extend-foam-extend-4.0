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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::expand
(
    const Field<Type>& pf
) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::expand"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << size()
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
        expandField[addr[i]] = pf[i];
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
    const Field<Type>& pf
) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != shadow().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::interpolate"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect slave patch field size.  Field size: "
            << pf.size() << " patch size: " << shadow().size()
            << abort(FatalError);
    }

    // Expand the field to zone size
    Field<Type> expandField = shadow().expand(pf);

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
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::interpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = interpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Type>
void Foam::ggiPolyPatch::bridge
(
    const Field<Type>& bridgeField,
    Field<Type>& ff
) const
{
    if (bridgeOverlap())
    {
        // Parallelisation

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
    }
    else
    {
        // HJ, temporary
        InfoIn
        (
            "void bridge(const Field<Type>& bridgeField, "
            "Field<Type>& ff) const"
        )   << "Bridging is switched off " << endl;
    }
}


// ************************************************************************* //
