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

#include "ggiAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > ggiAMGInterface::fastExpand(const UList<Type>& ff) const
{
    // Rewrite, 1/Jun/2016
    // To avoid creating zone-sized data and gather-scatter communication
    // to the master, the optimised map-distribute call is implemented.
    // The field is filled with local data which is then sent where needed
    // through map-distribute.
    // On return, the field is expanded to zone size but only filled with
    // the data which is needed for the shadow
    // HJ, 1/Jun/2016

    if (ff.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiAMGInterface::fastExpand"
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
    else
    {
        // Optimised mapDistribute

        // Execute init reduce to calculate addressing if not already done
        map();

        // Prepare for distribute.  Note: field will be expanded to zone size
        // during the distribute operation
        tmp<Field<Type> > tresult(new Field<Type>(ff));
        List<Type>& expand = tresult();

        map().distribute(expand);

        return tresult;
    }
}


template<class Type>
tmp<Field<Type> > ggiAMGInterface::fastReduce(const UList<Type>& ff) const
{
    // Old algorithm: OBOSLETE
    // Local processor contains faceCells part of the zone and requires
    // zoneAddressing part.
    // For fast communications, each processor will send the faceCells and
    // zoneAddressing to the master.  Master will assemble global zone
    // and send off messages to all processors containing only
    // the required data
    // HJ, 24/Jun/2011

    // Rewrite, 1/Jun/2016
    // To avoid creating zone-sized data and gather-scatter communication
    // to the master, the optimised map-distribute call is implemented.
    // The field is filled with local data which is then sent where needed
    // through map-distribute.
    // On return, the field is expanded to zone size but only filled with
    // the data which is needed for the shadow
    // Having received the zone data, shadow data is extracted from the
    // field size.  Note: this works only on coarse levels, where one-on-one
    // mapping applies
    // HJ, 1/Jun/2016

    if (ff.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiAMGInterface::fastReduce"
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
    else
    {
        // Optimised mapDistribute

        // Execute init reduce to calculate addressing if not already done
        map();

        // Prepare for distribute.  Note: field will be expanded to zone size
        // during the distribute operation
        List<Type> expand = ff;

        map().distribute(expand);

        const labelList& shadowZa = shadowInterface().zoneAddressing();

        // Prepare return field: zone size
        tmp<Field<Type> > tresult
        (
            new Field<Type>(shadowZa.size())
        );
        Field<Type>& result = tresult();

        // Filter from expanded field to zone size
        forAll (shadowZa, shadowZaI)
        {
            result[shadowZaI] = expand[shadowZa[shadowZaI]];
        }

        return tresult;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
