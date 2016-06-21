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

Class
    Foam::FieldMapper

Description
    Abstract base class to hold the Field mapping addressing and weights.

\*---------------------------------------------------------------------------*/

#ifndef FieldMapper_H
#define FieldMapper_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FieldMapper Declaration
\*---------------------------------------------------------------------------*/

class FieldMapper
{

public:

    // Constructors

        //- Null constructor
        FieldMapper()
        {}


    // Destructor

        virtual ~FieldMapper()
        {}


    // Member Functions

        virtual label size() const = 0;

        virtual label sizeBeforeMapping() const = 0;

        virtual bool direct() const = 0;

        virtual const unallocLabelList& directAddressing() const
        {
            FatalErrorIn("FieldMapper::directAddressing() const")
                << "attempt to access null direct addressing"
                << abort(FatalError);

            return unallocLabelList::null();
        }

        virtual const labelListList& addressing() const
        {
            FatalErrorIn("FieldMapper::addressing() const")
                << "attempt to access null interpolation addressing"
                << abort(FatalError);

            return labelListList::null();
        }

        virtual const scalarListList& weights() const
        {
            FatalErrorIn("FieldMapper::weights() const")
                << "attempt to access null interpolation weights"
                << abort(FatalError);

            return scalarListList::null();
        }


    // Member Operators

        template<class Type>
        tmp<Field<Type> > operator()(const Field<Type>& f) const
        {
            return tmp<Field<Type> >(new Field<Type>(f, *this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
