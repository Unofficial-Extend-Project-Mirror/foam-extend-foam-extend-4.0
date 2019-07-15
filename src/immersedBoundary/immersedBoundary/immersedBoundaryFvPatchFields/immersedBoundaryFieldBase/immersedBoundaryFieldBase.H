/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::immersedBoundaryFieldBase

Description

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    immersedBoundaryFieldBaseI.H
    immersedBoundaryFieldBase.C
    immersedBoundaryFieldBaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFieldBase_H
#define immersedBoundaryFieldBase_H

#include "immersedBoundaryFvPatch.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class immersedBoundaryFieldBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFieldBase
{
    // Private data

        //- Local reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch_;


        // Dead cell controls

            //- Set dead cell value
            Switch setDeadValue_;

            //- Dead cell value
            Type deadValue_;


public:

    // Constructors

        //- Construct from components
        immersedBoundaryFieldBase
        (
            const fvPatch& p,
            const bool setDeadValue,
            const Type deadValue
        );


    //- Destructor - default


    // Member Functions

        //- Return reference to immersed boundary patch
        inline const immersedBoundaryFvPatch& ibPatch() const
        {
            return ibPatch_;
        }


        //- Dead value handling

            //- Return set dead cell value
            inline bool setDeadValue() const
            {
                return setDeadValue_;
            }

            // Return access to set dead cell value
            inline Switch& setDeadValue()
            {
                return setDeadValue_;
            }

            //- Return dead cell value
            inline const Type& deadValue() const
            {
                return deadValue_;
            }

            // Return access to dead cell value
            inline Type& deadValue()
            {
                return deadValue_;
            }

            //- Set values in dead cells
            void setDeadValues(Field<Type>& psiI) const;

            //- Set matrix constraints in dead cells
            void setDeadValues(fvMatrix<Type>& matrix) const;


        // I/O

            //- Write
            void writeDeadData(Ostream& os) const;

            //- Write
            void writeField(const fvPatchField<Type>& f) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFieldBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
