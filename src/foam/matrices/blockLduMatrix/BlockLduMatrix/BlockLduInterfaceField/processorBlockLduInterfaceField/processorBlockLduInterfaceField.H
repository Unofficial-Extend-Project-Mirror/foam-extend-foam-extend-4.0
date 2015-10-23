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

Class
    Foam::processorBlockLduInterfaceField

Description
    Abstract base class for processor coupled interfaces.

Author
    Klas Jareteg, 2013-02-08

SourceFiles
    processorBlockLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef processorBlockLduInterfaceField_H
#define processorBlockLduInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "typeInfo.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class processorBlockLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class processorBlockLduInterfaceField
{

public:

    //- Runtime type information
    TypeName("processorBlockLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        processorBlockLduInterfaceField()
        {}


    // Destructor

        virtual ~processorBlockLduInterfaceField();


    // Member Functions

        // Access

            //- Return processor number
            virtual int myProcNo() const = 0;

            //- Return neigbour processor number
            virtual int neighbProcNo() const = 0;

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;


        //- Transform given patch component field
        void transformCoupleField
        (
            scalarField& f,
            const direction cmpt
        ) const;

        //- Transform given patch field
        void transformCoupleField
        (
            Field<Type>& f
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "processorBlockLduInterfaceField.C"
#endif

// ************************************************************************* //

#endif

// ************************************************************************* //
