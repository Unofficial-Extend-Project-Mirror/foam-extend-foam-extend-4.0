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
    PODEigenBase

Description
    Class which assembles scalar eigen-base given a set of scalar fields.
    This is performed by establishing a matrix of snapshots, calculating and
    sorting them and providing corresponding eigen-vectors.  Eigen-values are
    sorted in increasing order.  Snapshots are added one at a time and
    when the base is completed, the calculation is triggered

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    PODEigenBase.C

\*---------------------------------------------------------------------------*/

#ifndef PODEigenBase_H
#define PODEigenBase_H

#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "EigenSolver.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PODEigenBase Declaration
\*---------------------------------------------------------------------------*/

class PODEigenBase
{
    // Private data

        //- Eigen-values
        scalarField eigenValues_;

        //- Cumulative relative eigen-values.  Renormalised to sum to 1
        scalarField cumEigenValues_;

        //- Eigen-vectors
        FieldField<Field, scalar> eigenVectors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PODEigenBase(const PODEigenBase&);

        //- Disallow default bitwise assignment
        void operator=(const PODEigenBase&);


        //- Calculate eigen base
        void calcEigenBase(const scalarSquareMatrix& orthMatrix);


public:

    // Constructors

        //- Construct given a list of fields
        PODEigenBase(const PtrList<volScalarField>& snapshots);


    // Destructor - default


    // Member Functions

        //- Return eigen-values sorted in decreasing order
        const scalarField& eigenValues() const
        {
            return eigenValues_;
        }

        //- Return cumulative eigen-values
        const scalarField& cumulativeEigenValues() const
        {
            return cumEigenValues_;
        }

        //- Return eigen-vectors
        const FieldField<Field, scalar>& eigenVectors() const
        {
            return eigenVectors_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
