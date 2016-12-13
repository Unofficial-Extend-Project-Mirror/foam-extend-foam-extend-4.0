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
    pointConstraint

Description
    Accumulates point constraints through successive applications of the
    applyConstraint function.  After all the constraints have been entered the
    resulting transformation tensor is returned by the constraintTransformation
    function.

SourceFiles
    pointConstraintI.H

\*---------------------------------------------------------------------------*/

#ifndef pointConstraint_H
#define pointConstraint_H

#include "Tuple2.H"
#include "label.H"
#include "vector.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointConstraint Declaration
\*---------------------------------------------------------------------------*/

class pointConstraint
:
    public Tuple2<label, vector>
{

public:

    // Constructors

        //- Construct null
        inline pointConstraint();

        //- Construct from Istream
        inline pointConstraint(Istream&);


    // Member Functions

        //- Apply and accumulate the effect of the given constraint direction
        inline void applyConstraint(const vector& cd);

        //- Return the accumulated constraint transformation tensor
        inline tensor constraintTransformation() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointConstraintI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
