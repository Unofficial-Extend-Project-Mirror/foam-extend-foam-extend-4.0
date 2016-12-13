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
    Foam::geomDecomp

Description
    Geometrical domain decomposition

SourceFiles
    geomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef geomDecomp_H
#define geomDecomp_H

#include "decompositionMethod.H"
#include "labelVector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class geomDecomp Declaration
\*---------------------------------------------------------------------------*/

class geomDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        geomDecomp(const geomDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const geomDecomp&);


protected:

    // Protected data

        //- Geometric decomposition dictionary
        const dictionary& geomDecomDict_;

        //- Number of splits in the (x y z) direction
        labelVector n_;

        //- Small rotation to achieve smooth geometric decomposition
        scalar delta_;

        //- Small rotation tensor, calculated from delta
        tensor rotDelta_;


public:

    // Constructors

        //- Construct given the decomposition dictionary
        //  and the derived type name
        geomDecomp
        (
            const dictionary& decompositionDict,
            const word& derivedType
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
