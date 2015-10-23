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
    Foam::noEngineMesh

Description
    Dummy engine mesh

Author
    Tommaso Lucchini, Politecnico di Milano.

SourceFiles
    noEngineMesh.C

\*---------------------------------------------------------------------------*/

#ifndef noEngineMesh_H
#define noEngineMesh_H

#include "engineTopoChangerMesh.H"
#include "simpleEnginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noEngineMesh Declaration
\*---------------------------------------------------------------------------*/

class noEngineMesh
:
    public engineTopoChangerMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noEngineMesh(const noEngineMesh&);

        //- Disallow default bitwise assignment
        void operator=(const noEngineMesh&);


public:

    //- Runtime type information
    TypeName("noEngineMesh");


    // Constructors

        //- Construct from database
        explicit noEngineMesh(const IOobject& io);


    // Destructor

        virtual ~noEngineMesh();

    // Member Functions


        //- Move and morph
        virtual bool update()
        {
            return false;
        }

        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
