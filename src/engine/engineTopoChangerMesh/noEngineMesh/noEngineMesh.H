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
