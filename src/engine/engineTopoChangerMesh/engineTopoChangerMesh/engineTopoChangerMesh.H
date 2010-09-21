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
    engineTopoChangerMesh

Description
    Base class for engine mesh with topological changes

SourceFiles
    engineTopoChangerMesh.C
    newTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef engineTopoChangerMesh_H
#define engineTopoChangerMesh_H

#include "fvMesh.H"
#include "polyTopoChanger.H"
#include "engineTime.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class engineTopoChangerMesh Declaration
\*---------------------------------------------------------------------------*/

class engineTopoChangerMesh
:
    public fvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        engineTopoChangerMesh(const engineTopoChangerMesh&);

        //- Disallow default bitwise assignment
        void operator=(const engineTopoChangerMesh&);


protected:

    // Protected Data

        polyTopoChanger topoChanger_;

        //- Engine database
        const engineTime& engineTime_;


public:

    //- Runtime type information
    TypeName("engineTopoChangerMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            engineTopoChangerMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit engineTopoChangerMesh(const IOobject& io);


    // Selectors

        //- Select null constructed
        static autoPtr<engineTopoChangerMesh> New(const IOobject& io);

    // Destructor

        virtual ~engineTopoChangerMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;

    // Access

        const engineTime& engTime() const
        {
            return engineTime_;
        }


    // Member Functions

        //- Set Boundary Velocities
        virtual void setBoundaryVelocity(volVectorField& U) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
