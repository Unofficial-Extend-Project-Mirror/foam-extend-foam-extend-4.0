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
    Foam::magLongDelta

Description
    Long distance to the neighbour via face centre

SourceFiles
    magLongDelta.C

\*---------------------------------------------------------------------------*/

#ifndef magLongDelta_H
#define magLongDelta_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "primitiveFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class magLongDelta Declaration
\*---------------------------------------------------------------------------*/

class magLongDelta
:
    public MeshObject<fvMesh, magLongDelta>
{
    // Private data

        //- Magnitude of the normal distances to the face
        mutable surfaceScalarField* magLongDeltaPtr_;

        //- Magnitude of the normal distances on coupled patches
        //- Data in magLongDeltaPtr_ is reused, if it exists
        mutable List<scalarField*> magLongDeltaBnd_;


    // Private member functions

        //- Construct Long distance
        void makeMagLongDistance() const;

        //- Construct Long distance
        void makeMagLongDistance(label patchi) const;

        //- Clear Data
        void clearData() const;

public:

    // Declare name of the class and its debug switch
    TypeName("magLongDelta");


    // Constructors

        //- Construct given an fvMesh
        explicit magLongDelta(const fvMesh&);


    // Destructor

        virtual ~magLongDelta();


    // Member functions

        //- Return reference to the magnitude of the long distance
        const surfaceScalarField& magDelta() const;

        //- Return reference to the magnitude of the long distance
        //- for patch i
        const scalarField& magDelta(const label patchi) const;

        //- Update after mesh motion:
        //  Delete data when the mesh moves
        virtual bool movePoints() const;

        //- Update after topo change:
        //  Delete data when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
