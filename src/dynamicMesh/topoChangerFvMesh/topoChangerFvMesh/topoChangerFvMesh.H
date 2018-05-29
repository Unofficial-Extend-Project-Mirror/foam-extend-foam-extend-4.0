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
    Foam::topoChangerFvMesh

Description
    Abstract base class for a topology changing fvMesh

SourceFiles
    topoChangerFvMesh.C
    newTopoFvMesh.C
    topoChangerFvMeshLoadBalance.C
    topoChangerFvMeshLoadBalanceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef topoChangerFvMesh_H
#define topoChangerFvMesh_H

#include "dynamicFvMesh.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "fvMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class topoChangerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class topoChangerFvMesh
:
    public dynamicFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoChangerFvMesh(const topoChangerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const topoChangerFvMesh&);


        // Load balancing functions

            //- Send fields to other processor
            template<class GeoField, class Decomposer>
            void sendFields
            (
                const HashTable<const GeoField*>& geoFields,
                const Decomposer& decomposer,
                Ostream& toProc
            ) const;

            //- Insert local fields
            template<class GeoField, class Decomposer>
            void insertFields
            (
                const HashTable<const GeoField*>& geoFields,
                const Decomposer& decomposer,
                List<PtrList<GeoField> >& localFields
            ) const;

            //- Receive fields from other processor
            template<class GeoMesh, class GeoField>
            void receiveFields
            (
                const label procIndex,
                List<PtrList<GeoField> >& receivedFields,
                const GeoMesh& procMesh,
                Istream& fromProc
            ) const;

            //- Resize and rebuild fields
            template<class GeoField, class Reconstructor>
            void rebuildFields
            (
                const HashTable<const GeoField*>& geoFields,
                const Reconstructor& reconstructor,
                const List<PtrList<GeoField> >& receivedFields,
                const mapPolyMesh& meshMap
            ) const;


protected:

    // Protected Data

        //- Topo changer
        polyTopoChanger topoChanger_;


public:

    //- Runtime type information
    TypeName("topoChangerFvMesh");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit topoChangerFvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addPatches() member function
        topoChangerFvMesh
        (
            const IOobject& io,
            const Xfer<pointField>& points,
            const Xfer<faceList>& faces,
            const Xfer<labelList>& allOwner,
            const Xfer<labelList>& allNeighbour,
            const bool syncPar = true
        );

        //- Construct from components with cells rather than owner
        //  and neighbour without boundary.
        //  Boundary is added using addPatches() member function
        topoChangerFvMesh
        (
            const IOobject& io,
            const Xfer<pointField>& points,
            const Xfer<faceList>& faces,
            const Xfer<cellList>& cells,
            const bool syncPar = true
        );


    //- Destructor
    virtual ~topoChangerFvMesh();


    // Member Functions

        //- Load balance the mesh in parallel execution
        bool loadBalance(const dictionary& dict);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "topoChangerFvMeshLoadBalanceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
