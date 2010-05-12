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
    Foam::polyTopoChanger

Description
    List of mesh modifiers defining the mesh dynamics.

SourceFiles
    polyTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChanger_H
#define polyTopoChanger_H

#include "regIOobject.H"
#include "PtrList.H"
#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;
class mapPolyMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class polyTopoChanger Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChanger
:
    public PtrList<polyMeshModifier>,
    public regIOobject
{
    // Private data

        //- Current time index for mesh morph
        label morphIndex_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyTopoChanger(const polyTopoChanger&);

        //- Disallow default bitwise assignment
        void operator=(const polyTopoChanger&);

        void readModifiers();


protected:

    // Protected data

        //- Reference to mesh
        polyMesh& mesh_;


    // Protected Member Functions

        // Topology changes

            //- Rotate a face nShift positions in anticlockwise direction
            static face rotateFace(const face& f, const label nShift);

            //- Determine ordering of faces in coupled patches.
            //  Calculate mapping to shuffle faces inside coupled patches and
            //  rotation to make 0th vertex in faces align.
            static bool reorderCoupledPatches
            (
                const polyBoundaryMesh& boundary,
                const labelList& patchStarts,
                const labelList& patchSizes,
                const faceList& faces,
                const pointField& points,
                labelList& faceMap,
                labelList& rotation
            );

            //- Sync communications required for couple patch reordering when
            //  there is no local topological change
            void syncCoupledPatches();


public:

    //- Runtime type information
    TypeName("polyTopoChanger");


    // Constructors

        //- Read constructor given IOobject and a polyMesh
        polyTopoChanger(const IOobject&, polyMesh&);

        //- Read constructor for given polyMesh
        explicit polyTopoChanger(polyMesh&);


    // Destructor

        virtual ~polyTopoChanger()
        {}


    // Member functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of patch names
        wordList names() const;

        //- Is topology change required
        bool changeTopology() const;

        //- Return topology change request
        autoPtr<polyTopoChange> topoChangeRequest() const;

        //- Modify point motion
        void modifyMotionPoints(pointField&) const;

        static autoPtr<mapPolyMesh> changeMesh
        (
            polyMesh&,
            const polyTopoChange&
        );

        autoPtr<mapPolyMesh> changeMesh();

        //- Force recalculation of locally stored data on topological change
        void update(const mapPolyMesh& m);

        //- Add given set of topology modifiers to the topoChanger
        void addTopologyModifiers(const List<polyMeshModifier*>& tm);

        //- Find modifier given a name
        label findModifierID(const word& modName) const;

        //- Return morph index
        label morphIndex() const;


        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        bool operator!=(const polyTopoChanger&) const;
        bool operator==(const polyTopoChanger&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyTopoChanger&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
