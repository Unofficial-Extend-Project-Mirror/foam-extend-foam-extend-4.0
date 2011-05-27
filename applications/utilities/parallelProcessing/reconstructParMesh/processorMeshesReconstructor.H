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
    Foam::processorMeshesReconstructor

Description

SourceFiles
    processorMeshesReconstructor.C
    processorMeshesRebuild.C

\*---------------------------------------------------------------------------*/

#ifndef processorMeshesReconstructor_H
#define processorMeshesReconstructor_H

#include "PtrList.H"
#include "Time.H"
#include "fvMesh.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class processorPolyPatch;

/*---------------------------------------------------------------------------*\
                  Class processorMeshesReconstructor Declaration
\*---------------------------------------------------------------------------*/

class processorMeshesReconstructor
{
    // Private data

        //- Processor databases
        PtrList<Time>& databases_;

        const word meshName_;

        //- List of processor meshes
        PtrList<fvMesh> meshes_;

        //- List of processor point addressing lists
        PtrList<labelIOList> pointProcAddressing_;

        //- List of processor face addressing lists
        PtrList<labelIOList> faceProcAddressing_;

        //- List of processor cell addressing lists
        PtrList<labelIOList> cellProcAddressing_;

        //- List of processor boundary addressing lists
        PtrList<labelIOList> boundaryProcAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        processorMeshesReconstructor(const processorMeshesReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const processorMeshesReconstructor&);


        //- Read all meshes
        void readMeshes();

        //- Attempt to read mapping.  If not available, return false
        bool readMapping();

        //- Return neighbour processor patch
        const processorPolyPatch& neighbourProcPatch
        (
            const processorPolyPatch& procPatch
        ) const;

        //- Return neighbour processor patch
        labelList procPatchPointMapping
        (
            const processorPolyPatch& procPatch
        ) const;

        //- Clear reconstruction maps
        void clearMaps();


public:

    // Constructors

        //- Construct from components
        processorMeshesReconstructor
        (
            PtrList<Time>& databases,
            const word& meshName
        );


    // Destructor - default


    // Member Functions

        //- Update the meshes based on the mesh files saved in
        //  time directories
        fvMesh::readUpdateState readUpdate();

        //- Reconstruct mesh after topology change in parallel
        autoPtr<fvMesh> reconstructMesh(const Time& db);

        //- Reconstruct point position after motion in parallel
        void reconstructPoints(fvMesh& mesh) const;

        //- Return meshes
        const PtrList<fvMesh>& meshes()
        {
            return meshes_;
        }

        //- Return point-processor addressing arrays
        const PtrList<labelIOList>& pointProcAddressing() const;

        //- Return face-processor addressing arrays
        const PtrList<labelIOList>& faceProcAddressing() const;

        //- Return cell-processor addressing arrays
        const PtrList<labelIOList>& cellProcAddressing() const;

        //- Return boundary-processor addressing arrays
        const PtrList<labelIOList>& boundaryProcAddressing() const;


        //- Write addressing
        void writeAddressing();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
