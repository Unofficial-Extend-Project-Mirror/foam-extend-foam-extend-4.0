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
    Foam::processorMeshes

Description
    Container for processor mesh addressing.

SourceFiles
    processorMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef processorMeshes_H
#define processorMeshes_H

#include "PtrList.H"
#include "fvMesh.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorMeshes Declaration
\*---------------------------------------------------------------------------*/

class processorMeshes
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
        processorMeshes(const processorMeshes&);

        //- Disallow default bitwise assignment
        void operator=(const processorMeshes&);


        //- Read all meshes
        void readMeshes();

        //- Read addressing
        void readAddressing();


public:

    // Constructors

        //- Construct from components
        processorMeshes(PtrList<Time>& databases, const word& meshName);


    // Member Functions

        //- Update the meshes based on the mesh files saved in
        //  time directories
        fvMesh::readUpdateState readUpdate();

        //- Reconstruct point position after motion in parallel
        void reconstructPoints(fvMesh& mesh);

        PtrList<fvMesh>& meshes()
        {
            return meshes_;
        }

        //- Return point-processor addressing arrays
        const PtrList<labelIOList>& pointProcAddressing() const
        {
            return pointProcAddressing_;
        }

        //- Return face-processor addressing arrays
        const PtrList<labelIOList>& faceProcAddressing() const
        {
            return faceProcAddressing_;
        }

        //- Return cell-processor addressing arrays
        const PtrList<labelIOList>& cellProcAddressing() const
        {
            return cellProcAddressing_;
        }

        //- Return boundary-processor addressing arrays
        const PtrList<labelIOList>& boundaryProcAddressing() const
        {
            return boundaryProcAddressing_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
