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
    Foam::immersedBoundaryFvPatch

Description
    Immersed boundary FV patch

Notes
    Rewrite of the immersed boundary method for increased efficiency
    and robustness.  Main changes
    - using octree for fast search of eligible cells using indexedOctree
    - parallel communication using mapDistribute

    Complete rewrite of Immersed Boundary methodology: use cell cutting
    at fvMesh level.  Hrvoje Jasak.

Author
    Zeljko Tukovic
    Reorganisation and optimisation by Hrvoje Jasak

SourceFiles
    immersedBoundaryFvPatch.C
    immersedBoundaryFvPatchTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvPatch_H
#define immersedBoundaryFvPatch_H

#include "fvPatch.H"
#include "immersedBoundaryPolyPatch.H"
#include "volFieldsFwd.H"
#include "slicedVolFieldsFwd.H"
#include "slicedSurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class immersedBoundaryFvPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryFvPatch
:
    public fvPatch
{
    // Private data

        //- Reference to immersed boundary patch
        const immersedBoundaryPolyPatch& ibPolyPatch_;

        //- Finite volume mesh reference
        const fvMesh& mesh_;


    // Static data

        //- Live volume factor
        //  Fraction of cell volume for which the non-orthogonal correction
        //  is deactivated
        static const debug::tolerancesSwitch nonOrthogonalFactor_;


protected:

    // Protected Member Functions

        // Geometry correction functions

            //- Make mesh face centres
            virtual void makeCf(slicedSurfaceVectorField&) const;

            //- Make mesh face areas
            virtual void makeSf(slicedSurfaceVectorField&) const;

            //- Make mesh cell centres
            virtual void makeC(slicedVolVectorField&) const;

            //- Make mesh cell volumes
            virtual void makeV(scalarField&) const;


        // Discretisation correction functions

            //- Make patch face non-orthogonality correction vectors
            virtual void makeCorrVecs(fvsPatchVectorField&) const;


public:

    //- Runtime type information
    TypeName(immersedBoundaryPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        immersedBoundaryFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~immersedBoundaryFvPatch()
    {}


    // Member Functions

        // Access

            //- Return size equal to number of intersected cells
            virtual label size() const;

            //- Return faceCells next to the IB surface
            virtual const unallocLabelList& faceCells() const;


        // Access functions for geometrical data

            //- Return face normals
            //  Reconsider. HJ, 11/Dec/2017
            tmp<vectorField> nf() const;

            //- Return cell-centre to face-centre vector
            //  Reconsider. HJ, 11/Dec/2017
            virtual tmp<vectorField> delta() const;


        // Access to immersed boundary components

            //- Return reference to immersed boundary polyPatch
            const immersedBoundaryPolyPatch& ibPolyPatch() const
            {
                return ibPolyPatch_;
            }

            //- Return immersed boundary surface mesh
            const triSurfaceMesh& ibMesh() const
            {
                return ibPolyPatch_.ibMesh();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
