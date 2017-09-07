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
    Foam::dynamicRefinePolyFvMesh

Description
    A fvMesh with built-in refinement of arbitrary polyhedral cells.

    Determines which cells to refine/unrefine and does all in update().

SourceFiles
    dynamicRefinePolyFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefinePolyFvMesh_H
#define dynamicRefinePolyFvMesh_H

#include "dynamicFvMesh.H"
#include "polyRef.H"
#include "PackedBoolList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dynamicRefinePolyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefinePolyFvMesh
:
    public dynamicFvMesh
{

    // Private data

        // Helper variables to enable switching between a single and multiple
        // mesh motion updates within a time step (if update() is called more
        // than once in a single time step)

            //- Switch for single motion update (true by default)
            Switch singleMotionUpdate_;

            //- Helper varaible: current time index
            label curTimeIndex_;


protected:

        //- Mesh cutting engine
        polyRef meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        List<Pair<word> > correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;


    // Private Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Read the projection parameters from dictionary
        void readDict();


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to un/refine

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point min of connected cell
            scalarField minCellField(const volScalarField&) const;

            //- Simple (non-parallel) interpolation by averaging
            scalarField cellToPoint(const scalarField& vFld) const;

            //- Calculate error (= -1 by default or distance from inbetween
            //  levels
            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalarField& vFld,
                PackedBoolList& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell
            void extendMarkedCells(PackedBoolList& markedCell) const;


private:

        //- Disallow default bitwise copy construct
        dynamicRefinePolyFvMesh(const dynamicRefinePolyFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefinePolyFvMesh&);

public:

    //- Runtime type information
    TypeName("dynamicRefinePolyFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefinePolyFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicRefinePolyFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const polyRef& meshCutter() const
        {
            return meshCutter_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
