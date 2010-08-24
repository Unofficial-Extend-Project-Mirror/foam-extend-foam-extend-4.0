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

\*---------------------------------------------------------------------------*/

#ifndef pistonRefine_H
#define pistonRefine_H

#include "engineTopoChangerMesh.H"
#include "simpleEnginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"

#include "dynamicFvMesh.H"
#include "hexRef8.H"
#include "PackedList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class pistonRefine Declaration
\*---------------------------------------------------------------------------*/

class pistonRefine
:
    public engineTopoChangerMesh
{

protected:

         //- Piston patch
        simpleEnginePiston piston_;        

        //- Piston Position
        scalar pistonPosition_;

        //- deckHeight
        scalar deckHeight_;

       //- Mesh cutting engine
        hexRef8 meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        List<Pair<word> > correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedList<1> protectedCell_;


    // Private Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedList<1>&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedList<1>& unrefineableCell) const;

        //- Read the projection parameters from dictionary
        void readDict();


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point min of connected cell
            scalarField minCellField(const volScalarField&) const;

            scalarField cellToPoint(const scalarField& vFld) const;

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
                PackedList<1>& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedList<1>& candidateCell
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedList<1>& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedList<1>& markedCell) const;

private:
/*
        //- Disallow default bitwise copy construct
        dynamicRefineFvMesh(const dynamicRefineFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineFvMesh&);
*/

////////////////////////////////////////////////    // Private data

               
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pistonRefine(const pistonRefine&);

        //- Disallow default bitwise assignment
        void operator=(const pistonRefine&);
        
        //- Check if all patches exist, then calculate virtualPistonPosition, 
        //- pistonPosition and deckHeight for the first time
        void checkAndCalculate();
        
        //- Add valve and piston zones and modifiers
//        void addZonesAndModifiers(){}

public:

    //- Runtime type information
    TypeName("pistonRefine");


    // Constructors

        //- Construct from database
        explicit pistonRefine(const IOobject& io);


    // Destructor 

        virtual ~pistonRefine();

    // Member Functions

        //- Direct access to the refinement engine
        const hexRef8& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedList<1>& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedList<1>& protectedCell()
        {
            return protectedCell_;
        }

    // Member Functions

        
        const scalar& pistonPosition() const
        {
            return pistonPosition_;
        }

        scalar& pistonPosition() 
        {
            return pistonPosition_;
        }
        
        
        const scalar& deckHeight() const
        {
            return deckHeight_;
        }

        scalar& deckHeight() 
        {
            return deckHeight_;
        }
        

        const simpleEnginePiston& piston() const
        {
            return piston_;
        }


        //- Move and morph
        virtual bool update();
        
        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U);
        
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
