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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    adaptiveOverlapFringe

Description
    Adaptive overlap fringe algorithm with fall-back mechanism

    Holds and manages following data:

      - fringeHolesPtr_: holds fringe holes which are basically discarded
        acceptors from previous iterations
      - acceptorsPtr_: holds current set of acceptors to consider
      - finalDonorAcceptorsPtr_: final set of donor/acceptors which will be used
        for the overset interpolation
      - donorSuitability_: run-time selectable Donor Suitability Function
      - initPatchNames_: optional list of patches that we use to start the
        acceptor search from
      - minGlobalFraction_: minimum allowable global fraction of suitable pairs
      - cumulativeDonorAcceptorsPtr_: list of suitable donor/acceptor pairs
      - fringeIter_: iteration counter

    Member function updateIteration controls the iterative process:
      - Checks whether given donor/acceptors are suitable and filters them,
      - If the user-specified criterion controlled by minGlobalFraction_ is not
        reached (i.e. suitable overlap has not been found), then we loop through
        unsuitable acceptors and mark their neighbours (excluding existing
        acceptors and holes) as candidate acceptors for the next iteration,
      - Unsuitable acceptors become fringe holes.

Authors
    Josip Katavic, FMENA Zagreb.  All rights reserved.

SourceFiles
    adaptiveOverlapFringe.C

Notes
    Based on original overlapFringe assembly algorithm.

\*---------------------------------------------------------------------------*/

#ifndef adaptiveOverlapFringe_H
#define adaptiveOverlapFringe_H

#include "overlapFringe.H"
#include "donorSuitability.H"
#include "optimisationSwitch.H"
#include "volFields.H"
#include "FIFOStack.H"
#include "iterationData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class adaptiveOverlapFringe Declaration
\*---------------------------------------------------------------------------*/

class adaptiveOverlapFringe
:
    public oversetFringe
{
    // Private data

        //- Fringe hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;

        //- Final donor/acceptor pairs for this region (fringe)
        mutable donorAcceptorList* finalDonorAcceptorsPtr_;

        //- Iteration history - filtered donor/acceptor pairs, fringe holes and
        // average donor/acceptor suitability fraction
        mutable FIFOStack<iterationData> iterationDataHistory_;


        // Data for iteration control

            //- Run-time selectable donor/acceptor suitability handler
            autoPtr<donorSuitability::donorSuitability> donorSuitability_;

            //- List of patches to start the iterative fringe assembly process
            //  (empty list by default)
            const wordList initPatchNames_;

            //- Iteration counter
            mutable label fringeIter_;

            //- User defined number of iterations - by default 4 iterations will
            // be made
            const label specifiedIterationsNumber_;

            //- Relative iteration counter
            mutable label relativeCounter_;

            //- Additional iteration switch
            const Switch additionalIterations_;

            //- Penalty coefficient for orhpan donors, i.e. orhpan donors
            // suitability
            const scalar orphanSuitability_;

            //- Minimum local suitability - if donor/acceptor suitability
            //  of donor/acceptor pair is lower than minimum local
            //  fraction, neighbours of this acceptor are new acceptor
            //  candidates
            const scalar minLocalSuit_;

            // Suitable donor/acceptor pairs
            mutable donorAcceptorDynamicList suitableDAPairs_;

            //- Suitable pairs cumulative suitability
            mutable scalar suitablePairsSuit_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            adaptiveOverlapFringe(const adaptiveOverlapFringe&);

            //- Disallow default bitwise assignment
            void operator=(const adaptiveOverlapFringe&);


        // Addressing manipulation functions

            //- Calculate hole-acceptor addressing
            void calcAddressing() const;

            //- Clear hole-acceptor addressing
            void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("adaptiveOverlap");


    // Constructors

        //- Construct from dictionary
        adaptiveOverlapFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    //- Destructor
    virtual ~adaptiveOverlapFringe();


    // Member Functions

        //- Update iteration.
        virtual bool updateIteration
        (
            donorAcceptorList& donorAcceptorRegionData
        ) const;

        //- Return list of hole cells
        //  Fringe hole cells are collected in addition to geometric hole
        //  cells, which fall outside of all donor regions
        virtual const labelList& fringeHoles() const;

        //- Return list of acceptor cells
        virtual const labelList& candidateAcceptors() const;

        //- Return list of final donor acceptor pairs. Note: caller may
        //  invalidate finalDonorAcceptorsPtr_ for optimisation purposes
        virtual donorAcceptorList& finalDonorAcceptors() const;

        //- Update the fringe
        virtual void update() const;

        //- Sort function - third parameter for max_element function
        static bool sortBySuitabilityValue
        (
            const iterationData& lhsObject,
            const iterationData& rhsObject
        );

        //- Loop through all iteration data objects from iterationDataHistory_
        //  data member and calculate slope, i.e. suitability fraction values
        //  trend using linear least squares method.
        void suitabilityFractionSlope
        (
            FIFOStack<iterationData>,
            scalar& alpha,
            scalar& beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
