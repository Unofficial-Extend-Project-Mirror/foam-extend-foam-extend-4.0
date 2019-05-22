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
    cuttingPatchFringe

Description
    Fringe algorithm based on collecting all the cells cut by its donor region
    patch (assuming faceCells algorithm is used). Inside cells are marked and
    moved nLayers towards the interior to define the acceptors as robustly as
    possible.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    cuttingPatchFringe.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPatchFringe_H
#define cuttingPatchFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cuttingPatchFringe Declaration
\*---------------------------------------------------------------------------*/

class cuttingPatchFringe
:
    public oversetFringe
{
    // Private data

        //- Names of the cuttingPatches
        wordList cuttingPatchNames_;

        //- How many layers to move away from connected region donors to define
        //  acceptor (and holes)
        label nLayers_;

        //- Fringe hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;

        //- Final donor/acceptor pairs for this region (fringe)
        mutable donorAcceptorList* finalDonorAcceptorsPtr_;


    // Private Member Functions

        //- Calculate hole and acceptor addressing
        void calcAddressing() const;

        //- Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("cuttingPatch");


    // Constructors

        //- Construct from dictionary
        cuttingPatchFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construct
        cuttingPatchFringe
        (
            const cuttingPatchFringe&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const cuttingPatchFringe&) = delete;


    // Destructor

        virtual ~cuttingPatchFringe();


    // Member Functions

        //- Update iteration. Note: invalidates parameter
        virtual bool updateIteration
        (
            donorAcceptorList& donorAcceptorRegionData
        ) const;

        //- Return list of deactivated (hole) cells
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
