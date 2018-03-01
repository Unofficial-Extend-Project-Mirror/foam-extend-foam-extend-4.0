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
    donorBasedLayeredFringeFringe

Description
    Fringe algorithm based on layered approach, where the acceptors of the
    region are defined as neighbours of the donors toward the interior of the
    region.

    This is achieved by waiting for the acceptors of the donor region to be
    found and then defining acceptors.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    donorBasedLayeredFringeFringe.C

\*---------------------------------------------------------------------------*/

#ifndef donorBasedLayeredFringeFringe_H
#define donorBasedLayeredFringeFringe_H

#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class donorBasedLayeredFringeFringe Declaration
\*---------------------------------------------------------------------------*/

class donorBasedLayeredFringeFringe
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

        //- Number of layers from donors to acceptors of this region
        label nLayers_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        donorBasedLayeredFringeFringe(const donorBasedLayeredFringeFringe&);

        //- Disallow default bitwise assignment
        void operator=(const donorBasedLayeredFringeFringe&);


        // Calculate hole and acceptor addressing
        void calcAddressing() const;

        // Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("donorBasedLayeredFringe");


    // Constructors

        //- Construct from dictionary
        donorBasedLayeredFringeFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    // Destructor

        virtual ~donorBasedLayeredFringeFringe();


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
