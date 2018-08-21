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
    Foam::iterationData

Description
    Representation of the data type needed to track the history of the
    iterative process in which the most appropriate donor/acceptor pairs are
    sought.

    Holds and manages following data:

    - fringeHoles_: list of fringe holes from one iteration. In the first
      iteration fringe holes are hole cells generated from the hole cutting
      process. In the following iterations fringeHoles_ represents holes
      from the hole cutting process and discarded acceptors from previous
      iterations.
    - donorAcceptorPairs_: donorAcceptorList which holds donor/acceptor pairs
      from one iteration. This list is not filtered.
    - suitability_: average suitability fraction from one iteration. If this
      value is higher, donor/acceptor pairs are more suitable.
    - iteration_: iteration number which is needed for calculating slope using
      least square method.

SourceFiles
    iterationData.C

Author
    Josip Katavic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef iterationData_H
#define iterationData_H

#include "donorAcceptorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class iterationData Declaration
\*---------------------------------------------------------------------------*/

class iterationData
{
    // Private data

        //- Pre-final donor/acceptor pairs
        const donorAcceptorList donorAcceptorPairs_;

        //- List of fringe holes
        const labelList fringeHoles_;

        //- Average suitability fraction for all donor/acceptor pairs
        const scalar suitability_;

        //- Iteration number
        const label iteration_;


public:

    // Constructors

        //- Construct from Foam::adaptiveOverlap::updateIteration data by
        //  transfering the ownership from parameters to data members
        iterationData
        (
            const Xfer<donorAcceptorList>& donorAcceptorPairs,
            const Xfer<labelList>& fringeHoles,
            const scalar averageSuitability,
            const label iteration
        )
        :
            donorAcceptorPairs_(donorAcceptorPairs),
            fringeHoles_(fringeHoles),
            suitability_(averageSuitability),
            iteration_(iteration)
        {}

        // Note: There is no need for 'the big three' implementation because
        // this class does not hold pointers as data members, so use implicit
        // copy constructor, assignment operator and destructor.


    // Access

        //- Return reference to donor/acceptor pairs
        const donorAcceptorList& donorAcceptorPairs() const
        {
            return donorAcceptorPairs_;
        }

        //- Return reference to fringe holes
        const labelList& fringeHoles() const
        {
            return fringeHoles_;
        }

        //- Return reference to suitability fraction
        const scalar& suitability() const
        {
            return suitability_;
        }

        //- Return reference to iteration number
        const label& iteration() const
        {
            return iteration_;
        }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
