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
    Foam::repatchCoverage

Description
    Repatch boundary mesh modifier.  This modifier is used to move uncovered
    faces of a partially covered patch pair into a separate patch for
    master and slave side separately.

SourceFiles
    repatchCoverage.C

\*---------------------------------------------------------------------------*/

#ifndef repatchCoverage_H
#define repatchCoverage_H

#include "polyMeshModifier.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "ggiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class repatchCoverage Declaration
\*---------------------------------------------------------------------------*/

class repatchCoverage
:
    public polyMeshModifier
{
    // Private data

        //- Master covered patch ID
        polyPatchID masterCoveredPatchID_;

        //- Master uncovered patch ID
        polyPatchID masterUncoveredPatchID_;

        //- Slave covered patch ID
        polyPatchID slaveCoveredPatchID_;

        //- Slave uncovered patch ID
        polyPatchID slaveUncoveredPatchID_;

        //- Repatch threshold
        scalar repatchThreshold_;


        //- Demand-driven data

            //- Uncovered master face indicator
            mutable boolList uncMaster_;

            //- Uncovered slave face indicator
            mutable boolList uncSlave_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        repatchCoverage(const repatchCoverage&);

        //- Disallow default bitwise assignment
        void operator=(const repatchCoverage&);


        // Helper functions

            //- Check validity of construction data
            void checkDefinition();


public:

    //- Runtime type information
    TypeName("repatchCoverage");


    // Constructors

        //- Construct from components
        repatchCoverage
        (
            const word& name,
            const label index,
            const polyTopoChanger& mme,
            const word& masterCoveredPatchName,
            const word& masterUncoveredPatchName,
            const word& slaveCoveredPatchName,
            const word& slaveUncoveredPatchName,
            const scalar repatchThreshold
        );

        //- Construct from dictionary
        repatchCoverage
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    // Destructor

        virtual ~repatchCoverage();


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);


        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
