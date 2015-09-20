/*---------------------------------------------------------------------------* \
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
    Foam::uniformSet

Description

SourceFiles
    uniformSet.C

\*---------------------------------------------------------------------------*/

#ifndef uniformSet_H
#define uniformSet_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class passiveParticle;
template<class Type> class Particle;

/*---------------------------------------------------------------------------*\
                           Class uniformSet Declaration
\*---------------------------------------------------------------------------*/

class uniformSet
:
    public sampledSet
{
    // Private data

        //- Starting point
        point start_;

        //- End point
        point end_;

        //- Number of points
        label nPoints_;

    // Private Member Functions

        //- Calculates - starting at samplePt - the first sampling point
        //  on or after currentPt. smallDist is the tolerance used to compare
        //  positions. Returns false if end of samples reached.
        bool nextSample
        (
            const point& currentPt,
            const vector& offset,
            const scalar smallDist,
            point& samplePt,
            label& sampleI
        ) const;

        //- Samples from startTrackPt/CellI. Updates particle/samplePt/sampleI
        //  and puts
        //  samples in the DynamicLists. Returns false if end of all samples
        //  reached
        bool trackToBoundary
        (
            Particle<passiveParticle>& singleParticle,
            point& samplePt,
            label& sampleI,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Samples from start_ to end_. samplingSegments contains segmentNo
        //  for each sample.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("uniform");

    // Constructors

        //- Construct from components
        uniformSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& axis,
            const point& start,
            const point& end,
            const label nPoints
        );

        //- Construct from dictionary
        uniformSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~uniformSet();


    // Member Functions

        //- Get reference point
        virtual point getRefPoint(const List<point>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
