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
    lengthScaleEstimator

Description
    Utility class used to provide length-scale estimates at various
    mesh locations. These estimates are based on specified boundary
    conditions provided through dictionary entries.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    lengthScaleEstimator.C
    lengthScaleEstimatorI.H

\*---------------------------------------------------------------------------*/

#ifndef lengthScaleEstimator_H
#define lengthScaleEstimator_H

#include "polyMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class lengthScaleEstimator Declaration
\*---------------------------------------------------------------------------*/

class lengthScaleEstimator
{
    // Private data

        //- Const reference to polyMesh
        const polyMesh& mesh_;

        //- Edge bisection/collapse criteria
        scalar ratioMin_;
        scalar ratioMax_;
        scalar growthFactor_;
        scalar minLengthScale_;
        scalar maxLengthScale_;
        scalar curvatureDeviation_;

        //- Specific to proximity-based refinement
        label spatialRes_;
        boundBox proxBoundBox_;
        labelListList proximityBins_;

        //- Specific to mesh-slicing operations
        scalar sliceThreshold_;
        label sliceHoldOff_;
        List<boundBox> sliceBoxes_;

        // Buffers for parallel length-scale calculations
        labelListList sendLblBuffer_;
        labelListList recvLblBuffer_;
        labelListList sendLvlBuffer_;
        labelListList recvLvlBuffer_;
        scalarListList sendSclBuffer_;
        scalarListList recvSclBuffer_;

        //- Sub-dictionary which specifies
        //  fixed length-scales for patches
        dictionary fixedPatches_;

        //- Sub-dictionary which specifies
        //  floating length-scales for patches
        dictionary freePatches_;

        //- Sub-dictionary which specifies
        //  patches with curvature-based length-scale
        dictionary curvaturePatches_;

        //- Sub-dictionary which specifies
        //  patches with proximity-based length-scale
        dictionary proximityPatches_;

        //- Patches for which edge-refinements
        //  are to be avoided
        labelList noModPatchIDs_;

        //- Field-based refinement
        word field_;
        scalar fieldLength_;
        scalar lowerRefineLevel_;
        scalar upperRefineLevel_;

        //- Specify limits for refinement criteria
        scalar meanScale_;
        label maxRefineLevel_;

    // Private Member Functions

        // Check for legitimacy of patches
        void checkPatches(const wordList& patchList) const;

        // Prepare for proximity-based refinement, if necessary
        void prepareProximityPatches();

        // Perform spatial hashing on a set of points
        static void spatialHash
        (
            const pointField& pointLocations,
            const labelList& pointIndices,
            const boundBox& box,
            const label resolution,
            labelListList& bins
        );

        // Send length-scale info across processors
        void writeLengthScaleInfo
        (
            const labelList& cellLevels,
            const scalarList& lengthScale
        );

        // Receive length-scale info across processors
        void readLengthScaleInfo
        (
            const label level,
            label& visitedCells,
            labelList& cellLevels,
            UList<scalar>& lengthScale,
            labelHashSet& levelCells
        ) const;

public:

    // Declare the name of the class and its debug switch
        TypeName("lengthScaleEstimator");

    // Constructors

        //- Construct from polyMesh and dictionary
        explicit lengthScaleEstimator
        (
            const polyMesh& mesh
        );

    // Destructor

        virtual ~lengthScaleEstimator();

    // Member Functions

        //- Read edge refinement options from the dictionary
        void readRefinementOptions
        (
            const dictionary& refineDict,
            bool reRead = false,
            bool mandatory = false
        );

        //- Set explicitly coupled patch information
        void setCoupledPatches
        (
            const dictionary& coupledPatches
        );

        //- Calculate the length scale field
        void calculateLengthScale(UList<scalar>& lengthScale);

        //- Return refinement criteria
        inline scalar ratioMin() const;
        inline scalar ratioMax() const;
        inline scalar growthFactor() const;

        //- Limit length scale for surface-edges
        inline void limitScale(scalar& scale) const;

        //- Check if a particular patch is free-floating
        //  (i.e., gets its length-scale from the interior)
        inline bool isFreePatch(const label pIndex) const;

        //- Check if a particular patch is flagged
        //  for proximity-based refinement
        inline bool isProximityPatch(const label pIndex) const;

        //- Check if a particular patch is flagged
        //  for curvature-based refinement
        inline bool isCurvaturePatch(const label pIndex) const;

        //- Return reference curvature deviation
        inline scalar curvatureDeviation() const
        {
            return curvatureDeviation_;
        }

        //- Check whether a particular point is too close
        //  to a previous mesh slice location
        inline bool checkOldSlices(const vector& gCentre) const;

        //- Add a boundBox to the existing set of sliceBoxes
        inline void appendBox(const boundBox& bBox);

        //- Clear the list of sliceBoxes
        inline void clearBoxes();

        //- Set an initial hold-off value
        inline void setHoldOff(const label hVal)
        {
            sliceHoldOff_ = hVal;
        }

        //- Return the current holdOff value
        inline label holdOff() const
        {
            return sliceHoldOff_;
        }

        //- Decrement the current holdOff value
        inline void decrementHoldOff()
        {
            sliceHoldOff_--;
        }

        //- Check whether a particular patch permits refinement
        inline bool checkRefinementPatch
        (
            const label pIndex
        ) const;

        //- Return the appropriate length-scale for boundary face
        inline scalar fixedLengthScale
        (
            const label fIndex,
            const label pIndex,
            bool usePolyMesh = false
        ) const;

        //- Test for proximity to patch faces
        inline bool testProximity
        (
            const vector& gCentre,
            const vector& gNormal,
            const scalar testStep,
            label& proxFace,
            scalar& proxDistance
        ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lengthScaleEstimatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
