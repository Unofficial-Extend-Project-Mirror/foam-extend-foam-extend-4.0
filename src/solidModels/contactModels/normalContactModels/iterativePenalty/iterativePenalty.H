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
    iterativePenalty

Description
    Iterative penalty method for normal contact law.
    Using an inital penalty factor, force is incrementally added to
    the contact patches until they are within the contact tolerance.

SourceFiles
    iterativePenalty.C

Author
    Philip Cardiff UCD
    Cardiff et al. 2012. Development of a finite volume contact solver based
    on the penalty method. Computational Materials Science.

\*---------------------------------------------------------------------------*/

#ifndef iterativePenalty_H
#define iterativePenalty_H

#include "normalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "patchToPatchInterpolation.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class iterativePenalty Declaration
\*---------------------------------------------------------------------------*/

class iterativePenalty
:
    public normalContactModel
{
    // Private data

         //- Contact law coefficients
        dictionary normalContactModelDict_;

         //- reference to mesh
        const fvMesh& mesh_;

       //- slave normal displacement - not used by iterative penalty
       vectorField slaveDisp_;

       //- slave pressure
       vectorField slavePressure_;

       // slave valueFraction - always zero for iterative penalty
       symmTensorField slaveValueFrac_;

       // global and local slave point penetrations
       scalarField globalSlavePointPenetration_;
       scalarField slavePointPenetration_;

       // if the penetration is limited
       Switch limitPenetration_;
       scalar penetrationLimit_;

       // correct pointDistToInter misses
       // pointPoints are needed
       Switch correctMissedVertices_;
       labelListList* slavePointPointsPtr_;

       // penalty factor and returning penalty factor
       scalar* penaltyFactorPtr_;
       scalar* returnPenaltyFactorPtr_;
       scalar penaltyScale_;
       scalar returnScale_;

       // current slave point tractions
       // IO to allow restarting
       vectorIOField totalSlavePointTrac_;

       // contact gap tolerance
       scalar contactGapTol_;

       // contact iteration number
       label contactIterNum_;

       // oscillation correction factor
       // 1 -> no correction
       // 0 -> fully smoothed
       // values as high as 0.999 can help eliminate oscillations
       Switch oscillationCorr_;
       scalar oscillationCorrFac_;

       // contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow copy construct
        iterativePenalty(const iterativePenalty&);

        //- Disallow default bitwise assignment
        void operator=(const iterativePenalty&);



public:

    //- Runtime type information
    TypeName("iterativePenalty");


    // Constructors

        //- Construct from dictionary
        iterativePenalty
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
        const label masterPatchID,
        const label slavePatchID,
        const label masterFaceZoneID,
        const label slaveFaceZoneID,
        const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
        const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
        );


    // Destructor

        virtual ~iterativePenalty()
        {}


    // Member Functions


        //- Correct contact model
        virtual void correct
    (
     const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
     const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
     const intersection::algorithm alg,
     const intersection::direction dir,
     word fieldName,
     const Switch orthotropic,
     const word nonLinear,
     vectorField& slaveFaceNormals,
     GGIInterpolation< PrimitivePatch< face, List, pointField >,
               PrimitivePatch< face, List, pointField >
               >* ggiInterpolatorPtr
     );


        //- Return slave displacement
        virtual const vectorField& slaveDisp() const
        {
      // not used for iterativePenalty
      return slaveDisp_;
    }

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
      return slavePressure_;
    }

        //- Return slave valueFraction
        virtual const symmTensorField& slaveValueFrac() const
        {
      // always zero for iterativePenalty
      return slaveValueFrac_;
    }

        //- calculate slave point points
        void calcSlavePointPoints(const label slavePatchID);

        //- calculate penalty factor
        void calcPenaltyFactor();

        //- Return penalty factor
        virtual scalar penaltyFactor()
        {
            if (!penaltyFactorPtr_)
            {
                calcPenaltyFactor();
            }

            return *penaltyFactorPtr_;
        }

        //- Return penalty factor
        // we reduce the traction at a different rate (probably more quickly)
        // on points which come out of contact
        virtual scalar returnPenaltyFactor()
        {
            if (!returnPenaltyFactorPtr_)
            {
                returnPenaltyFactorPtr_ =
                    new scalar(returnScale_*penaltyFactor());
            }

            return *returnPenaltyFactorPtr_;
    }

        //- Correct oscillations in the contact pressure
        virtual void correctOscillations
        (
            vectorField& slavePressure,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
            );

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
