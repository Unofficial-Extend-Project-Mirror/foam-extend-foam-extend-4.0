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
    standardPenalty

Description
    standard penalty method for normal contact law.
    Using a fixed penalty factor, force is applied to
    the contact faces based on the gap.

SourceFiles
    standardPenalty.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef standardPenalty_H
#define standardPenalty_H

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
                         Class standardPenalty Declaration
\*---------------------------------------------------------------------------*/

class standardPenalty
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
       scalar penaltyScale_;
       scalar relaxFac_;

       // current slave point tractions
       // IO to allow restarting
       //vectorIOField totalSlavePointTrac_;
       scalarField totalSlavePointTrac_;

       // contact iteration number
       label contactIterNum_;

       // oscillation correction factor
       // 1 -> no correction
       // 0 -> fully smoothed
       // values as high as 0.999 can help eliminate oscillations
       Switch oscillationCorr_;
       //scalar oscillationCorrFac_;
       int smoothingSteps_;

       // info frequency - how often to print to contact file
       int infoFreq_;

       // contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow copy construct
        standardPenalty(const standardPenalty&);

        //- Disallow default bitwise assignment
        void operator=(const standardPenalty&);



public:

    //- Runtime type information
    TypeName("standardPenalty");


    // Constructors

        //- Construct from dictionary
        standardPenalty
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

        virtual ~standardPenalty()
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
            // not used for standardPenalty
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
            // always zero for standardPenalty
            return slaveValueFrac_;
        }

        //- calculate slave point points
        void calcSlavePointPoints(const label slavePatchID);

        //- calculate penalty factor
        void calcPenaltyFactor();

        //- penalty factor
        virtual scalar penaltyFactor()
        {
            if (!penaltyFactorPtr_)
            {
                calcPenaltyFactor();
            }

            return *penaltyFactorPtr_;
        }

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
