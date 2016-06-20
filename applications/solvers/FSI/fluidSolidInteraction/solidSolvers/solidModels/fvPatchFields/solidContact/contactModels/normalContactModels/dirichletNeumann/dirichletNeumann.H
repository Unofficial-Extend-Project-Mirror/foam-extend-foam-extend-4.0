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
    dirichletNeumann

Description
    Dirichlet-Neumann approach for normal contact model.
    The slave is treated as a Dirichlet condition (displacement)
    and the master is treated as a Neumann condition (traction).

    The slave is displaced back to the master surface and the
    traction induced on the slave is passed to the master.
    The values are under-relaxed.
    This is the same approach as employed in many partitioned FSI
    procedures.

SourceFiles
    dirichletNeumann.C

Author
    Philip Cardiff UCD
    Based on Jasak and Weller 2000 contacting linear elastic bodies
    and contactStressedFoam, newContactStressedFoam and
    Vlado Tropsa FOAM contact solvers.

\*---------------------------------------------------------------------------*/

#ifndef dirichletNeumann_H
#define dirichletNeumann_H

#include "normalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "patchToPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "IOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dirichletNeumann Declaration
\*---------------------------------------------------------------------------*/

class dirichletNeumann
:
    public normalContactModel
{
    // Private data

         //- Contact model coefficients
        dictionary normalContactModelDict_;

         //- reference to mesh
        const fvMesh& mesh_;

       //- slave normal displacement
       vectorField slaveDisp_;

       //- slave pressure
       vectorField slavePressure_;
       vectorField oldSlavePressure_;
       vectorField zeroSlavePressure_;

       // Contact residual
       scalar residual_;
       scalar maxNorm_;
       scalar averagePenetration_;

       // Nonlinear type
       word nonLinearType_;

       //- if the face is currently in contact
       // fraction of face in contact
       //scalarField touchFraction_;
       // slave patch area in contact
       // 0.0 = face not in contact
       // 1.0 = face is fully in contact
       scalarField areaInContact_;

       // slave valueFraction
       symmTensorField slaveValueFrac_;
       symmTensorField oldSlaveValueFrac_;

       // if the penetration is limited
       Switch limitPenetration_;
       scalar penetrationLimit_;

       // if the slave pressure is limited
       // to help convergence
       Switch limitPressure_;
       scalar pressureLimit_;

       // correct pointDistToInter misses
       // pointPoints are needed
       Switch correctMissedVertices_;
       labelListList* slavePointPointsPtr_;

       // contact gap tolerance
       scalar contactGapTol_;

       // contact iteration number
       label contactIterNum_;

       //- Under-relaxation factors for displacement and traction
       scalar urDisp_;
       scalar urTrac_;

       //- Method used to calculate penetration distances (face or point)
       word distanceMethod_;

       // Current time index and outer iteration
       label curTimeIndex_;
       label iCorr_;
       //vectorField slaveDispPrevIter_;

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

        // //- Calculate slave face normals
        // void calculateSlaveFaceNormals
        // (
        //     vectorField& slaveFaceNormals,
        //     const nonLinearGeometry::nonLinearType nonLinear,
        //   const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
        //     const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
        //     ExtendedGGIInterpolation<
        //         PrimitivePatch< face, Foam::List, pointField >,
        //         PrimitivePatch< face, Foam::List, pointField >
        //         >* ggiInterpolatorPtr
        // );

        //- Disallow copy construct
        //dirichletNeumann(const dirichletNeumann&);

        //- Disallow default bitwise assignment
        void operator=(const dirichletNeumann&);



public:

    //- Runtime type information
    TypeName("dirichletNeumann");


    // Constructors

        dirichletNeumann
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID,
            const primitiveFacePatch& masterFaceZonePatch,
            const primitiveFacePatch& slaveFaceZonePatch
        );

        //- Construct as copy
        dirichletNeumann(const dirichletNeumann&);

        //- Construct and return a clone
        virtual autoPtr<normalContactModel> clone() const
        {
            return autoPtr<normalContactModel>(new dirichletNeumann(*this));
        }

    // Destructor

        virtual ~dirichletNeumann()
        {}


    // Member Functions


        //- Correct contact law
        virtual void correct
        (
            const primitiveFacePatch& masterFaceZonePatch,
            const primitiveFacePatch& slaveFaceZonePatch,
            const intersection::algorithm alg,
            const intersection::direction dir,
            word fieldName,
            const Switch orthotropic,
            const nonLinearGeometry::nonLinearType nonLinear,
            vectorField& slaveFaceNormals,
            const extendedGgiZoneInterpolation& zoneToZone
         );


        //- Return slave displacement
        virtual const vectorField& slaveDisp() const
        {
            return slaveDisp_;
        }

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            // always zero for dirichletNeumann
            return zeroSlavePressure_;
        }

        //- Return slave valueFraction
        virtual const symmTensorField& slaveValueFrac() const
        {
            return slaveValueFrac_;
        }

        //- calculate slave point points
        void calcSlavePointPoints(const label slavePatchID);

        //- Return contact residual
        virtual scalar residual();

        //- Return area in contact of slave patch
        virtual const scalarField& areaInContact() const
        {
            return areaInContact_;
        };

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
