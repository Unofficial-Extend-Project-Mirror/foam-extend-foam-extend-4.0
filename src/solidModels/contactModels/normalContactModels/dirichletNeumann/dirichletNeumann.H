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

       //- old slaveDispMag for under-relaxation
       scalarField oldSlaveDispMag_;

       //- slave pressure
       vectorField slavePressure_;
       vectorField oldSlavePressure_;

       //- if the face is currently in contact
       // IOList<bool> faceInContact_;
       //List<bool> faceInContact_;
       // fraction of face in contact
       scalarField touchFraction_;

       // slave valueFraction
       symmTensorField slaveValueFrac_;
       symmTensorField oldSlaveValueFrac_;

       // global and local slave point penetrations
       scalarField globalSlavePointPenetration_;
       scalarField slavePointPenetration_;

       // if the penetration is limited
       Switch limitPenetration_;
       scalar penetrationLimit_;

       // if the slave pressure is limited
       // to help convergence
       Switch limitPressure_;
       scalar pressureLimit_;

       // sometimes the contact has essentially converged
       // but slightly oscillates between two states indefinitely
       // so we can turn off contact corrections after the specified
       // number of corrections
       Switch settleContact_;
       label settleIterationNumber_;

       // correct pointDistToInter misses
       // pointPoints are needed
       Switch correctMissedVertices_;
       labelListList* slavePointPointsPtr_;

       // contact gap tolerance
       scalar contactGapTol_;

       // contact iteration number
       label contactIterNum_;

       // under-relaxation factor for slave displacement
       scalar relaxFactor_;

       // method used to calculate penetration distances (face or point)
       word distanceMethod_;

       // if Aitken's under-relaxation method is used for slave displacement
       Switch aitkenRelaxation_;
       label curTimeIndex_;
       label iCorr_;
       scalar aitkenRes0_;
       scalar aitkenTheta_;
       vectorField aitkenDelta_;
       vectorField aitkenDeltaPrevIter_;
       vectorField slaveDispPrevIter_;

       // oscillation correction factor
       // 1 -> no correction
       // 0 -> fully smoothed
       // values as high as 0.999 can help eliminate oscillations
       Switch oscillationCorr_;
       //scalar oscillationCorrFac_;
       int smoothingSteps_;

       // info frequency - how often to print to contact file
       int infoFreq_;

       // debug file
       //OFstream* debugFilePtr_;

       // contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow copy construct
        dirichletNeumann(const dirichletNeumann&);

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
        const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
        const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
        );


    // Destructor

        virtual ~dirichletNeumann()
        {}


    // Member Functions


        //- Correct contact law
        virtual void correct
    (
     //const PrimitivePatchInterpolation
     //< PrimitivePatch<face, List, pointField> >& slaveInterpolator,
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
      return slaveDisp_;
    }

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
      // always zero for dirichletNeumann
      return slavePressure_;
    }

        //- Return slave valueFraction
        virtual const symmTensorField& slaveValueFrac() const
        {
      return slaveValueFrac_;
    }

        //- calculate slave point points
        void calcSlavePointPoints(const label slavePatchID);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
