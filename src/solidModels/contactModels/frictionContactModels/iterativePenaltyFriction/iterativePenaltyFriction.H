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
    iterativePenaltyFriction

Description
    iterativePenaltyFriction contact friction law
    Essentially Coulomb friction
    Increments of friction traction are added to the slave faces until the
    the increment of relative motion is zero or until the friction reaches
    fricCoeff*pressure

SourceFiles
    iterativePenaltyFriction.C

Author
    Philip Cardiff UCD
    Matthias Schnaithmann Universitat Stuttgart

\*---------------------------------------------------------------------------*/

#ifndef iterativePenaltyFriction_H
#define iterativePenaltyFriction_H

#include "frictionContactModel.H"
#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class iterativePenaltyFriction Declaration
\*---------------------------------------------------------------------------*/

class iterativePenaltyFriction
:
    public frictionContactModel
{
    // Private data

         //- the contact model contains the contact procedure
         //- whereas the contact law contains is the physcial law eg Coulomb
        dictionary frictionContactModelDict_;

         //- friction law
        frictionLaw* frictionLawPtr_;

         //- reference to mesh
        const fvMesh& mesh_;

       //- slave traction - not used for iterativePenaltyFriction
       vectorField slaveDisp_;

       //- slave traction
       vectorField slaveTraction_;

       //- slave friction value fraction
       symmTensorField slaveValueFrac_;

       // frictional penalty factor and scaling factor
       scalar* frictionPenaltyFactorPtr_;
       scalar frictionPenaltyScale_;

       // contact iteration number
       label contactIterNum_;

       // info frequency - how often to print to contact file
       int infoFreq_;

       // oscillation correction factor
       // 1 -> no correction
       // 0 -> fully smoothed
       // values as high as 0.999 can help eliminate oscillations
       Switch oscillationCorr_;
       scalar oscillationCorrFac_;

       //- Contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        iterativePenaltyFriction(const iterativePenaltyFriction&);

        //- Disallow default bitwise assignment
        void operator=(const iterativePenaltyFriction&);

public:

    //- Runtime type information
    TypeName("iterativePenalty");

    // Static data members


    // Constructors

        //- Construct from dictionary
        iterativePenaltyFriction
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID
        );


    // Destructor

        virtual ~iterativePenaltyFriction()
        {}


    // Member Functions

        //- Correct the rheological model
        virtual void correct
        (
            const vectorField& slavePressure,
            const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
            const intersection::algorithm alg,
            const intersection::direction dir,
            const word interpolationMethod,
            const word fieldName,
            const Switch orthotropic,
            const word nonLinear,
            const vectorField& slaveFaceNormals
            );

        //- Return slave friction displacement
        virtual const vectorField& slaveDisp() const
        {
            return slaveDisp_;
        }

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slave friction value fraction
        virtual const symmTensorField& slaveValueFrac() const
        {
            return slaveValueFrac_;
        }

        //- calculate penalty factor
        void calcFrictionPenaltyFactor();

        //- Return penalty factor
        virtual scalar frictionPenaltyFactor()
        {
            if (!frictionPenaltyFactorPtr_)
            {
                calcFrictionPenaltyFactor();
            }

            return *frictionPenaltyFactorPtr_;
        }

        //- Correct oscillations in the contact shear tractions
        virtual void correctOscillations
        (
            vectorField& slaveTraction,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
            );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
