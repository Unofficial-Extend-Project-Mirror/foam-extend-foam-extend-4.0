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
    dirichletNeumannFriction

Description
    dirichletNeumannFriction contact friction law
    Essentially Coulomb friction.
    A shear displacement is applied to the slave and the generated
    shear traction is passed to the master
    i.e. the slave is a Dirichlet condition and the master is a
    Neumann condition.

SourceFiles
    dirichletNeumannFriction.C

Author
    Philip Cardiff UCD
    similar to approach of newContactStressedFoam

\*---------------------------------------------------------------------------*/

#ifndef dirichletNeumannFriction_H
#define dirichletNeumannFriction_H

#include "frictionContactModel.H"
#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dirichletNeumannFriction Declaration
\*---------------------------------------------------------------------------*/

class dirichletNeumannFriction
:
    public frictionContactModel
{
    // Private data

         //- Contact model coefficients
         //- the contact model contains the contact procedure
         //- whereas the contact law contains is the physcial law eg Coulomb
        dictionary frictionContactModelDict_;

         //- friction law
        frictionLaw* frictionLawPtr_;

         //- reference to mesh
        const fvMesh& mesh_;

       //- slave shear displacement
       vectorField slaveDisp_;
       vectorField oldSlaveDisp_;

       //- old slip
       vectorField oldSlip_;

       //- slave traction - always zero for dirichletNeumannFriction
       vectorField slaveTraction_;
       vectorField oldSlaveTraction_;

       // slave valueFraction
       symmTensorField slaveValueFrac_;
       symmTensorField oldSlaveValueFrac_;

       // under-relaxation factor
       scalar relaxationFactor_;

       // contact iteration number
       label contactIterNum_;

       // info frequency - how often to print to contact file
       int infoFreq_;

       // oscillation correction factor
       // 1 -> no correction
       // 0 -> fully smoothed
       // values as high as 0.999 can help eliminate oscillations
       Switch oscillationCorr_;
       //scalar oscillationCorrFac_;
       int smoothingSteps_;

       // old stickSlip field so we can see how it is changing
       scalarField oldStickSlip_;

       // contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dirichletNeumannFriction(const dirichletNeumannFriction&);

        //- Disallow default bitwise assignment
        void operator=(const dirichletNeumannFriction&);

public:

    //- Runtime type information
    TypeName("dirichletNeumann");

    // Static data members


    // Constructors

        //- Construct from dictionary
        dirichletNeumannFriction
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

        virtual ~dirichletNeumannFriction()
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

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slave friction displacement
        virtual const vectorField& slaveDisp() const
        {
            return slaveDisp_;
        }

        //- Return slave valueFraction
        virtual const symmTensorField& slaveValueFrac() const
        {
            return slaveValueFrac_;
        }

        //- Correct oscillations in the contact shear tractions/displacements
        // virtual void correctOscillations
    // (
    //  const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
    //  );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
