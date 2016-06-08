/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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
         //- whereas the contact law contains is the physcial law e.g. Coulomb
        dictionary frictionContactModelDict_;

         //- Friction law
        frictionLaw* frictionLawPtr_;

        //- Const reference to mesh
        const fvMesh& mesh_;

       //- slave shear displacement
       vectorField slaveDisp_;
       vectorField oldSlaveDisp_;

       //- Slip
       vectorField slip_;
       vectorField oldSlip_;
       const scalar slipTol_;

       //- Slave traction - always zero for dirichletNeumannFriction
       vectorField slaveTraction_;
       vectorField oldSlaveTraction_;
       vectorField slaveTractionForMaster_;
       vectorField oldSlaveTractionForMaster_;

       // Slave valueFraction
       symmTensorField slaveValueFrac_;
       symmTensorField oldSlaveValueFrac_;

       // Displacement under-relaxation factor
       scalar urDisp_;

       // Traction under-relaxation factor
       scalar urTrac_;

       // contact iteration number
       label contactIterNum_;

       // info frequency - how often to print to contact file
       int infoFreq_;

       // Oscillation correction switch
       Switch oscillationCorr_;

       // Number of oscillation correction smoothing steps
       int smoothingSteps_;

       // Old stickSlip field so we can see how it is changing
       scalarField oldStickSlip_;

       // Contact info file
       OFstream* contactFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //dirichletNeumannFriction(const dirichletNeumannFriction&);

        //- Disallow default bitwise assignment
        void operator=(const dirichletNeumannFriction&);

public:

    //- Runtime type information
    TypeName("dirichletNeumann");


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

        //- Construct as copy
        dirichletNeumannFriction(const dirichletNeumannFriction&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const
        {
            return
                autoPtr<frictionContactModel>
                (
                    new dirichletNeumannFriction(*this)
                );
        }

    // Destructor

        virtual ~dirichletNeumannFriction()
        {}


    // Member Functions

        //- Correct the rheological model
        virtual void correct
        (
            const vectorField& slavePressure,
            const primitiveFacePatch& masterFaceZonePatch,
            const primitiveFacePatch& slaveFaceZonePatch,
            const intersection::algorithm alg,
            const intersection::direction dir,
            const word interpolationMethod,
            const word fieldName,
            const Switch orthotropic,
            const nonLinearGeometry::nonLinearType nonLinear,
            const vectorField& slaveFaceNormals
        );

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slave friction traction for master. This has the same size as
        // the slave patch
        virtual const vectorField& slaveTractionForMaster() const
        {
            return slaveTractionForMaster_;
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

        //- Return slip field
        virtual const vectorField& slip() const
        {
            return slip_;
        }

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
