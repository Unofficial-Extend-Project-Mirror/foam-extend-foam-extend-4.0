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
    constitutiveModel

Description
    Class in charge of constitutive laws and mechanical properties

SourceFiles
    constitutiveModel.C

Author
    Philip Cardiff UCD
    philip.cardiff@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef constitutiveModel_H
#define constitutiveModel_H

#include "IOdictionary.H"
#include "plasticityStressReturn.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "rheologyLaw.H"
#include "cohesiveLaw.H"
#include "Switch.H"
#include "fvc.H"
#include "IOReferencer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declaration
class solidInterface;

/*---------------------------------------------------------------------------*\
                        Class constitutiveModel Declaration
\*---------------------------------------------------------------------------*/

class constitutiveModel
:
  public IOdictionary
{
    // Private data

        //- Reference to stress field
        const volSymmTensorField& sigma_;

        //- Rheology law
        autoPtr<rheologyLaw> rheologyLawPtr_;

        //- Cohesive dictionary and law
        IOdictionary* cohesiveDictPtr_;
        autoPtr<cohesiveLaw> cohesiveLawPtr_;

        //- Plane stress
        Switch planeStress_;

        //- Run-time selectable solidInterface method for correct
        //  discretisation on bi-material interfaces
        IOReferencer<solidInterface>* solidInterfacePtr_;

        // if solidInterface is on/off
        bool solidInterfaceActive_;

        //- Run-time selectable method to calculate DEpsilonP
        //- eg radial return
        autoPtr<plasticityStressReturn> plasticityStressReturnPtr_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveModel(const constitutiveModel&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveModel&);


public:

    //- Runtime type information
    TypeName("constitutiveModel");


    // Constructors

        //- Construct from components
        constitutiveModel
        (
            const volSymmTensorField& sigma,
            const volVectorField& U
        );


    // Destructor

        virtual ~constitutiveModel();


    // Member Functions

        // Access

        //- If plasticity is active
        bool plasticityActive() const
        {
            if (plasticityStressReturnPtr_.valid())
            {
                return plasticityStressReturnPtr_->plasticityActive();
            }

            return false;
        }

        //- If visco effects are active
        bool viscoActive() const
        {
            return rheologyLawPtr_->viscoActive();
        }

        //- Return reference to stress field
        const volSymmTensorField& sigma() const
        {
            return sigma_;
        }

        //- Return true for plane stress
        const Switch& planeStress() const
        {
            return planeStress_;
        }

        //- Return rheology law
        const rheologyLaw& law() const
        {
            return rheologyLawPtr_();
        }

        //- Return cohesive law
        const cohesiveLaw& cohLaw() const
        {
            return cohesiveLawPtr_();
        }

        //- Return density
        tmp<volScalarField> rho() const
        {
            return rheologyLawPtr_->rho();
        }

        tmp<volScalarField> rho(scalar t) const
        {
            return rheologyLawPtr_->rho(t);
        }

        //- Return first Lame's coefficient
        tmp<volScalarField> mu() const;
        tmp<volScalarField> mu(scalar t) const;
        tmp<volScalarField> mu(const volScalarField& epsilonEq) const;
        tmp<surfaceScalarField> muf() const;

        //- Return second Lame's coefficient
        tmp<volScalarField> lambda() const;
        tmp<volScalarField> lambda(scalar t) const;
        tmp<volScalarField> lambda(const volScalarField& epsilonEq) const;
        tmp<surfaceScalarField> lambdaf() const;

        //- Return threeK
        tmp<volScalarField> threeK() const;
        tmp<surfaceScalarField> threeKf() const;

        //- Return yield stress
        tmp<volScalarField> sigmaY() const
        {
            return rheologyLawPtr_->sigmaY();
        }

        //- Return yield stress
        scalar sigmaY(const scalar epsilonPEq, const label cellID) const
        {
            return rheologyLawPtr_->sigmaY(epsilonPEq, cellID);
        }

        //- Return plastic modulus
        tmp<volScalarField> Ep() const
        {
            return rheologyLawPtr_->Ep();
        }
        tmp<volScalarField> Ep(const volScalarField& epsilonEq) const
        {
            return rheologyLawPtr_->Ep(epsilonEq);
        }

        //- Return plastic strain increment
        const volSymmTensorField& DEpsilonP() const;

        //- orthotropic tensorial diffusivity
        //- Note: this is not the bulk modulus
        tmp<volDiagTensorField> K() const;
        tmp<surfaceDiagTensorField> Kf() const;

        //- fourth order elastic stiffness tensor
        tmp<volSymmTensor4thOrderField> C() const;
        tmp<surfaceSymmTensor4thOrderField> Cf() const;

        //- Correct plastic strain increment
        virtual void correct();

        //- Return reference to solidInterface
        virtual const solidInterface& solInterface() const
        {
            return solidInterfacePtr_->operator()();
        }

        //- Return access to solidInterface
        virtual solidInterface& solInterface()
        {
            return solidInterfacePtr_->operator()();
        }

        //- Return true if solidInterface is active
        virtual bool solidInterfaceActive()
        {
            return solidInterfaceActive_;
        }

        //- Update yield stress
        void updateYieldStress();

        //- Read plasticityProperties dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
