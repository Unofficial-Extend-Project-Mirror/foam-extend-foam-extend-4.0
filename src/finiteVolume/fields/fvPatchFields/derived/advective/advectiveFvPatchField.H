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
    Foam::advectiveFvPatchField

Description
    Advective outflow boundary condition based on solving DDt(psi, U) = 0
    at the boundary.

    The standard (Euler, backward, CrankNicolson) time schemes are
    supported.  Additionally an optional mechanism to relax the value at
    the boundary to a specified far-field value is provided which is
    switched on by specifying the relaxation length-scale lInf and the
    far-field value fieldInf.

    The flow/wave speed at the outlet is provided by the virtual function
    advectionSpeed() the default implementation of which requires the name of
    flux field a the outlet (phi) and optionally the density (rho) if the
    mass-flux rather than the volumetric-flux is given.
    @verbatim
        outlet
        {
            type                    advective;
            phi                     phi;
            // rho                  rho; // Not needed, phi volumetric
            // fieldInf             1e5; // Optional
            // lInf                 0.1; // Optional

            inletOutlet             true; // Correct for back-flow
            correctSupercritical    true;
        }
    @endverbatim

    The flow/wave speed at the outlet can be changed by deriving a specialised
    BC from this class and overriding advectionSpeed() e.g. in
    waveTransmissiveFvPatchField the advectionSpeed() calculates and returns
    the flow-speed plus the acoustic wave speed creating an acoustic wave
    transmissive boundary condition.

    Steady-state handling has been added by Hrvoje Jasak.
    Supercritical outlet handling has been added by Hrvoje Jasak.
    HJ, 27/Oct/2009

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    advectiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef advectiveFvPatchField_H
#define advectiveFvPatchField_H

#include "mixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class advectiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class advectiveFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Field value of the far-field
        Type fieldInf_;

        //- Relaxation length-scale
        scalar lInf_;

        //- Inlet-outlet treatment
        Switch inletOutlet_;

        //- Supercritical correction treatment
        Switch correctSupercritical_;


public:

    //- Runtime type information
    TypeName("advective");


    // Constructors

        //- Construct from patch and internal field
        advectiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        advectiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given advectiveFvPatchField
        //  onto a new patch
        advectiveFvPatchField
        (
            const advectiveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        advectiveFvPatchField
        (
            const advectiveFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new advectiveFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        advectiveFvPatchField
        (
            const advectiveFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new advectiveFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the field at infinity
            const Type& fieldInf() const
            {
                return fieldInf_;
            }

            //- Return reference to the field at infinity to allow adjustment
            Type& fieldInf()
            {
                return fieldInf_;
            }

            //- Return the relaxation length-scale
            scalar lInf() const
            {
                return lInf_;
            }

            //- Return reference to the relaxation length-scale
            //  to allow adjustment
            scalar& lInf()
            {
                return lInf_;
            }

            //- Return inlet-outlet treatment switch
            bool inletOutlet() const
            {
                return inletOutlet_;
            }

            //- Return reference to inlet-outlet treatment switch
            //  to allow adjustment
            bool& inletOutlet()
            {
                return inletOutlet_;
            }

            //- Return supercritical outlet treatment switch
            bool correctSupercritical() const
            {
                return correctSupercritical_;
            }

            //- Return reference to supercritical treatment switch
            //  to allow adjustment
            bool& correctSupercritical()
            {
                return correctSupercritical_;
            }


        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;

            //- Calculate and return the supercritical switch at the boundary
            //  Supercritical = 1 converts the outlet to zeroGradient
            //  Supercritical = 0 no correction
            virtual tmp<scalarField> supercritical() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "advectiveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
