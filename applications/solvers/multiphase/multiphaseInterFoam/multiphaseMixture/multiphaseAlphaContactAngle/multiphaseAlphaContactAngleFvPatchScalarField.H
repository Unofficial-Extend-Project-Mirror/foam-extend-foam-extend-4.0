/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::multiphaseAlphaContactAngleFvPatchScalarField

Description
    Contact-angle boundary condition for multi-phase interface-capturing
    simulations.  Used in conjuction with multiphaseMixture.

SourceFiles
    multiphaseAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseAlphaContactAngleFvPatchScalarField_H
#define multiphaseAlphaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"
#include "multiphaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class multiphaseAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class multiphaseAlphaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    class interfaceThetaProps
    {
        //- Equilibrium contact angle
        scalar theta0_;

        //- Dynamic contact angle velocity scale
        scalar uTheta_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receeding contact angle
        scalar thetaR_;


    public:

        // Constructors
            interfaceThetaProps()
            {}

            interfaceThetaProps(Istream&);


        // Member functions

            //- Return the equilibrium contact angle theta0
            scalar theta0(bool matched=true) const
            {
                if (matched) return theta0_;
                else return 180.0 - theta0_;
            }

            //- Return the dynamic contact angle velocity scale
            scalar uTheta() const
            {
                return uTheta_;
            }

            //- Return the limiting advancing contact angle
            scalar thetaA(bool matched=true) const
            {
                if (matched) return thetaA_;
                else return 180.0 - thetaA_;
            }

            //- Return the limiting receeding contact angle
            scalar thetaR(bool matched=true) const
            {
                if (matched) return thetaR_;
                else return 180.0 - thetaR_;
            }


        // IO functions

            friend Istream& operator>>(Istream&, interfaceThetaProps&);
            friend Ostream& operator<<(Ostream&, const interfaceThetaProps&);
    };

    typedef HashTable
    <
        interfaceThetaProps,
        multiphaseMixture::interfacePair,
        multiphaseMixture::interfacePair::hash
    > thetaPropsTable;


private:

    // Private data

        thetaPropsTable thetaProps_;


public:

    //- Runtime type information
    TypeName("multiphaseAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        multiphaseAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        multiphaseAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  multiphaseAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        multiphaseAlphaContactAngleFvPatchScalarField
        (
            const multiphaseAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new multiphaseAlphaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        multiphaseAlphaContactAngleFvPatchScalarField
        (
            const multiphaseAlphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new multiphaseAlphaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the contact angle properties
        const thetaPropsTable& thetaProps() const
        {
            return thetaProps_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
