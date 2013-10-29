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
    rheologyModel

Description
    Material rheology for solids.

SourceFiles
    rheologyModel.C

\*---------------------------------------------------------------------------*/

#ifndef rheologyModel_H
#define rheologyModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "rheologyLaw.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rheologyModel Declaration
\*---------------------------------------------------------------------------*/

class rheologyModel
:
    public IOdictionary
{
    // Private data

        //- Reference to stress field
        const volSymmTensorField& sigma_;

        //- Plane stress
        Switch planeStress_;

        //- Rheology law
        autoPtr<rheologyLaw> lawPtr_;


    // Private Member Functions

        //- Disallow copy construct
        rheologyModel(const rheologyModel&);

        //- Disallow default bitwise assignment
        void operator=(const rheologyModel&);


public:

    //- Runtime type information
    TypeName("rheologyModel");

    // Constructors

        //- Construct from dictionary
        rheologyModel
        (
            const volSymmTensorField& sigma
        );


    // Destructor

        virtual ~rheologyModel()
        {}


    // Member Functions

        //- Return true for plane stress
        const Switch& planeStress() const
        {
            return planeStress_;
        }

        //- Return rheology law
        const rheologyLaw& law() const
        {
            return lawPtr_();
        }

        //- Return density
        tmp<volScalarField> rho() const
        {
            return lawPtr_->rho();
        }

        //- Return first Lame's coefficient
        tmp<volScalarField> mu() const;

        //- Return second Lame's coefficient
        tmp<volScalarField> lambda() const;

        //- Return threeK
        tmp<volScalarField> threeK() const;

        //- Return density
        tmp<volScalarField> rho(scalar t) const
        {
            return lawPtr_->rho(t);
        }

        //- Return first Lame's coefficient
        tmp<volScalarField> mu(scalar t) const;

        //- Return second Lame's coefficient
        tmp<volScalarField> lambda(scalar t) const;

        //- Correct the rheological model
        void correct();

        //- Read rheologyProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
