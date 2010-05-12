/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    viscoelasticModel

Description
    Abstract base class for viscoelastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoelasticModel.C
    newViscoelasticModel.C

\*---------------------------------------------------------------------------*/

#ifndef viscoelasticModel_H
#define viscoelasticModel_H

#include "viscoelasticLaw.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoelasticModel Declaration
\*---------------------------------------------------------------------------*/

class viscoelasticModel
:
    public IOdictionary
{
private:

    // Private data

        //- Viscoelastic law
        autoPtr<viscoelasticLaw> lawPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoelasticModel(const viscoelasticModel&);

        //- Disallow default bitwise assignment
        void operator=(const viscoelasticModel&);


public:

    //- Runtime type information
    TypeName("viscoelasticModel");


    // Constructors

        //- Construct from components
        viscoelasticModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~viscoelasticModel()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoelastic stress
        virtual void correct();

        //- Read viscoelasticProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
