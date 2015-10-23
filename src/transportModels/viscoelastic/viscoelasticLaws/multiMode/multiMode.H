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
    multiMode

Description
    Multi-mode viscoelastic fluid model.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    multiMode.C

\*---------------------------------------------------------------------------*/

#ifndef multiMode_H
#define multiMode_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiMode Declaration
\*---------------------------------------------------------------------------*/

class multiMode
:
    public viscoelasticLaw
{
    // Private data

        //- Transported viscoelastic stress
        mutable volSymmTensorField tau_;

        //- List of models
        PtrList<viscoelasticLaw> models_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMode(const multiMode&);

        //- Disallow default bitwise assignment
        void operator=(const multiMode&);


public:

    //- Runtime type information
    TypeName("multiMode");

    // Constructors

        //- Construct from components
        multiMode
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMode()
        {}


    // Member Functions

        virtual tmp<volSymmTensorField> tau() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoelastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
