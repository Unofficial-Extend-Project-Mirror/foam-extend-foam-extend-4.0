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
    PODODE

Description
    Virtual base class for ODEs derived by Proper Orthogonal Decomposition

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    PODODE.C

\*---------------------------------------------------------------------------*/

#ifndef PODODE_H
#define PODODE_H

#include "ODE.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PODODE Declaration
\*---------------------------------------------------------------------------*/

class PODODE
:
    public ODE
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to dictionary
        const dictionary& dict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PODODE(const PODODE&);

        //- Disallow default bitwise assignment
        void operator=(const PODODE&);


public:

    //- Runtime type information
    TypeName("PODODE");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            PODODE,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected POD model
        static autoPtr<PODODE> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        PODODE
        (
            const fvMesh& mesh,
            const dictionary& dict
        )
        :
            mesh_(mesh),
            dict_(dict)
        {}


    // Destructor

        virtual ~PODODE()
        {}


    // Member Functions

        //- Return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }


        // Solution variables

            //- Return reference to interpolation coefficients
            virtual scalarField& coeffs() = 0;

            //- Return reference to interpolation coefficients
            virtual const scalarField& coeffs() const = 0;

            //- Clear ortho-normal base
            virtual void clearBase() const = 0;

            //- Update reconstructed fields
            virtual void updateFields() const = 0;

            //- Clear reconstructed field
            virtual void clearFields() const = 0;


        // Update

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta)
            {}


        //- Write reconstructed fields
        virtual void write() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
