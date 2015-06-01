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
    Foam::basicRhoThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    basicRhoThermo.C
    newBasicRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicRhoThermo_H
#define basicRhoThermo_H

#include "basicThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class basicRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class basicRhoThermo
:
    public basicThermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rhoThermo' to avoid (potential) conflict with solver density
        volScalarField rho_;


    // Protected member functions

        //- Construct as copy (not implemented)
        basicRhoThermo(const basicRhoThermo&);


public:

    //- Runtime type information
    TypeName("basicRhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicRhoThermo,
        fvMesh,
        (const fvMesh& mesh, const objectRegistry& obj),
        (mesh, obj)
    );


    // Constructors

        //- Construct from mesh and object registry
        basicRhoThermo(const fvMesh&, const objectRegistry&);


    //- Selector
    static autoPtr<basicRhoThermo> New(const fvMesh&, const objectRegistry&);

    static autoPtr<basicRhoThermo> New(const fvMesh& mesh)
    {
        return New(mesh, mesh);
    }

    //- Destructor
    virtual ~basicRhoThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
