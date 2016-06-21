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
    Foam::thermoCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    thermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef thermoCloud_H
#define thermoCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class thermoCloud Declaration
\*---------------------------------------------------------------------------*/

class thermoCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoCloud(const thermoCloud&);

        //- Disallow default bitwise assignment
        void operator=(const thermoCloud&);


public:

    //- Runtime type information
    TypeName("thermoCloud");

    // Constructors

        //- Null constructor
        thermoCloud();


    //- Destructor
    virtual ~thermoCloud();


    // Member Functions

        // Access

            // Radiation

                //- Return tmp equivalent particulate emission
                virtual tmp<volScalarField> Ep() const = 0;

                //- Return tmp equivalent particulate absorption
                virtual tmp<volScalarField> ap() const = 0;

                //- Return tmp equivalent particulate scattering factor
                virtual tmp<volScalarField> sigmap() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
