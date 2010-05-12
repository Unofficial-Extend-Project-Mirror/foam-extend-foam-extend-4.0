/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::forcesFiltered

Description
    Derived from the forces function object, creates a specialisation to
    calculate filtered lift and drag forces.

SourceFiles
    forcesFiltered.C
    IOforcesFiltered.H

\*---------------------------------------------------------------------------*/

#ifndef forcesFiltered_H
#define forcesFiltered_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class forcesFiltered Declaration
\*---------------------------------------------------------------------------*/

class forcesFiltered
:
    public forces
{
    // Private data

        //- Filter field name
        word filterFieldName_;


    // Private member functions

        //- Disallow default bitwise copy construct
        forcesFiltered(const forcesFiltered&);

        //- Disallow default bitwise assignment
        void operator=(const forcesFiltered&);


public:

    //- Runtime type information
    TypeName("forcesFiltered");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        forcesFiltered
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~forcesFiltered();


    // Member Functions

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Calculate and return forces and moments from the calculated the
        //  stress filtered by the filter field provided
        virtual forcesMoments calcForces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
