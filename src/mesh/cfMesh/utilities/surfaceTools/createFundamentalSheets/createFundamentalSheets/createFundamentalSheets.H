/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    createFundamentalSheets

Description
    A base class for various method to generate fundamental sheets
    necessary to capture feature edges

SourceFiles
    createFundamentalSheets.C

\*---------------------------------------------------------------------------*/

#ifndef createFundamentalSheets_H
#define createFundamentalSheets_H

#include "polyMeshGenModifier.H"
#include "labelLongList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class createFundamentalSheets Declaration
\*---------------------------------------------------------------------------*/

class createFundamentalSheets
{
protected:
    // Protected data

        //- type of sheet creation
        word name_;

        //- reference to mesh
        polyMeshGen& mesh_;

        //- shall the procedure create the intial wrapper sheet
        const bool createWrapperSheet_;

    // Protected member functions

        //- Disallow default bitwise copy construct
        createFundamentalSheets(const createFundamentalSheets&);

        //- Disallow default bitwise assignment
        void operator=(const createFundamentalSheets&);

public:

    // Runtime type information
        TypeName("createFundamentalSheets");

    // Constructors

        //- Construct from mesh data
        createFundamentalSheets
        (
            polyMeshGen& mesh,
            const bool createWrapperSheet = true
        );

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            createFundamentalSheets,
            polyMeshGen,
            (
                polyMeshGen& mesh,
                const bool createWrapperSheet
            ),
            (mesh, createWrapperSheet)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<createFundamentalSheets> New
        (
            polyMeshGen& mesh, const bool createWrapperSheet = true
        );

    // Destructor

        virtual ~createFundamentalSheets();

    // Member Functions
        //- return the method name
        const word& name() const
        {
            return name_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
