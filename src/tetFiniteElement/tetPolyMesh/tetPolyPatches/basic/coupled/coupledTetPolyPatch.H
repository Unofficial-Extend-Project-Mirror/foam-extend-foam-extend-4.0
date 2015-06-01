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
    coupledTetPolyPatch

Description
    Coupled patch base class.

SourceFiles
    coupledTetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTetPolyPatch_H
#define coupledTetPolyPatch_H

#include "coupledPolyPatch.H"
#include "tetPolyPatch.H"
#include "lduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
               Class coupledTetPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledTetPolyPatch
:
    public lduInterface,
    public tetPolyPatch
{
    // Private Member Functions

        //- Disallow default construct as copy
        coupledTetPolyPatch(const coupledTetPolyPatch&);

        //- Disallow default assignment
        void operator=(const coupledTetPolyPatch&);

public:

    typedef tetPolyBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledTetPolyPatch(const tetPolyBoundaryMesh& bm);


    // Destructor

        virtual ~coupledTetPolyPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
