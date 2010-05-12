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
    coupledTetPolyPatchCellDecomp

Description
    Coupled patch base class.

SourceFiles
    coupledTetPolyPatchCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTetPolyPatchCellDecomp_H
#define coupledTetPolyPatchCellDecomp_H

#include "lduInterface.H"
#include "coupledPolyPatch.H"
#include "tetPolyPatchCellDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetPolyBoundaryMeshCellDecomp;

/*---------------------------------------------------------------------------*\
               Class coupledTetPolyPatchCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class coupledTetPolyPatchCellDecomp
:
    public lduInterface,
    public tetPolyPatchCellDecomp
{
    // Private Member Functions

        //- Disallow default construct as copy
        coupledTetPolyPatchCellDecomp(const coupledTetPolyPatchCellDecomp&);

        //- Disallow default assignment
        void operator=(const coupledTetPolyPatchCellDecomp&);

public:

    typedef tetPolyBoundaryMeshCellDecomp BoundaryMesh;


    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledTetPolyPatchCellDecomp(const tetPolyBoundaryMeshCellDecomp& bm);


    // Destructor

        virtual ~coupledTetPolyPatchCellDecomp();


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
