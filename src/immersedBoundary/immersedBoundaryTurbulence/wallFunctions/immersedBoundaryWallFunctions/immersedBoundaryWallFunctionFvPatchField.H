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
    Foam::incompressible::RASModels::immersedBoundaryWallFunctionFvPatchField

Description
    Boundary condition for passive turbulence variables (U, k, q, R, nut)
    when using wall functions on an immersed boundary patch

    Apart form standard immersed boundary capability, the patch field allows
    only some values to be fixed, as defined by the wall function calculation.

    This is handled using the ibWallValue and ibWallMask arrays, corresponding
    to the number of ibCells

    The implementation of wall functions on an immersed boundary will fix
    the values of k, epsilon, nut and tangential velocity in cells within
    the log-law layer.  For the cells in the laminar sublayer, k and epsilon
    (omega) will be calculated using a zero gradient condition,
    nut will be set to zero and the tangential component of the
    velocity will be corrected.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    immersedBoundaryWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryWallFunctionFvPatchField_H
#define immersedBoundaryWallFunctionFvPatchField_H

#include "fvPatchFields.H"
#include "immersedBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
         Class immersedBoundaryWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryWallFunctionFvPatchField
:
    public immersedBoundaryFvPatchField<Type>
{
    // Private data

        //- Value to fix in IB cell
        mutable Field<Type> wallValue_;

        //- Indicator on values to fix
        mutable boolList wallMask_;


protected:

    // Protected Member Functions

        //- Execute immersed boundary patch field motion update
        virtual void motionUpdate() const;

        //- Set IB cell values: contains data manipulation
        virtual void setIbCellValues(const Field<Type>&) const;


public:

    //- Runtime type information
    TypeName("immersedBoundaryWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // immersedBoundaryWallFunctionFvPatchField
        //  onto a new patch
        immersedBoundaryWallFunctionFvPatchField
        (
            const immersedBoundaryWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryWallFunctionFvPatchField
        (
            const immersedBoundaryWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryWallFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryWallFunctionFvPatchField
        (
            const immersedBoundaryWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryWallFunctionFvPatchField(*this, iF)
            );
        }


    //- Destructor
    virtual ~immersedBoundaryWallFunctionFvPatchField()
    {}


    // Member functions

        // Access

            //- Access to value to fix in IB cell.  Note non-const access
            Field<Type>& wallValue() const
            {
                // Note: on a moving mesh, the intersection has changed and
                // wallValue fields should be cleared and recalculated.
                // This should happen only once, but I cannot see the mechanism
                // HJ, 17/Oct/2012
                if (wallValue_.empty())
                {
                    wallValue_.setSize
                    (
                        this->ibPatch().ibCells().size(),
                        pTraits<Type>::zero
                    );
                }

                return wallValue_;
            }

            //- Access to indicator on fixed values.  Note non-const access
            boolList& wallMask() const
            {
                if (wallMask_.empty())
                {
                    wallMask_.setSize
                    (
                        this->ibPatch().ibCells().size(),
                        false
                    );
                }

                return wallMask_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryWallFunctionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
