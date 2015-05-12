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
    Foam::incompressible::RASModels::
    immersedBoundaryVelocityWallFunctionFvPatchVectorField

Description
    Boundary condition for velocity when using wall functions
    - uses tangential velocity as prescribed by the epsilon boundary condition
      to enforce into the fit

SourceFiles
    immersedBoundaryVelocityWallFunctionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryVelocityWallFunctionFvPatchVectorField_H
#define immersedBoundaryVelocityWallFunctionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "immersedBoundaryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
    Class immersedBoundaryVelocityWallFunctionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryVelocityWallFunctionFvPatchVectorField
:
    public immersedBoundaryFvPatchVectorField
{
    // Private data

        //- Tangential velocity value to fix in IB cell
        mutable scalarField wallTangentialValue_;

        //- Wall shear stress
        mutable vectorField tauWall_;

        //- Indicator on values to fix
        mutable boolList wallMask_;


protected:

    // Protected Member Functions

        //- Set IB cell values: contains data manipulation
        virtual void setIbCellValues(const vectorField&) const;


public:

    //- Runtime type information
    TypeName("immersedBoundaryVelocityWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryVelocityWallFunctionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryVelocityWallFunctionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  immersedBoundaryVelocityWallFunctionFvPatchVectorField
        //  onto a new patch
        immersedBoundaryVelocityWallFunctionFvPatchVectorField
        (
            const immersedBoundaryVelocityWallFunctionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryVelocityWallFunctionFvPatchVectorField
        (
            const immersedBoundaryVelocityWallFunctionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new immersedBoundaryVelocityWallFunctionFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryVelocityWallFunctionFvPatchVectorField
        (
            const immersedBoundaryVelocityWallFunctionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new immersedBoundaryVelocityWallFunctionFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~immersedBoundaryVelocityWallFunctionFvPatchVectorField()
    {}


    // Member functions

        // Access

            //- Access to tangential velocity value to fix in IB cell
            //  Note non-const access
            scalarField& wallTangentialValue() const
            {
                if (wallTangentialValue_.empty())
                {
                    wallTangentialValue_.setSize
                    (
                        this->ibPatch().ibCells().size(),
                        0
                    );
                }

                return wallTangentialValue_;
            }

            //- Return wall shear stress
            const vectorField& wallShearStress() const;

            //- Access to wall shear stress in IB cell
            //  Note non-const access
            vectorField& tauWall() const
            {
                if (tauWall_.empty())
                {
                    tauWall_.setSize
                    (
                        this->ibPatch().ibCells().size(),
                        vector::zero
                    );
                }

                return tauWall_;
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

#endif

// ************************************************************************* //
