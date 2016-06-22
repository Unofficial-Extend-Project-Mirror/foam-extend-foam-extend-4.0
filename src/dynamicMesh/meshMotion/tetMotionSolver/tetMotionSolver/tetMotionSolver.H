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
    tetMotionSolver

Description
    Virtual base class for tetDecomposition mesh motion solvers.

SourceFiles
    tetMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef tetMotionSolver_H
#define tetMotionSolver_H

#include "motionSolver.H"
#include "polyMesh.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "tetFemMatrices.H"
#include "tetPolyMeshMapper.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class tetMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class tetMotionSolver
:
    public motionSolver
{
    // Private data

        //- Motion mesh
        tetPolyMesh tetMesh_;

        //- Motion field
        tetPointVectorField motionU_;

        //- List of constrained points
        dynamicLabelList fixedPoints_;

        //- List of constrained velocities
        DynamicList<vector> fixedVelocity_;

        //- Total displacement field
        mutable tetPointVectorField* totDisplacementPtr_;


protected:

    // Protected member functions

        //- Apply constraints to the matrix
        void applyConstraints(tetFemVectorMatrix& matrix);


public:

    //- Runtime type information
    TypeName("tetMotionSolver");


    // Constructors

        //- Construct from components
        tetMotionSolver(const polyMesh& mesh);


    // Destructor

        virtual ~tetMotionSolver();


    // Member Functions

        // Access

            //- Return reference to motion mesh
            const tetPolyMesh& tetMesh() const
            {
                return tetMesh_;
            }

            //- Return reference to motion field
            tetPointVectorField& motionU()
            {
                return motionU_;
            }

            void storeTotDisplacement() const
            {
                if (!totDisplacementPtr_)
                {
                    totDisplacementPtr_ =
                        new tetPointVectorField
                        (
                            IOobject
                            (
                                "totalMotionU",
                                mesh().time().timeName(),
                                mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            tetMesh(),
                            dimensionedVector("zero", dimLength, vector::zero)
                        );
                }
            }

            bool needTotDisplacement() const
            {
                return totDisplacementPtr_;
            }


            //- Return reference to total displacement field
            tetPointVectorField& totDisplacement()
            {
                return *totDisplacementPtr_;
            }

            const tetPointVectorField& totDisplacement() const
            {
                return *totDisplacementPtr_;
            }

            //- Provide current points for motion.  Uses current motion field
            virtual tmp<pointField> curPoints() const;


        // Solving for motion

            //- Set point motion constraints
            void setConstraint
            (
                const label pointID,
                const vector& fixedVel
            );

            //- Clear constraints
            virtual void clearConstraints();

            //- Solve for motion
            virtual void solve() = 0;

            //- Update topology
            virtual void updateMesh(const mapPolyMesh&);


        // Mesh quality measures

            //- Distortion energy density
            tmp<elementScalarField> distortionEnergy() const;

            //- Deformation energy density
            tmp<elementScalarField> deformationEnergy() const;

            //- Total distortion energy density
            tmp<elementScalarField> totDistortionEnergy() const;

            //- Total deformation energy density
            tmp<elementScalarField> totDeformationEnergy() const;

            //- Return mean cell non-orthogonality
            tmp<scalarField> nonOrthogonality() const;

            //- Return mean cell skewness
            tmp<Foam::scalarField> cellSkewness() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
