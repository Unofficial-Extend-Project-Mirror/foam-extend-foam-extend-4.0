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

Namespace
    Foam::sixDoFRigidBodyMotionRestraints

Description
    Namespace for six DoF motion restraints


Class
    Foam::sixDoFRigidBodyMotionRestraint

Description
    Base class for defining restraints for sixDoF motions

SourceFiles
    sixDoFRigidBodyMotionRestraint.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionRestraint_H
#define sixDoFRigidBodyMotionRestraint_H

#include "objectRegistry.H"
#include "foamTime.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDoFRigidBodyMotion;


/*---------------------------------------------------------------------------*\
                Class sixDoFRigidBodyMotionRestraint Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionRestraint
{

protected:

    // Protected data

        //- Restraint model specific coefficient dictionary
        dictionary sDoFRBMRCoeffs_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFRigidBodyMotionRestraint,
            dictionary,
            (const dictionary& sDoFRBMRDict),
            (sDoFRBMRDict)
        );


    // Constructors

        //- Construct from the sDoFRBMRDict dictionary and Time
        sixDoFRigidBodyMotionRestraint
        (
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMRDict dictionary and Time
        static autoPtr<sixDoFRigidBodyMotionRestraint> New
        (
            const dictionary& sDoFRBMRDict
        );


    // Destructor

        virtual ~sixDoFRigidBodyMotionRestraint();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRDict);

        // Access

            // Return access to sDoFRBMRCoeffs
            inline const dictionary& coeffDict() const
            {
                return sDoFRBMRCoeffs_;
            }

        //- Write
        virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
