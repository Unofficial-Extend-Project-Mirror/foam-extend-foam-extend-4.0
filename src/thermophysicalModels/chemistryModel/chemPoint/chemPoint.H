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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Point of the binary tree.

\*---------------------------------------------------------------------------*/

#ifndef chemPoint_H
#define chemPoint_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"

namespace Foam
{

class binaryNode;

/*---------------------------------------------------------------------------*\
                           Class chemPoint Declaration
\*---------------------------------------------------------------------------*/

class chemPoint
{
public:

        //- Reference to the node in the binary tree
        binaryNode* node_;


private:


        //- Number of times the element has been used
        label nUsed_;

        //- Vector storing the composition, temperature and pressure
        scalarField v0_;

        //- Reaction mapping for v0
        scalarField r_;

        //- Ellipsoid of accuracy
        scalarField EOA_;

        //- Ellipsoid of accuracy for the solution
        scalarField solutionsEOA_;

        //- tolerance for the Ellipsoid of accuracy
        scalar absErr_;

        //- Use logarithm of temperature
        Switch logT_;

        //- Time step, stored but not actually used
        scalar deltaT_;

        //- Ellipsoid of accuracy for the time step
        scalar timeEOA_;


public:

    // Constructors

        //- Construct from components
        chemPoint
        (
            const scalarField& v0,
            const scalarField& r,
            const scalarField& tolerances,
            const scalarField& tolerancesSolutions,
            const scalar& absErr,
            const Switch& logT,
            const scalar& deltaT
        );

        //- Construct from components and reference to a binary node
        chemPoint
        (
            const scalarField& v0,
            const scalarField& r,
            const scalarField& tolerances,
            const scalarField& tolerancesSolutions,
            const scalar& absErr,
            const Switch& logT,
            const scalar& deltaT,
            binaryNode* node
        );

        //- Construct from another chemPoint and reference to a binary node
        chemPoint
        (
            const chemPoint& p,
            binaryNode* node
        );

        //- Construct from another chemPoint
        chemPoint
        (
            const chemPoint& p
        );


        //- Access

        inline label& nUsed()
        {
            return nUsed_;
        }

        inline const label& nUsed() const
        {
            return nUsed_;
        }

        inline scalarField& v0()
        {
            return v0_;
        }

        inline const scalarField& v0() const
        {
            return v0_;
        }

        inline scalarField& r()
        {
            return r_;
        }

        inline const scalarField& r() const
        {
            return r_;
        }

        inline scalarField& EOA()
        {
            return EOA_;
        }

        inline const scalarField& EOA() const
        {
            return EOA_;
        }

        inline scalarField& solutionsEOA()
        {
            return solutionsEOA_;
        }

        inline const scalarField& solutionsEOA() const
        {
            return solutionsEOA_;
        }

        inline const scalar& absErr() const
        {
            return absErr_;
        }

        inline Switch& logT()
        {
            return logT_;
        }

        inline const Switch& logT() const
        {
            return logT_;
        }

        inline binaryNode* node()
        {
            return node_;
        }

        inline scalar& deltaT()
        {
            return deltaT_;
        }

        inline const scalar& deltaT() const
        {
            return deltaT_;
        }

        inline scalar& timeEOA()
        {
            return timeEOA_;
        }

        inline const scalar& timeEOA() const
        {
            return timeEOA_;
        }


//         scalar difference
//         (
//             const scalarField& point,
//             const scalarField& scaleFactor,
//             const bool& logT
//         );

        // Is the point in the ellipsoid of accuracy?
        bool inEOA
        (
            const scalarField& point,
            const scalarField& Wi,
            const scalar& rhoi,
            const scalar& deltaT,
            const scalarField& scaleFactor
        );

        // Grow the ellipsoid of accuracy?
        void grow
        (
            const scalarField& v0,
            const scalarField& Wi,
            const scalar& rhoi,
            const scalarField& v,
            const scalar& deltaT
        );

        // Check if the new solution is in the ellipsoid of accuracy?
        bool checkSolution
        (
            const scalarField& v0,
            const scalarField& v,
            const scalarField& Wi,
            const scalar& rhoi,
            const scalar& T,
            const scalar& p,
            const scalar& deltaT,
            const scalarField& tolerance
        );

        //- Set free the point from its node, used for replacing
        //  in the binary tree
        void setFree();

        // clear all the stored data
        void clearData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
