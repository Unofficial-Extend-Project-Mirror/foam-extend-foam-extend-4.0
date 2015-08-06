/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coordinateModifier.H"
#include "plane.H"

namespace Foam
{

void coordinateModifier::checkForValidInverse() const
{
    if( modifiers_.size() > 1 )
    {
        //- the if the modifiers allow combinations
        forAll(modifiers_, modI)
            if( !modifiers_[modI].combiningPossible() )
            {
                FatalErrorIn
                (
                    "void coordinateModifier::checkForValidInverse() const"
                ) << modifiers_[modI].name() << " cannot be combined with"
                  << " other anisotropic sources. The operation"
                  << " cannot be reverted!" << exit(FatalError);
            }

        //- check if the modifications overlap
        forAll(modifiers_, modI)
        {
            PtrList<plane> bndPlanes;
            modifiers_[modI].boundingPlanes(bndPlanes);

            # ifdef DEBUGCoordinateModifier
            Info << "Checking planes for object " << modifiers_[modI].name()
                 << " which are " << bndPlanes << endl;
            # endif

            for(label modJ=modI+1;modJ<modifiers_.size();++modJ)
            {
                PtrList<plane> otherBndPlanes;
                modifiers_[modJ].boundingPlanes(otherBndPlanes);

                # ifdef DEBUGCoordinateModifier
                Info << "Bnd planes planes for " << modifiers_[modJ].name()
                     << " are " << otherBndPlanes << endl;
                # endif

                for(label i=0;i<bndPlanes.size();i+=2)
                {
                    const plane& pl = bndPlanes[i];

                    for(label j=0;j<otherBndPlanes.size();j+=2)
                    {
                        const plane& opl = otherBndPlanes[j];

                        const scalar dn = mag(pl.normal() & opl.normal());
                        if( dn > SMALL )
                        {
                            if( dn < (1.0 - SMALL) )
                            {
                                FatalErrorIn
                                (
                                    "void coordinateModifier::"
                                    "checkForValidInverse() const"
                                ) << "Bounding planes of the objects "
                                  << modifiers_[modI].name()
                                  << " and " << modifiers_[modJ].name()
                                  << " are not parallel. This combination of"
                                  << " modifications cannot be reverted!"
                                  << exit(FatalError);
                            }
                            else
                            {
                                //- check if the scaling regions overlap
                                const scalar tMax =
                                    (
                                        bndPlanes[i+1].refPoint() -
                                        pl.refPoint()
                                    ) & pl.normal();

                                const scalar t0 =
                                    (
                                        otherBndPlanes[j].refPoint() -
                                        pl.refPoint()
                                    ) & pl.normal();

                                const scalar t1 =
                                    (
                                        otherBndPlanes[j+1].refPoint() -
                                        pl.refPoint()
                                    ) & pl.normal();

                                # ifdef DEBUGCoordinateModifier
                                Info << "tMax " << tMax << endl;
                                Info << "t0 " << t0 << endl;
                                Info << "t1 " << t1 << endl;
                                # endif

                                //- check if the intervals overlap
                                if( (t1 >= 0) && (t0 < tMax) )
                                {
                                    FatalErrorIn
                                    (
                                        "void coordinateModifier::"
                                        "checkForValidInverse() const"
                                    ) << "Scaling regions of objects "
                                      << modifiers_[modI].name()
                                      << " and " << modifiers_[modJ].name()
                                      << " are overlapping each other."
                                      << " This combination of"
                                      << " modifications cannot be reverted!"
                                      << exit(FatalError);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coordinateModifier::coordinateModifier(const dictionary& geomModDict)
:
    modificationDict_(geomModDict),
    modifiers_(),
    backwardModifiers_()
{
    const wordList modifiers = modificationDict_.toc();

    //- setup modification
    modifiers_.setSize(modifiers.size());
    backwardModifiers_.setSize(modifiers.size());
    forAll(modifiers, modI)
    {
        const word& mName = modifiers[modI];
        const dictionary& modDict = modificationDict_.subDict(mName);
        modifiers_.set(modI, coordinateModification::New(mName, modDict));

        backwardModifiers_.set
        (
            modI,
            coordinateModification::New(mName, modDict)
        );
    }

    //- setup backward modification
    forAll(backwardModifiers_, modI)
    {
        vector disp(vector::zero);
        const point pOrigin = backwardModifiers_[modI].origin();

        forAll(modifiers_, i)
            disp += modifiers_[i].displacement(pOrigin);

        backwardModifiers_[modI].translateAndModifyObject(disp);
    }

    checkForValidInverse();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coordinateModifier::~coordinateModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

point coordinateModifier::modifiedPoint(const point& p) const
{
    point pNew = p;

    forAll(modifiers_, modI)
    {
        pNew += modifiers_[modI].displacement(p);
    }

    return pNew;
}

point coordinateModifier::backwardModifiedPoint(const point& p) const
{
    point pNew = p;

    forAll(backwardModifiers_, modI)
        pNew += backwardModifiers_[modI].backwardDisplacement(p);

    return pNew;
}

void coordinateModifier::printObjects() const
{
    Info << "Modification objects " << modifiers_ << endl;

    Info << "Backward modification objects " << backwardModifiers_ << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
