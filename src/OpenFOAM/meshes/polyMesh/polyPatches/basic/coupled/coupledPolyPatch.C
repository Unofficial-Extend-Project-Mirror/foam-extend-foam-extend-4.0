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

\*---------------------------------------------------------------------------*/

#include "coupledPolyPatch.H"
#include "ListOps.H"
#include "transform.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::coupledPolyPatch, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::coupledPolyPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nr,
    const scalarField& smallDist,
    const scalar absTol
) const
{
    if (debug)
    {
        Pout<< "coupledPolyPatch::calcTransformTensors : " << name() << endl
            << "    (half)size:" << Cf.size() << nl
            << "    absTol:" << absTol << nl
            << "    sum(mag(nf & nr)):" << sum(mag(nf & nr)) << endl;
    }

    // Tolerance calculation.
    // - normal calculation: assume absTol is the absolute error in a
    // single normal/transformation calculation. Consists both of numerical
    // precision (on the order of SMALL and of writing precision
    // (from e.g. decomposition)
    // Then the overall error of summing the normals is sqrt(size())*absTol
    // - separation calculation: pass in from the outside an allowable error.

    if (size() == 0)
    {
        // Dummy geometry.
        separation_.setSize(0);
        forwardT_ = I;
        reverseT_ = I;
    }
    else
    {
        scalar error = absTol*Foam::sqrt(1.0*Cf.size());

        if (debug)
        {
            Pout<< "    error:" << error << endl;
        }

        if (sum(mag(nf & nr)) < Cf.size() - error)
        {
            // Rotation, no separation

            separation_.setSize(0);

            forwardT_.setSize(Cf.size());
            reverseT_.setSize(Cf.size());

            forAll (forwardT_, facei)
            {
                forwardT_[facei] = rotationTensor(-nr[facei], nf[facei]);
                reverseT_[facei] = rotationTensor(nf[facei], -nr[facei]);
            }

            if (debug)
            {
                Pout<< "    sum(mag(forwardT_ - forwardT_[0])):"
                    << sum(mag(forwardT_ - forwardT_[0]))
                    << endl;
            }

            if (sum(mag(forwardT_ - forwardT_[0])) < error)
            {
                forwardT_.setSize(1);
                reverseT_.setSize(1);

                if (debug)
                {
                    Pout<< "    difference in rotation less than"
                        << " local tolerance "
                        << error << ". Assuming uniform rotation." << endl;
                }
            }
        }
        else
        {
            forwardT_.setSize(0);
            reverseT_.setSize(0);

            separation_ = (nf&(Cr - Cf))*nf;

            // Three situations:
            // - separation is zero. No separation.
            // - separation is same. Single separation vector.
            // - separation differs per face. Separation vectorField.

            // Check for different separation per face
            bool sameSeparation = true;

            forAll(separation_, facei)
            {
                scalar smallSqr = sqr(smallDist[facei]);

                if (magSqr(separation_[facei] - separation_[0]) > smallSqr)
                {
                    if (debug)
                    {
                        Pout<< "    separation " << separation_[facei]
                            << " at " << facei
                            << " differs from separation[0] " << separation_[0]
                            << " by more than local tolerance "
                            << smallDist[facei]
                            << ". Assuming non-uniform separation." << endl;
                    }
                    sameSeparation = false;
                    break;
                }
            }

            if (sameSeparation)
            {
                // Check for zero separation (at 0 so everywhere)
                if (magSqr(separation_[0]) < sqr(smallDist[0]))
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation_[0])
                            << " less than local tolerance " << smallDist[0]
                            << ". Assuming zero separation." << endl;
                    }

                    separation_.setSize(0);
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation_[0])
                            << " more than local tolerance " << smallDist[0]
                            << ". Assuming uniform separation." << endl;
                    }

                    separation_.setSize(1);
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "    separation_:" << separation_.size() << nl
            << "    forwardT size:" << forwardT_.size() << endl;
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::~coupledPolyPatch()
{}


// ************************************************************************* //
