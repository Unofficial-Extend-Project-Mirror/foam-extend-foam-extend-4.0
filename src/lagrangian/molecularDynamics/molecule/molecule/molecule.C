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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "Random.H"
#include "Time.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool molecule::move(molecule::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    scalar deltaT = cloud().pMesh().time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    moleculeCloud::integrationMethods method
            = td.molCloud().integrationMethod();

    if (method == moleculeCloud::imVerletLeapfrog)
    {
        if (td.part() == 1)  // Leapfrog 1st Part
        {
            if (stepFraction() < VSMALL)
            {
                U_ += 0.5*deltaT*A_;
            }

            while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
            {
                // set the lagrangian time-step
                scalar dt = min(dtMax, tEnd);

                dt *= trackToFace(position() + dt*U_, td);

                tEnd -= dt;
                stepFraction() = 1.0 - tEnd/deltaT;
            }
        }
        else if (td.part() == 2)  // Leapfrog 2nd Part
        {
            U_ += 0.5*deltaT*A_;
        }
        else
        {
            FatalErrorIn("molecule::move(molecule::trackData& td)") << nl
                << td.part()
                << " is an invalid part of integration method: "
                << method << nl
                << abort(FatalError);
        }
    }
    else if (method == moleculeCloud::imPredictorCorrector)
    {
        if (td.part() == 1) // Predictor Part
        {

        }
        else if (td.part() == 2) // Corrector Part
        {

        }
        else
        {
            FatalErrorIn("molecule::move(molecule::trackData& td)") << nl
                << td.part() << " is an invalid part of integration method: "
                << method
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("molecule::move(molecule::trackData& td)") << nl
                << "Unknown integration method: "
                << method
                << abort(FatalError);
    }

    return td.keepParticle;
}


void molecule::transformProperties(const tensor& T)
{}


void molecule::transformProperties(const vector& separation)
{
    if (tethered_)
    {
        tetherPosition_ += separation;
    }
}


void molecule::hitProcessorPatch
(
    const processorPolyPatch&,
    molecule::trackData& td
)
{
    td.switchProcessor = true;
}


void molecule::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void molecule::hitWallPatch
(
    const wallPolyPatch& wpp,
    molecule::trackData& td
)
{
    vector nw = wpp.faceAreas()[wpp.whichFace(face())];
    nw /= mag(nw);

    scalar Un = U_ & nw;
//     vector Ut = U_ - Un*nw;

//     Random rand(clock::getTime());

//     scalar tmac = 0.8;

//     scalar wallTemp = 2.5;

//     if (rand.scalar01() < tmac)
//     {
//         // Diffuse reflection
//
//         vector tw1 = Ut/mag(Ut);
//
//         vector tw2 = nw ^ tw1;
//
//         U_ = sqrt(wallTemp/mass_)*rand.GaussNormal()*tw1
//                 + sqrt(wallTemp/mass_)*rand.GaussNormal()*tw2
//                 - mag(sqrt(wallTemp/mass_)*rand.GaussNormal())*nw;
//     }

//     else
//     {
        // Specular reflection

        if (Un > 0)
        {
            U_ -= 2*Un*nw;
        }

//     }

}


void molecule::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


void molecule::hitPatch
(
    const polyPatch&,
    molecule::trackData& td
)
{
    td.keepParticle = false;
}


void molecule::hitPatch
(
    const polyPatch&,
    int&
)
{}

} // End namespace Foam


// ************************************************************************* //
