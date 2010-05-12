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

Global
    meanMomentumEnergyAndNMols.H

Description
    Calculates and prints the mean momentum and energy in the system
    and the number of molecules.

\*---------------------------------------------------------------------------*/


vector singleStepTotalMomentum(vector::zero);

scalar singleStepMaxVelocityMag = 0.0;

scalar singleStepTotalMass = 0.0;

scalar singleStepTotalKE = 0.0;

scalar singleStepTotalPE = 0.0;

scalar singleStepTotalrDotf = 0.0;

{
    IDLList<molecule>::iterator mol(molecules.begin());

    for
    (
        mol = molecules.begin();
        mol != molecules.end();
        ++mol
    )
    {
        const scalar molM(mol().mass());

        const vector& molU(mol().U());

        singleStepTotalMomentum += molU * molM;

        singleStepTotalMass += molM;

        if(mag(molU) > singleStepMaxVelocityMag)
        {
            singleStepMaxVelocityMag = mag(molU);
        }

        singleStepTotalKE += 0.5*molM*magSqr(molU);

        singleStepTotalPE += mol().potentialEnergy();

        singleStepTotalrDotf += tr(mol().rf());
    }
}

label singleStepNMols = molecules.size();

if (Pstream::parRun())
{
    reduce(singleStepTotalMomentum, sumOp<vector>());

    reduce(singleStepMaxVelocityMag, maxOp<scalar>());

    reduce(singleStepTotalMass, sumOp<scalar>());

    reduce(singleStepTotalKE, sumOp<scalar>());

    reduce(singleStepTotalPE, sumOp<scalar>());

    reduce(singleStepTotalrDotf, sumOp<scalar>());

    reduce(singleStepNMols, sumOp<label>());
}

if (singleStepNMols)
{
    Info << "Number of mols in system = "
            << singleStepNMols << nl
        << "Overall number density = "
            << singleStepNMols/meshVolume << " m^-3" << nl
        << "Overall mass density = "
            << singleStepTotalMass/meshVolume << " kg/m^3" << nl
        << "Average velocity per mol = "
            << singleStepTotalMomentum/singleStepTotalMass << " m/s" << nl
        << "Maximum |velocity| = "
            << singleStepMaxVelocityMag << " m/s" << nl
        << "Average KE per mol = "
            << singleStepTotalKE/singleStepNMols << " J" << nl
        << "Average PE per mol = "
            << singleStepTotalPE/singleStepNMols << " J" << nl
        << "Average TE per mol = "
            << (singleStepTotalKE + singleStepTotalPE)/singleStepNMols << " J"
        << endl;

//     Info << singleStepNMols << " "
// //     << singleStepTotalMomentum/singleStepTotalMass << " "
//     << singleStepMaxVelocityMag << " "
//     << singleStepTotalKE/singleStepNMols << " "
//     << singleStepTotalPE/singleStepNMols << " "
//     << (singleStepTotalKE + singleStepTotalPE)/singleStepNMols << endl;
}
else
{
    Info << "No molecules in system" << endl;
}


// ************************************************************************* //
