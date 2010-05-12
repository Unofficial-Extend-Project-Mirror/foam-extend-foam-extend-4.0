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

#include "molecule.H"
#include "IOstreams.H"

#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecule::molecule
(
    const Cloud<molecule>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<molecule>(cloud, is)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            id_ = readLabel(is);
            mass_ = readScalar(is);
            is >> U_;
            is >> A_;
            is >> potentialEnergy_;
            is >> rf_;
            is >> tethered_;
            is >> tetherPosition_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&mass_),
                sizeof(mass_)
                + sizeof(U_)
                + sizeof(A_)
                + sizeof(tetherPosition_)
                + sizeof(potentialEnergy_)
                + sizeof(rf_)
                + sizeof(tethered_)
                + sizeof(id_)
            );
        }
    }

    // Check state of Istream
    is.check("Foam::molecule::molecule(Foam::Istream&)");
}


namespace Foam
{

void molecule::readFields(moleculeCloud& mC)
{
    if (!mC.size())
    {
        return;
    }

    IOField<label> id(mC.fieldIOobject("id"));
    mC.checkFieldIOobject(mC, id);

    IOField<scalar> mass(mC.fieldIOobject("mass"));
    mC.checkFieldIOobject(mC, mass);

    IOField<vector> U(mC.fieldIOobject("U"));
    mC.checkFieldIOobject(mC, U);

    IOField<vector> A(mC.fieldIOobject("A"));
    mC.checkFieldIOobject(mC, A);

    IOField<label> tethered(mC.fieldIOobject("tethered"));
    mC.checkFieldIOobject(mC, tethered);

    IOField<vector> tetherPositions(mC.fieldIOobject("tetherPositions"));
    mC.checkFieldIOobject(mC, tetherPositions);

    label i = 0;
    forAllIter(moleculeCloud, mC, iter)
    {
        molecule& mol = iter();

        mol.id_ = id[i];
        mol.mass_ = mass[i];
        mol.U_ = U[i];
        mol.A_ = A[i];
        mol.potentialEnergy_ = 0.0;
        mol.rf_ = tensor::zero;
        mol.tethered_ = tethered[i];
        mol.tetherPosition_ = tetherPositions[i];
        i++;
    }
}


void molecule::writeFields(const moleculeCloud& mC)
{
    Particle<molecule>::writeFields(mC);

    label np =  mC.size();

    IOField<label> id(mC.fieldIOobject("id"), np);
    IOField<scalar> mass(mC.fieldIOobject("mass"), np);
    IOField<vector> U(mC.fieldIOobject("U"), np);
    IOField<vector> A(mC.fieldIOobject("A"), np);
    IOField<label> tethered(mC.fieldIOobject("tethered"), np);
    IOField<vector> tetherPositions(mC.fieldIOobject("tetherPositions"), np);

    label i = 0;
    forAllConstIter(moleculeCloud, mC, iter)
    {
        const molecule& mol = iter();

        id[i] = mol.id_;
        mass[i] = mol.mass_;
        U[i] = mol.U_;
        A[i] = mol.A_;
        tethered[i] = mol.tethered_;
        tetherPositions[i] = mol.tetherPosition_;
        i++;
    }

    id.write();
    mass.write();
    U.write();
    A.write();
    tethered.write();
    tetherPositions.write();
}

};  // end of namespace Foam


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const molecule& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << mol.id_
                << token::SPACE << mol.mass_
                << token::SPACE << static_cast<const Particle<molecule>&>(mol)
                << token::SPACE << mol.face()
                << token::SPACE << mol.stepFraction()
                << token::SPACE << mol.U_
                << token::SPACE << mol.A_
                << token::SPACE << mol.potentialEnergy_
                << token::SPACE << mol.rf_
                << token::SPACE << mol.tethered_
                << token::SPACE << mol.tetherPosition_;
    }
    else
    {
        os  << static_cast<const Particle<molecule>&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.mass_),
            sizeof(mol.mass_)
            + sizeof(mol.U_)
            + sizeof(mol.A_)
            + sizeof(mol.tetherPosition_)
            + sizeof(mol.potentialEnergy_)
            + sizeof(mol.rf_)
            + sizeof(mol.tethered_)
            + sizeof(mol.id_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::molecule&)"
    );

    return os;
}


// ************************************************************************* //
