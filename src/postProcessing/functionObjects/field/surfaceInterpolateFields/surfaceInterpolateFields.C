/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

#include "surfaceInterpolateFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceInterpolateFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceInterpolateFields::surfaceInterpolateFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fieldSet_()
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "surfaceInterpolateFields::surfaceInterpolateFields"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceInterpolateFields::~surfaceInterpolateFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceInterpolateFields::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::surfaceInterpolateFields::execute()
{
    if (active_)
    {
        Info<< type() << " " << name_ << " output:" << nl;

        // Clear out any previously loaded fields
        ssf_.clear();
        svf_.clear();
        sSpheretf_.clear();
        sSymmtf_.clear();
        stf_.clear();

        interpolateFields<scalar>(ssf_);
        interpolateFields<vector>(svf_);
        interpolateFields<sphericalTensor>(sSpheretf_);
        interpolateFields<symmTensor>(sSymmtf_);
        interpolateFields<tensor>(stf_);

        Info<< endl;
    }
}


void Foam::surfaceInterpolateFields::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::surfaceInterpolateFields::timeSet()
{
    // Do nothing
}


void Foam::surfaceInterpolateFields::write()
{
    if (active_)
    {
        Info<< type() << " " << name_ << " output:" << nl;

        Info<< "    Writing interpolated surface fields to "
            << obr_.time().timeName() << endl;

        forAll(ssf_, i)
        {
            ssf_[i].write();
        }
        forAll(svf_, i)
        {
            svf_[i].write();
        }
        forAll(sSpheretf_, i)
        {
            sSpheretf_[i].write();
        }
        forAll(sSymmtf_, i)
        {
            sSymmtf_[i].write();
        }
        forAll(stf_, i)
        {
            stf_[i].write();
        }
    }
}


// ************************************************************************* //
