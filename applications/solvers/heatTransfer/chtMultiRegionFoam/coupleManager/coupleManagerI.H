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

inline const Foam::fvPatch& Foam::coupleManager::patch() const
{
    return patch_;
};


inline const Foam::word& Foam::coupleManager::neighbourRegionName() const
{
    return neighbourRegionName_;
};


inline const Foam::word& Foam::coupleManager::neighbourPatchName() const
{
    return neighbourPatchName_;
};


inline const Foam::word& Foam::coupleManager::neighbourFieldName() const
{
     return neighbourFieldName_;
};


inline const Foam::fvMesh& Foam::coupleManager::neighbourRegion() const
{
    return localRegion_.objectRegistry::parent()
        .lookupObject<fvMesh>(neighbourRegionName_);
}


inline Foam::label Foam::coupleManager::neighbourPatchID() const
{
    return neighbourRegion().boundaryMesh().findPatchID(neighbourPatchName_);
}


inline const Foam::fvPatch& Foam::coupleManager::neighbourPatch() const
{
    return neighbourRegion().boundary()[neighbourPatchID()];
}


template<class Type>
inline const Foam::fvPatchField<Type>&
Foam::coupleManager::neighbourPatchField() const
{
    return neighbourPatch().lookupPatchField
        <GeometricField<Type, fvPatchField, volMesh>, Type>
            (neighbourFieldName_);
}


inline Foam::word& Foam::coupleManager::neighbourRegionName()
{
    return neighbourRegionName_;
};


inline Foam::word& Foam::coupleManager::neighbourPatchName()
{
    return neighbourPatchName_;
};


inline Foam::word& Foam::coupleManager::neighbourFieldName()
{
    return neighbourFieldName_;
};


// ************************************************************************* //
