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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const PtrList<equation>& equationReader::eqns() const
{
    return eqns_;
}


inline const PtrList<scalar>&
    equationReader::internalScalars() const
{
    return internalScalars_;
}


inline const UPtrList<const dictionary>&
    equationReader::dictSources() const
{
    return dictSources_;
}


inline const PtrList<word>& equationReader::dictLookups() const
{
    return dictLookups_;
}


inline const UPtrList<const equationVariable>&
    equationReader::activeSources() const
{
    return activeSources_;
}


inline const wordList& equationReader::activeSourceNames() const
{
    return activeSourceNames_;
}


inline const equationSource<scalar>&
    equationReader::scalarSources() const
{
    return scalarSources_;
}


inline equationSource<scalar>& equationReader::scalarSources()
{
    return scalarSources_;
}


inline const equationSource<vector>&
    equationReader::vectorSources() const
{
    return vectorSources_;
}


inline equationSource<vector>& equationReader::vectorSources()
{
    return vectorSources_;
}


inline const equationSource<tensor>&
    equationReader::tensorSources() const
{
    return tensorSources_;
}


inline equationSource<tensor>& equationReader::tensorSources()
{
    return tensorSources_;
}


inline const equationSource<diagTensor>&
    equationReader::diagTensorSources() const
{
    return diagTensorSources_;
}


inline equationSource<diagTensor>&
    equationReader::diagTensorSources()
{
    return diagTensorSources_;
}


inline const equationSource<symmTensor>&
    equationReader::symmTensorSources() const
{
    return symmTensorSources_;
}


inline equationSource<symmTensor>&
    equationReader::symmTensorSources()
{
    return symmTensorSources_;
}


inline const equationSource<sphericalTensor>&
    equationReader::sphericalTensorSources() const
{
    return sphericalTensorSources_;
}


inline equationSource<sphericalTensor>&
    equationReader::sphericalTensorSources()
{
    return sphericalTensorSources_;
}


/*inline const label& equationReader::geoIndex() const
{
    return geoIndex_;
}


inline void equationReader::setGeoIndex(label newIndex)
{
    geoIndex_ = newIndex;
}


inline const label& equationReader::cellIndex() const
{
    return cellIndex_;
}


inline void equationReader::setCellIndex(label newIndex)
{
    cellIndex_ = newIndex;
}*/

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
