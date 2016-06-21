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

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Face>
Foam::autoPtr< Foam::MeshedSurface<Face> >
Foam::MeshedSurface<Face>::New(const fileName& name, const word& ext)
{
    if (debug)
    {
        Info<< "MeshedSurface::New(const fileName&, const word&) : "
            "constructing MeshedSurface"
            << endl;
    }

    typename fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(ext);

    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        // no direct reader, delegate if possible
        wordHashSet supported = FriendType::readTypes();
        if (supported.found(ext))
        {
            // create indirectly
            autoPtr< MeshedSurface<Face> > surf(new MeshedSurface<Face>);
            surf().transfer(FriendType::New(name, ext)());

            return surf;
        }

        // nothing left to try, issue error
        supported += readTypes();

        FatalErrorIn
        (
            "MeshedSurface<Face>::New(const fileName&, const word&) : "
            "constructing MeshedSurface"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << nl
            << supported
            << exit(FatalError);
    }

    return autoPtr< MeshedSurface<Face> >(cstrIter()(name));
}


template<class Face>
Foam::autoPtr< Foam::MeshedSurface<Face> >
Foam::MeshedSurface<Face>::New(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return New(name, ext);
}

// ************************************************************************* //
