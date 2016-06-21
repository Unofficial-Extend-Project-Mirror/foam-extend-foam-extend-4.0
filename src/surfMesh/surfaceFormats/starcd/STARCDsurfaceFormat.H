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

Class
    Foam::fileFormats::STARCDsurfaceFormat

Description
    Read/write the surface shells from pro-STAR vrt/cel files.

Note
    Uses the extension @a .inp (input) to denote the format.

See Also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDsurfaceFormat_H
#define STARCDsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "STARCDsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class STARCDsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class STARCDsurfaceFormat
:
    public MeshedSurface<Face>,
    public STARCDsurfaceFormatCore
{
    // Private Data

        //- STAR-CD identifier for shell shapes (2d elements)
        static const int starcdShellShape_ = 3;

        //- STAR-CD identifier for shell type (shells vs. baffles)
        static const int starcdShellType_  = 4;


    // Private Member Functions

        static inline void writeShell
        (
            Ostream&,
            const Face&,
            const label cellId,
            const label cellTableId
        );

        //- Disallow default bitwise copy construct
        STARCDsurfaceFormat(const STARCDsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const STARCDsurfaceFormat<Face>&);


public:

    // Constructors

        //- Construct from file name
        STARCDsurfaceFormat(const fileName&);


    // Selectors

        //- Read file and return surface
        static autoPtr<MeshedSurface<Face> > New(const fileName& name)
        {
            return autoPtr<MeshedSurface<Face> >
            (
                new STARCDsurfaceFormat<Face>(name)
            );
        }


    //- Destructor
    virtual ~STARCDsurfaceFormat()
    {}


    // Member Functions

        //- Write surface mesh components by proxy
        static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, MeshedSurfaceProxy<Face>(*this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "STARCDsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
