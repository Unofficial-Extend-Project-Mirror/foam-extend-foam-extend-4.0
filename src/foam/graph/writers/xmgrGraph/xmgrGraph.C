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

#include "xmgrGraph.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::xmgrGraph, 0);
const Foam::word Foam::xmgrGraph::ext_("agr");

namespace Foam
{
    typedef graph::writer graphWriter;
    addToRunTimeSelectionTable(graphWriter, xmgrGraph, word);
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::xmgrGraph::write(const graph& g, Ostream& os) const
{
    os  << "@title " << g.title() << endl
        << "@xaxis label " << g.xName() << endl
        << "@yaxis label " << g.yName() << endl;

    label fieldI = 0;

    for (graph::const_iterator iter = g.begin(); iter != g.end(); ++iter)
    {
        os  << "@s" << fieldI << " legend "
            << iter()->name() << endl
            << "@target G0.S" << fieldI << endl
            << "@type xy" << endl;

        writeXY(g.x(), *iter(), os);

        os << endl;

        fieldI++;
    }
}


// ************************************************************************* //
