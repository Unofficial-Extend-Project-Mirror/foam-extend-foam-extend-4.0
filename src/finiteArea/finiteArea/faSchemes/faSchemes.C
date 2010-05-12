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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "faSchemes.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

int faSchemes::debug(Foam::debug::debugSwitch("faSchemes", false));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faSchemes::faSchemes(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "faSchemes",
            obr.time().system(),
            obr,
//             IOobject::MUST_READ,
            IOobject::READ_IF_PRESENT,  // Allow default dictionary creation
            IOobject::NO_WRITE
        )
    ),
    ddtSchemes_(ITstream("ddtSchemes", tokenList())()),
    defaultDdtScheme_("default", tokenList()),
    d2dt2Schemes_(ITstream("d2dt2Schemes", tokenList())()),
    defaultD2dt2Scheme_("default", tokenList()),
    interpolationSchemes_(ITstream("interpolationSchemes", tokenList())()),
    defaultInterpolationScheme_("default", tokenList()),
    divSchemes_(ITstream("divSchemes", tokenList())()),
    defaultDivScheme_("default", tokenList()),
    gradSchemes_(ITstream("gradSchemes", tokenList())()),
    defaultGradScheme_("default", tokenList()),
    lnGradSchemes_(ITstream("lnGradSchemes", tokenList())()),
    defaultLnGradScheme_("default", tokenList()),
    laplacianSchemes_(ITstream("laplacianSchemes", tokenList())()),
    defaultLaplacianScheme_("default", tokenList()),
    fluxRequired_(ITstream("fluxRequired", tokenList())())
{
    if (!headerOk())
    {
        if (debug)
        {
            InfoIn
            (
                "faSchemes::faSchemes(const objectRegistry& obr)"
            )   << "faSchemes dictionary not found.  Creating default."
                << endl;
        }
    }

    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool faSchemes::read()
{
    bool readOk = false;

    if (headerOk())
    {
        readOk = regIOobject::read();
    }

    if (readOk)
    {
        const dictionary& dict = schemesDict();

        // ddt Schemes
        if (dict.found("ddtSchemes"))
        {
            ddtSchemes_ = dict.subDict("ddtSchemes");
        }
        else
        {
            ddtSchemes_.add("default", "none");
        }

        if
        (
            ddtSchemes_.found("default")
         && word(ddtSchemes_.lookup("default")) != "none"
        )
        {
            defaultDdtScheme_ = ddtSchemes_.lookup("default");
        }


        // d2dt2 schemes
        if (dict.found("d2dt2Schemes"))
        {
            d2dt2Schemes_ = dict.subDict("d2dt2Schemes");
        }
        else if (dict.found("timeScheme"))
        {
            // For backward compatibility.
            // The timeScheme will be deprecated with warning or removed in 2.4.

            word timeSchemeName(dict.lookup("timeScheme"));

            if (timeSchemeName == "EulerImplicit")
            {
                timeSchemeName = "Euler";
            }
            else if (timeSchemeName == "SteadyState")
            {
                timeSchemeName = "steadyState";
            }

            if (d2dt2Schemes_.found("default"))
            {
                d2dt2Schemes_.remove("default");
            }

            d2dt2Schemes_.add("default", timeSchemeName);
        }
        else
        {
            d2dt2Schemes_.add("default", "none");
        }

        if
        (
            d2dt2Schemes_.found("default")
         && word(d2dt2Schemes_.lookup("default")) != "none"
        )
        {
            defaultD2dt2Scheme_ = d2dt2Schemes_.lookup("default");
        }


        // Interpolation schemes
        if (dict.found("interpolationSchemes"))
        {
            interpolationSchemes_ = dict.subDict("interpolationSchemes");
        }
        else
        {
            interpolationSchemes_.add("default", "linear");
        }

        if
        (
            interpolationSchemes_.found("default")
         && word(interpolationSchemes_.lookup("default")) != "none"
        )
        {
            defaultInterpolationScheme_ =
                interpolationSchemes_.lookup("default");
        }


        // Div schemes
        if (dict.found("divSchemes"))
        {
            divSchemes_ = dict.subDict("divSchemes");

        }
        else
        {
            divSchemes_.add("default", "none");
        }

        if
        (
            divSchemes_.found("default")
         && word(divSchemes_.lookup("default")) != "none"
        )
        {
            defaultDivScheme_ = divSchemes_.lookup("default");
        }


        // Grad schemes
        if (dict.found("gradSchemes"))
        {
            gradSchemes_ = dict.subDict("gradSchemes");

        }
        else
        {
            gradSchemes_.add("default", "none");
        }

        if
        (
            gradSchemes_.found("default")
         && word(gradSchemes_.lookup("default")) != "none"
        )
        {
            defaultGradScheme_ = gradSchemes_.lookup("default");
        }


        // lnGrad schemes
        if (dict.found("lnGradSchemes"))
        {
            lnGradSchemes_ = dict.subDict("lnGradSchemes");
        }
        else
        {
            lnGradSchemes_.add("default", "corrected");
        }

        if
        (
            lnGradSchemes_.found("default")
         && word(lnGradSchemes_.lookup("default")) != "none"
        )
        {
            defaultLnGradScheme_ = lnGradSchemes_.lookup("default");
        }


        // laplacian schemes
        if (dict.found("laplacianSchemes"))
        {
            laplacianSchemes_ = dict.subDict("laplacianSchemes");

        }
        else
        {
            laplacianSchemes_.add("default", "none");
        }

        if
        (
            laplacianSchemes_.found("default")
         && word(laplacianSchemes_.lookup("default")) != "none"
        )
        {
            defaultLaplacianScheme_ = laplacianSchemes_.lookup("default");
        }


        // Flux required
        if (dict.found("fluxRequired"))
        {
            fluxRequired_ = dict.subDict("fluxRequired");
        }
    }

    return readOk;
}


const dictionary& faSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


ITstream& faSchemes::ddtScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup ddtScheme for " << name << endl;
    }

    if
    (
        ddtSchemes_.found(name)
     || !defaultDdtScheme_.size()
    )
    {
        return ddtSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDdtScheme_).rewind();
        return const_cast<ITstream&>(defaultDdtScheme_);
    }
}


ITstream& faSchemes::d2dt2Scheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup d2dt2Scheme for " << name << endl;
    }

    if
    (
        d2dt2Schemes_.found(name)
     || !defaultD2dt2Scheme_.size()
    )
    {
        return d2dt2Schemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultD2dt2Scheme_).rewind();
        return const_cast<ITstream&>(defaultD2dt2Scheme_);
    }
}


ITstream& faSchemes::interpolationScheme(const word& name) const
{
    if
    (
        interpolationSchemes_.found(name)
     || !defaultInterpolationScheme_.size()
    )
    {
        return interpolationSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultInterpolationScheme_).rewind();
        return const_cast<ITstream&>(defaultInterpolationScheme_);
    }
}


ITstream& faSchemes::divScheme(const word& name) const
{
    if (divSchemes_.found(name) || !defaultDivScheme_.size())
    {
        return divSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDivScheme_).rewind();
        return const_cast<ITstream&>(defaultDivScheme_);
    }
}


ITstream& faSchemes::gradScheme(const word& name) const
{
    if (gradSchemes_.found(name) || !defaultGradScheme_.size())
    {
        return gradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultGradScheme_).rewind();
        return const_cast<ITstream&>(defaultGradScheme_);
    }
}


ITstream& faSchemes::lnGradScheme(const word& name) const
{
    if (lnGradSchemes_.found(name) || !defaultLnGradScheme_.size())
    {
        return lnGradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultLnGradScheme_).rewind();
        return const_cast<ITstream&>(defaultLnGradScheme_);
    }
}


ITstream& faSchemes::laplacianScheme(const word& name) const
{
    if (laplacianSchemes_.found(name) || !defaultLaplacianScheme_.size())
    {
        return laplacianSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultLaplacianScheme_).rewind();
        return const_cast<ITstream&>(defaultLaplacianScheme_);
    }
}


bool faSchemes::fluxRequired(const word& name) const
{
    return fluxRequired_.found(name);
}


bool faSchemes::writeData(Ostream& os) const
{
    // Write dictionaries
    os << nl << "ddtSchemes";
    ddtSchemes_.write(os, true);

    os << nl << "d2dt2Schemes";
    d2dt2Schemes_.write(os, true);

    os << nl << "interpolationSchemes";
    interpolationSchemes_.write(os, true);

    os << nl << "divSchemes";
    divSchemes_.write(os, true);

    os << nl << "gradSchemes";
    gradSchemes_.write(os, true);

    os << nl << "lnGradSchemes";
    lnGradSchemes_.write(os, true);

    os << nl << "laplacianSchemes";
    laplacianSchemes_.write(os, true);

    os << nl << "fluxRequired";
    fluxRequired_.write(os, true);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
