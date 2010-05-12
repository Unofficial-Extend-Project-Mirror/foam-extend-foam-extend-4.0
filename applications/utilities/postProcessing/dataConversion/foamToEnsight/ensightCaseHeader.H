if (Pstream::master())
{
    ensightCaseFile << "FORMAT" << nl;
    ensightCaseFile << "type: ensight gold" << nl << nl;
}
