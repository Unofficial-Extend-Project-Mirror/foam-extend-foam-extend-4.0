#include "dictionary.H"
#include "fileNameList.H"
#include "IFstream.H"
#include "OSspecific.H"

using namespace Foam;

int main()
{
    Info << "\nReading Roots" << endl;

    IFstream rootsFile(home()/".foam/apps/openDX/roots");
    fileNameList rootsList(dictionary(rootsFile).lookup("roots"));

    char** rootsStrings = new char*[rootsList.size() + 1];
    rootsStrings[rootsList.size()] = 0;

    if (rootsList.size())
    {
        for (int i=0; i<rootsList.size(); i++)
        {
            rootsStrings[i] = new char[rootsList[i].size() + 1];
            strcpy(rootsStrings[i], rootsList[i].c_str());

            Info<< rootsStrings[i] << endl;
        }
    }

    return 0;
}
