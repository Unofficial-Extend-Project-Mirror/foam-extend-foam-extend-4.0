// write time values to case file
{
    scalar timeCorrection = 0;
    if (timeDirs[0].value() < 0)
    {
        timeCorrection = - timeDirs[0].value();
        Info<< "Correcting time values. Adding " << timeCorrection << endl;
    }

    caseFile.setf(ios_base::scientific, ios_base::floatfield);
    caseFile.precision(5);

    // time set 1 - geometry and volume fields
    if (fieldFileNumbers.size())
    {
        caseFile
            << "time set:        " << 1 << nl
            << "number of steps: " << fieldFileNumbers.size() << nl
            << "filename numbers:" << nl;

        label count = 0;
        forAll (fieldFileNumbers, i)
        {
            caseFile
                << " " << setw(12) << fieldFileNumbers[i];

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }

        caseFile
            << nl << "time values:" << nl;

        count = 0;
        forAll (fieldFileNumbers, i)
        {
            caseFile
                << " " << setw(12)
                << timeIndices[fieldFileNumbers[i]] + timeCorrection;

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }
        caseFile << nl << nl;
    }

    // time set 2 - lagrangian fields
    if (hasLagrangian && sprayFileNumbers.size())
    {
        caseFile
            << "time set:        " << 2 << nl
            << "number of steps: " << sprayFileNumbers.size() << nl
            << "filename numbers:" << nl;

        label count = 0;
        forAll (sprayFileNumbers, i)
        {
            caseFile
                << " " << setw(12) << sprayFileNumbers[i];

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }

        caseFile
            << nl << "time values:" << nl;

        count = 0;
        forAll (sprayFileNumbers, i)
        {
            caseFile
                << " " << setw(12)
                << timeIndices[sprayFileNumbers[i]] + timeCorrection;

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }
        caseFile << nl << nl;
    }


    caseFile << "# end" << nl;
}
