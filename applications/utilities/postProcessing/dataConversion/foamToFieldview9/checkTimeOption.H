    if (args.options().found("time"))
    {
        scalar time(readScalar(IStringStream(args.options()["time"])()));

        int nearestIndex = -1;
        scalar nearestDiff = Foam::GREAT;

        forAll(Times, timeIndex)
        {
            scalar diff = fabs(Times[timeIndex].value() - time);
            if (diff < nearestDiff)
            {
                nearestDiff = diff;
                nearestIndex = timeIndex;
            }
        }

        startTime = nearestIndex;
        endTime = nearestIndex + 1;
    }
    else
    {
        startTime = 0;
        endTime = Times.size();
    }

