    if (args.options().found("time"))
    {
        scalar timeValue(readScalar(IStringStream(args.options()["time"])()));

        startTime = Time::findClosestTimeIndex(Times, timeValue);
        endTime = startTime + 1;
    }
