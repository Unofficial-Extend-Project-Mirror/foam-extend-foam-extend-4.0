    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();
    if (args.options().found("sourceTime"))
    {
        if ((args.options()["sourceTime"]) == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            scalar sourceTime =
                readScalar(IStringStream(args.options()["sourceTime"])());
            sourceTimeIndex = getTimeIndex(sourceTimes, sourceTime);
        }
    }
    else
    {
        scalar sourceTime = runTimeTarget.time().value();
        sourceTimeIndex = getTimeIndex(sourceTimes, sourceTime);
    }
