int USERD_get_timeset_description
(
    int timeset_number,
    char timeset_description[Z_BUFL]
)
{
   
#ifdef ENSIGHTDEBUG
    Info << "Entering: USERD_get_timeset_description" << endl 
        << flush;
#endif

    if (TimeList[1].value() < 0)
    {
        strncpy(timeset_description, "CAD", Z_BUFL);
    }
    else
    {
        strncpy(timeset_description, "seconds", Z_BUFL);
    }

#ifdef ENSIGHTDEBUG
    Info << "Leaving: USERD_get_timeset_description" << endl 
        << flush;
#endif

    return Z_OK;
}
