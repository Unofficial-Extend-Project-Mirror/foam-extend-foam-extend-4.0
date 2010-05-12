if(Pstream::master())
{
    pathFile << runTime.value() << tab;

    pathFile << XF.value().x() << tab
        << XF.value().y() << tab
        << XF.value().z() << tab
        << mag(UF.value()) << tab
        << mag(aF.value()) << endl;
}
