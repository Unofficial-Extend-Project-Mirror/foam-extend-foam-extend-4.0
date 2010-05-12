
        bubbleProperties.remove("Vb");
        Vb = gSum((1.0 - interface.fluidIndicator().internalField())*mesh.V());
        bubbleProperties.add("Vb", Vb);

        if(Pstream::master())
        {
            Info << "Bubble volume reduction: " << Vb/Vbf << endl;
        }
