/*
    chemPointList_[n].node_->vT() = chemPointList_[n].node_->elementRight_->v0() - chemPointList_[n].node_->elementLeft_->v0();

    scalarField fiMedium = 0.5*(chemPointList_[n].node_->elementLeft_->v0() + chemPointList_[n].node_->elementRight_->v0());
    chemPointList_[n].node_->a() = sum(fiMedium*chemPointList_[n].node_->vT());  
*/

    scalarField vsR = scaleComposition(chemPointList_[n]->node_->elementRight_->v0());    
    scalarField vsL = scaleComposition(chemPointList_[n]->node_->elementLeft_->v0());
    
    chemPointList_[n]->node_->vT() = vsR - vsL;

    scalarField fiMedium = 0.5*(vsL + vsR);
    chemPointList_[n]->node_->a() = sum(fiMedium*chemPointList_[n]->node_->vT());  

