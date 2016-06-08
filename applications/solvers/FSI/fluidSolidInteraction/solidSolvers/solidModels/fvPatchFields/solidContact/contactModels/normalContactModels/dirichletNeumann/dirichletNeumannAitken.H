       // Aitken's method for under-relaxation
   if (curTimeIndex_ != mesh.time().timeIndex())
     {
       //Info << "aitkenRes0_ is " << aitkenRes0_ << endl;
       aitkenRes0_ = 1.0;
       aitkenTheta_ = relaxFactor_;
     }

       if (iCorr_ == 1)
     {
       //const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
       //const fvPatchField<vector>& dispField =
       //slavePatch.lookupPatchField<volVectorField, vector>(fieldName);
       const volVectorField& dispField =
         mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
       //all residuals will be normalised to the max mag of the DU/U field
       aitkenRes0_ = gMax(mag(dispField.internalField()));
       //Info << "updating aitkenRes0_ to " << aitkenRes0_ << endl;
     }

       // update delta
       // normalised with repsect of aitkenRes0 to avoid very small numbers
       slaveDisp_ = slaveFaceNormals*slaveDispMag;
       aitkenDelta_ = (slaveDisp_ - slaveDispPrevIter_) / aitkenRes0_;

       // first iteration set to relaxFactor
       if (iCorr_ > 0)
     {
       //Info << "slaveDisp_ " << slaveDisp_ << endl;
       vectorField b = aitkenDelta_ - aitkenDeltaPrevIter_;
       // scalar sumMagB = gSum(mag(b));
       scalar sumMagB = gSum(magSqr(b));
       if (sumMagB < SMALL)
         {
           Warning << "Dirichlet Normal Contatc Aitken under-relaxation: denominator less then SMALL"
               << endl;
           sumMagB += SMALL;
         }
       aitkenTheta_ = -aitkenTheta_*
         gSum(aitkenDeltaPrevIter_ & b)
         /
         sumMagB;
     }

       // update disp with Aitken correction
       slaveDisp_ += aitkenTheta_*aitkenDelta_*aitkenRes0_;

       // update previous iteration value
       slaveDispPrevIter_ = slaveDisp_;
       aitkenDeltaPrevIter_ = aitkenDelta_;
