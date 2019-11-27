fLightPath.CalcByPositionPartial( startPos, DU::Utility::Get()->GetPMTInfo().GetPosition( pmt.GetID() ) );
 double distInUpperTarget = fLightPath.GetDistInUpperTarget();
 double distInLowerTarget = fLightPath.GetDistInLowerTarget();
 double distInAV = fLightPath.GetDistInAV();
 double distInWater = fLightPath.GetDistInWater();

 const double transitTime = DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );
