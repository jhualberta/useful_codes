import ROOT
from ROOT import TGeoManager

geoManager = ROOT.TGeoManager("toto","toto");
material = ROOT.TGeoMaterial("Vaccum",0,0,0);
medium = ROOT.TGeoMedium("Vaccum",1,material);
# worldShape = ROOT.TGeoSphere(0.,10.);
worldShape = ROOT.TGeoCone(10.,0.,0.,0.,10.);
world = ROOT.TGeoVolume("top",worldShape,medium);
geoManager.SetTopVolume(world);
geoManager.CloseGeometry();
world.SetLineColor(8);
geoManager.SetTopVisible();
world.Draw()
raw_input()
