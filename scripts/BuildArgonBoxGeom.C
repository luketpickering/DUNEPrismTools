#include "TGeoManager.h"

void BuildArgonBoxGeom(){

  TGeoManager *gman = new TGeoManager("world", "DUNEPrism LAr Box");
  TGeoMaterial *mat_vac = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMaterial *mat_lar = new TGeoMaterial("LAr",40,18,1.3954);

  TGeoMedium *med_vac = new TGeoMedium("Vacuum",1,mat_vac);
  TGeoMedium *med_lar = new TGeoMedium("LAr",1,mat_lar);

  const double m2cm = 100;

  TGeoVolume *vol_top = gGeoManager->MakeBox("Top",med_vac,40*m2cm,40*m2cm,40*m2cm);
  gGeoManager->SetTopVolume(vol_top);

  TGeoVolume *vol_larbox = gGeoManager->MakeBox("LAr",med_lar,18*m2cm,1.5*m2cm,2.5*m2cm);
  TGeoTranslation *tr_latoffset = new TGeoTranslation(-15*m2cm,0,0);

  TGeoTranslation *tr_latoffset_z = new TGeoTranslation(0,0,15*m2cm);

  vol_top->AddNode(vol_larbox,1,tr_latoffset);
  vol_top->AddNode(vol_larbox,1,tr_latoffset_z);

  gman->CloseGeometry();

  vol_top->SetLineColor(kMagenta);
  vol_larbox->SetLineColor(kBlue);
  vol_larbox->SetVisibility(true);
  gman->SetTopVisible();

  gman->Export("DUNEPrismLArBox.geom.root");
  gman->Export("DUNEPrismLArBox.geom.gdml");

  vol_top->Draw();
}
