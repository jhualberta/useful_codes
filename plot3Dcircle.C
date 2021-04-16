{
   c1 = new TCanvas("c1","PolyLine3D & PolyMarker3D Window",200,10,500,500);

   // Creating a view
   TView3D *view = TView::CreateView(1);
   view->SetRange(5,5,5,25,25,25);

   // Create a first PolyLine3D
   TPolyLine3D *pl3d1 = new TPolyLine3D(5);
   pl3d1->SetPoint(0, 10, 10, 10);
   pl3d1->SetPoint(1, 15, 15, 10);
   pl3d1->SetPoint(2, 20, 15, 15);
   pl3d1->SetPoint(3, 20, 20, 20);
   pl3d1->SetPoint(4, 10, 10, 20);

   // Create a first PolyMarker3D
   TPolyMarker3D *pm3d1 = new TPolyMarker3D(12);
   pm3d1->SetPoint(0, 10, 10, 10);
   pm3d1->SetPoint(1, 11, 15, 11);
   pm3d1->SetPoint(2, 12, 15, 9);
   pm3d1->SetPoint(3, 13, 17, 20);
   pm3d1->SetPoint(4, 14, 16, 15);
   pm3d1->SetPoint(5, 15, 20, 15);
   pm3d1->SetPoint(6, 16, 18, 10);
   pm3d1->SetPoint(7, 17, 15, 10);
   pm3d1->SetPoint(8, 18, 22, 15);
   pm3d1->SetPoint(9, 19, 28, 25);
   pm3d1->SetPoint(10, 20, 12, 15);
   pm3d1->SetPoint(11, 21, 12, 15);
   pm3d1->SetMarkerSize(2);
   pm3d1->SetMarkerColor(4);
   pm3d1->SetMarkerStyle(2);

   // Draw
   pl3d1->Draw();
   pm3d1->Draw();
}
