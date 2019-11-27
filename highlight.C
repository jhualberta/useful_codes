{
   TH2D h;
   h.SetBins(10,0,10,10,0,10);
   h.Fill(1,1,1);
   h.Fill(3,3,3);
   h.Fill(5,5,5);
   h.Draw("TEXT");

   // Highlight the maximum
   Int_t bin = h.GetMaximumBin();
   Int_t binx, biny, binz;
   h.GetBinXYZ(bin, binx, biny, binz);
   Double_t x1 = h.GetXaxis()->GetBinLowEdge(binx);
   Double_t x2 = h.GetXaxis()->GetBinUpEdge(binx);
   Double_t y1 = h.GetYaxis()->GetBinLowEdge(biny);
   Double_t y2 = h.GetYaxis()->GetBinUpEdge(biny);
   TBox b(x1, y1, x2, y2);
   b.SetFillStyle(0);
   b.SetLineWidth(4);
   b.SetLineColor(kRed);
   b.Draw();
}
