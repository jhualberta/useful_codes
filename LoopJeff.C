/*Jeff Tseng*/

void CountFits(const char* filename) {
  int naf = 0;
  RAT::DU::DSReader ds(filename);
  for (size_t ientry = 0; ientry < ds.GetEntryCount(); ientry++) {
    const RAT::DS::Entry& entry = ds.GetEntry(ientry);
    for (size_t iev = 0; iev < entry.GetEVCount(); iev++) {
      const RAT::DS::EV& ev = entry.GetEV(iev);
      vector<string> s = ev.GetFitNames();
      if (ientry < 100) {
        cout << ientry << "(" << ev.GetNhits() << "," << ev.GetNhitsCleaned() << "):  ";
        for (size_t i = 0; i < s.size(); i++) cout << s[i] << " ";
        cout << endl;
      }
      if (ev.FitResultExists("albertaFitter")) {
        naf++;
        const RAT::DS::FitResult& result = ev.GetFitResult("albertaFitter");
        const RAT::DS::FitVertex& vertex = result.GetVertex(0);
        TVector3 v = vertex.GetPosition();
        cout << "alberta:  (" << v.X() << "," << v.Y() << "," << v.Z() << ") " << vertex.GetTime() << endl;
      }
      if (ev.FitResultExists("MultiPathProcessor")) {
        naf++;
        const RAT::DS::FitResult& result = ev.GetFitResult("MultiPathProcessor");
        const RAT::DS::FitVertex& vertex = result.GetVertex(0);
        TVector3 v = vertex.GetPosition();
        cout << "MPWProc:  (" << v.X() << "," << v.Y() << "," << v.Z() << ") " << vertex.GetTime() << endl;
      }
    }
    if (naf > 100) break;
  }
  cout << "Number of events = " << naf;
}
