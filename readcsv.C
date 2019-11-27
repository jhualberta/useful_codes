#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
void readcsv(void) {
  const int histosize = 252;

  ULong64_t timeStamp[4096]; //l, 64 bit unsigned integers
  UShort_t energy[1024]; //s, 16 bit unsigned integers
  UShort_t energyShort[1024]; //s, 16 bit unsigned integers
  UShort_t flag[1024]; //s, 16 bit unsigned integers
  //UInt_t nWave[2048]; //i, 32 bit unsigned integers
  vector<UShort_t> wave[histosize]; //s, 16 bit unsigned integers

//  FILE *fin = fopen("ch0_2-5-1616_Data_test_dump.bin", "r");
  FILE *fin = fopen("line.csv", "r");

  if (!fin) {
    printf("Error : data not found!\n");
    return;
  }
  TFile *fout = TFile::Open("dump_ch0.root", "recreate");
  TTree *tree = new TTree("tree", "TTree ch0");
  // "==> Case A” in http://root.cern.ch/root/html/TTree.html 17
  //tree->Branch("Board", board, "board[1024]/s");
  //tree->Branch("Channel", channel, "channel[1024]/s");
  tree->Branch("TimeStamp", timeStamp, "timeStamp[4096]/s");
  tree->Branch("Energy", energy, "energy[1024]/s");
  tree->Branch("EnergyShort", energyShort, "energyShort[1024]/s");
  tree->Branch("Flag", flag, "flag[1024]/s");
  //tree->Branch("Nwave", nWave, "nWave[2048]/s");
  tree->Branch("Wave", wave, "wave[1024]/s");
  int eventID = 0;
  while ( fin )  {





  }
  tree->Fill();
  }
  fclose(fin); // no longer needed
  tree->Write();
  delete fout; // automatically deletes the "tree”, too
  return;
}
