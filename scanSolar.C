{

TFile f("FilteredSolar_oct2018_rat6174.root");
TTree* T = (TTree*)f->Get("T");
//T->Scan("runNumber:eventGTID:Nhits:energy:posX:posY:(posZ-108):dirX:dirY:dirZ:sunDirX:sunDirY:sunDirZ","sqrt(posX**2+posY**2+(posZ-108)**2)<5500 && itr>0.55 && iso>-0.12 && iso<0.95 && Nhits>30");

T->Scan("energy","sqrt(posX**2+posY**2+(posZ-108)**2)<5500 && itr>0.55 && iso>-0.12 && iso<0.95 && Nhits>30");

}
