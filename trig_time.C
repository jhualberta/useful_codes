{
float trigger_time[10000];
int flag;
float trig_counts;
int evtNum = 10000;
int recordLength = 252;
float attenuation = 0.75;//0.25;
float delay = 1;
float threshold[10000];
float baseline[10000];
float baseline2[10000];

float energy[10000];

TFile *f0 = new TFile("../data/feb22/dump_ch0_Data_TeSOP_241Am_coin248ns_2200V_3_0.root");

    vector<TH1F*> hlistProcessWf_ch0;
    //vector<TH1F*> hlistWaveform_ch1;
    for(int i = 0;i<evtNum;i++)
    {
     TH1F *hprocess = new TH1F("hprocess","", recordLength, 0, recordLength);
     //histProcessWf_ch0.push_back((TH1F*)f0->Get(Form("hwf%u",i)));
     hlistProcessWf_ch0.push_back(hprocess);
     hlistProcessWf_ch0[i]->SetName(Form("hpwf%u",i));
     //delete htemp;
    }

    
    //vector<TH1F*> hlistSsignal_ch1;
//trigger

vector<TH1F*> hlistSsignal_ch0;


TH1F *trig = new TH1F("triger_time","", 100, 0, 100);
TH1F *hEnergy = new TH1F("energy spectrum","", 1000, 0, 12000);


for(int i = 0; i<evtNum; i++){	
	flag = 0;
		//calculate baseline[i]	
		baseline[i] = ((TH1F*)f0->Get(Form("hwf%u",i)))->Integral(0,32)/32;
		baseline2[i] = baseline[i]*(attenuation-1);
		//cout<<"baseline2[i] = "<<baseline2[i]<<endl;
	for(int q = 0;q<recordLength;q++){
		//invert and delay and fill histogram


        hlistProcessWf_ch0[i]->SetBinContent(q,(((TH1F*)f0->Get(Form("hwf%u",i)))->GetBinContent(q+1+delay))*-1  +  attenuation * (((TH1F*)f0->Get(Form("hwf%u",i)))->GetBinContent(q+1))); 
		hlistProcessWf_ch0[i]->SetBinContent(252,baseline2[i]);  
	}


	threshold[i] = (hlistProcessWf_ch0[i]->GetMaximum()- hlistProcessWf_ch0[i]->GetMinimum())/3;
	//cout<<"threshold[i] = "<<threshold[i]<<endl;

	for(int q = 0;q<recordLength;q++){
		//find trigger time
 		trig_counts = hlistProcessWf_ch0[i]->GetBinContent(q);
		if((trig_counts < (baseline2[i] - threshold[i])) and flag == 0){flag = 1;}
		if((flag == 1) and (trig_counts > baseline2[i])){trigger_time[i] = q-1;break;}  
	}

	if((trigger_time[i] == 0) or trigger_time[i] >50.){
		//cout<<"i "<<i<<" trigger_time "<<trigger_time[i]<<" baseline2[i] "<<baseline2[i]<<" threshold[i] "<<threshold[i]<<endl;
	}
	//cout<<"i = "<<i<<" trigger_time "<<trigger_time[i]<<endl;
	trig->Fill(trigger_time[i]);
	energy[i] = ((TH1F*)f0->Get(Form("hwf%u",i)))->Integral(trigger_time[i]-8,trigger_time[i]+7) - baseline[i]*15;
	hEnergy->Fill(energy[i]);

}


TFile *fprocess = new TFile("processedWf.root","recreate");
fprocess->cd();
for(int i = 0;i<hlistProcessWf_ch0.size();i++)
{
  hlistProcessWf_ch0[i]->Write();
}
fprocess->Close();


hEnergy->Draw();

//+8ns, short gate 15ns
}
