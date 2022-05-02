#include"libs.hh"
#include <TMathBase.h>

using namespace std;

const double NSIGMA = 5;
int DET_NUM = 4;

bool is_good_strip(UInt_t det, UInt_t strip)
{
    if(det==1)
    {
        if(        strip==507 || strip==509 || strip==510  ||
                strip==508 || strip==512 || strip==94   ||
                strip==93  || strip==629 || strip==630  ||
                strip==627  || strip==632 || strip==630  ||
                strip==450 || strip==455 
          )
            return false;
    }
    else if(det==7)
    {
        if(        strip==343
          )
            return false;
    }
    //if((i%64)>62 || (i%64)<3) return false;
    return true;
}

void FOOT_ana(int firstEvent, int max_events)
{



    TApplication* theApp = new TApplication("App", 0, 0);

    TCanvas * canvas = new TCanvas("canvas","canvas",2000,1000);
    canvas->Divide(3,2);

    TH2D * h2d_raw = new TH2D("h2d_raw","Raw pedestals",640,1,641,700,0,700);

    TH1D * h1d_coarse= new TH1D("h1d_coarse","global_baseline",100,-50,50);
    TH2D * h2d_coarse= new TH2D("h2d_coarse","Only pedestal subtraction",640,1,641,100,-50,50);

    TH1D * h1d_fine = new TH1D("h1d_fine","Global baseline",100,-50,50);
    TH2D * h2d_fine = new TH2D("h2d_fine","Fine correction of baseline",640,1,641,200,-50,50);

    TH1D * h1d_mul= new TH1D("h1d_mul","Strip multiplicity per trigger",21,-0.5,20.5);
    //============== Drawing raw data ===============
    TChain * ch = new TChain("h101");
    //TString filename = Form("/lustre/land/vpanin/FOOT/pedestals_10FOOT_biason_02052022.root",DET_NUM);
    TString filename = Form("/lustre/land/vpanin/FOOT/pedestals_10FOOT_biason_02052022.root");

    //ch->Add("../../FOOT_test_08042022/test_08_04_2022_FOOT8_biason.root");
    ch->Add(filename.Data());
    UInt_t  FOOT10I[640];
    UInt_t  FOOT10E[640];
    ch->SetBranchAddress("FOOT10I",FOOT10I);
    ch->SetBranchAddress("FOOT10E",FOOT10E);
    Int_t Nevents = ch->GetEntries();
    if(max_events>0) Nevents = max_events; 

    Int_t     mul_good_hits = 0;
    Int_t     mul_good_events = 0;
    Double_t  mean_ssd=0 ;
    Double_t  asic_offset[10];
    Double_t  asic_offset_fine[10];
    Double_t  high_mul_events =0 ;
    Double_t  signal = 0;
    int       counter_asic =0;
    Int_t     stat=0;
    Int_t     i=0;

    //Start from drawing raw pedestals
    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        ch->GetEntry(ev);
        for(i=0; i<640; i++)
        {
            if(is_good_strip(DET_NUM,FOOT10I[i]))
            {
                h2d_raw->Fill(FOOT10I[i],FOOT10E[i]);
            }
            else
                h2d_raw->Fill(FOOT10I[i],0);

        }
    }    

    ////ch->Draw("FOOT10E:FOOT10I>>h(640,1,641,500,0,500)","","colz");
    canvas->cd(1);
    gPad->SetLogz();
    h2d_raw->GetXaxis()->SetTitle("Strip No.");
    h2d_raw->GetYaxis()->SetTitle("ADC channel");
    h2d_raw->Draw("colz");
    h2d_raw->GetYaxis()->SetRangeUser(0,800);
    for(int i_asic=1; i_asic<10; i_asic++)
    {
        TLine* l = new TLine(64.5*i_asic,0,64.5*i_asic,800);
        l->Draw("same");
        l->SetLineStyle(7);
            l->SetLineWidth(1);
            l->SetLineColor(13);
        }

        ////============== Slicing, fitting, drawing params ===============
        TF1* foo = new TF1("foo","gaus",0,1000);
        h2d_raw->FitSlicesY(foo,1,640,0,"QNR",0);
        TH1D * h2d_raw_1 = (TH1D*)gDirectory->Get("h2d_raw_1")->Clone("h2d_raw_1");
        TH1D * h2d_raw_2 = (TH1D*)gDirectory->Get("h2d_raw_2")->Clone("h2d_raw_2");
        h2d_raw_1->SetMarkerStyle(kFullCircle);
        h2d_raw_1->SetMarkerSize(0.2);
        h2d_raw_1->SetMarkerColor(kBlack);
        h2d_raw_1->SetLineColor(kBlack);
        h2d_raw_1->Draw("same");

        canvas->cd(2);
        h2d_raw_2->SetMarkerStyle(kFullCircle);
        h2d_raw_2->SetMarkerSize(0.3);
        h2d_raw_2->SetMarkerColor(kBlue);
        h2d_raw_2->SetLineColor(kBlue);
        h2d_raw_2->GetYaxis()->SetRangeUser(-2,10);
        h2d_raw_2->GetXaxis()->SetTitle("Strip No.");
        h2d_raw_2->GetYaxis()->SetTitle("ADC sigma");
        h2d_raw_2->SetTitle("Sigmas before baseline correction");
        h2d_raw_2->Draw();
        for(int i_asic=1; i_asic<10; i_asic++)
        {
            TLine* l = new TLine(64.5*i_asic,-2,64.5*i_asic,10);
            l->Draw("same");
            l->SetLineStyle(7);
            l->SetLineWidth(1);
            l->SetLineColor(13);
        }

        TText *t_before = new TText(50,8,"Only pedetsal subtraction");
        t_before->SetTextColor(kBlue);
        t_before->Draw("same");
        TText *t_after = new TText(50,7,"After baseline correciton");
        t_after->SetTextColor(kRed);
        t_after->Draw("same");

        ////============== Writing pedestals ===============
        auto pedfilename = "pedestal.dat";
        std::ofstream fout;
        fout.open(pedfilename);
        for(Int_t i = 0; i<640; i++){
            fout << "1  " << i+1 << "  " <<  h2d_raw_1->GetBinContent(i+1) << "  " <<  h2d_raw_2->GetBinContent(i+1) << endl;
            cout << "1  " << i+1 << "  " <<  h2d_raw_1->GetBinContent(i+1) << "  " <<  h2d_raw_2->GetBinContent(i+1) << endl;
        }
        fout.close();

        //============== Reading pedestals ===============
        ifstream pedfile(pedfilename,ifstream::in);
        if ( !pedfile.is_open() ){
            cout << "Cannot open Ped file" << endl;
        }
        int Nlines = 0;
        std::string line;
        while( std::getline(pedfile, line) )
            Nlines++;
        pedfile.clear();
        pedfile.seekg( 0, std::ios::beg );
        Int_t DetId_[Nlines];
        Int_t StripId_[Nlines];
        Double_t Ped_[Nlines];
        Double_t Sig_[Nlines];
        for(Int_t i=0 ; i<Nlines ; i++){
            pedfile >> DetId_[i] >> StripId_[i] >> Ped_[i] >> Sig_[i];
        }
        pedfile.close();

        //============== Applying pedestals to root data ===============
        for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
        {
            cout << "\r-- Event # : " << ev << flush;
            ch->GetEntry(ev);

            //======== Global base line correction for this event 
            mean_ssd=0;
            stat=0;
            for(i=0; i<640; i++)
            {
                if(!is_good_strip(DET_NUM,FOOT10I[i])) continue; 
                signal = FOOT10E[i] - Ped_[i];
                stat++;
                mean_ssd += signal;

            }
            mean_ssd = mean_ssd/stat;

            if(mean_ssd>10){
                cout << "\nMean ssd = " << mean_ssd << endl;
                //   continue;
            }

            //======== Fine base line correction for individual asics
            stat=0;            
            counter_asic=0;            
            for(i=0; i<10; i++){  asic_offset_fine[i]=0; }//reset asic baselines
            for(i=0; i<640; i++)
            {
                if(is_good_strip(DET_NUM,FOOT10I[i]))
                {
                    signal = FOOT10E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic];
                    if(fabs(signal) < (NSIGMA * Sig_[i]) )
                    {
                        stat++;
                        asic_offset_fine[counter_asic] += signal;
                    }
                }
                if((FOOT10I[i]%64)==0) 
                {
                    asic_offset_fine[counter_asic] /= stat;
                    //cout << "\n Calculated asic_offset = " << asic_offset[counter_asic] << " in " << counter_asic << " asic" <<  endl; 
                    counter_asic++;
                    stat=0;
                }
            }

            //======== Get number of "good" hits
            counter_asic=0;
            mul_good_hits=0;
            double _asic;
            for(i=0; i<640; i++)
            {
                if(FOOT10I[i]%64 == 1 && FOOT10I[i]>1) counter_asic++;
                if(!is_good_strip(DET_NUM,FOOT10I[i])) continue;
                signal = FOOT10E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic] - asic_offset_fine[counter_asic];
                //    cout << "\n Signal " << signal << " in strip " << FOOT10I[i];
                //    cout << "\n Raw data = " <<  FOOT10E[i];
                //    cout << "\n Ped_[i] = " << Ped_[i];
                //    cout << "\n Mean ssd = " << mean_ssd;
                //    cout << "\n asic_offset = " << asic_offset[counter_asic] << " in " << counter_asic << " asic" <<  endl; 
                if(signal > (NSIGMA * Sig_[i]) ) mul_good_hits++;
                //if(signal > 2.5 * NSIGMA ) mul_good_hits++;
            }
            h1d_mul->Fill(mul_good_hits);
            if(mul_good_hits>0) mul_good_events++;
            if(mul_good_hits>1)
            {
                cout << "\n--High multiplitcity: " << mul_good_hits << endl;
                high_mul_events++;
            }

            //========= Filling final histograms
            counter_asic=0;
            for(i=0; i<640; i++)
            {
                if((FOOT10I[i]%64) == 1 && FOOT10I[i]>1) counter_asic++;
                if(!is_good_strip(DET_NUM,FOOT10I[i])) continue;
                signal = FOOT10E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic] - asic_offset_fine[counter_asic];

                h1d_coarse->Fill(FOOT10E[i] - Ped_[i]);
                h2d_coarse->Fill(FOOT10I[i], (FOOT10E[i] - Ped_[i]));

                h1d_fine  ->Fill(signal);
            h2d_fine  ->Fill(FOOT10I[i], signal);
        }

    }//end of eventloop
    cout << "\n--Identified " << mul_good_events << " good events" << " and " << high_mul_events<< "chrismas-tree events" <<  endl;
    
    canvas->cd(3);
    h2d_coarse->Draw("colz");
    gPad->SetLogz();
    h2d_coarse->GetXaxis()->SetTitle("Strip No.");
    h2d_coarse->GetYaxis()->SetTitle("ADC channel");



    for(int i_asic=1; i_asic<10; i_asic++)
    {
        TLine* l = new TLine(65*i_asic,-50,65*i_asic,50);
        l->Draw("same");
        l->SetLineStyle(7);
        l->SetLineWidth(1);
        l->SetLineColor(13);
    }


    canvas->cd(6);
    h2d_fine->Draw("colz");
    gPad->SetLogz();
    h2d_fine->GetXaxis()->SetTitle("Strip No.");
    h2d_fine->GetYaxis()->SetTitle("ADC channel");

    h2d_fine->FitSlicesY(foo,1,640,0,"QNR",0);
    TH1D * h2d_fine_1 = (TH1D*)gDirectory->Get("h2d_fine_1")->Clone("h2d_fine_1");
    TH1D * h2d_fine_2 = (TH1D*)gDirectory->Get("h2d_fine_2")->Clone("h2d_fine_2");
    for(int i_asic=1; i_asic<10; i_asic++)
    {
        TLine* l = new TLine(64*i_asic,-50,64*i_asic,50);
        l->Draw("same");
        l->SetLineStyle(7);
        l->SetLineWidth(1);
        l->SetLineColor(13);
    }

    canvas->cd(2);
    h2d_fine_2->GetYaxis()->SetRangeUser(-2,10);
    h2d_fine_2->SetMarkerStyle(kFullCircle);
    h2d_fine_2->SetMarkerSize(0.3);
    h2d_fine_2->SetMarkerColor(kRed);
    h2d_fine_2->SetLineColor(kRed);
    h2d_fine_2->SetTitle("Pedestal sigmas");
    h2d_fine_2->Draw("same");

    canvas->cd(4);
    h1d_coarse->SetLineColor(kBlue);
    h1d_coarse->SetLineWidth(2);

    h1d_fine->SetLineColor(kRed);
    h1d_fine->SetLineWidth(2);

    h1d_fine->Draw();
    h1d_fine->Fit("gaus");

    h1d_coarse->Draw("same");
    h1d_coarse->Fit("gaus");

    TF1 *fit_coarse = h1d_coarse->GetFunction("gaus");
    fit_coarse->SetLineColor(kBlue);

    h1d_fine->GetXaxis()->SetTitle("ADC channel");
    h1d_fine->GetYaxis()->SetTitle("Counts");

    canvas->cd(5);
    gPad->SetLogy();
    h1d_mul->GetXaxis()->SetTitle("Strip multiplicity");
    h1d_mul->GetYaxis()->SetTitle("Probability %");
    h1d_mul->Sumw2();
    h1d_mul->Scale(100./h1d_mul->Integral());
    h1d_mul->SetBarWidth(0.5);
    h1d_mul->SetBarOffset(0.25);
    h1d_mul->SetFillColor(kRed);
    h1d_mul->Draw("HISTO B");

    theApp->Run();
    return;
}

int main(Int_t argc, Char_t* argv[])
{
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gROOT->Macro("rootlogon.C");
    gStyle->SetPalette(kRainBow);
    FOOT_ana(0,-1);
    return 0;
}
