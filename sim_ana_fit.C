
/*
  sim_ana_fit.C

  Simulation analogue of Marco's ana_fit.C. Reads the event tree produced
  by sim_event_builder.C and, per crystal (0-54) x trigger (0-3):
    1. Builds charge histograms (hQ_range for range-finding, hQ for fitting).
    2. Fits a Landau+exponential (noise tail) to seed initial parameters.
    3. Fits the full Landau-Gauss convolution + exponential (langauExpo) to
       extract the MPV (most probable charge deposit).
    4. Uses the MPV vs longitudinal position (4 trigger points) to compute
       light yield (MeV deposited / MeV expected from Ediff * path length)
       and attenuation profiles per crystal, grouped by crystal "type".

  Output: multi-page PDF (outname) + light-yield summary text file
  (out/LY_batch_<batch>.dat).

  langau()/langauExpo(): shared fit function definitions (Landau convolved
  with a Gaussian resolution term, optionally + an exponential noise tail).
*/

gErrorIgnoreLevel = kWarning;

double langau(double *x, double *par) {

    //Fit parameters:
       //par[0]=Width (scale) parameter of Landau density
       //par[1]=Most Probable (MP, location) parameter of Landau density
       //par[2]=Total area (integral -inf to inf, normalization constant)
       //par[3]=Width (sigma) of convoluted Gaussian function
      
    // Numeric constants
          Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
          Double_t mpshift  = -0.22278298;       // Landau maximum location

          // Control constants
          Double_t np = 1000.0;      // number of convolution steps
          Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

          // Variables
          Double_t xx;
          Double_t mpc;
          Double_t fland;
          Double_t sum = 0.0;
          Double_t xlow,xupp;
          Double_t step;
          Double_t i;

          // MP shift correction
          mpc = par[1] - mpshift * par[0];

          // Range of convolution integral
          xlow = x[0] - sc * par[3];
          xupp = x[0] + sc * par[3];

          step = (xupp-xlow) / np;

          // Convolution integral of Landau and Gaussian by sum
          for(i=1.0; i<=np/2; i++) {
             xx = xlow + (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);

             xx = xupp - (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);
          }

          return (par[2]*step*sum*invsq2pi / par[3]);
}

double langauExpo(double *x, double *par)
{
    // par[0] = Landau width
    // par[1] = MPV
    // par[2] = area
    // par[3] = Gaussian sigma
    // par[4] = exponential amplitude
    // par[5] = exponential decay constant

    double lang = langau(x, par);

    double expo = par[4] * TMath::Exp(-x[0] / par[5]);

    return lang + expo;
}

void sim_ana_fit(){

    gStyle->SetOptFit(1111); 
    int batch = 5;

    // ── Input/output files ───────────────────────────────────────────────
    string fname = "~/Desktop/BDX/data/sim_reco/out/evt/out_typ8_evt.root";
    string outname  = "./pdf/batches/out_typ8_fit.pdf";

    // ── Load crystal-serial lookup map (channel index -> serial number) ──
    // Batch map file: each line is a row of the stack, 5 crystal ids per
    // row, column 0 = leftmost, column 4 = rightmost. The FIRST line in the
    // file is the BOTTOM layer of the stack (iy=0), increasing downward in
    // the file as iy increases.
    // Channel numbering must match the same convention used when building
    // the simulation geometry (geometry_crs.pl):
    //   ch = iy*5 + (4-ti) + 1;   // col4=rightmost=offset0, col0=leftmost=offset4
    // i.e. within a row, channels count right-to-left, so we can't just use
    // a straight running counter over the file's (left-to-right) token order.
    ifstream crs_map_file(Form("data/batch_%i_map.dat", batch));
    map<int, int> crs_map; 

    int crs_id, iy = 0, ti = 0; 
    while(crs_map_file >> crs_id){
        int ch = iy*5 + (4-ti) + 1; 
        crs_map[ch] = crs_id; 

        ti++; 
        if(ti == 5){ ti = 0; iy++; }
    }
    for(int ii = 56; ii < 60; ii++){
        crs_map[ii] = 8000+ii; 
    }
    for(int ii = 60; ii < 64; ii++){
        crs_map[ii] = 9000+ii; 
    }

    TFile *f = new TFile(fname.c_str(), "read"); 
    TTree *t = (TTree*) f->Get("t"); 

    int ch, seed; 
    double q; 

    t->SetBranchAddress("ch", &ch); 
    t->SetBranchAddress("seed", &seed); 
    t->SetBranchAddress("q", &q); 
    
  // ── hQ_range: wide-range histogram per crystal/trigger, used only to ────
  // estimate a sensible [xmin,xmax] quantile window for the real hQ below.
  TH1D* hQ_range[55][4]; 
  for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        hQ_range[crs][pl] = new TH1D(Form("hQ_range_%i_%i", crs, pl), 
        Form("Charge crs %i (ch %i) - trg %i; Q (MeV); A.U.", crs_map[crs+1], crs, pl), 100,  0, 120); 
    }
  }

  for(int ii = 0; ii < t->GetEntries(); ii++){
    t->GetEntry(ii); 
    hQ_range[ch][seed]->Fill(q/1000.); 
  }

 
// ── hQ: final per-crystal/trigger charge histograms used for fitting ────
// Range is set from the 2%-99% quantiles of hQ_range (crude noise cut on the
// low side, while keeping most of the long Landau tail on the high side so
// the fit isn't squeezed/truncated against the upper edge).
TH1D* hQ[55][4]; 
  for(int crs = 0; crs < 55; crs++){
    double xmin = 200e3, xmax = 0; 

    for(int pl = 0; pl < 4; pl++){
      double probs[2] = {0.02, 0.99};
      double q[2];

      hQ_range[crs][pl]->GetQuantiles(2, q, probs);
      xmin = q[0];
      xmax = q[1];

      // Guard against degenerate/empty hQ_range (e.g. phantom channels with
      // no crystal installed, or very low statistics): GetQuantiles can
      // return xmin==xmax in that case, which creates a zero-width TH1D and
      // later crashes ROOT (TPad::Range "illegal world coordinates range").
      if(xmax <= xmin) xmax = xmin + 1;

      hQ[crs][pl] = new TH1D(Form("hQ_%i_%i", crs, pl), 
      Form("Charge crs %i (ch %i) - trg %i; Q (MeV); A.U.", crs_map[crs+1], crs, pl), 25 , xmin, xmax); 
      hQ[crs][pl]->SetLineColor(pl+1); 
    }
  }

  // ── hcount: 2D hit-count maps (col vs row) per trigger, for hit-map canvas
  TH2D* hcount[4]; 
  for(int ii = 0; ii < 4; ii++){
    hcount[ii] = new TH2D(Form("hcount_%i", ii), Form("Count trg %i", ii), 5, 0, 5, 12, 0, 12); 
  }

  // ── Fill hQ and hcount from the event tree ───────────────────────────
  for(int ii = 0; ii < t->GetEntries(); ii++){
    t->GetEntry(ii); 
    
    int x = 4 - ch % 5, y = (int) ch / 5; 
    hcount[seed]->Fill(x, y); 
    
    hQ[ch][seed]->Fill(q/1000.); 
  }

// Trigger/seed colors and labels: seed index pl (0-3) corresponds to the
// top+bottom plastic pair that fired: pl=0 -> scint 56(top)+60(bottom), etc.
int trig_colors[4]  = {kBlack, kRed, kGreen+2, kBlue};
const char* trig_labels[4] = {
    "Trig 0 (scint 56+60)",
    "Trig 1 (scint 57+61)",
    "Trig 2 (scint 58+62)",
    "Trig 3 (scint 59+63)"
};

// ══════════════════════════════════════════════════════════════════════
// CANVAS 1 (page 1): overview grid, 1 pad/crystal, all 4 triggers overlaid
// and normalized to their own peak height, for a quick by-eye shape check.
// ══════════════════════════════════════════════════════════════════════
TCanvas* c = new TCanvas(); 
c->Divide(7, 8); 
for(int crs = 0; crs < 55; crs ++){
    c->cd(crs+1); 
    hQ_range[crs][0]->SetLineColor(trig_colors[0]); 
    hQ_range[crs][0]->Scale(1./hQ_range[crs][0]->GetBinContent(hQ_range[crs][0]->GetMaximumBin()));
    hQ_range[crs][0]->Draw("histo"); 
    for(int pl = 1; pl < 4; pl++){
        hQ_range[crs][pl]->SetLineColor(trig_colors[pl]); 
        hQ_range[crs][pl]->Scale(1./hQ_range[crs][pl]->GetBinContent(hQ_range[crs][pl]->GetMaximumBin()));
        hQ_range[crs][pl]->Draw("samehisto"); 
    }
}

// 7x8 = 56 pads, only 55 crystals used -> pad 56 is empty, put legend there
c->cd(56);
TLegend* leg = new TLegend(0.05, 0.05, 0.95, 0.95);
leg->SetBorderSize(0);
leg->SetTextSize(0.09);
for(int pl = 0; pl < 4; pl++){
    leg->AddEntry(hQ_range[0][pl], trig_labels[pl], "l");
}
leg->Draw();

c->Print(Form("%s(", outname.c_str())); 

// ── Page(s) 2+: one page per crystal, all 4 triggers overlaid (un-normalized)
for(int crs = 0; crs < 55; crs ++){
TCanvas* cqi = new TCanvas(); 
    for(int pl = 0; pl < 4; pl++){
        cqi->cd(); 
        hQ_range[crs][pl]->Draw("samehisto");
    }
    cqi->Print(outname.c_str()); 
}

// ── Pre-allocate summary canvases (one per trigger, 7x8 grid of crystals) ──
TCanvas* cfitlandau[4]; 
for(int ii = 0; ii < 4; ii++){
    cfitlandau[ii] = new TCanvas(); 
    cfitlandau[ii]->Divide(7, 8); 
}
TCanvas* cfitlangau[4]; 
for(int ii = 0; ii < 4; ii++){
    cfitlangau[ii] = new TCanvas(); 
    cfitlangau[ii]->Divide(7, 8); 
}


TCanvas* cfit = new TCanvas(); 
cfit->Divide(2,1); 

TF1* flandau[55][4]; 
TF1* flangau[55][4]; 

// ══════════════════════════════════════════════════════════════════════
// MAIN FIT LOOP: for every crystal x trigger, fit Landau+expo first (to
// seed mu/sigma), then the full Landau-Gauss-convolution+expo fit (langauExpo)
// ══════════════════════════════════════════════════════════════════════
for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        cfitlandau[pl]->cd(1+crs);

TH1D *hist = hQ[crs][pl];
// ── Find Landau peak via smoothed (3-bin) histogram maximum ────────────
// Note: hQ is already built from the 2%-99% quantile range of hQ_range, so
// the low-energy noise bump is cut at histogram-creation time. skipBins is
// no longer needed to avoid noise and is kept at 0 so the full fit range
// is searched (a nonzero value here would eat into the real Landau rise).
int skipBins = 0;
int peakBin = skipBins;
double maxSmooth = -1;
for(int ib = skipBins; ib < hist->GetNbinsX()-1; ib++){
    double smooth =
        hist->GetBinContent(ib-1) +
        hist->GetBinContent(ib)   +
        hist->GetBinContent(ib+1);
    if(smooth > maxSmooth){
        maxSmooth = smooth;
        peakBin = ib;
    }
}

double peak = hist->GetBinCenter(peakBin);

double rms = hQ[crs][pl]->GetRMS();

// Range più largo
double fitLow  = std::max(0.0, peak - 1.5*rms);
double fitHigh = peak + 2.5*rms;


// ============================================================
// LANDAU + EXPONENTIAL
// ============================================================

flandau[crs][pl] = new TF1(
    Form("flandau_%i_%i", crs, pl),
    "landau(0) + expo(3)",
    fitLow,
    fitHigh
);

// Landau parameters
flandau[crs][pl]->SetParameter(0, 0.5*hQ[crs][pl]->GetEntries());
flandau[crs][pl]->SetParameter(1, peak);
flandau[crs][pl]->SetParameter(2, 0.3*rms);

flandau[crs][pl]->SetParLimits(0, 0, 10*hQ[crs][pl]->GetEntries());
flandau[crs][pl]->SetParLimits(1, peak-rms, peak+rms);
flandau[crs][pl]->SetParLimits(2, 0.02*rms, rms);

// Exponential initial parameters
// expo(3) = exp(p3 + p4*x)
flandau[crs][pl]->SetParameter(3, TMath::Log(std::max(1.0, hQ[crs][pl]->GetBinContent(hQ[crs][pl]->FindBin(fitLow)))));

flandau[crs][pl]->SetParameter(4, -1.0/rms);
flandau[crs][pl]->SetParLimits(4, -0.1, -1e-6);

hQ[crs][pl]->Fit(flandau[crs][pl], "QBL");

hQ[crs][pl]->Draw();
flandau[crs][pl]->Draw("same");

// pad 1 of cfit: show the Landau+expo (pre-fit) result
cfit->cd(1);

hQ[crs][pl]->Draw();
flandau[crs][pl]->Draw("same");

cfitlangau[pl]->cd(1+crs);

// Seed the full langau fit from the Landau+expo fit results
double Afit = flandau[crs][pl]->GetParameter(0);

double mu =  flandau[crs][pl]->GetParameter(1);

double sigma = flandau[crs][pl]->GetParameter(2);


// ============================================================
// LANDAU-GAUSSIAN + EXPONENTIAL
// ============================================================

// Range ancora più largo rispetto a prima
double langauLow  = std::max(0.0, mu - 3.0*sigma);
double langauHigh = mu + 6.0*sigma;     // mu and sigma from Landauexpo fit

flangau[crs][pl] = new TF1(Form("flangau_%i_%i", crs, pl),
    langauExpo,
    langauLow,
    langauHigh,
    6
);

flangau[crs][pl]->SetNpx(1000);


// Langau initial parameters
flangau[crs][pl]->SetParameter(0, sigma);
flangau[crs][pl]->SetParameter(1, mu);
flangau[crs][pl]->SetParLimits(1, mu - 10.0, mu + 10.0);
flangau[crs][pl]->SetParameter(2, Afit);
flangau[crs][pl]->SetParameter(3, 0.1);


// Exponential initial parameters
double expAmplitude =
    hQ[crs][pl]->GetBinContent(
        hQ[crs][pl]->FindBin(langauLow)
    );

flangau[crs][pl]->SetParameter(
    4,
    std::max(1.0, expAmplitude)
);

flangau[crs][pl]->SetParameter(
    5,
    -1.0/flandau[crs][pl]->GetParameter(4)
);


// Parameter limits
flangau[crs][pl]->SetParLimits(0, 2.0, 8.0);

flangau[crs][pl]->SetParLimits(3, 0.0, 1.0);

// Exponential amplitude
flangau[crs][pl]->SetParLimits(
    4,
    0,
    1e6
);

// Exponential decay scale
flangau[crs][pl]->SetParLimits(
    5,
    0.02,
    100.0
);


hQ[crs][pl]->Fit(flangau[crs][pl], "QBL");

hQ[crs][pl]->Draw();
flangau[crs][pl]->Draw("same");

// pad 2 of cfit: show the final langauExpo result, side by side with pad 1
cfit->cd(2);

hQ[crs][pl]->Draw();
flangau[crs][pl]->Draw("same");

cfit->Print(outname.c_str());

    }
  }

  // print the per-trigger 7x8 summary canvases (Landau+expo canvases disabled)
  for(int ii = 0; ii < 4; ii++){
    //cfitlandau[ii]->Print(outname.c_str()); 
    cfitlangau[ii]->Print(outname.c_str()); 
  }




 
  // ══════════════════════════════════════════════════════════════════
  // LIGHT YIELD / ATTENUATION SECTION: for each crystal, plot fitted MPV
  // (gr_att) and MPV/(Ediff*dE/dx) light yield (gr_ly) vs longitudinal
  // position x, using the 4 trigger positions as the 4 x-points.
  // ══════════════════════════════════════════════════════════════════
  TCanvas* cg = new TCanvas(); 
  cg->Divide(2,1); 

  TH1D* hly_distrib = new TH1D("hly_distrib", "Distribution of light yield; MeV; A.U.", 100, 0, 0);

  TGraphErrors* gr_att[55]; 
  TGraphErrors* gr_ly[55]; 
  TGraphErrors* gr_ly_norm[55]; 

  // H/h: crystal geometry (long/short face heights) per crystal "type" 1-8,
  // used to interpolate the muon path length y(x) through the trapezoidal
  // crystal cross-section at longitudinal position x.
  double H[8] = {4.96, 5.31, 5.27, 6.65, 6.28, 6.08, 5.80, 5.75};
  double h[8] = {2.44, 2.38, 2.13, 2.45, 2.08, 1.88, 1.74, 1.69};
  double Ediff = 2*7;  // energy loss per cm for a MIP muon (supervisor-provided, do not change)

  TCanvas* catt = new TCanvas(); 

  vector<double> ly_map(55); 
  double ly_map_pl[4][55] = {0}; // per-trigger (pl=0..3) light yield, for the 4 separate 2D maps

  ofstream ofile_ly_err(Form("out/LY/batch_%i_simu.dat", batch));

  catt->Divide(10, 6); 
  for(int crs = 0; crs < 55; crs++){

    // Skip phantom channels with no crystal installed (crs_map defaults to 0
    // for unset keys). Without this, channels still receiving simulated hits
    // get plotted/titled as "Crs 0" instead of being excluded entirely.
    if(crs_map[crs+1] == 0) continue;

    ofile_ly_err << crs_map[crs+1] << endl;

    gr_att[crs] = new TGraphErrors(); 
    gr_att[crs]->SetTitle(Form("Crs %i; x (cm); MeV", crs_map[crs+1])); 

    gr_ly[crs] = new TGraphErrors();
    gr_ly[crs]->SetTitle(Form("Crs %i; x (cm); MeV/MeV", crs_map[crs+1]));

    gr_ly_norm[crs] = new TGraphErrors();
    gr_ly_norm[crs]->SetTitle(Form("Crs %i; x (cm); MeV/MeV (normalized to average)", crs_map[crs+1]));

    double ly_average = 0; 
    for(int pl = 0; pl < 4; pl++){

        // x: longitudinal position along crystal (cm) for this trigger index,
        // direction depends on crystal orientation in the 5x11 grid.
        double x = (crs < 5 || 
                    (crs > 9 && crs < 15) ||
                     (crs > 19 && crs < 25) || 
                     (crs > 29 && crs < 35) || 
                     (crs > 39 && crs < 45) || 
                     (crs > 49 && crs < 55) ) ? 20-pl*6 : 2+pl*6;

        ofile_ly_err << x << " " << flandau[crs][pl]->GetParameter(1) << " " << flandau[crs][pl]->GetParError(1) << endl;
      
        if(hQ[crs][pl]->GetEntries() < 20) continue;

        int type = (int) crs_map[crs+1]/100;
        if(type > 8) type = 16-type; 

        double  y = H[type-1] + (h[type-1]-H[type-1])/24*x;   // path length (cm) at x

        gr_att[crs]->SetPoint(pl, x, flandau[crs][pl]->GetParameter(1)); 
        gr_att[crs]->SetPointError(pl, 1, flandau[crs][pl]->GetParError(1));

        // Light yield = fitted MPV / (energy loss per cm * path length)
        gr_ly[crs]->SetPoint(pl, x, flandau[crs][pl]->GetParameter(1)/(Ediff*y));
        gr_ly[crs]->SetPointError(pl, 1, flandau[crs][pl]->GetParError(1)/(Ediff*y));

        gr_ly_norm[crs]->SetPoint(pl, x, flandau[crs][pl]->GetParameter(1)/(Ediff*y));
        gr_ly_norm[crs]->SetPointError(pl, 1, flandau[crs][pl]->GetParError(1)/(Ediff*y));

        ly_average += flandau[crs][pl]->GetParameter(1)/(Ediff*y);

        // per-trigger LY value (MeV/MeV), kept separate from ly_average/ly_map
        ly_map_pl[pl][crs] = flandau[crs][pl]->GetParameter(1)/(Ediff*y);
    }

    hly_distrib->Fill(ly_average/4.);

    ly_map[crs] = ly_average/4.;

    // Some channels have no crystal installed (new setup), so their graphs
    // stay empty (GetN()==0). Drawing an empty TGraphErrors with "AP" makes
    // ROOT compute an illegal axis range (xmin==xmax) and can crash/segfault.
    // Skip drawing (and the per-crystal page in cg) for those empty channels.
    if(gr_att[crs]->GetN() == 0) continue;

    catt->cd(crs+1); 
    gr_att[crs]->SetMarkerColor(20); 
    gr_att[crs]->Draw("AP"); 

    cg->cd(1);
    gr_att[crs]->Draw("ap"); 
    cg->cd(2);
    gr_ly[crs]->Draw("ap");
    cg->Print(outname.c_str()); 
  }


  TCanvas* cly_distrib = new TCanvas(); 
  hly_distrib->Draw("histo");
  cly_distrib->Print(outname.c_str());

  // ══════════════════════════════════════════════════════════════════════
  // LIGHT-YIELD-BY-POSITION: 2D map of light yield (col vs row)
  // One map per trigger position (pl=0..3), with the crystal serial number
  // overlaid. col = 4 - crs%5 matches the hcount convention used in the
  // working hit-count maps below, so orientation is consistent between the
  // two sets of plots. Lets us check whether light yield correlates with
  // position in the stack, at each individual trigger point.
  // ══════════════════════════════════════════════════════════════════════

  TH2D* hly_map2d[4];
  TCanvas* cly_map2d = new TCanvas();
  cly_map2d->Divide(2,2);
  TLatex* ltx_ly = new TLatex();
  ltx_ly->SetTextAlign(22);
  ltx_ly->SetTextSize(0.03);
  ltx_ly->SetTextColor(kWhite);

  for(int pl = 0; pl < 4; pl++){
    hly_map2d[pl] = new TH2D(Form("hly_map2d_%i", pl), Form("Light yield map - trg %i; col; row", pl), 5, 0, 5, 12, 0, 12);
    for(int crs = 0; crs < 55; crs++){
      if(crs_map[crs+1] == 0) continue; // no crystal installed at this position
      int col = 4 - crs % 5, row = crs / 5;
      hly_map2d[pl]->SetBinContent(col+1, row+1, ly_map_pl[pl][crs]);
    }
    cly_map2d->cd(pl+1);
    hly_map2d[pl]->SetStats(0);
    hly_map2d[pl]->Draw("COLZ");
    for(int crs = 0; crs < 55; crs++){
      if(crs_map[crs+1] == 0) continue; // no crystal installed at this position
      int col = 4 - crs % 5, row = crs / 5;
      ltx_ly->DrawLatex(col+0.5, row+0.5, Form("%i", crs_map[crs+1]));
    }
  }
  cly_map2d->Print(outname.c_str());


  gStyle->SetOptStat(0); 

// ── Per-trigger hit-count maps (col vs row, 2x2 grid) ──────────────────
TCanvas* ccount = new TCanvas(); 
ccount->Divide(2,2); 
for(int ii = 0; ii < 4; ii++){
  ccount->cd(ii+1); 
  //ccount->cd(ii+1)->SetLogz(); 
  hcount[ii]->Draw("histo"); 
}
ccount->Print(outname.c_str()); 

// ══════════════════════════════════════════════════════════════════════
// LIGHT-YIELD-BY-TYPE SECTION: group crystals by "type" (1-8, encoded in
// the hundreds digit of the serial number) and compare normalized light
// yield profiles + distribution per type.
// ══════════════════════════════════════════════════════════════════════
TCanvas* cly_type[8]; 
for(int ii = 0; ii < 8; ii++){
  cly_type[ii] = new TCanvas(); 
  cly_type[ii]->Divide(2,1); 
}

TH1D* hly_type[8];
for(int ii = 0; ii < 8; ii++){
  hly_type[ii] = new TH1D(Form("hly_type_%i", ii), Form("Distribution of light yield type %i; MeV/MeV; A.U.", ii+1), 100, 0, 0);
}

for(int ii = 0; ii < 55; ii++){
  // Skip empty channels (no crystal installed): crs_map defaults to 0 for
  // unset keys, and ly_map stays 0 since the fit loop never filled it,
  // which would otherwise give an invalid type index (hly_type[-1]) or a
  // divide-by-zero in the normalization below.
  if(crs_map[ii+1] == 0 || ly_map[ii] == 0) continue;

  int type = (int) crs_map[ii+1]/100;
  if(type > 8) type = 16-type; 

  hly_type[type-1]->Fill(ly_map[ii]);

  cly_type[type-1]->cd(1);
  gr_ly_norm[ii]->SetMarkerStyle(20);
  gr_ly_norm[ii]->SetMarkerColor(hly_type[type-1]->GetEntries());
  // normalize each crystal's 4-point profile to its own average light yield
  for(int jj = 0; jj < 4; jj++){
    gr_ly_norm[ii]->SetPoint(jj, jj, gr_ly_norm[ii]->GetPointY(jj)/ly_map[ii]);
    gr_ly_norm[ii]->SetPointError(jj, 0, gr_ly_norm[ii]->GetErrorY(jj)/ly_map[ii]);
  }

  if(hly_type[type-1]->GetEntries() == 1) gr_ly_norm[ii]->Draw("AP");
  else gr_ly_norm[ii]->Draw("P same");
}

for(int ii = 0; ii < 8; ii++){
  cly_type[ii]->cd(2);
  hly_type[ii]->Draw("histo");
  cly_type[ii]->Print(outname.c_str());
}

// ── Write per-crystal light yield (MPV-derived) to a text summary file ──
ofstream ofile_ly(Form("out/LY_batch_%i.dat", batch));

for(int ii = 0; ii < 55; ii++){
    ofile_ly << crs_map[ii+1] << "\t" ; 
  gr_ly[ii]->Sort();

    for(int jj = 0; jj < 4; jj ++){
      ofile_ly << gr_ly[ii]->GetPointY(jj) << "\t";
    } 
    ofile_ly << ly_map[ii] << endl; 
}

  // Close out the multi-page PDF (catt was opened across the loop above)
  catt->Print((outname + ")").c_str()); 

}