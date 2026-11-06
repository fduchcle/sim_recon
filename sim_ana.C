/*
  sim_ana_fit.C

  Simple analysis script for simulation event-builder output
  (from sim_event_builder.C).

  Before the main 55-crystal charge-distribution scan, this script first
  looks for HIGH-ANGLE coincidence events (theta > 45 deg from vertical) --
  these were flagged by inspection as outliers in the theta distribution --
  and produces one diagnostic canvas per such event: a single-pad 2D hit
  map of the BGO crystal array (5 col x 11 row) coloured by energy
  deposited in that single event, with theta shown in the pad title.

  It then does the same for a sample of "typical" events with theta ~ 10
  deg (first 10 found), for comparison.

  It also prints to the terminal, for each selected event, the plastic
  scintillator signal (sector + energy in phe) for that event.

  Usage:
    root -l -q 'sim_ana.C("../sim_data/out_final_10.22.root")'

  Output:
    ./pdf/<basename>_highTheta.pdf   (one page per theta>45 deg event)
    ./pdf/<basename>_lowTheta.pdf    (one page per theta~10 deg event, first 10)
    ./pdf/<basename>_ana.pdf         (55 pages, one per crystal)
*/

#include <string>
#include <vector>
using namespace std;

static const int N_COL = 5;
static const int N_ROW = 11;

void sim_ana(const char* filename = "../sim_data/out_typ5_typ2.root") {

    gErrorIgnoreLevel = kWarning;  // suppress "Current canvas added to pdf file" info messages
    gStyle->SetOptStat(1111);

    string fname = filename;
    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');
    const string basename = fname.substr(start, end - start);
    const string inname  = "./out/" + basename + "_evt.root";
    const string outname = "./pdf/batches/" + basename + "_ana.pdf";
    const string hthetaname = "./pdf/batches/" + basename + "_highTheta.pdf";

    gSystem->mkdir("./pdf/batches", true);

    TFile *f = new TFile(inname.c_str(), "read");
    if (!f || f->IsZombie()) {
        cout << "Error: cannot open " << inname
             << " -- did you run sim_event_builder.C first?" << endl;
        return;
    }

    // ══════════════════════════════════════════════════════════════════
    // Section 1: high-theta (> 45 deg) event diagnostics
    // ══════════════════════════════════════════════════════════════════
    TTree *tevt = (TTree*) f->Get("tevt");
    if (!tevt) {
        cout << "Warning: tree 'tevt' not found in " << inname
             << " -- skipping high-theta event diagnostics"
             << " (re-run the updated sim_event_builder.C to produce it)." << endl;
    } else {
        int evtId, seed;
        double theta;
        double tpx, tpy, tpz, tvx, tvy, tvz;
        vector<int>    *ch_v = nullptr, *scint_sector_v = nullptr;
        vector<double> *q_v = nullptr,  *scint_A_v = nullptr;

        tevt->SetBranchAddress("evt",          &evtId);
        tevt->SetBranchAddress("seed",         &seed);
        tevt->SetBranchAddress("theta",        &theta);
        tevt->SetBranchAddress("ch",           &ch_v);
        tevt->SetBranchAddress("q",            &q_v);
        tevt->SetBranchAddress("scint_sector", &scint_sector_v);
        tevt->SetBranchAddress("scint_A",      &scint_A_v);
        tevt->SetBranchAddress("px", &tpx);
        tevt->SetBranchAddress("py", &tpy);
        tevt->SetBranchAddress("pz", &tpz);
        tevt->SetBranchAddress("vx", &tvx);
        tevt->SetBranchAddress("vy", &tvy);
        tevt->SetBranchAddress("vz", &tvz);

        const double THETA_CUT = 40.0;  // deg

        Long64_t Nevt = tevt->GetEntries();
        vector<Long64_t> highThetaIdx;
        for (Long64_t ii = 0; ii < Nevt; ii++) {
            tevt->GetEntry(ii);
            if (theta > THETA_CUT) highThetaIdx.push_back(ii);
        }

        cout << "==============================================" << endl;
        cout << "High-theta (> " << THETA_CUT << " deg) coincidence events: "
             << highThetaIdx.size() << " / " << Nevt << endl;
        cout << "==============================================" << endl;

        if (!highThetaIdx.empty()) {
            TCanvas *ch2 = new TCanvas("cHighTheta", "High-theta event diagnostics", 700, 700);

            gStyle->SetPalette(kBird);
            ch2->Print((hthetaname + "[").c_str());

            for (Long64_t idx : highThetaIdx) {
                tevt->GetEntry(idx);

                cout << "----------------------------------------------" << endl;
                cout << "Event " << evtId << "  (theta = " << theta << " deg, seed = " << seed << ")" << endl;
                cout << "  Plastic scintillator signal:" << endl;
                for (size_t jj = 0; jj < scint_sector_v->size(); jj++) {
                    cout << "    sector=" << (*scint_sector_v)[jj]
                         << "  E = " << (*scint_A_v)[jj] << " phe" << endl;
                }
                cout << "  Primary truth: p=(" << tpx << ", " << tpy << ", " << tpz
                     << ") keV,  vertex=(" << tvx << ", " << tvy << ", " << tvz << ") cm" << endl;

                // ── Converted values ready for a GEMC event-display gcard ──
                double p_keV = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);
                double p_GeV = p_keV / 1.0e6;
                double thetaDir = acos(tpz / p_keV) * 180.0 / TMath::Pi();
                double phiDir   = atan2(tpy, tpx) * 180.0 / TMath::Pi();
                cout << "  GEMC BEAM_P: mu-," << p_GeV << "*GeV," << thetaDir << "*deg," << phiDir << "*deg" << endl;
                cout << "  GEMC BEAM_V: (" << tvx << "," << tvy << "," << tvz << ")cm" << endl;

                // ── Single pad: 2D BGO hit map (energy-weighted) for this event ──
                ch2->cd();
                gPad->Clear();
                gPad->SetRightMargin(0.15);
                gPad->SetLeftMargin(0.10);
                TH2D *hmap = new TH2D(Form("hmap_%lld", idx),
                    Form("Event %i (#theta=%.1f deg); z_{ch} (col); y_{ch} (row)", evtId, theta),
                    N_COL, -0.5, N_COL - 0.5,
                    N_ROW, -0.5, N_ROW - 0.5);
                for (size_t jj = 0; jj < ch_v->size(); jj++) {
                    int ch0 = (*ch_v)[jj];        // already 0-54
                    int ch  = ch0 + 1;            // back to 1-55 to invert row/col
                    int row = (ch - 1) / 5;
                    int col = 4 - (ch - 1) % 5;
                    hmap->Fill(col, row, (*q_v)[jj]);
                }
                hmap->SetStats(0);
                hmap->Draw("COLZ");

                for (int col = 0; col <= N_COL; col++) {
                    TLine *l = new TLine(col - 0.5, -0.5, col - 0.5, N_ROW - 0.5);
                    l->SetLineColor(kGray+1); l->SetLineWidth(1); l->Draw();
                }
                for (int row = 0; row <= N_ROW; row++) {
                    TLine *l = new TLine(-0.5, row - 0.5, N_COL - 0.5, row - 0.5);
                    l->SetLineColor(kGray+1); l->SetLineWidth(1); l->Draw();
                }

                ch2->Print(hthetaname.c_str());
            }

            ch2->Print((hthetaname + "]").c_str());
        }

        // ══════════════════════════════════════════════════════════════
        // Section 1b: sample of "typical" events with theta ~ 0 deg
        // ══════════════════════════════════════════════════════════════
        const double LOW_THETA  = 0.0;  // deg, target angle
        const double LOW_WINDOW = 2.0;   // deg, +/- window around target
        const int    LOW_MAXEVT = 100;    // first N matching events

        vector<Long64_t> lowThetaIdx;
        for (Long64_t ii = 0; ii < Nevt && (int)lowThetaIdx.size() < LOW_MAXEVT; ii++) {
            tevt->GetEntry(ii);
            if (fabs(theta - LOW_THETA) <= LOW_WINDOW) lowThetaIdx.push_back(ii);
        }

        if (!lowThetaIdx.empty()) {
            const string lowthetaname = "./pdf/" + basename + "_lowTheta.pdf";

            TCanvas *ch3 = new TCanvas("cLowTheta", "theta~10 deg event diagnostics", 700, 700);

            gStyle->SetPalette(kBird);
            ch3->Print((lowthetaname + "[").c_str());

            for (Long64_t idx : lowThetaIdx) {
                tevt->GetEntry(idx);

                cout << "----------------------------------------------" << endl;
                cout << "Event " << evtId << "  (theta = " << theta << " deg, seed = " << seed << ")" << endl;
                cout << "  Plastic scintillator signal:" << endl;
                for (size_t jj = 0; jj < scint_sector_v->size(); jj++) {
                    cout << "    sector=" << (*scint_sector_v)[jj]
                         << "  E = " << (*scint_A_v)[jj] << " phe" << endl;
                }
                cout << "  Primary truth: p=(" << tpx << ", " << tpy << ", " << tpz
                     << ") keV,  vertex=(" << tvx << ", " << tvy << ", " << tvz << ") cm" << endl;

                // ── Converted values ready for a GEMC event-display gcard ──
                double p_keV_low = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);
                double p_GeV_low = p_keV_low / 1.0e6;
                double thetaDir_low = acos(tpz / p_keV_low) * 180.0 / TMath::Pi();
                double phiDir_low   = atan2(tpy, tpx) * 180.0 / TMath::Pi();
                cout << "  GEMC BEAM_P: mu-," << p_GeV_low << "*GeV," << thetaDir_low << "*deg," << phiDir_low << "*deg" << endl;
                cout << "  GEMC BEAM_V: (" << tvx << "," << tvy << "," << tvz << ")cm" << endl;

                // ── Single pad: 2D BGO hit map (energy-weighted) for this event ──
                ch3->cd();
                gPad->Clear();
                gPad->SetRightMargin(0.15);
                gPad->SetLeftMargin(0.10);
                TH2D *hmap = new TH2D(Form("hmap_low_%lld", idx),
                    Form("Event %i (#theta=%.1f deg); z_{ch} (col); y_{ch} (row)", evtId, theta),
                    N_COL, -0.5, N_COL - 0.5,
                    N_ROW, -0.5, N_ROW - 0.5);
                for (size_t jj = 0; jj < ch_v->size(); jj++) {
                    int ch0 = (*ch_v)[jj];
                    int ch  = ch0 + 1;
                    int row = (ch - 1) / 5;
                    int col = 4 - (ch - 1) % 5;
                    hmap->Fill(col, row, (*q_v)[jj]);
                }
                hmap->SetStats(0);
                hmap->Draw("COLZ");

                for (int col = 0; col <= N_COL; col++) {
                    TLine *l = new TLine(col - 0.5, -0.5, col - 0.5, N_ROW - 0.5);
                    l->SetLineColor(kGray+1); l->SetLineWidth(1); l->Draw();
                }
                for (int row = 0; row <= N_ROW; row++) {
                    TLine *l = new TLine(-0.5, row - 0.5, N_COL - 0.5, row - 0.5);
                    l->SetLineColor(kGray+1); l->SetLineWidth(1); l->Draw();
                }

                ch3->Print(lowthetaname.c_str());
            }

            ch3->Print((lowthetaname + "]").c_str());
        }

        // ══════════════════════════════════════════════════════════════
        // Section 1c: angular distribution vs truth energy, all events
        // ══════════════════════════════════════════════════════════════
        const string angEname = "./pdf/" + basename + "_angVsE.pdf";

        // Truth kinetic/total energy proxy: |p| (keV), as already used above
        // for the GEMC conversion (tpx,tpy,tpz are truth momentum components).
        // ── First pass: find the truth-energy range to set sensible bins ──
        double EtruthMin = 1e300, EtruthMax = -1e300;
        for (Long64_t ii = 0; ii < Nevt; ii++) {
            tevt->GetEntry(ii);
            double Etruth = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);  // keV
            if (Etruth < EtruthMin) EtruthMin = Etruth;
            if (Etruth > EtruthMax) EtruthMax = Etruth;
        }
        if (EtruthMax <= EtruthMin) EtruthMax = EtruthMin + 1.0;
        const double Epad = 0.05 * (EtruthMax - EtruthMin);

        TH2D *hThetaE = new TH2D("hThetaE",
            "Angle vs truth energy (all events); E_{truth} (keV); #theta (deg)",
            100, EtruthMin - Epad, EtruthMax + Epad, 90, 0, 90);

        TH1D *hEtruth = new TH1D("hEtruth",
            "Truth energy, all events; E_{truth} (keV); Entries",
            100, EtruthMin - Epad, EtruthMax + Epad);

        for (Long64_t ii = 0; ii < Nevt; ii++) {
            tevt->GetEntry(ii);
            double Etruth = sqrt(tpx*tpx + tpy*tpy + tpz*tpz);  // keV
            hThetaE->Fill(Etruth, theta);
            hEtruth->Fill(Etruth);
        }

        TCanvas *c4 = new TCanvas("cAngVsE", "Angular distribution vs truth energy", 1200, 600);
        c4->Divide(2, 1);

        c4->cd(1);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        hThetaE->SetStats(0);
        hThetaE->Draw("COLZ");

        c4->cd(2);
        gPad->SetLogy();
        hEtruth->SetLineColor(kBlue+1);
        hEtruth->Draw("histo");

        c4->Print(angEname.c_str());

        cout << "==============================================" << endl;
        cout << "Angular-vs-energy plot written to " << angEname << endl;
        cout << "==============================================" << endl;
    }

    // ══════════════════════════════════════════════════════════════════
    // Section 2: per-crystal, per-seed charge distributions (all events)
    // ══════════════════════════════════════════════════════════════════
    TTree *t = (TTree*) f->Get("t");
    if (!t) { cout << "Error: tree 't' not found in " << inname << endl; return; }

    int ch, seed;
    double q, theta;
    t->SetBranchAddress("ch",    &ch);
    t->SetBranchAddress("seed",  &seed);
    t->SetBranchAddress("q",     &q);
    t->SetBranchAddress("theta", &theta);

    // ── Book one histogram per crystal (0-54) x seed (0-3) ──────────────
    TH1D* hQ[55][4];
    for (int crs = 0; crs < 55; crs++) {
        for (int pl = 0; pl < 4; pl++) {
            hQ[crs][pl] = new TH1D(
                Form("hQ_%i_%i", crs, pl),
                Form("Channel %i - seed %i; Q (keV); Entries", crs, pl),
                100, 0, 2e5);
        }
    }

    Long64_t N = t->GetEntries();
    for (Long64_t ii = 0; ii < N; ii++) {
        t->GetEntry(ii);
        if (ch < 0 || ch >= 55 || seed < 0 || seed >= 4) continue;
        hQ[ch][seed]->Fill(q);
    }

    // ── One 1x4 canvas per crystal, saved as pages of a single PDF ──────
    TCanvas *c = new TCanvas("c", "c", 1600, 500);
    c->Divide(4, 1);

    c->Print((outname + "[").c_str());  // open multi-page PDF

    for (int ch = 0; ch < 55; ch++) {
        for (int pl = 0; pl < 4; pl++) {
            c->cd(pl + 1);
            hQ[ch][pl]->SetLineColor(pl + 1);
            hQ[ch][pl]->Draw("histo");
        }
        c->Print(outname.c_str());
    }

    c->Print((outname + "]").c_str());  // close multi-page PDF

    f->Close();
}
