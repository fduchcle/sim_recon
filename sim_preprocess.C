/*
  sim_preprocess.C

  Simulation analogue of Marco's preprocess.C (BGO_testing-main).

  Reads a raw simulation ROOT file (trees: scint, crs, generated) where each
  tree entry already corresponds to one simulated event, and produces a
  simplified "preprocessed" file with ONE ENTRY PER EVENT containing:

    - the summed BGO energy deposition per crystal channel (1-55) that fired
      above threshold (defensive sum: a crystal is expected to appear at most
      once per event in the raw tree, but we sum just in case),
    - the plastic scintillator hits (sector 56-63) above threshold,
    - the primary particle truth momentum (px,py,pz) from the 'generated'
      tree, needed later to cut on the muon trajectory angle.

  This mirrors preprocess.C's role of turning raw per-hit data into a clean,
  per-event structure before event building. Unlike the experimental data,
  the simulation has no waveforms/timestamps, so there is no baseline
  subtraction or time-slicing step here.

  Usage:
    root -l -q 'sim_preprocess.C("../sim_data/out_60crs_8.7.root")'

  Output:
    ./out/<basename>_preprocessed.root   (tree "tout")
*/

#include <vector>
#include <map>
#include <string>
using namespace std;

void sim_preprocess(const char* filename = "../sim_data/out_typ6_typ3.root") {

    const double crs_th   = 5000.0;  // keV — BGO energy threshold
    const double scint_th = 10000;   // photoelectrons — plastic scintillator threshold, equivalent to 10 phe

    string fname = filename;
    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');
    const string basename = fname.substr(start, end - start);
    const string outname  = "./out/" + basename + "_preprocessed.root";

    gSystem->Unlink(outname.c_str());  // remove any stale file to avoid ROOT cycle buildup

    TFile *f = new TFile(fname.c_str(), "read");
    if (!f || f->IsZombie()) { cout << "Error: cannot open " << fname << endl; return; }

    TTree *tscint = (TTree*) f->Get("scint");
    TTree *tcrs   = (TTree*) f->Get("crs");
    TTree *tgen   = (TTree*) f->Get("generated");
    if (!tscint) { cout << "Error: tree 'scint' not found" << endl; return; }
    if (!tcrs)   { cout << "Error: tree 'crs' not found" << endl; return; }
    if (!tgen)   { cout << "Error: tree 'generated' not found" << endl; return; }

    Long64_t N = tscint->GetEntries();
    if (tcrs->GetEntries() != N || tgen->GetEntries() != N) {
        cout << "Error: scint/crs/generated trees have different entry counts" << endl;
        return;
    }

    vector<double> *sA = nullptr, *sSector = nullptr;
    tscint->SetBranchAddress("ADC_ADC", &sA);
    tscint->SetBranchAddress("sector",  &sSector);

    vector<double> *cA = nullptr, *cSector = nullptr;
    tcrs->SetBranchAddress("ADC_ADC", &cA);
    tcrs->SetBranchAddress("sector",  &cSector);

    vector<int>    *gpid = nullptr;
    vector<double> *gpx = nullptr, *gpy = nullptr, *gpz = nullptr;
    vector<double> *gvx = nullptr, *gvy = nullptr, *gvz = nullptr;
    tgen->SetBranchAddress("pid", &gpid);
    tgen->SetBranchAddress("px",  &gpx);
    tgen->SetBranchAddress("py",  &gpy);
    tgen->SetBranchAddress("pz",  &gpz);
    tgen->SetBranchAddress("vx",  &gvx);
    tgen->SetBranchAddress("vy",  &gvy);
    tgen->SetBranchAddress("vz",  &gvz);

    TFile *fout = new TFile(outname.c_str(), "recreate");
    TTree *tout = new TTree("tout", "Preprocessed per-event BGO + scint hits (sim)");

    vector<int>    ch_v, scint_sector_v;
    vector<double> q_v, scint_A_v;
    double px, py, pz;
    double vx, vy, vz;  // primary generation vertex (cm, GEMC "generated" tree convention)

    tout->Branch("ch",           &ch_v);           // BGO channel, 1-55
    tout->Branch("q",            &q_v);            // summed energy per channel (keV)
    tout->Branch("scint_sector", &scint_sector_v); // plastic sector, 56-63
    tout->Branch("scint_A",      &scint_A_v);      // plastic amplitude (keV, same ADC_ADC bank as BGO)
    tout->Branch("px", &px);
    tout->Branch("py", &py);
    tout->Branch("pz", &pz);
    tout->Branch("vx", &vx);
    tout->Branch("vy", &vy);
    tout->Branch("vz", &vz);

    long n_crs_hits_raw = 0, n_crs_hits_summed = 0;

    for (Long64_t ii = 0; ii < N; ii++) {
        tscint->GetEntry(ii);
        tcrs->GetEntry(ii);
        tgen->GetEntry(ii);

        ch_v.clear(); q_v.clear();
        scint_sector_v.clear(); scint_A_v.clear();

        // ── Sum BGO energy per crystal channel (defensive: normally 1 hit/ch/event) ──
        map<int, double> qsum;
        for (size_t jj = 0; jj < cA->size(); jj++) {
            double a = (*cA)[jj];
            n_crs_hits_raw++;
            if (a < crs_th) continue;
            int sec = (int)(*cSector)[jj];
            if (sec < 1 || sec > 55) continue;   // sector IS the 1-55 channel id
            qsum[sec] += a;
        }
        for (auto &p : qsum) {
            ch_v.push_back(p.first);
            q_v.push_back(p.second);
            n_crs_hits_summed++;
        }

        // ── Plastic scintillator hits above threshold ──────────────────────
        for (size_t jj = 0; jj < sA->size(); jj++) {
            double a = (*sA)[jj];
            if (a < scint_th) continue;
            int sec = (int)(*sSector)[jj];
            if (sec < 56 || sec > 63) continue;
            scint_sector_v.push_back(sec);
            scint_A_v.push_back(a);
        }

        // ── Primary particle truth momentum + vertex (single primary expected) ─────
        px = py = pz = 0;
        vx = vy = vz = 0;
        if (gpid->size() > 0) {
            px = (*gpx)[0]; py = (*gpy)[0]; pz = (*gpz)[0];
            vx = (*gvx)[0]; vy = (*gvy)[0]; vz = (*gvz)[0];
        }

        tout->Fill();
    }

    fout->cd();
    tout->Write("", TObject::kOverwrite);

    cout << "==============================================" << endl;
    cout << "sim_preprocess: " << fname << endl;
    cout << "Events processed:          " << N << endl;
    cout << "Raw BGO hits (all sectors, any amplitude): " << n_crs_hits_raw << endl;
    cout << "BGO channel-hits after threshold+sum:       " << n_crs_hits_summed << endl;
    cout << "Output written to: " << outname << endl;
    cout << "==============================================" << endl;

    fout->Close();  // avoid stale ROOT cycles (tout;1, tout;2, ...) on re-run
    f->Close();
}
