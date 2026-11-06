/*
  sim_event_builder.C

  Simulation analogue of Marco's event_builder.C (BGO_testing-main).

  Reads the file produced by sim_preprocess.C and:
    1. Requires a top+bottom plastic coincidence at the SAME Z position
       (sector pair: 56&60, 57&61, 58&62, 59&63) -> defines "seed" 0-3,
       exactly like the experimental trigger pairs. The coincidence must
       be EXCLUSIVE: exactly one top plastic and exactly one bottom
       plastic fired in the event, and nothing else - any additional
       plastic hit (top or bottom, any Z) disqualifies the event.
    2. Records (but does NOT cut on) the primary muon truth direction angle,
       measured from the VERTICAL axis, which is Y in this geometry
       (plastic_up at y=+8cm, plastic_down at y=-8cm; see geometry_crs.pl).
       theta = acos(|py| / |p|). This lets the later analysis stage show
       the angular distribution of coincidence events and decide on a cut
       there, rather than baking one in here.
    3. Writes the final flat event tree in Marco's schema, extended with
       theta so it survives into the analysis stage:
         t{ ch (0-54), seed (0-3), q (energy, keV), theta (deg) }
       (ch/seed/q match the branches read by ana_fit.C; theta is additive.)
    4. ALSO writes a per-event tree "tevt" (one row per coincidence event,
       vectors preserved) so downstream analysis can inspect individual
       events (e.g. high-theta outliers) in full: which BGO crystals fired
       with what energy, and which plastic sectors fired with what
       amplitude.

  Usage:
    root -l -q 'sim_event_builder.C("../sim_data/out_60crs_8.7.root", 1)'
    (2nd arg = batch number, selects data/batch_<batch>_map.dat for the
     channel -> crystal serial number lookup; default 1.)

  Output:
    ./out/<basename>_evt.root
      tree "t"      : ch, seed, q, theta          (one row per BGO channel hit)
      tree "tevt"   : evt, seed, theta,
                      ch (vector<int>), crs (vector<int>, mapped serial #),
                      q (vector<double>),
                      scint_sector (vector<int>), scint_A (vector<double>)
                      (one row per coincidence EVENT)
      hist "htheta" : theta for ALL coincidence events (one entry/event)
*/

#include <vector>
#include <string>
using namespace std;

void sim_event_builder(const char* filename = "../sim_data/out_typ8.root", int batch = 1) {

    string fname = filename;
    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');
    const string basename = fname.substr(start, end - start);
    const string inname  = "./out/preprocessed/" + basename + "_preprocessed.root";
    const string outname = "./out/evt/" + basename + "_evt.root";

    gSystem->Unlink(outname.c_str());  // remove any stale file to avoid ROOT cycle buildup

    // ── Load crystal-serial lookup map (channel index -> serial number) ──
    // Same convention as sim_ana_fit.C: batch_<batch>_map.dat lists crystals
    // starting from the BOTTOM row of the stack (ch=0..4), direct 1:1 mapping.
    ifstream crs_map_file(Form("data/batch_%i_map.dat", batch));
    map<int, int> crs_map;
    int crs_id, crs_counter = 0;
    while (crs_map_file >> crs_id) {
        crs_counter++;
        crs_map[crs_counter] = crs_id;
    }
    for (int ii = 56; ii < 60; ii++) crs_map[ii] = 8000 + ii;
    for (int ii = 60; ii < 64; ii++) crs_map[ii] = 9000 + ii;

    TFile *fin = new TFile(inname.c_str(), "read");
    if (!fin || fin->IsZombie()) {
        cout << "Error: cannot open " << inname
             << " -- did you run sim_preprocess.C first?" << endl;
        return;
    }
    TTree *tin = (TTree*) fin->Get("tout");
    if (!tin) { cout << "Error: tree 'tout' not found in " << inname << endl; return; }

    vector<int>    *ch_v = nullptr, *scint_sector_v = nullptr;
    vector<double> *q_v = nullptr,  *scint_A_v = nullptr;
    double px, py, pz;
    double vx, vy, vz;

    tin->SetBranchAddress("ch",           &ch_v);
    tin->SetBranchAddress("q",            &q_v);
    tin->SetBranchAddress("scint_sector", &scint_sector_v);
    tin->SetBranchAddress("scint_A",      &scint_A_v);
    tin->SetBranchAddress("px", &px);
    tin->SetBranchAddress("py", &py);
    tin->SetBranchAddress("pz", &pz);
    tin->SetBranchAddress("vx", &vx);
    tin->SetBranchAddress("vy", &vy);
    tin->SetBranchAddress("vz", &vz);

    TFile *fout = new TFile(outname.c_str(), "recreate");
    TTree *tout = new TTree("t", "t");
    int chout, seedout;
    double qout, thetaout;
    tout->Branch("ch",    &chout);
    tout->Branch("seed",  &seedout);
    tout->Branch("q",     &qout);
    tout->Branch("theta", &thetaout);  // deg from vertical (Y) - NOT cut on here

    // ── Per-event tree (one row per coincidence event, vectors preserved) ───
    TTree *tevt = new TTree("tevt", "Per-event coincidence data (vectors preserved)");
    int evt_evt, evt_seed;
    double evt_theta;
    double evt_px, evt_py, evt_pz;
    double evt_vx, evt_vy, evt_vz;
    vector<int>    evt_ch_v, evt_crs_v, evt_scint_sector_v;
    vector<double> evt_q_v,  evt_scint_A_v;
    tevt->Branch("evt",          &evt_evt);
    tevt->Branch("seed",         &evt_seed);
    tevt->Branch("theta",        &evt_theta);
    tevt->Branch("px", &evt_px);
    tevt->Branch("py", &evt_py);
    tevt->Branch("pz", &evt_pz);
    tevt->Branch("vx", &evt_vx);
    tevt->Branch("vy", &evt_vy);
    tevt->Branch("vz", &evt_vz);
    tevt->Branch("ch",           &evt_ch_v);           // BGO channel, 0-54
    tevt->Branch("crs",          &evt_crs_v);          // BGO crystal serial number (crs_map lookup)
    tevt->Branch("q",            &evt_q_v);            // BGO energy, keV
    tevt->Branch("scint_sector", &evt_scint_sector_v); // plastic sector, 56-63
    tevt->Branch("scint_A",      &evt_scint_A_v);      // plastic amplitude, keV

    TH1D *htheta = new TH1D("htheta",
        "Primary muon angle from vertical (Y) - coincidence events; #theta (deg); Events", 90, 0, 90);

    Long64_t N = tin->GetEntries();
    long nCoinc = 0, nHitsWritten = 0;

    for (Long64_t ii = 0; ii < N; ii++) {
        tin->GetEntry(ii);

        // ── Require EXCLUSIVE same-Z top+bottom plastic coincidence ─────────
        // Exactly one top plastic AND exactly one bottom plastic fired, and
        // they must be at the same Z position (seed pair). Any extra plastic
        // hit (top or bottom, any Z) disqualifies the event.
        int nTop = 0, nBot = 0;
        int topSeed = -1, botSeed = -1;
        for (size_t jj = 0; jj < scint_sector_v->size(); jj++) {
            int sec = (*scint_sector_v)[jj];
            if (sec >= 56 && sec <= 59) { nTop++; topSeed = sec - 56; }
            else if (sec >= 60 && sec <= 63) { nBot++; botSeed = sec - 60; }
        }
        if (nTop != 1 || nBot != 1) continue;      // more/fewer than 1+1 hits
        if (topSeed != botSeed) continue;          // not the same Z position
        int pair = topSeed;
        nCoinc++;

        // ── Record (do not cut on) the primary muon truth direction angle ───
        double p = sqrt(px*px + py*py + pz*pz);
        double theta = (p > 0) ? acos(fabs(py) / p) * 180.0 / TMath::Pi() : 999;
        htheta->Fill(theta);

        // ── Write one output row per BGO channel hit in this event ──────────
        int seed = pair;  // 0-3, matches Marco's trigger-position index
        for (size_t jj = 0; jj < ch_v->size(); jj++) {
            chout    = (*ch_v)[jj] - 1;   // convert 1-55 -> 0-54 (Marco's convention)
            seedout  = seed;
            qout     = (*q_v)[jj];
            thetaout = theta;
            tout->Fill();
            nHitsWritten++;
        }

        // ── Write one per-event row (vectors preserved) ──────────────────────
        evt_evt   = (int) ii;
        evt_seed  = seed;
        evt_theta = theta;
        evt_px = px; evt_py = py; evt_pz = pz;
        evt_vx = vx; evt_vy = vy; evt_vz = vz;
        evt_ch_v.clear(); evt_crs_v.clear(); evt_q_v.clear();
        evt_scint_sector_v.clear(); evt_scint_A_v.clear();
        for (size_t jj = 0; jj < ch_v->size(); jj++) {
            int chMapped = (*ch_v)[jj] - 1;  // 0-54
            evt_ch_v.push_back(chMapped);
            evt_crs_v.push_back(crs_map[chMapped + 1]);  // crs_map keyed 1-55
            evt_q_v.push_back((*q_v)[jj]);
        }
        for (size_t jj = 0; jj < scint_sector_v->size(); jj++) {
            evt_scint_sector_v.push_back((*scint_sector_v)[jj]);
            evt_scint_A_v.push_back((*scint_A_v)[jj]);
        }
        tevt->Fill();
    }

    fout->cd();
    tout->Write("", TObject::kOverwrite);
    tevt->Write("", TObject::kOverwrite);
    htheta->Write("", TObject::kOverwrite);

    cout << "==============================================" << endl;
    cout << "sim_event_builder: " << inname << endl;
    cout << "Preprocessed events read:                " << N << endl;
    cout << "Top+bottom same-pair coincidence: " << nCoinc << endl;
    cout << "Final BGO channel-hits written:            " << nHitsWritten << endl;
    cout << "Output written to: " << outname << endl;
    cout << "(theta stored per-hit in tree 't' and per-event in 'htheta')" << endl;
    cout << "==============================================" << endl;

    fout->Close();  // avoid stale ROOT cycles (t;1, t;2, ...) on re-run
    fin->Close();
}
