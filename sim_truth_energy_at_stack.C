/*
  sim_truth_energy_at_stack.C

  Finds, event-by-event, which generated (thrown) muons actually PUNCHED
  THROUGH the 3 concrete floors and deposited energy along a genuine
  straight-line track in the BGO crystal stack ("crs" tree), and extracts
  their GENERATOR-LEVEL truth momentum/energy (from the "generated" tree,
  i.e. BEFORE any energy loss in the concrete floors). This is the number
  you want to use to set a minimum momentum/energy cutoff in the
  COSMICRAYS gcard option, so you stop wasting statistics on muons that
  can never reach the stack.

  IMPORTANT -- rejecting secondary-particle contamination:
  A naive "any crs hit" selection also picks up events where the primary
  muon stopped/showered in the concrete and only a LOW-ENERGY SECONDARY
  (knock-on electron, scattered gamma, etc.) leaked into the stack and
  deposited a small, isolated amount of energy in an off-axis crystal.
  These are NOT muons that reached the stack.

  The crystal stack is a N_COL x N_ROW (5 x 11) array. A muon that truly
  travels straight down through it deposits a near-MIP energy (tens of
  MeV, given the dense BGO) in SEVERAL CONSECUTIVE CRYSTALS OF THE SAME
  COLUMN (one column = one vertical row of crystals the muon track stays
  in as it descends). Secondaries instead show up as isolated, low-energy,
  off-axis hits scattered into other columns.

  So here we require, per event:
    1) hits above a MIP-like energy threshold (TRACK_E_MIN, keV),
    2) at least MIN_TRACK_HITS of them falling in the SAME dominant
       column (the "track column"),
  before accepting the event as a genuine muon that reached the stack.
  Stray hits in other columns (secondaries) are simply ignored/discarded,
  they do not spoil the event if a real track is also present, but an
  event with ONLY scattered off-axis hits (no real column track) is
  correctly rejected.

  Usage:
    root -l -q 'sim_truth_energy_at_stack.C("../sim_data/out_planegen_10.7.root")'

  Output:
    - prints: events processed, events with a genuine muon track reaching
      the stack, and min/max/mean generator-level momentum and kinetic
      energy for those events
    - writes ./out/<basename>_truthE.root containing:
        TH1D "hP_all"    - generator momentum, ALL events
        TH1D "hP_hit"    - generator momentum, events with a genuine muon track
        TH1D "hKE_hit"   - generator kinetic energy, events with a genuine muon track
*/

#include <vector>
#include <string>
#include <map>

using namespace std;

static const int N_COL = 5;
static const int N_ROW = 11;

void sim_truth_energy_at_stack(const char* filename = "../sim_data/out_planegen_10.7.root") {

    const double muMass = 0.105658; // GeV, muon mass

    string fname = filename;
    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');
    const string basename = fname.substr(start, end - start);
    const string outname  = "./out/" + basename + "_truthE.root";

    gSystem->mkdir("./out", true);
    gSystem->Unlink(outname.c_str());

    TFile *f = new TFile(fname.c_str(), "read");
    if (!f || f->IsZombie()) { cout << "Error: cannot open " << fname << endl; return; }

    TTree *tcrs = (TTree*) f->Get("crs");
    TTree *tgen = (TTree*) f->Get("generated");
    if (!tcrs) { cout << "Error: tree 'crs' not found" << endl; return; }
    if (!tgen) { cout << "Error: tree 'generated' not found" << endl; return; }

    Long64_t N = tgen->GetEntries();
    if (tcrs->GetEntries() != N) {
        cout << "Error: crs/generated trees have different entry counts ("
             << tcrs->GetEntries() << " vs " << N << ")" << endl;
        return;
    }

    vector<double> *cSector = nullptr;
    vector<double> *cADC    = nullptr;
    tcrs->SetBranchAddress("sector",  &cSector);
    tcrs->SetBranchAddress("ADC_ADC", &cADC);

    vector<int>    *gpid = nullptr;
    vector<double> *gpx = nullptr, *gpy = nullptr, *gpz = nullptr;
    tgen->SetBranchAddress("pid", &gpid);
    tgen->SetBranchAddress("px",  &gpx);
    tgen->SetBranchAddress("py",  &gpy);
    tgen->SetBranchAddress("pz",  &gpz);

    // ── Track-selection tuning (rejects secondary-only hits) ──────────────
    const double TRACK_E_MIN    = 5000.0; // keV, MIP-like single-crystal threshold
    const int    MIN_TRACK_HITS = 3;      // min # of same-column MIP hits required

    TFile *fout = new TFile(outname.c_str(), "recreate");
    TH1D *hP_all  = new TH1D("hP_all",  "Generator momentum, ALL events;p (GeV);Events",      100, 0, 100);
    TH1D *hP_hit  = new TH1D("hP_hit",  "Generator momentum, events with genuine muon track;p (GeV);Events", 100, 0, 100);
    TH1D *hKE_hit = new TH1D("hKE_hit", "Generator kinetic energy, events with genuine muon track;KE (GeV);Events", 100, 0, 100);

    long nHitEvents = 0;
    double pMin = 1e9, pMax = -1e9, pSum = 0;
    double keMin = 1e9, keMax = -1e9, keSum = 0;

    for (Long64_t ii = 0; ii < N; ii++) {
        tgen->GetEntry(ii);
        tcrs->GetEntry(ii);

        if (gpid->size() == 0) continue; // no primary recorded (shouldn't happen)

        double px = (*gpx)[0], py = (*gpy)[0], pz = (*gpz)[0];
        double p  = sqrt(px*px + py*py + pz*pz) / 1000.0; // MeV -> GeV

        hP_all->Fill(p);

        // ── Count MIP-like hits per column, keep only genuine straight tracks ──
        map<int,int> colCount; // column -> # of hits above TRACK_E_MIN in that column
        for (size_t jj = 0; jj < cSector->size(); jj++) {
            if ((*cADC)[jj] < TRACK_E_MIN) continue;  // reject low-E secondary-like hits
            int ch0 = (int)(*cSector)[jj] - 1;        // 0-54
            int col = 4 - (ch0 % N_COL);
            colCount[col]++;
        }
        int bestColHits = 0;
        for (auto &kv : colCount) if (kv.second > bestColHits) bestColHits = kv.second;

        bool hasTrack = (bestColHits >= MIN_TRACK_HITS);
        if (!hasTrack) continue;

        nHitEvents++;
        double ke = sqrt(p*p + muMass*muMass) - muMass;

        hP_hit->Fill(p);
        hKE_hit->Fill(ke);

        if (p < pMin) pMin = p;
        if (p > pMax) pMax = p;
        pSum += p;

        if (ke < keMin) keMin = ke;
        if (ke > keMax) keMax = ke;
        keSum += ke;
    }

    fout->cd();
    hP_all->Write();
    hP_hit->Write();
    hKE_hit->Write();

    cout << "==============================================" << endl;
    cout << "sim_truth_energy_at_stack: " << fname << endl;
    cout << "Events processed:                      " << N << endl;
    cout << "Track selection: >=" << MIN_TRACK_HITS << " hits in same column, each >= "
         << TRACK_E_MIN << " keV" << endl;
    cout << "Events with genuine muon track in stack: " << nHitEvents << endl;
    if (nHitEvents > 0) {
        cout << "--- generator-level momentum (GeV) for stack-reaching muons ---" << endl;
        cout << "  min  = " << pMin << endl;
        cout << "  max  = " << pMax << endl;
        cout << "  mean = " << pSum / nHitEvents << endl;
        cout << "--- generator-level kinetic energy (GeV), muMass=" << muMass << " GeV ---" << endl;
        cout << "  min  = " << keMin << endl;
        cout << "  max  = " << keMax << endl;
        cout << "  mean = " << keSum / nHitEvents << endl;
        cout << endl;
        cout << "Suggested COSMICRAYS cminp ~ " << keMin
             << " GeV (kinetic) i.e. momentum ~ " << pMin << " GeV" << endl;
    } else {
        cout << "No events with a crs hit found - nothing to report." << endl;
    }
    cout << "Output written to: " << outname << endl;
    cout << "==============================================" << endl;

    fout->Close();
    f->Close();
}
