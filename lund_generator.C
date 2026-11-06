/*
  lund_generator.C

  Simplified LUND generator for this analysis: reads the per-event tree
  "tevt" written by sim_event_builder.C (e.g. out/final_11.4_evt.root),
  selects events with theta > 45 deg, and writes a LUND file containing
  each selected event's primary muon px, py, pz, vx, vy, vz.

  Header/body line format follows Marco's new_lund_generator.C
  (see lines 138-139 there).
*/

void lund_generator(){
    const char* filename = "out/final_11.4_evt.root";

    double theta_min = 35.0; 

    TFile *f = new TFile(filename, "read");
    TTree *tevt = (TTree*) f->Get("tevt");

    double theta, px, py, pz, vx, vy, vz;
    tevt->SetBranchAddress("theta", &theta);
    tevt->SetBranchAddress("px",    &px);
    tevt->SetBranchAddress("py",    &py);
    tevt->SetBranchAddress("pz",    &pz);
    tevt->SetBranchAddress("vx",    &vx);
    tevt->SetBranchAddress("vy",    &vy);
    tevt->SetBranchAddress("vz",    &vz);

    vector <int> *scint_sector = nullptr;
    vector <double> *scint_A = nullptr;
    tevt->SetBranchAddress("scint_sector", &scint_sector);
    tevt->SetBranchAddress("scint_A", &scint_A);

    ofstream ofile("highangle.lund");

    int PID = 13;   // mu-
    int nSelected = 0;

    Long64_t N = tevt->GetEntries();
    for(Long64_t ii = 0; ii < N; ii++){
        tevt->GetEntry(ii);

        if(theta <= theta_min) continue;
        nSelected++;

        string header = "1 0 0 0 0 0 0 0 0 0";
        string body = "1 0 1 "+to_string(PID)+" 0 0 "+to_string(px/1000.)+" "+to_string(py/1000.)+" "+to_string(pz/1000.)+" 0 0 "+to_string(vx/10.)+" "+to_string(vy/10.)+" "+to_string(vz/10.);
        ofile << header << endl;
        ofile << body << endl;
        
        cout << endl ;

        // get the scint_sector edep for each event
        double q = 0;
        for (size_t i = 0; i < scint_sector->size(); i++) {
            cout << "scint sector: " << (*scint_sector)[i] << ", edep: " << (*scint_A)[i] << endl;
        }
    }

    ofile.close();

    cout << "Events read:      " << N << endl;
    cout << "Events selected (theta > " << theta_min << " deg): " << nSelected << endl;
    cout << "LUND file written to: hightheta_events.lund" << endl;
}
