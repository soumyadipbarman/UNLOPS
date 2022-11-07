class VetoSecondEmissionDyn : public UserHooks {

public:

  VetoSecondEmissionDyn(Pythia* pythiaPtrIn, int nFinalOffsetIn) : nISR(0), nFSR(0),
  npNLO(0), pythiaPtr(pythiaPtrIn), nFinalOffset(nFinalOffsetIn) {}

  bool canVetoISREmission() { return true; }
  bool canVetoFSREmission() { return true; }

  bool doVetoISREmission(int, const Event&, int/*, bool acceptEmission*/) {
    bool doVeto = false;
    // acceptEmission needed if variations are on, since veto comes later!!!!
    // Assume off for now
    bool acceptEmission = true;
    // Veto as if emission has been accepted
    if (nISR + nFSR >= 1 || isReal) doVeto = true;
    // Keep accepted emissions only
    if (acceptEmission) nISR++;
    // Done.
    return doVeto;
  }

  bool doVetoFSREmission(int, const Event&, int, bool/* acceptEmission */) {
    // acceptEmission needed if variations are on, since veto comes later!!!!
    // Assume off for now
    bool acceptEmission = true;
    bool doVeto = false;
    // Veto as if emission has been accepted
    if (nISR + nFSR >= 1 || isReal) doVeto = true;
    // Keep accepted emissions only
    if (acceptEmission) nFSR++;
    // Done.
    return doVeto;
  }

  bool canVetoProcessLevel() { return true; }
  bool doVetoProcessLevel(Event& process) {
    // Initailize and store resonance decay products.
    nISR = nFSR = 0;

    string nps_nlo = infoPtr->getEventAttribute("npNLO",true);
    npNLO     = (nps_nlo != "") ? atoi((char*)nps_nlo.c_str()) : -1;

    nFinal  = 0 ;
    for (int i=0; i < process.size(); ++i) if (process[i].isFinal()) nFinal++;
    const int nFinalProcess = 2;
    int npNow = nFinal-nFinalProcess;
    if (npNLO < 0 || npNow > npNLO) isReal=true;
    else isReal = false;

    //cout << endl << endl << "Parsing the following process: " << endl;
    //process.list();
    //cout << "Found npNow, npNLO = " << npNow << ", " << npNLO << endl;
    //cout << "Setting isReal = " << isReal << endl;
    return false;
  }

  // Members.
  int nISR, nFSR, npNLO, nFinal;
  bool isReal;
  Pythia* pythiaPtr;
  int nFinalOffset;
};
