// This is a simple test program. It performs the first emission on 
// aMC@NLO like S events to prepare for further processing. 

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <sstream>

using namespace Pythia8;

#include "VetoSecondEmissionDyn.h"

double MASSTOL=1e-4;
double MOMTOL=1e-7;

double val ( double in, double tol=1e-7){
  if (abs(in) < tol) return 0.0;
  return in;
}

void print_end_tag( ogzstream & outfile ) { 
  outfile << "</LesHouchesEvents>" << "\n";
}

void print_begin_tag( ogzstream & outfile, const Info* infoPtr, 
  igzstream & infile, ifstream & settingsfile ){

  outfile << "<LesHouchesEvents version=\"3.0\">\n";
  outfile << "<header>\n";
  string line, tag;
  // Print input file header.
  do {
    if (!getline(infile, line)) break;
    if (line.find_first_not_of(" \n\t\v\b\r\f\a") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) break;
    }
    outfile << "# " << line << "\n";
  } while (tag != "<init>" && tag != "<init" && tag != "</header>");
  // Print Pythia settings
  if (settingsfile.is_open()) {
    while ( getline(settingsfile,line) )
      outfile << "# " << line << "\n";
    settingsfile.close();
  }
  infoPtr->initrwgt->list( outfile );
  outfile << "</header>\n";
  outfile << "<init>\n"
          << "  " << infoPtr->idA() << "  " << infoPtr->idB()
          << "  " << infoPtr->eA() << "  " << infoPtr->eB()
          << "  0  0  0  0  -4  1\n"
          << "  "
          << scientific << setprecision(12)
          << infoPtr->sigmaGen()*1e9 << "  0.000000e+00  1.000000e+00   1\n"
          << "</init>\n";
}

void print_event_begin_tag( int npLO, int npNLO, ogzstream & outfile ) { 
  outfile << "<event npLO=\" " << npLO << " \" npNLO=\" " << npNLO << " \">\n";
}

void print_event_begin_tag( int npLO, int npNLO, int zdir, double emissionscale,  ogzstream & outfile ) { 
  outfile << "<event npLO=\" " << npLO << " \" npNLO=\" " << npNLO << " \" zsgn=\" " <<zdir << " \" lastpt=\" " << emissionscale <<" \">\n";
}

void print_event_end_tag( ogzstream & outfile ){ outfile << "</event>\n"; }

void print_event_info( ogzstream & outfile, const Event& event, const Info* infoPtr ){
    int nParticles = 2;
    // Get number of intermediate bosons.
    for (int i = 1; i < event.size(); ++i) {
      if ( ( event[i].status() == -22
          || event[i].status() == -44
          || event[i].status() == -51
          || event[i].status() == -52)
        && event[i].canDecay()
        && event[i].mayDecay()
        && event[i].isResonance()
        && event[i].daughter1() != event[i].daughter2()){
        int iD1 = event[i].daughter1();
        int iD2 = event[i].daughter2();
        if ( iD1 != 0 && iD2 != 0 && event[iD1].canDecay()  
          && event[iD1].mayDecay() && event[iD1].isResonance()
          && event[iD1].id() == event[i].id()) continue;
        if ( iD1 != 0 && iD2 != 0 && event[iD2].canDecay()  
          && event[iD2].mayDecay() && event[iD2].isResonance()
          && event[iD2].id() == event[i].id()) continue;
        nParticles++;
      }
    }

    for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal()) nParticles++;
    outfile << "  " << nParticles << "  1  " << scientific << setprecision(12)
            << infoPtr->weight()
            << "  "  << infoPtr->QFac() << "  " << infoPtr->alphaEM()
            << "  " << infoPtr->alphaS() << "\n";
}

void print_event( ogzstream & outfile, const Event& event, const Info* ,
   ParticleData* particleDataPtr){

  // Get index of first incoming
  int in1 = 0;
  for (int i=event.size()-1; i > 0; --i)
    if(abs(event[i].status()) == 41 ){
      in1 = i;
      break;
    }

  // Get index of second incoming
  int in2 = 0;
  for (int i=event.size()-1; i > 0; --i)
    if(abs(event[i].status()) == 42 ){
      in2 = i;
      break;
    }

  // If no incoming of the cascade are found, try incoming
  if(in1 == 0){
    for(int i=0; i < int(event.size()); ++i)
      if(event[i].mother1() == 1) in1 = i;
  }
  if(in2 == 0){
    for(int i=0; i < int(event.size()); ++i)
      if(event[i].mother1() == 2) in2 = i;
  }

  if (event[in1].pz() < 0.) {int temp=in1; in1=in2; in2=temp;}

  // Construct the output event
  Event NewEvent = Event();
  NewEvent.init("(dummy)", particleDataPtr);
  NewEvent.clear();

  int i1 = NewEvent.append(event[in1]);
  NewEvent[i1].mothers(-1,0);
  int i2 = NewEvent.append(event[in2]);
  NewEvent[i2].mothers(-1,0);
 
  // Attach resonances and decay products.
  vector<int> iAttached;
  // Loop through event and attach remaining decays
  int iDec = 0;
  do {

    if ( ( event[iDec].status() == -22
        || event[iDec].status() == -44
        || event[iDec].status() == -51
        || event[iDec].status() == -52)
      && event[iDec].canDecay()
      && event[iDec].mayDecay()
      && event[iDec].isResonance()
      && event[iDec].daughter1() != event[iDec].daughter2() ) {

      bool skip = false;
      for (int j = 0; j < int(iAttached.size()); ++j)
        if (iDec == iAttached[j]) {skip = true; break;}

      int iD1 = event[iDec].daughter1();
      int iD2 = event[iDec].daughter2();

      if ( iD1 != 0 && iD2 != 0 && event[iD1].canDecay()
        && event[iD1].mayDecay() && event[iD1].isResonance()
        && event[iD1].id() == event[iDec].id()) continue;
      if ( iD1 != 0 && iD2 != 0 && event[iD2].canDecay()
        && event[iD2].mayDecay() && event[iD2].isResonance()
        && event[iD2].id() == event[iDec].id()) continue;

      int iRes = 0;
      if (!skip){
        iRes = NewEvent.append(event[iDec]);
        NewEvent[iRes].mothers(0,1);
        NewEvent[iRes].status(-22);
        iAttached.push_back(iDec);  
      }

      // Done if no daughters exist.
      if ( iD1 == 0 || iD2 == 0 ) continue;

      // If the resonance had already been attached earlier, find it in NewEvent, in order
      // to get the correct "mother" index for the daughters to be attached.
      if (iRes == 0) {
        for (int i = 0; i < NewEvent.size(); ++i) {
          Particle res = event[iDec];
          Particle now = NewEvent[i];
          if ( res.id()   == now.id()
            && res.col()  == now.col()
            && res.acol() == now.acol()
            && abs(res.px()-now.px()) < 1e-8
            && abs(res.py()-now.py()) < 1e-8
            && abs(res.pz()-now.pz()) < 1e-8
            && abs(res.e() -now.e())  < 1e-8){
            iRes = i;
            break;
          }
        }
      }

      bool skipD = false;
      for (int j = 0; j < int(iAttached.size()); ++j)
        if (iD1 == iAttached[j] || iD2 == iAttached[j]) {skipD = true; break;}
      if (skipD) continue;

      // Attach daughters.
      for ( int k = iD1; k <= iD2; ++k ) {
        int iNow = NewEvent.append(event[k]);
        NewEvent[iNow].mothers(iRes,0);
        if ( NewEvent[iNow].canDecay() && NewEvent[iNow].mayDecay()
          && NewEvent[iNow].isResonance()) NewEvent[iNow].status(-22);
        iAttached.push_back(k);  
      }

    }

  } while (++iDec < event.size());

  // Attach all remaining particles.
  for (int i = 1; i < event.size(); ++i){

    if (!event[i].isFinal()) continue;
    bool skip = false;
    for (int j = 0; j < int(iAttached.size()); ++j)   if (i == iAttached[j]) {skip = true; break;}
    if (skip) continue;

    int iNow = NewEvent.append(event[i]);
    NewEvent[iNow].mothers(0,1);

  }

  // Write LHEF.
  for (int i = 0; i < NewEvent.size(); ++i){
    int m1 = NewEvent[i].mother1(), m2 = 0;
    int status = (NewEvent[i].status() == -22) ? 2 : (NewEvent[i].mother1() >= 0) ? 1 : 0;
    if (NewEvent[i].mother1() < 0)  {m1 = m2 = 0;}
    if (NewEvent[i].mother1() > 0)  {m1++;}
    if (NewEvent[i].mother1() == 0) {m1 = 1; m2 = 2;}
    outfile << "  " << NewEvent[i].id() << "   " << status << "  " << m1 << "  " << m2 << "  "
        << NewEvent[i].col() << "  " << NewEvent[i].acol()
        << scientific << setprecision(12)
        << "  " << val(NewEvent[i].px()) << "  " << val(NewEvent[i].py())
        << "  " << val(NewEvent[i].pz()) << "  " << val(NewEvent[i].e())
        << "  " << val(NewEvent[i].mCalc(),MASSTOL)
        << "  0  " << NewEvent[i].pol() << "\n";
  }

}

void print_event_weights( ogzstream & outfile, const Info* infoPtr ){
  infoPtr->weights->list( outfile );
  infoPtr->rwgt->list( outfile );
  infoPtr->scales->list( outfile );
}

int main(int argc, char* argv[]) {

  // Check for correct number of arguments
  if (argc < 5) {
    cout << R"(  Please provide the following arguments:
        - config file\n
        - final state multiplicity for np = 0\n
        - input lhe file\n
        - output lhe file)" << endl;
    return 1;
  }
  // Get command-line arguments
  vector<string> arguments;
  for (int i = 1; i < argc; ++i) { 
    arguments.push_back(string(argv[i]));
  }

  Pythia* pythiaPtr = new Pythia();
  shared_ptr<UserHooks> myUserHooksPtr;

  // Read input files.
  string input_file = arguments[2];
  myUserHooksPtr = make_shared<VetoSecondEmissionDyn>(pythiaPtr, atoi(argv[2]));

  // Output lhe file
  string output=string(arguments[3]);
  string outputLHE = output;
  ogzstream outfile; // ogzstream für gzip aus pythia::Streams
  const char* cstring = outputLHE.c_str();
  outfile.open(cstring);

  igzstream infile(input_file.c_str());
  string settingscard = string(arguments[0]);
  ifstream settingsfile(settingscard.c_str());

  pythiaPtr->readFile(arguments[0]);
  pythiaPtr->settings.word("Beams:LHEF", input_file.c_str());
  cout << "Reading events from " << input_file << endl;
  pythiaPtr->settings.mode("Beams:frameType", 4);
  pythiaPtr->setUserHooksPtr( myUserHooksPtr );

  pythiaPtr->init();


  bool printedBeginTag = false;
  
  while (true) {
    if ( !pythiaPtr->next() ) {
      if ( pythiaPtr->info.atEndOfFile() ) { break;}
      else continue;
    }
    if ( !printedBeginTag ) {
      print_begin_tag( outfile, &(pythiaPtr->info), infile, settingsfile);
      printedBeginTag = true;
    }
    string nps_nlo = pythiaPtr->info.getEventAttribute("npNLO",true);
    int npNLO     = (nps_nlo != "") ? atoi((char*)nps_nlo.c_str()) : -1;
    string nps_lo  = pythiaPtr->info.getEventAttribute("npLO",true);
    int npLO      = (nps_lo != "") ? atoi((char*)nps_lo.c_str()) : -1;
    if (npNLO == -1 && npLO == -1) {
      cout << "Error: missing npNLO and npLO attributes for individual events!"
        << endl;
      abort();
    }
    print_event_begin_tag(npLO,npNLO, outfile);
    print_event_info(outfile, pythiaPtr->event, &(pythiaPtr->info));
    print_event(outfile, pythiaPtr->event, &(pythiaPtr->info), &(pythiaPtr->particleData));
    print_event_weights( outfile, &(pythiaPtr->info) );
    print_event_end_tag(outfile);
  }

  pythiaPtr->stat();

  print_end_tag( outfile );
  outfile.close();

  // Done.
  return 0;

}
