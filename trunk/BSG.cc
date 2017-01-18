#include <iostream>
#include "spectrum.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <argp.h>

using namespace std;

int main(int argc, char** argv) {
  /** how much information do we write to the screen?
   * 0 - nothing
   * 1 - normal
   * 2 - not implemented
   */
  int verbose = 1;

  /// default constructor, we will initialize later
  Spectrum* mySpectrum = new Spectrum();

  /// step with which the SpectrumShape is generated, in keV
  double StepSize = 0.1;

  /// file name where to write the resulting spectrum, will be overwritten!
  char* ResultsFileName = NULL;

  /// file name from where we read the decay information
  char* IniFileName = NULL;

  // stream for the resulst
  ofstream ofile;

  // stream for the decay information output
  ostream* infostream = &cout;

  // parse program options
  char c;
  while ((c = getopt(argc, argv, "QqDapnhfbwcClRrxmuvi:o:s:")) != -1) {
    // cout <<"Parameter " << c << endl;
    switch (c) {
      case 'p':
        mySpectrum->SetPhaseSpace(false);
        break;
      case 'f':
        mySpectrum->SetFermiFunction(false);
        break;
      case 'c':
        mySpectrum->SetCCorrection(false);
        break;
      case 'C':
        mySpectrum->SetCICorrection(false);
        break;
      case 'R':
        mySpectrum->SetRelativisticCorrection(false);
        break;
      case 'D':
        mySpectrum->SetDeformationCorrection(false);
        break;
      case 'l':
        mySpectrum->SetL0Correction(false);
        break;
      case 'u':
        mySpectrum->SetUCorrection(false);
        break;
      case 'Q':
        mySpectrum->SetQCorrection(false);
        break;
      case 'r':
        mySpectrum->SetRadiativeCorrection(false);
        break;
      case 'n':
        mySpectrum->SetRecoilCorrection(false);
        break;
      case 'a':
        mySpectrum->SetAtomicScreeningCorrection(false);
        break;
      case 'x':
        mySpectrum->SetAtomicExchangeCorrection(false);
        break;
      case 'm':
        mySpectrum->SetAtomicMismatchCorrection(false);
        break;
      case 'v':
        mySpectrum->SetCalculateNeutrinoSpectrum(true);
        break;
      case 'b':
        mySpectrum->SetWeakMagnetism(atof(optarg));
        break;
      case 'o':
        ResultsFileName = optarg;
        break;
      case 'i':
        IniFileName = optarg;
        break;
      case 'q':  // quiet mode, don't write initialization information
        verbose = 0;
        break;
      case 's':
        StepSize = atof(optarg);
        cout << "StepSize: " << StepSize << endl;
        break;
      case 'w':  // write initialization information into output file
        infostream = &ofile;
        break;
      case 'h':
      default:
        cout << "Usage: betaspectgen [OPTIONS] -i INIFILE -o OUTPUTFILE"
             << endl;
        cout << "Options:" << endl;
        cout << "\t -p \t turn off phase space factors" << endl;
        cout << "\t -f \t turn off Fermi function" << endl;
        cout << "\t -l \t turn off L0 corrections" << endl;
        cout << "\t -c \t turn off C corrections" << endl;
        cout << "\t -C \t turn off CI corrections" << endl;
        cout << "\t -R \t turn off relativistic corrections" << endl;
        cout << "\t -D \t turn off deformation corrections" << endl;
        cout << "\t -u \t turn off U corrections" << endl;
        cout << "\t -Q \t turn off Q corrections" << endl;
        cout << "\t -r \t turn off radiative corrections" << endl;
        cout << "\t -n \t turn off recoiling nucleus corrections" << endl;
        cout << "\t -a \t turn off atomic screening corrections" << endl;
        cout << "\t -x \t turn off atomic exchange corrections" << endl;
        cout << "\t -m \t turn off atomic mismatch corrections" << endl;
        cout << "\t -b \t shortcut for specifying the b/Ac value" << endl;
        cout << "\t -s STEPSIZE \t set energy steps for the spectrum generator"
             << endl;
        cout << "\t -q \t quiet mode, don't write anything to the screen"
             << endl;
        /*cout << "\t -e FILENAME \t enable recoil corrections, read them from "
                "FILENAME" << endl;*/
        cout << "\t -w \t write decay information and other options into "
                "OUTPUTFILE" << endl;
        cout << "\t -v \t calculate accompanying neutrino spectrum" << endl;
        return 0;
    }
  }

  if (IniFileName == NULL) {
    cout << "Usage: betaspectgen [OPTIONS] -i INIFILE -o OUTPUTFILE" << endl;
    return 0;
  }
  mySpectrum->ReadFromFile(IniFileName);
  if (verbose > 0) cout << "Initializing from: " << IniFileName << endl;
  mySpectrum->Initialize();

  if (ResultsFileName == NULL) {
    // how to allocate memory for this unfortunate variable?
    char tempfname[100];
    sprintf(tempfname, "%s.spectrum", IniFileName);
    if (verbose > 0) cout << "Output to " << tempfname << endl;
    ResultsFileName = &tempfname[0];
  }
  if (verbose > 0) cout << "Writing into: " << ResultsFileName << endl;

  ofile.open(ResultsFileName, std::ofstream::out);
  if (!ofile.is_open()) {
    cerr << "Error opening " << ResultsFileName << endl;
    return 1;
  }

  if (verbose > 0) mySpectrum->PrintStatus(infostream);
  mySpectrum->WriteSpectrumToStream(&ofile, StepSize);

  ofile.close();

  return 0;
}
