#include "marlin/Processor.h"
#include "EVENT/Track.h"
#include "lcio.h"
#include "TFile.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ParticleIDImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include <marlin/Global.h>
#include "gear/BField.h"
#include "TLorentzVector.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace lcio;

	/** LCIOToFile:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class LCIOToFile : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new LCIOToFile ; }

  LCIOToFile(const LCIOToFile&) = delete ;
  LCIOToFile& operator=(const LCIOToFile&) = delete ;

  LCIOToFile() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindTracks(LCEvent* evt);
 
  std::vector<double> getTrackPxPyPz(Track* t);
  std::vector<double> getTrackXYZ(Track* t);
  
  //printing utility
  void printTrack(Track* t);




  protected:
  int nEvt{};
  
  //vector to hold the tracks for the event
  std::vector<Track*> _trackvec{};
  int   _printing{};
  std::string _outFilename{};
  int _RW{};

   //need to no BField to calculate stuff
   double BField{};

//filestreams
	std::ofstream file;
 
  
// _inputTrackCollectionName 
  std::string _outputParticleCollectionName{};
  std::string _inputTrackCollectionName{};
//  std::string m_rootFile{};
};
