#include "LCIOToFile.h"

LCIOToFile aLCIOToFile;


LCIOToFile::LCIOToFile() : Processor("LCIOToFile") {


  // register steering parameters: name, description, class-variable, default value

	registerProcessorParameter( "Printing" ,
	                            "Print certain messages"  ,
	                             _printing,
	                             (int)5 ) ;

	
   	std::string outFilename = "x.trks";
	registerProcessorParameter( "outFilename" ,
	                            "output flat file name"  ,
	                             _outFilename,
	                             outFilename ) ;

	registerProcessorParameter( "RW" ,
	                            "Read from flat file = 1 or Write to flat file = 2"  ,
	                             _RW,
	                            0 ) ;
	

   	std::string inputTrackCollectionName = "x";
  	registerInputCollection( LCIO::TRACK,
                                 "InputTrackCollectionName" ,
                                 "Input Track Collection Name " ,
                                 _inputTrackCollectionName,
                                 inputTrackCollectionName);  

  	std::string outputParticleCollectionName = "x";
  	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             	"OutputParticleCollectionName" ,
			     	"Output Particle Collection Name "  ,
                             	_outputParticleCollectionName,
                             	outputParticleCollectionName);

}

void LCIOToFile::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  if(_RW= 2){
 	
  	file.open (_outFilename);
  }

 


  // usually a good idea to
  printParameters() ;
  nEvt = 0;
 
  BField = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();

}

void LCIOToFile::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ;
}
void LCIOToFile::printTrack(Track* t){
	std::cout<<"Track: (d0,phi,ome,z0,tanL) "<< 
		t->getD0()<<" "<<
		t->getPhi()<<" "<<
		t->getOmega()<<" "<<
		t->getZ0()<<" "<<
		t->getTanLambda()<<std::endl;
}

bool LCIOToFile::FindTracks( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _trackvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
    if(*itname==_inputTrackCollectionName){
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
        Track* track = dynamic_cast<Track*>(col->getElementAt(i));
        _trackvec.push_back(track);
      }
    }
  }

  if(_printing>1)std::cout << "FindTracks : " << tf << std::endl;

  return tf;
}
std::vector<double> LCIOToFile::getTrackPxPyPz(Track* t){
	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = BField*c*mm2m*eV2GeV;
 
 	
	/*double cosLambda = 1 / sqrt(1 + t->getTanLambda()*t->getTanLambda() );
	double P = (eB/fabs(t->getOmega()))/cosLambda;
	double sinLambda = t->getTanLambda()*cosLambda;
	double cosPhi = cos(t->getPhi());
	double sinPhi = sin(t->getPhi());
	double px = P*cosLambda*cosPhi;
	double py = P*cosLambda*sinPhi;
	double pz = P*sinLambda;
	*/
	double omega = t->getOmega();
	
	double q = omega/fabs(omega);
	
	double cosPhi = cos(t->getPhi());
	double sinPhi = sin(t->getPhi());	

	double px = q*eB/omega  * cosPhi;
	double py = q*eB/omega  * sinPhi;
	double pz = q*eB/omega * t->getTanLambda();

	std::vector<double> txtytz;
	txtytz.push_back(px);
	txtytz.push_back(py);
	txtytz.push_back(pz);
	return txtytz;
}


std::vector<double> getTrackXYZ(Track* t){
	std::vector<double> xyz(3);
	const float* ref = t->getReferencePoint();
	double phi = t->getPhi();
	double d0 = t->getD0();
	double z0 = t->getZ0();
	xyz.at(0) = -d0*sin(phi) + ref[0];
	xyz.at(1) = d0*cos(phi) + ref[1];
	xyz.at(2) = z0;
	//do things
	return xyz;
}

void LCIOToFile::processEvent( LCEvent * evt ) {
 //FindMCParticles(evt);
// = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 FindTracks(evt);
  
 // LCCollectionVec * partCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 //EVENT::LCCollection* partCollection = evt->getCollection("NewPfoCol");

  streamlog_out(MESSAGE) << " start processing event " << std::endl;
 

   //write to file stuff
    if(_RW == 2){
 	file<<nEvt<<" "<<_trackvec.size()<<std::endl;

	//loop over tracks
	for(unsigned int i =0; i<_trackvec.size(); i++){
		Track* t = _trackvec.at(i);		
		file<<t->getD0()<<" "<<t->getPhi()<<" "<<t->getOmega()<<" "<<t->getZ0()<<" "<<t->getTanLambda()<<" ";
		std::vector<double> xyz = getTrackXYZ(t);
		file<<xyz.at(0)<<" "<<xyz.at(1)<<" "<<xyz.at(2)<<std::endl;
		std::vector<float> cov = t->getCovMatrix();
		file<< sqrt(cov.at(0)) <<" "<< sqrt(cov.at(2)) <<" "<<sqrt(cov.at(5))<<" "<<sqrt(cov.at(9))<<" "<<sqrt(cov.at(14))<< " ";
		file<< 0.0 <<" "<<0.0<<" "<<0.0<<std::endl;
		
	}
    }

  

  nEvt++;

  // Add new collection to event
//comment this next line when appending to collection
 // evt->addCollection(partCollection , _outputParticleCollectionName.c_str() ); 
 std::cout << "======================================== event " << nEvt << std::endl ;

}
void LCIOToFile::end(){
	if(_RW == 2){
		file.close();
	}
}

