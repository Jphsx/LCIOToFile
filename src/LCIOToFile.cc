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

	std::string inputMcParticleCollectionName = "x";
	registerInputCollection( LCIO::MCPARTICLE,
				"McParticleCollectionName" ,
				"Name of the MCParticle input collection" ,
				_inputMcParticleCollectionName,
				inputMcParticleCollectionName);

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
bool LCIOToFile::FindMCParticles( LCEvent* evt ){
   
	bool collectionFound = false;

  	// clear old global MCParticle vector
  	_mcpartvec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;

	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    
		//if found print name and number of elements 
		if(*itname==_inputMcParticleCollectionName){
      			LCCollection* collection = evt->getCollection(*itname);
     			std::cout<< "Located MC Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
      			collectionFound = true;
      
			//add the collection elements to the global vector
			for(int i=0;i<collection->getNumberOfElements();i++){
				MCParticle* mcPart = dynamic_cast<MCParticle*>(collection->getElementAt(i));
				_mcpartvec.push_back(mcPart);

       
      			}
    		}
  	}

  	if(!collectionFound){
		std::cout<<"MC Collection "<< _inputMcParticleCollectionName << "not found"<<std::endl;
	}
  
  	return collectionFound;
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


std::vector<double> LCIOToFile::getTrackXYZ(Track* t){
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
std::vector<double> LCIOToFile::findMCTrack(Track* t){
	
	//returns the track parameters and vertex of the mcparticle
	std::vector<float> cov = t->getCovMatrix();
	std::vector<double> errors{};
	errors.push_back(sqrt(cov.at(0)));
	errors.push_back(sqrt(cov.at(2)));
	errors.push_back(sqrt(cov.at(5)));
	errors.push_back(sqrt(cov.at(9)));
	errors.push_back(sqrt(cov.at(14)));

	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = BField*c*mm2m*eV2GeV;

	//try to match
	//int indexOfMatch = -1;
	std::vector<double> mcvec{};//5 track params and xyz	

	for(unsigned int i =0; i<_mcpartvec.size(); i++){
		MCParticle* mcp = _mcpartvec.at(i);
		if(mcp->getCharge() == 0) continue;
		const double* mcvtx = mcp->getVertex();
		const double* mcpp = mcp->getMomentum();
		//use same ref as track
		const float* ref = t->getReferencePoint();
		double pt = sqrt(mcpp[1]*mcpp[1] + mcpp[0]*mcpp[0]);
		double q = (double)mcp->getCharge();

		double phimc = acos(mcpp[0]/pt);
		double d0mc = -(mcvtx[0] - (double)ref[0])*sin(phimc) + (mcvtx[1] - (double)ref[1])*cos(phimc);
		
		double ommc = q*eB/pt;
		double z0mc = mcvtx[2];
		double tlmc = mcpp[2]/pt;
		
		double d0 = t->getD0();
		double phi = t->getPhi();
		double om = t->getOmega();
		double z0 = t->getZ0();
		double tl = t->getTanLambda();

		std::cout<<"comparing tracks: "<<std::endl;
		std::cout<<d0mc<<" "<<phimc<<" "<<ommc<<" "<<z0mc<<" "<<tlmc<<std::endl;
		std::cout<<d0<<" "<<phi<<" "<<om<<" "<<z0<<" "<<tl<<std::endl;

		//allow all matching with 1.5sigma for every parameter
		double factor = 5.0;
		if(// (fabs(d0mc-d0) < factor*errors.at(0)) &&
		   // (fabs(phimc-phi) < factor*errors.at(1)) &&
		    (fabs(ommc-om) < factor*errors.at(2)) ){
		   // (fabs(z0mc-z0) < factor*errors.at(3)) &&
		   // (fabs(tlmc-tl) < factor*errors.at(4)) ){
		
			std::cout<<"found match at index "<<i<<std::endl;
			//indexOfMatch = i;
			mcvec.push_back(d0mc);
			mcvec.push_back(phimc);
			mcvec.push_back(ommc);
			mcvec.push_back(z0mc);
			mcvec.push_back(tlmc);
			mcvec.push_back(mcvtx[0]);
			mcvec.push_back(mcvtx[1]);
			mcvec.push_back(mcvtx[2]);
		}
		
	}
	//MCParticle* mcp = _mcpartvec.at(indexOfMatch);
	//return indexOfMatch;
	//look for -0 and get rid of it
	for(unsigned int i=0; i<mcvec.size(); i++){
		if(mcvec.at(i) == -0.0){
			mcvec.at(i) = 0.0;
		}
	}
	return mcvec;
	

}
void LCIOToFile::processEvent( LCEvent * evt ) {
 //FindMCParticles(evt);
// = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  
  
 // LCCollectionVec * partCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 //EVENT::LCCollection* partCollection = evt->getCollection("NewPfoCol");

  streamlog_out(MESSAGE) << " start processing event " << std::endl;
 
	file.precision(8);
	std::fixed;
   //write to file stuff
    if(_RW == 2){
	FindTracks(evt);
	FindMCParticles( evt );
 	file<<nEvt<<" "<<_trackvec.size()<<std::endl;

	//loop over tracks
	for(unsigned int i =0; i<_trackvec.size(); i++){
		Track* t = _trackvec.at(i);		
		file<<t->getD0()<<" "<<t->getPhi()<<" "<<t->getOmega()<<" "<<t->getZ0()<<" "<<t->getTanLambda()<<" ";
		std::vector<double> xyz = getTrackXYZ(t);
		file<<xyz.at(0)<<" "<<xyz.at(1)<<" "<<xyz.at(2)<<std::endl;
		std::vector<float> cov = t->getCovMatrix();
		int covrow = 0;
		int inc = 2;
		for(int i=0; i<15; i++){
			file<<cov[i]<< " ";
			if(i == covrow){
				file<<std::endl;
				covrow += inc;
				inc++;
				
			}
		}
		std::vector<double> mcvec = findMCTrack(t);
		
		if(mcvec.size() == 0){
			file<< -1 <<" "<<-1<<" "<<-1<<" "<<-1<<" "<<-1<<" "<<-1<<" "<<-1<<" "<<-1<<std::endl;
		}
		else{
			file<<mcvec.at(0)<<" "<<mcvec.at(1)<<" "<<mcvec.at(2)<<" "<<mcvec.at(3)<<" "<<mcvec.at(4)<<" "<<mcvec.at(5)<<" "<<mcvec.at(6)<<" "<<mcvec.at(7)<<std::endl;
		}
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

