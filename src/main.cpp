#include <cmath>
#include <iostream>
#include <cstdlib>
#include <HepMC3/ReaderAscii.h>

#include "../include/analysis/AnalysisGeneral.h"
#include "../include/analysis/AnalysisGeneralRC.h"
#include "../include/analysis/AnalysisALU.h"
#include "../include/analysis/AnalysisTSlope.h"

#ifdef __APPLE__
    #include <filesystem>
    namespace fs = std::filesystem;
#endif

#ifdef __linux__
    #if __GNUC__ < 9
    	    #if __GNUC__ < 6

    			#define __USE_BOOST__

		 		#include <boost/filesystem.hpp>
           	 	namespace fs = boost::filesystem; 
	    #else
            	 #include <experimental/filesystem>
            	 namespace fs = std::experimental::filesystem;
    	    #endif
    #else
            #include <filesystem>
            namespace fs = std::filesystem;
    #endif
#endif

int main(int argc, char* argv[]){

	//check arguments
	if(argc == 1){
		std::cout << __func__ << " error: usage: " << argv[0] << " dir1 dir2 ..." << std::endl;
		exit(0);
	}

	//analysis objects
	AnalysisGeneral analysisGeneral;
	AnalysisGeneralRC analysisGeneralRC;
	AnalysisALU analysisALU;
	AnalysisTSlope analysisTSlope;

	//loop over directories
	for(size_t i = 1; i < argc; i++){

		//check if exists
		if(! fs::exists(fs::path(argv[i]))){
			
			std::cout << __func__ << " warning: directory: " << argv[i] << " does not exist" << std::endl;
			continue;
		}

		//loop over files
		#ifdef __USE_BOOST__
			const fs::recursive_directory_iterator end;
			for(fs::recursive_directory_iterator dirEntry(fs::path(argv[i])); dirEntry != end; dirEntry++){
		#else
			for(const auto& dirEntry : fs::recursive_directory_iterator(fs::path(argv[i]))){
		#endif
		
			//skip directories
			#ifdef __USE_BOOST__
				 if(! fs::is_regular_file(dirEntry->status())) continue;
			#else
				 if(! fs::is_regular_file(dirEntry.status())) continue;
			#endif

			//txt
			#ifdef __USE_BOOST__
				if(dirEntry->path().extension() == ".txt"){
			#else
				if(dirEntry.path().extension() == ".txt"){
			#endif

				//get path
				#ifdef __USE_BOOST__
					std::string path = dirEntry->path().string();
				#else
					std::string path = dirEntry.path().string();
				#endif
				
			
				//print
				std::cout << __func__ << 
					" info: reading: " << path << std::endl;

				//variables
				std::pair<double, double> crossSection;
				size_t nEvents;
				int beamPolarisation;
				bool isRCSample = false;

				//read to collect atributes ===============
				{
					HepMC3::ReaderAscii inputFile(path);

					//to check if RC sample
					size_t lastParticleSize = 0;

					//loop over events
					for(;;){

						//event
	                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
	               
	               		//read
	          			inputFile.read_event(evt);

	          			//if the number of particles is not fixed, we have RC sample
	          			if(evt.particles().size() != 0 && evt.particles().size() != lastParticleSize){

	          				if(lastParticleSize == 0){
	          					lastParticleSize = evt.particles().size();
	          				}else{
	          					isRCSample = true;
	          				}
	          			}

	                	//if reading failed - exit loop
	                	if(inputFile.failed() ) break;
					}

					//run info
					std::shared_ptr<HepMC3::GenRunInfo> runInfo = inputFile.run_info();

					crossSection.first = 
						std::stod((runInfo->attributes().find("integrated_cross_section_value")->second)->unparsed_string());
					crossSection.second = 
						std::stod((runInfo->attributes().find("integrated_cross_section_uncertainty")->second)->unparsed_string());
					nEvents = 
						std::stoul((runInfo->attributes().find("generated_events_number")->second)->unparsed_string());
					beamPolarisation = 
						std::stoi((runInfo->attributes().find("beam_polarisation")->second)->unparsed_string());

					//close
    				inputFile.close();
				}

    			//read to process events ===============
    			{
					HepMC3::ReaderAscii inputFile(path);

					//loop over events
					for(;;){

						//event
	                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
	               
	               		//read
	          			inputFile.read_event(evt);

	                	//if reading failed - exit loop
	                	if(inputFile.failed() ) break;

	                	//DVCS event 
	                	//TODO add beam charge and target polarisation
	               	 	DVCSEvent dvcsEvent(evt, beamPolarisation, -1, TVector3(0., 0., 0.), isRCSample);

	               	 	//fill
	               	 	//TODO add weight 
	               	 	analysisGeneral.fill(dvcsEvent, 1.);
	               	 	analysisGeneralRC.fill(dvcsEvent, 1.);
	               	 	analysisALU.fill(dvcsEvent, 1.);
						analysisTSlope.fill(dvcsEvent, 1.);
					}

					//close file
	    			inputFile.close();
    			}
			}
		}
	}
   
	//analyse
	analysisGeneral.analyse();
	analysisGeneralRC.analyse();
	analysisALU.analyse();
	analysisTSlope.analyse();

	//print
	analysisGeneral.plot("analysisGeneral.pdf");
	analysisGeneralRC.plot("analysisGeneralRC.pdf");
	analysisALU.plot("analysisALU.pdf");
	analysisTSlope.plot("analysisTSlope.pdf");

	return 0;
}
