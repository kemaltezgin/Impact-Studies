#include "../../include/bin/BinTSlope.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"

BinTSlope::BinTSlope(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	size_t nTBins,
	const std::pair<double, double>& rangeT
) : Bin("BinTSlope"), m_lumiALL(0.), m_lumiBH(0.) {

	//reset
	reset();

	//set ranges
	m_rangeXB = checkRange(rangeXB);
	m_rangeQ2 = checkRange(rangeQ2);
	m_rangeT = checkRange(rangeT);

	//make histograms
 	for(size_t i = 0; i < 4; i++){

 		//labels
		std::stringstream ss;

		if(i == 0) ss << "ALL rec lumi: ";
		if(i == 1) ss << "BH gen: ";
		if(i == 2) ss << "ALL gen: ";
		if(i == 3) ss << "ALL rec: ";

		ss << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
		
		//make histograms
		m_hDistributions.push_back(new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
						nTBins, rangeT.first, rangeT.second));

		//set sumw2
		m_hDistributions.back()->Sumw2();
	}

	//function for fitting
	m_fFit = new TF1((HashManager::getInstance()->getHash()).c_str(), "[0]*exp(-1*[1]*x)", 0., 2.);

	m_fFit->SetParameter(0, m_sumWeights);
	m_fFit->SetParameter(1, 5.);

}

BinTSlope::~BinTSlope(){
}

void BinTSlope::reset(){

	//run for parent class
	Bin::reset();

	//reset ranges
	m_rangeXB = std::make_pair(0., 0.);
	m_rangeQ2 = std::make_pair(0., 0.);
	m_rangeT = std::make_pair(0., 0.);

	//reset sums
	m_sumXB = 0.;
	m_sumQ2 = 0.;
	m_sumT = 0.;

	//reset number of events
	m_nEvents = 0;

	//reset histograms
	m_hDistributions.clear();
	m_hTSlope = nullptr;
	m_hTAcceptance = nullptr;
}

void BinTSlope::fill(DVCSEvent& event, double weight){

	//if BH+INT+DVCS
	if( event.checkSubProcessType(SubProcessType::BH) && 
		event.checkSubProcessType(SubProcessType::INT) && 
		event.checkSubProcessType(SubProcessType::DVCS)
	){

		if(weight > 0.){

			if(event.isReconstructed()) m_hDistributions.at(0)->Fill(-1 * event.getT());
			m_lumiALL += weight;

			//kinematics
			Bin::fill(event, weight);

			m_sumXB += weight * event.getXB();
			m_sumQ2 += weight * event.getQ2();
			m_sumT += weight * event.getT();
		}

		m_hDistributions.at(2)->Fill(-1 * event.getT(KinematicsType::True));
		if(event.isReconstructed()) m_hDistributions.at(3)->Fill(-1 * event.getT());
	}
	
	//if BH
	if( event.checkSubProcessType(SubProcessType::BH) && 
		(! event.checkSubProcessType(SubProcessType::INT)) && 
		(! event.checkSubProcessType(SubProcessType::DVCS))
	){

		m_hDistributions.at(1)->Fill(-1 * event.getT(KinematicsType::True));
		m_lumiBH += weight;
	}
}

void BinTSlope::analyse(){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry
	if(m_nEvents < 2000){
		return;
	}

	//make t-slope histogram
	m_hTSlope = static_cast<TH1*>(m_hDistributions.at(0)->Clone());
	m_hTSlope->SetTitle("DVCS corrected");

	//calculate acceptance
	m_hTAcceptance = static_cast<TH1*>(m_hDistributions.at(3)->Clone());
	m_hTAcceptance->SetTitle("ALL acceptance");
	m_hTAcceptance->Divide(m_hDistributions.at(2));

	m_hTSlope->Divide(m_hTAcceptance);

	//subtract BH
	m_hTSlope->Add(m_hDistributions.at(1), -1 * m_lumiALL/m_lumiBH);

	//fit
	//fit options: 
	//0: do not attempt to draw function

	int statusCode = int(m_hTSlope->Fit(m_fFit, "0B"));
	//int statusCode = int(m_hTSlope->Fit(m_fFit, "I"));

	//store results
	if(m_fitResult){

		delete m_fitResult;
		m_fitResult = nullptr;
	}

	m_fitResult = new FitResult();

	m_fitResult->setStatusCode(statusCode);
	m_fitResult->setNPoints(m_fFit->GetNumberFitPoints());
	m_fitResult->setChi2(m_fFit->GetChisquare());

	for(size_t i = 0; i < m_fFit->GetNpar(); i++){
		m_fitResult->addParameter(std::make_pair(m_fFit->GetParameter(i), m_fFit->GetParError(i)));
	}
}

void BinTSlope::print() const{

	//run for parent class
	Bin::print();

	//print
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (from events): " << getMeanXB() << std::endl;
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (from events): " << getMeanQ2() << std::endl;
}

const std::pair<double, double>& BinTSlope::getRangeXB() const{
	return m_rangeXB;
}

const std::pair<double, double>& BinTSlope::getRangeQ2() const{
	return m_rangeQ2;
}

const std::pair<double, double>& BinTSlope::getRangeT() const{
	return m_rangeT;
}

double BinTSlope::getMeanXB() const{
	return getMean(m_sumXB, m_sumWeights);
}

double BinTSlope::getMeanQ2() const{
	return getMean(m_sumQ2, m_sumWeights);
}

double BinTSlope::getMeanT() const{
	return getMean(m_sumT, m_sumWeights);
}

const std::vector<TH1*>& BinTSlope::getHDistributions() const{
	return m_hDistributions;
}

TH1* BinTSlope::getHTSlope() const{
	return m_hTSlope;
}

TH1* BinTSlope::getHTAcceptance() const{
	return m_hTAcceptance;
}
