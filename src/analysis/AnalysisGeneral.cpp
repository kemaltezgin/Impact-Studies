#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TLegend.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral() : Analysis("AnalysisGeneral"){

	//number of bins 
	int nbin = 100;

	//reconstruction probabilities
	for(size_t i = 0; i < 2; i++){

		m_resProbPOut[i] = 0;
		m_resProbEOut[i] = 0;
		m_resProbGOut[i] = 0;

		m_resProbExcl[i] = 0;
	}

	//set 2D histograms
	m_hXBvsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), "xB vs. Q2", 
		nbin, -4., 0., nbin, 0., 3.);

	//set 1D and 2D histograms
	for(size_t i = 0; i < 2; i++){
			
		std::string title;
	
		if(i == 0){
			title = " (Generated)";
		}else{
			title = " (Full reconstructed)";
		}
	
		m_hxBvsQ2[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 2.2);
		m_hXBvsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 1.6);
		m_hXBvsY[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 0.7);
		m_hQ2vsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 2.2, nbin, 0., 1.6);
		m_hYvsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 0.7, nbin, 0., 1.6);
		m_hQ2vsY[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 2.2, nbin, 0., 0.7);

		m_hXB[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -4., 0.);
		m_hQ2[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 3.);
		m_hT[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 2.);
		m_hPhi[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 2.*TMath::Pi());
		m_hPhiS[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 2.*TMath::Pi());
		m_hY[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",	nbin, -0.1, 1.1);
		m_hEtaEOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -4., 2.);
		m_hEtaPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 10.);
		m_hEtaGOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -10., 4.);
		m_hPPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 50., 110.);
		m_hPPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 1.4);
		m_hPThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 25.);
		m_hPPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hEPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 20.);
		m_hEPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0.,10.);
		m_hEThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 5.);
		m_hEPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hGPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 40.);
		m_hGPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 10.);
		m_hGThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, 0., 5.);
		m_hGPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
	}
}

AnalysisGeneral::~AnalysisGeneral(){
}

void AnalysisGeneral::fill(DVCSEvent& event, double weight){

	//reset weight
	weight = 1.;

	//fill 1D histograms
	for(size_t i = 0; i < 2; i++){

		//type
		KinematicsType::Type kinematicsType = (i == 0)?(KinematicsType::True):(KinematicsType::Observed);

		if (event.getPOut(kinematicsType).E() > 0.) {

			m_resProbPOut[i]++;

			m_hPPOut[i]->Fill(event.getPOut(kinematicsType).P(), weight);
			m_hPPtOut[i]->Fill(event.getPOut(kinematicsType).Pt(), weight);
			m_hPThOut[i]->Fill(event.getPOut(kinematicsType).Theta()*1000., weight);
			m_hPPhOut[i]->Fill(event.getPOut(kinematicsType).Phi(), weight);
			m_hEtaPOut[i]->Fill(event.getEtaPOut(kinematicsType), weight);
		}
		if (event.getEOut(kinematicsType).E() > 0.) { 

			m_resProbEOut[i]++;

			m_hEPOut[i]->Fill(event.getEOut(kinematicsType).P(), weight);
			m_hEPtOut[i]->Fill(event.getEOut(kinematicsType).Pt(), weight);
			m_hEThOut[i]->Fill(event.getEOut(kinematicsType).Theta(), weight);
			m_hEPhOut[i]->Fill(event.getEOut(kinematicsType).Phi(), weight);
			m_hEtaEOut[i]->Fill(event.getEtaEOut(kinematicsType), weight);
		}
		if (event.getGammaOut(kinematicsType).E() > 0.) { 

			m_resProbGOut[i]++;

			m_hGPOut[i]->Fill(event.getGammaOut(kinematicsType).P(), weight);
			m_hGPtOut[i]->Fill(event.getGammaOut(kinematicsType).Pt(), weight);
			m_hGThOut[i]->Fill(event.getGammaOut(kinematicsType).Theta(), weight);
			m_hGPhOut[i]->Fill(event.getGammaOut(kinematicsType).Phi(), weight);
			m_hEtaGOut[i]->Fill(event.getEtaGOut(kinematicsType), weight);
		}

		if (event.getPOut(kinematicsType).E() < 0. || event.getEOut(kinematicsType).E() < 0. || event.getGammaOut(kinematicsType).E() < 0.) continue; 

			m_resProbExcl[i]++;

			m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()), weight);

			m_hXB[i]->Fill(log10(event.getXB(kinematicsType)), weight);
			m_hQ2[i]->Fill(log10(event.getQ2(kinematicsType)), weight);
			m_hT[i]->Fill(fabs(event.getT(kinematicsType)), weight);
			m_hPhi[i]->Fill(event.getPhi(kinematicsType), weight);
			m_hPhiS[i]->Fill(event.getPhiS(kinematicsType), weight);
			m_hY[i]->Fill(event.getY(kinematicsType), weight);

			m_hxBvsQ2[i]->Fill(log10(event.getXB(kinematicsType)), log10(event.getQ2(kinematicsType)), weight);
			m_hXBvsT[i]->Fill(log10(event.getXB(kinematicsType)), fabs(event.getT(kinematicsType)), weight);
			m_hXBvsY[i]->Fill(log10(event.getXB(kinematicsType)), event.getY(kinematicsType), weight);
			m_hQ2vsT[i]->Fill(log10(event.getQ2(kinematicsType)), fabs(event.getT(kinematicsType)), weight);
			m_hYvsT[i]->Fill(event.getY(kinematicsType), fabs(event.getT(kinematicsType)), weight);
			m_hQ2vsY[i]->Fill(log10(event.getQ2(kinematicsType)), event.getY(kinematicsType), weight);
	}
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//print reconstruction probabilities
	std::cout << "info: " << __func__ << ":  p': generated: " << m_resProbPOut[0] << "\treconstructed: " << 
		m_resProbPOut[1] << "\tratio: " << m_resProbPOut[1]/double(m_resProbPOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ":  e': generated: " << m_resProbEOut[0] << "\treconstructed: " << 
		m_resProbEOut[1] << "\tratio: " << m_resProbEOut[1]/double(m_resProbEOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ": gam: generated: " << m_resProbGOut[0] << "\treconstructed: " << 
		m_resProbGOut[1] << "\tratio: " << m_resProbGOut[1]/double(m_resProbGOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ": all: generated: " << m_resProbExcl[0] << "\treconstructed: " << 
		m_resProbExcl[1] << "\tratio: " << m_resProbExcl[1]/double(m_resProbExcl[0]) << std::endl;

	//Clone and make ratio of histograms
	m_hRatio[0] = (TH1*)m_hPPOut[1]->Clone();
	m_hRatio[0]->Divide(m_hPPOut[0]); // PPOut ratio
	m_hRatio[1] = (TH1*)m_hPPtOut[1]->Clone();
	m_hRatio[1]->Divide(m_hPPtOut[0]); // PPtOut ratio
	m_hRatio[2] = (TH1*)m_hPThOut[1]->Clone();
	m_hRatio[2]->Divide(m_hPThOut[0]); // PThOut ratio
	m_hRatio[3] = (TH1*)m_hPPhOut[1]->Clone();
	m_hRatio[3]->Divide(m_hPPhOut[0]); // PPhOut ratio
	m_hRatio[4] = (TH1*)m_hEtaPOut[1]->Clone();
	m_hRatio[4]->Divide(m_hEtaPOut[0]); // EtaPOut ratio

	m_hRatio[5] = (TH1*)m_hEPOut[1]->Clone();
	m_hRatio[5]->Divide(m_hEPOut[0]); // EPOut ratio
	m_hRatio[6] = (TH1*)m_hEPtOut[1]->Clone();
	m_hRatio[6]->Divide(m_hEPtOut[0]); // EPtOut ratio
	m_hRatio[7] = (TH1*)m_hEThOut[1]->Clone();
	m_hRatio[7]->Divide(m_hEThOut[0]); // EThOut ratio
	m_hRatio[8] = (TH1*)m_hEPhOut[1]->Clone();
	m_hRatio[8]->Divide(m_hEPhOut[0]); // EPhOut ratio
	m_hRatio[9] = (TH1*)m_hEtaEOut[1]->Clone();
	m_hRatio[9]->Divide(m_hEtaEOut[0]); // EtaEOut ratio
	
	m_hRatio[10] = (TH1*)m_hGPOut[1]->Clone();
	m_hRatio[10]->Divide(m_hGPOut[0]); // GPOut ratio
	m_hRatio[11] = (TH1*)m_hGPtOut[1]->Clone();
	m_hRatio[11]->Divide(m_hGPtOut[0]); // GPtOut ratio
	m_hRatio[12] = (TH1*)m_hGThOut[1]->Clone();
	m_hRatio[12]->Divide(m_hGThOut[0]); // GThOut ratio
	m_hRatio[13] = (TH1*)m_hGPhOut[1]->Clone();
	m_hRatio[13]->Divide(m_hGPhOut[0]); // GPhOut ratio
	m_hRatio[14] = (TH1*)m_hEtaGOut[1]->Clone();
	m_hRatio[14]->Divide(m_hEtaGOut[0]); // EtaGOut ratio
	
	m_hRatio[15] = (TH1*)m_hXB[1]->Clone();
	m_hRatio[15]->Divide(m_hXB[0]); // XB ratio
	m_hRatio[16] = (TH1*)m_hQ2[1]->Clone();
	m_hRatio[16]->Divide(m_hQ2[0]); // Q2 ratio
	m_hRatio[17] = (TH1*)m_hT[1]->Clone();
	m_hRatio[17]->Divide(m_hT[0]); // |t| ratio
	m_hRatio[18] = (TH1*)m_hY[1]->Clone();
	m_hRatio[18]->Divide(m_hY[0]); // y ratio
	m_hRatio[19] = (TH1*)m_hPhi[1]->Clone();
	m_hRatio[19]->Divide(m_hPhi[0]); // phi ratio
	m_hRatio[20] = (TH1*)m_hPhiS[1]->Clone();
	m_hRatio[20]->Divide(m_hPhiS[0]); // phi_s ratio

	//canvases
	std::vector<TCanvas*> cans;
	
	leg[0] = new TLegend(0.4,0.65,0.75,0.8);
	leg[1] = new TLegend(0.6,0.55,0.88,0.7);
	leg[2] = new TLegend(0.5,0.65,0.85,0.8);
	leg[3] = new TLegend(0.15,0.65,0.5,0.8);
	leg[4] = new TLegend(0.5,0.55,0.85,0.7);
	leg[5] = new TLegend(0.15,0.65,0.5,0.8);

	L[0] = new TLine(-4,1,0,1);
	L[1] = new TLine(0,1,3,1);
	L[2] = new TLine(0,1,2,1);
	L[3] = new TLine(-0.1,1,1.1,1);
	L[4] = new TLine(0,1,2*TMath::Pi(),1);
	L[5] = new TLine(0,1,2*TMath::Pi(),1);
	L[6] = new TLine(50,1,110,1);
	L[7] = new TLine(0,1,20,1);
	L[8] = new TLine(0,1,40,1);
	L[9] = new TLine(0,1,1.4,1);
	L[10] = new TLine(0,1,10,1);
	L[11] = new TLine(0,1,10,1);
	L[12] = new TLine(0,1,25,1);
	L[13] = new TLine(0,1,5,1);
	L[14] = new TLine(0,1,5,1);
	L[15] = new TLine(-TMath::Pi()-0.2,1,TMath::Pi()+0.2,1);
	L[16] = new TLine(-TMath::Pi()-0.2,1,TMath::Pi()+0.2,1);
	L[17] = new TLine(-TMath::Pi()-0.2,1,TMath::Pi()+0.2,1);
	L[18] = new TLine(0,1,10,1);
	L[19] = new TLine(-4,1,2,1);
	L[20] = new TLine(-10,1,4,1);
	//loop over canvases
	for(size_t i = 0; i < 6; i++){

		//add canvas
		cans.push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), "")
		);
		leg[i]->SetFillColor(10);
		leg[i]->SetLineColor(10);
		leg[i]->SetShadowColor(10);
		leg[i]->AddEntry(m_hXB[0],"Generated");
		leg[i]->AddEntry(m_hXB[1],"Full reconstruction");
			
		for(size_t k = 0; k <21; k++){
			m_hRatio[k]->SetStats(0);
		//	m_hRatio[k]->SetTitle("Reconstructed / Generated");
			m_hRatio[k]->SetMinimum(0.2);
			m_hRatio[k]->SetMaximum(1.8);
			m_hRatio[k]->SetLineColor(4);
			m_hRatio[k]->SetLineWidth(1);
			m_hRatio[k]->GetYaxis()->SetTitle("Ratio Rec/Gen");
			m_hRatio[k]->SetTitleSize(0.15,"x");
			m_hRatio[k]->SetTitleSize(0.12,"y");
			m_hRatio[k]->SetLabelSize(0.1,"x");
			m_hRatio[k]->SetLabelSize(0.1,"y");
			m_hRatio[k]->SetTitleOffset(0.8,"x");
			m_hRatio[k]->SetTitleOffset(0.4,"y");
		
			L[k]->SetLineColor(15); 
			L[k]->SetLineStyle(2); 
		}

		if (i == 0){

			cans.back()->cd();
			m_hXBvsQ2->Draw("colz");
			m_hXBvsQ2->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsQ2->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();

		} if (i == 1){

			cans.back()->Divide(3,4);
 			padHisto[0] = (TPad*) cans.back()->cd(1);
			padPull[0] = (TPad*) cans.back()->cd(4);
 			padHisto[1] = (TPad*) cans.back()->cd(2);
			padPull[1] = (TPad*) cans.back()->cd(5);
 			padHisto[2] = (TPad*) cans.back()->cd(3);
			padPull[2] = (TPad*) cans.back()->cd(6);
 			padHisto[3] = (TPad*) cans.back()->cd(7);
			padPull[3] = (TPad*) cans.back()->cd(10);
 			padHisto[4] = (TPad*) cans.back()->cd(8);
			padPull[4] = (TPad*) cans.back()->cd(11);
 			padHisto[5] = (TPad*) cans.back()->cd(9);
			padPull[5] = (TPad*) cans.back()->cd(12);
			for (int k = 0; k < 6; k++){
				padHisto[k]->SetBottomMargin(0.01);
				padPull[k]->SetBottomMargin(0.3);
				padPull[k]->SetTopMargin(0.01);
			}
			double y1 = 0.15; double y2 = 0.5; double y3 = 0.65; 
			double x1 = 0.33; double x2 = 0.66;

			std::string options;
			int color;

			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}

				padHisto[0]->SetPad(0., y3, x1, 1.);
				padPull[0]->SetPad(0., y2, x1, y3);
				padHisto[0]->cd();
				cans.back()->cd(1);
				m_hXB[j]->SetLineColor(color);
				m_hXB[j]->SetTitle("log(X_{B})");
				m_hXB[j]->Draw(options.c_str());
				leg[1]->Draw();
				
				cans.back()->cd(4);
				m_hRatio[15]->GetXaxis()->SetTitle("log(X_{B})");
				m_hRatio[15]->Draw();
				L[0]->Draw("sames");
				
				padHisto[1]->SetPad(x1, y3, x2, 1.);
				padPull[1]->SetPad(x1, y2, x2, y3);
				padHisto[1]->cd();
				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hQ2[j]->SetLineColor(color);
				m_hQ2[j]->SetTitle("log(Q^{2}) [(GeV/c)^{2}]");
				m_hQ2[j]->Draw(options.c_str());
				
				cans.back()->cd(5);
				m_hRatio[16]->GetXaxis()->SetTitle("log(Q^{2}) [(GeV/c)^{2}]");
				m_hRatio[16]->Draw();
				L[1]->Draw("sames");

				padHisto[2]->SetPad(x2, y3, 1., 1.);
				padPull[2]->SetPad(x2, y2, 1., y3);
				padHisto[2]->cd();
				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hT[j]->SetLineColor(color);
				m_hT[j]->SetTitle("|t| [(GeV/c)^{2}]");
				m_hT[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hRatio[17]->GetXaxis()->SetTitle("|t| [(GeV/c)^{2}]");
				m_hRatio[17]->Draw();
				L[2]->Draw("sames");

				padHisto[3]->SetPad(0., y1, x1, y2);
				padPull[3]->SetPad(0., 0, x1, y1);
				padHisto[3]->cd();
				cans.back()->cd(7);
				m_hY[j]->SetLineColor(color);
				m_hY[j]->SetTitle("y");
				m_hY[j]->Draw(options.c_str());
				
				cans.back()->cd(10);
				m_hRatio[18]->GetXaxis()->SetTitle("y");
				m_hRatio[18]->Draw();
				L[3]->Draw("sames");
				
				padHisto[4]->SetPad(x1, y1, x2, y2);
				padPull[4]->SetPad(x1, 0, x2, y1);
				padHisto[4]->cd();
				cans.back()->cd(8);
				m_hPhi[j]->SetLineColor(color);
				m_hPhi[j]->SetTitle("#varphi [rad]");
				m_hPhi[j]->SetMinimum(0.);
				m_hPhi[j]->Draw(options.c_str());
				
				cans.back()->cd(11);
				m_hRatio[19]->GetXaxis()->SetTitle("#varphi [rad]");
				m_hRatio[19]->Draw();
				L[4]->Draw("sames");

				padHisto[5]->SetPad(x2, y1, 1., y2);
				padPull[5]->SetPad(x2, 0, 1., y1);
				padHisto[5]->cd();
				cans.back()->cd(9);
				m_hPhiS[j]->SetLineColor(color);
				m_hPhiS[j]->SetTitle("#varphi_{S} [rad]");
				m_hPhiS[j]->SetMinimum(0.);
				m_hPhiS[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hRatio[20]->GetXaxis()->SetTitle("#varphi_{S} [rad]");
				m_hRatio[20]->Draw();
				L[5]->Draw("sames");
			}
		}if (i == 2){
		
			for(size_t k = 0; k < 2; k++){
			
				m_hxBvsQ2[k]->SetTitleSize(0.06,"x");
				m_hxBvsQ2[k]->SetTitleSize(0.06,"y");
				m_hxBvsQ2[k]->SetLabelSize(0.05,"x");
				m_hxBvsQ2[k]->SetLabelSize(0.05,"y");
				m_hxBvsQ2[k]->SetLabelSize(0.05,"z");
				m_hxBvsQ2[k]->SetTitleOffset(0.7,"x");
				m_hxBvsQ2[k]->SetTitleOffset(0.8,"y");
				
				m_hXBvsT[k]->SetTitleSize(0.06,"x");
				m_hXBvsT[k]->SetTitleSize(0.06,"y");
				m_hXBvsT[k]->SetLabelSize(0.05,"x");
				m_hXBvsT[k]->SetLabelSize(0.05,"y");
				m_hXBvsT[k]->SetLabelSize(0.05,"z");
				m_hXBvsT[k]->SetTitleOffset(0.7,"x");
				m_hXBvsT[k]->SetTitleOffset(0.8,"y");
				
				m_hQ2vsT[k]->SetTitleSize(0.06,"x");
				m_hQ2vsT[k]->SetTitleSize(0.06,"y");
				m_hQ2vsT[k]->SetLabelSize(0.05,"x");
				m_hQ2vsT[k]->SetLabelSize(0.05,"y");
				m_hQ2vsT[k]->SetLabelSize(0.05,"z");
				m_hQ2vsT[k]->SetTitleOffset(0.8,"x");
				m_hQ2vsT[k]->SetTitleOffset(0.8,"y");
				
				m_hYvsT[k]->SetTitleSize(0.06,"x");
				m_hYvsT[k]->SetTitleSize(0.06,"y");
				m_hYvsT[k]->SetLabelSize(0.05,"x");
				m_hYvsT[k]->SetLabelSize(0.05,"y");
				m_hYvsT[k]->SetLabelSize(0.05,"z");
				m_hYvsT[k]->SetTitleOffset(0.8,"x");
				m_hYvsT[k]->SetTitleOffset(0.8,"y");
		
				m_hQ2vsY[k]->SetTitleSize(0.06,"x");
				m_hQ2vsY[k]->SetTitleSize(0.06,"y");
				m_hQ2vsY[k]->SetLabelSize(0.05,"x");
				m_hQ2vsY[k]->SetLabelSize(0.05,"y");
				m_hQ2vsY[k]->SetLabelSize(0.05,"z");
				m_hQ2vsY[k]->SetTitleOffset(0.8,"x");
				m_hQ2vsY[k]->SetTitleOffset(0.8,"y");
		
				m_hXBvsY[k]->SetTitleSize(0.06,"x");
				m_hXBvsY[k]->SetTitleSize(0.06,"y");
				m_hXBvsY[k]->SetLabelSize(0.05,"x");
				m_hXBvsY[k]->SetLabelSize(0.05,"y");
				m_hXBvsY[k]->SetLabelSize(0.05,"z");
				m_hXBvsY[k]->SetTitleOffset(0.7,"x");
				m_hXBvsY[k]->SetTitleOffset(0.8,"y");
		
			}
			
			cans.back()->Divide(4,4);
			
			cans.back()->cd(1);
			m_hxBvsQ2[0]->Draw("colz");
			m_hxBvsQ2[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[0]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();
			
			cans.back()->cd(2);
			m_hxBvsQ2[1]->Draw("colz");
			m_hxBvsQ2[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[1]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();

			cans.back()->cd(3);
			m_hXBvsT[0]->Draw("colz");
			m_hXBvsT[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(4);
			m_hXBvsT[1]->Draw("colz");
			m_hXBvsT[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(5);
			m_hQ2vsT[0]->Draw("colz");
			m_hQ2vsT[0]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(6);
			m_hQ2vsT[1]->Draw("colz");
			m_hQ2vsT[1]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();

			cans.back()->cd(7);
			m_hYvsT[0]->Draw("colz");
			m_hYvsT[0]->GetXaxis()->SetTitle("y");
			m_hYvsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(8);
			m_hYvsT[1]->Draw("colz");
			m_hYvsT[1]->GetXaxis()->SetTitle("y");
			m_hYvsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(9);
			m_hQ2vsY[0]->Draw("colz");
			m_hQ2vsY[0]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsY[0]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(10);
			m_hQ2vsY[1]->Draw("colz");
			m_hQ2vsY[1]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsY[1]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(11);
			m_hXBvsY[0]->Draw("colz");
			m_hXBvsY[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsY[0]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(12);
			m_hXBvsY[1]->Draw("colz");
			m_hXBvsY[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsY[1]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();

		} if (i == 3){

			cans.back()->Divide(3,4);
 			padHisto[6] = (TPad*) cans.back()->cd(1);
			padPull[6] = (TPad*) cans.back()->cd(4);
 			padHisto[7] = (TPad*) cans.back()->cd(2);
			padPull[7] = (TPad*) cans.back()->cd(5);
 			padHisto[8] = (TPad*) cans.back()->cd(3);
			padPull[8] = (TPad*) cans.back()->cd(6);
 			padHisto[9] = (TPad*) cans.back()->cd(7);
			padPull[9] = (TPad*) cans.back()->cd(10);
 			padHisto[10] = (TPad*) cans.back()->cd(8);
			padPull[10] = (TPad*) cans.back()->cd(11);
 			padHisto[11] = (TPad*) cans.back()->cd(9);
			padPull[11] = (TPad*) cans.back()->cd(12);
			for (int k = 6; k < 12; k++){
				padHisto[k]->SetBottomMargin(0.01);
				padPull[k]->SetBottomMargin(0.3);
				padPull[k]->SetTopMargin(0.01);
			}
			double y1 = 0.15; double y2 = 0.5; double y3 = 0.65; 
			double x1 = 0.33; double x2 = 0.66;

			std::string options;
			int color;
	
			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				padHisto[6]->SetPad(0., y3, x1, 1.);
				padPull[6]->SetPad(0., y2, x1, y3);
				padHisto[6]->cd();
				cans.back()->cd(1);
				cans.back()->cd(1)->SetLogy();
				m_hPPOut[j]->SetLineColor(color);
				m_hPPOut[j]->SetTitle("momentum p' [GeV/c]");
				m_hPPOut[j]->Draw(options.c_str());
				leg[3]->Draw();
				
				cans.back()->cd(4);
				m_hRatio[0]->GetXaxis()->SetTitle("p p' [GeV/c]");
				m_hRatio[0]->Draw();
				L[6]->Draw("sames");
			
				padHisto[7]->SetPad(x1, y3, x2, 1.);
				padPull[7]->SetPad(x1, y2, x2, y3);
				padHisto[7]->cd();
				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hEPOut[j]->SetLineColor(color);
				m_hEPOut[j]->SetTitle("momentum e' [GeV/c]");
				m_hEPOut[j]->Draw(options.c_str());
				
				cans.back()->cd(5);
				m_hRatio[5]->GetXaxis()->SetTitle("p e' [GeV/c]");
				m_hRatio[5]->Draw();
				L[7]->Draw("sames");

				padHisto[8]->SetPad(x2, y3, 1., 1.);
				padPull[8]->SetPad(x2, y2, 1., y3);
				padHisto[8]->cd();
				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hGPOut[j]->SetLineColor(color);
				m_hGPOut[j]->SetTitle("momentum #gamma' [GeV/c]");
				m_hGPOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hRatio[10]->GetXaxis()->SetTitle("p #gamma' [GeV/c]");
				m_hRatio[10]->Draw();
				L[8]->Draw("sames");

				padHisto[9]->SetPad(0., y1, x1, y2);
				padPull[9]->SetPad(0., 0, x1, y1);
				padHisto[9]->cd();
				cans.back()->cd(7);
				cans.back()->cd(7)->SetLogy();
				m_hPPtOut[j]->SetLineColor(color);
				m_hPPtOut[j]->SetTitle("transverse momentum p' [GeV/c]");
				m_hPPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(10);
				m_hRatio[1]->GetXaxis()->SetTitle("p_{T} p' [GeV/c]");
				m_hRatio[1]->Draw();
				L[9]->Draw("sames");

				padHisto[10]->SetPad(x1, y1, x2, y2);
				padPull[10]->SetPad(x1, 0, x2, y1);
				padHisto[10]->cd();
				cans.back()->cd(8);
				cans.back()->cd(8)->SetLogy();
				m_hEPtOut[j]->SetLineColor(color);
				m_hEPtOut[j]->SetTitle("transverse momentum e' [GeV/c]");
				m_hEPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(11);
				m_hRatio[6]->GetXaxis()->SetTitle("p_{T} e' [GeV/c]");
				m_hRatio[6]->Draw();
				L[10]->Draw("sames");

				padHisto[11]->SetPad(x2, y1, 1., y2);
				padPull[11]->SetPad(x2, 0, 1., y1);
				padHisto[11]->cd();
				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hGPtOut[j]->SetLineColor(color);
				m_hGPtOut[j]->SetTitle("transverse momentum #gamma' [GeV/c]");
				m_hGPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hRatio[11]->GetXaxis()->SetTitle("p_{T} #gamma' [GeV/c]");
				m_hRatio[11]->Draw();
				L[11]->Draw("sames");

			}
		} if (i == 4){

			cans.back()->Divide(3,4);
 			padHisto[12] = (TPad*) cans.back()->cd(1);
			padPull[12] = (TPad*) cans.back()->cd(4);
 			padHisto[13] = (TPad*) cans.back()->cd(2);
			padPull[13] = (TPad*) cans.back()->cd(5);
 			padHisto[14] = (TPad*) cans.back()->cd(3);
			padPull[14] = (TPad*) cans.back()->cd(6);
 			padHisto[15] = (TPad*) cans.back()->cd(7);
			padPull[15] = (TPad*) cans.back()->cd(10);
 			padHisto[16] = (TPad*) cans.back()->cd(8);
			padPull[16] = (TPad*) cans.back()->cd(11);
 			padHisto[17] = (TPad*) cans.back()->cd(9);
			padPull[17] = (TPad*) cans.back()->cd(12);
			for (int k = 12; k < 18; k++){
				padHisto[k]->SetBottomMargin(0.01);
				padPull[k]->SetBottomMargin(0.3);
				padPull[k]->SetTopMargin(0.01);
			}
			double y1 = 0.15; double y2 = 0.5; double y3 = 0.65; 
			double x1 = 0.33; double x2 = 0.66;

			std::string options;
			int color;
	
			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				padHisto[12]->SetPad(0., y3, x1, 1.);
				padPull[12]->SetPad(0., y2, x1, y3);
				padHisto[12]->cd();
				cans.back()->cd(1);
				cans.back()->cd(1)->SetLogy();
				m_hPThOut[j]->SetLineColor(color);
				m_hPThOut[j]->SetTitle("Polar angle #theta p' [mrad]");
				m_hPThOut[j]->Draw(options.c_str());
				leg[4]->Draw();
				
				cans.back()->cd(4);
				m_hRatio[2]->GetXaxis()->SetTitle("#theta p' [mrad]");
				m_hRatio[2]->Draw();
				L[12]->Draw("sames");

				padHisto[13]->SetPad(x1, y3, x2, 1.);
				padPull[13]->SetPad(x1, y2, x2, y3);
				padHisto[13]->cd();
				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hEThOut[j]->SetLineColor(color);
				m_hEThOut[j]->SetTitle("Polar angle #theta e' [#murad]");
				m_hEThOut[j]->Draw(options.c_str());
				
				cans.back()->cd(5);
				m_hRatio[7]->GetXaxis()->SetTitle("#theta e' [#murad]");
				m_hRatio[7]->Draw();
				L[13]->Draw("sames");

				padHisto[14]->SetPad(x2, y3, 1., 1.);
				padPull[14]->SetPad(x2, y2, 1., y3);
				padHisto[14]->cd();
				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hGThOut[j]->SetLineColor(color);
				m_hGThOut[j]->SetTitle("Polar angle #theta #gamma' [#murad]");
				m_hGThOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hRatio[12]->GetXaxis()->SetTitle("#theta #gamma' [#murad]");
				m_hRatio[12]->Draw();
				L[14]->Draw("sames");

				padHisto[15]->SetPad(0., y1, x1, y2);
				padPull[15]->SetPad(0., 0, x1, y1);
				padHisto[15]->cd();
				cans.back()->cd(7);
				m_hPPhOut[j]->SetLineColor(color);
				m_hPPhOut[j]->SetTitle("Azimuthal angle #phi p' [rad]");
				m_hPPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(10);
				m_hRatio[3]->GetXaxis()->SetTitle("#phi p' [rad]");
				m_hRatio[3]->Draw();
				L[15]->Draw("sames");

				padHisto[16]->SetPad(x1, y1, x2, y2);
				padPull[16]->SetPad(x1, 0, x2, y1);
				padHisto[16]->cd();
				cans.back()->cd(8);
				m_hEPhOut[j]->SetLineColor(color);
				m_hEPhOut[j]->SetTitle("Azimuthal angle #phi e' [rad]");
				m_hEPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(11);
				m_hRatio[8]->GetXaxis()->SetTitle("#phi e' [rad]");
				m_hRatio[8]->Draw();
				L[16]->Draw("sames");

				padHisto[17]->SetPad(x2, y1, 1., y2);
				padPull[17]->SetPad(x2, 0, 1., y1);
				padHisto[17]->cd();
				cans.back()->cd(9);
				m_hGPhOut[j]->SetLineColor(color);
				m_hGPhOut[j]->SetTitle("Azimuthal angle #phi #gamma' [rad]");
				m_hGPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hRatio[13]->GetXaxis()->SetTitle("#phi #gamma' [rad]");
				m_hRatio[13]->Draw();
				L[17]->Draw("sames");

			}
		} if (i == 5){
			double y1 = 0.65; double x1 = 0.33; double x2 = 0.66;
			cans.back()->Divide(3,2);
 			padHisto[18] = (TPad*) cans.back()->cd(1);
			padPull[18] = (TPad*) cans.back()->cd(4);
 			padHisto[19] = (TPad*) cans.back()->cd(2);
			padPull[19] = (TPad*) cans.back()->cd(5);
 			padHisto[20] = (TPad*) cans.back()->cd(3);
			padPull[20] = (TPad*) cans.back()->cd(6);
			for (int k = 18; k < 21; k++){
				padHisto[k]->SetBottomMargin(0.01);
				padPull[k]->SetBottomMargin(0.3);
				padPull[k]->SetTopMargin(0.01);
			}
			
			std::string options;
			int color;

			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				
	
				padHisto[18]->SetPad(0., y1, x1, 1.);
				padPull[18]->SetPad(0., 0.5, x1, y1);
				padHisto[18]->cd();
				cans.back()->cd(1);
				m_hEtaPOut[j]->SetLineColor(color);
				m_hEtaPOut[j]->SetTitle("#eta p'");
				m_hEtaPOut[j]->Draw(options.c_str());
				leg[5]->Draw();
				
				cans.back()->cd(4);
				m_hRatio[4]->GetXaxis()->SetTitle("#eta p'");
				m_hRatio[4]->Draw();
				L[18]->Draw("sames");
				
				padHisto[19]->SetPad(x1, y1, x2, 1.);
				padPull[19]->SetPad(x1, 0.5, x2, y1);
				padHisto[19]->cd();
				cans.back()->cd(2);
				m_hEtaEOut[j]->SetLineColor(color);
				m_hEtaEOut[j]->SetTitle("#eta e'");
				m_hEtaEOut[j]->Draw(options.c_str());
				
				cans.back()->cd(5);
				m_hRatio[9]->GetXaxis()->SetTitle("#eta e'");
				m_hRatio[9]->Draw();
				L[19]->Draw("sames");
				
				padHisto[20]->SetPad(x2, y1, 1., 1.);
				padPull[20]->SetPad(x2, 0.5, 1., y1);
				padHisto[20]->cd();
				cans.back()->cd(3);
				m_hEtaGOut[j]->SetLineColor(color);
				m_hEtaGOut[j]->SetTitle("#eta #gamma'");
				m_hEtaGOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hRatio[14]->GetXaxis()->SetTitle("#eta #gamma'");
				m_hRatio[14]->Draw();
				L[20]->Draw("sames");
			}
		}
	}

	//print
	for(size_t i = 0; i < cans.size(); i++){

		if(cans.size() > 1 && i == 0){
			cans[i]->Print((path+"(").c_str(), "pdf");
		}
		else if(cans.size() > 1 && i == cans.size() - 1){
			cans[i]->Print((path+")").c_str(), "pdf");
		}
		else{
			cans[i]->Print(path.c_str(), "pdf");
		}
	}
}
