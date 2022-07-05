#include <evt/Event.h>
#include <evt/SimEvent.h>
#include <evt/RecEvent.h>
#include <io/EventFileChain.h>
#include <io/IoCodes.h>
#include <utl/PDGParticleIds.h>

#include <TH2D.h>
#include <TH2I.h>
#include <TH1D.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <vector> 
#include <algorithm>

//#include "fileLoad.h"

using namespace std;
using namespace io;
using namespace utl;
using namespace evt;
using namespace evt::sim;


int main(int argc, char* argv[])
{
	if (argc < 2) {
		cerr << "usage: " << argv[0] << " <filenames>" << endl;
		return 1;
	}

	// --- write histogram of correction factors to ROOT file
	TFile outFile("simulationAnalysis.root", "RECREATE");

	TH1D* wszystko = new TH1D("wszystko", "wszystko", 20, 0, 20);

	wszystko-> GetXaxis () -> SetBinLabel(20, "11C");
 	wszystko-> GetXaxis () -> SetBinLabel(19, "10C");
 	wszystko-> GetXaxis () -> SetBinLabel(18, "9C");
	wszystko-> GetXaxis () -> SetBinLabel(17, "12B");
	wszystko-> GetXaxis () -> SetBinLabel(16, "11B");
	wszystko-> GetXaxis () -> SetBinLabel(15, "10B");
	wszystko-> GetXaxis () -> SetBinLabel(14, "8B");
	wszystko-> GetXaxis () -> SetBinLabel(13, "10Be");
	wszystko-> GetXaxis () -> SetBinLabel(12, "9Be");
	wszystko-> GetXaxis () -> SetBinLabel(11, "7Be");
	wszystko-> GetXaxis () -> SetBinLabel(10, "9Li"); 
	wszystko-> GetXaxis () -> SetBinLabel(9, "8Li"); 
	wszystko-> GetXaxis () -> SetBinLabel(8, "7Li");
	wszystko-> GetXaxis () -> SetBinLabel(7, "6Li");
	wszystko-> GetXaxis () -> SetBinLabel(6, "6He");
	wszystko-> GetXaxis () -> SetBinLabel(5, "4He");
	wszystko-> GetXaxis () -> SetBinLabel(4, "3He");
	wszystko-> GetXaxis () -> SetBinLabel(3, "3H");
	wszystko-> GetXaxis () -> SetBinLabel(2, "2H");
	wszystko-> GetXaxis () -> SetBinLabel(1, "1H");

	// ---  set chain of input files
	const vector<string> fileNames(argv + 1, argv + argc);
	EventFileChain eventFileChain(fileNames);
	
	// --- loop over events
	Event event;

	//vector<int> pidlist=loadPidList(argv[1]);
	vector<int>::iterator it;

	int loop=1;	

	while (eventFileChain.Read(event) == eSuccess) 
	{
	  	loop++;
		//cout<<"------------------------------------"<<endl;
	  	//cout << "loop" << loop << endl;

    	using namespace evt::sim;		
		const SimEvent& simEvent = event.GetSimEvent();

		//int krotnosc = 0;
		//int Ztotal = 0;
		//int Atotal = 0;
		double zinteraction = 0.;		
		double xinteraction = 0.;		
		double yinteraction = 0.;		
		
		for (list<VertexTrack>::const_iterator 
		       simTrackIter = simEvent.Begin<VertexTrack>(), 
		       simTrackEnd = simEvent.End<VertexTrack>(); 
		     simTrackIter != simTrackEnd; ++simTrackIter)
		{
		  
		  const VertexTrack& simTrack = *simTrackIter;        
		  //const Vector& momentum = simTrack.GetMomentum();
		  
		  
		  int pid = simTrack.GetParticleId();
		  int Z = simTrack.GetCharge();			
		  int A = (pid - 1000000000 - Z*10000)/10; 
		  if( pid == 2112){ //neutron
		    //Z = 0;
		    A = 1;
		  }
		  if( pid == 2212){ //proton
		    //Z = 1;
		    A = 1;
		  }
 		  
		  //int hits = 0;
		  
		  //double mommag = momentum.GetMag();
		  //double momperp = momentum.GetPerp();
		  //cout<<momperp<<endl;
		  //double number = 47.02;               
		  
		  double zpositionstop = -1000;				         
		  double xpositionstop = 0;
		  double ypositionstop = 0;

		  if (simTrack.GetStartVertexIndex().IsValid())
		    {
		      const auto& startVertex = simEvent.Get(simTrack.GetStartVertexIndex());
		      const double xpositionstart = startVertex.GetPosition().GetX()/cm;
		      const double ypositionstart = startVertex.GetPosition().GetY()/cm;
		      const double zpositionstart = startVertex.GetPosition().GetZ()/cm;
		      if (simTrack.GetStopVertexIndex().IsValid())
			{
			  const auto& stopVertex = simEvent.Get(simTrack.GetStopVertexIndex());							
			  zpositionstop = stopVertex.GetPosition().GetZ()/cm;
			  xpositionstop = stopVertex.GetPosition().GetX()/cm;
			  ypositionstop = stopVertex.GetPosition().GetY()/cm;
			} 
		      
		      if (!simTrack.GetStopVertexIndex().IsValid())
			{
			  zpositionstop = 1000;
			  xpositionstop = 1000;
			  ypositionstop = 1000;
			}
		      
		      if(zpositionstart == -595) 
		     {
			zinteraction = zpositionstop;
			xinteraction = xpositionstop;
			yinteraction = ypositionstop;
		      }

		    if(zpositionstart == zinteraction && xpositionstart == xinteraction && ypositionstart == yinteraction && A >= 0 && zpositionstart <= -580.03 && zpositionstart >= -582.03 && xpositionstart <= 1.25 && xpositionstart >= -1.25 && ypositionstart <= 1.25 && ypositionstart >= -1.25 && zpositionstop > -580.03)
		      {

				if (simTrack.GetParticleId() == 1000060110)
				{
					wszystko->AddBinContent(20, 1);
				}
				if (simTrack.GetParticleId() == 1000060100)
				{
					wszystko->AddBinContent(19, 1);
				}	
				if (simTrack.GetParticleId() == 1000060090)
				{
					wszystko->AddBinContent(18, 1);
				}						
				if (simTrack.GetParticleId() == 1000050120)
				{
					wszystko->AddBinContent(17, 1);
				}
				if (simTrack.GetParticleId() == 1000050110)
				{
					wszystko->AddBinContent(16, 1);
				}
				if (simTrack.GetParticleId() == 1000050100)
				{
					wszystko->AddBinContent(15, 1);
				}
				if (simTrack.GetParticleId() == 1000050080)
				{
					wszystko->AddBinContent(14, 1);
				}
				if (simTrack.GetParticleId() == 1000040100)
				{
					wszystko->AddBinContent(13, 1);
				}
				if (simTrack.GetParticleId() == 1000040090)
				{
					wszystko->AddBinContent(12, 1);
				}
				if (simTrack.GetParticleId() == 1000040070)
				{
					wszystko->AddBinContent(11, 1);
				}
				if (simTrack.GetParticleId() == 1000030090)
				{
					wszystko->AddBinContent(10, 1);
				}
				if (simTrack.GetParticleId() == 1000030080)
				{
					wszystko->AddBinContent(9, 1);
				}
				if (simTrack.GetParticleId() == 1000030070)
				{
					wszystko->AddBinContent(8, 1);
				}
				if (simTrack.GetParticleId() == 1000030060)
				{
					wszystko->AddBinContent(7, 1);
				}
				if (simTrack.GetParticleId() == 1000020060)
				{
					wszystko->AddBinContent(6, 1);
				}	
				if (simTrack.GetParticleId() == 1000020040)
				{
					wszystko->AddBinContent(5, 1);
				}	
				if (simTrack.GetParticleId() == 1000020030)
				{
					wszystko->AddBinContent(4, 1);
				}	
				if (simTrack.GetParticleId() == 1000010030)
				{
					wszystko->AddBinContent(3, 1);
				}	
				if (simTrack.GetParticleId() == 1000010020)
				{
					wszystko->AddBinContent(2, 1);
				}	
				if (simTrack.GetParticleId() == 2212)
				{
					wszystko->AddBinContent(1, 1);
				}
		      }
		      						
		    }						
		}
   	}

	outFile.Write();
	outFile.Close();

	return 0;
}
