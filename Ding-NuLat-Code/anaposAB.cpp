#include <iostream>
#include <fstream>
#include <cstdlib> // for exit()
#include <string.h>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "anafunction.h"
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;
typedef vector<double> RowPsd;
typedef vector<RowPsd> MatrixPsd;


int main(int argc, char* argv[])
{
	/*initalizae the variable using to reading different files*/
	int lineMaximum=5000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*initialize the option this analysis needs to do
		
		T: Just for test, will creat 100 event 2D energy map
		C: Threshold cut and then output 2D historgram of energy map and 1D spectrum for each channel
		RC: Recalibration the energy just be reasonable number (needs to modified incase the number out of range)	
	*/
	string option;
	cout << " Choose which command you need : AB pick (AB) ; B event analysis(B); Energy Mapping (m) ; Timing check for AB (T); AB pick out with timing check (ABT) " << endl;
	cin >> option;
	/* open different file ready to be input
		
		energy	psd	peak	timing	cubeID	event
		
	*/

	string analysisname=argv[1];
	ifstream fincondition;
	fincondition.open ("condition.txt");
	ifstream finenergyA;
	finenergyA.open (analysisname+" energyA.txt");
	ifstream finpsdA;
	finpsdA.open (analysisname+" psdA.txt");
	ifstream finpeakA;
	finpeakA.open (analysisname+" peakA.txt");
	ifstream fintimingA;
	fintimingA.open (analysisname+" timingA.txt");
	ifstream finenergyB;
	finenergyB.open (analysisname+" energyB.txt");
	ifstream finpsdB;
	finpsdB.open (analysisname+" psdB.txt");
	ifstream finpeakB;
	finpeakB.open (analysisname+" peakB.txt");
	ifstream fintimingB;
	fintimingB.open (analysisname+" timingB.txt");	
	ifstream fincube;
	fincube.open (analysisname+" cubeID.txt");
	ifstream finevent;
	finevent.open (analysisname+" event.txt");
	ifstream fintimingAB;
	fintimingAB.open (analysisname+" timeAB.txt");
	/* reading event file and create a event vector to be used and event cateory in the future, also deine the number of line needs to be read from other file by counting event number*/	
	vector<int> eventID;
	int n=0; 
	while (finevent.getline(str,lineMaximum))
	{
		p=strtok (str,d);		
		eventID.push_back (atoi(p));
		//cout << p << endl;
		n++;
	}
	/*Fill the map matrix*/
	Matrix map(10,Row(10));
	ifstream finmap;
	finmap.open ("map.txt");
	int a=0;
	int b=0;

	while (finmap.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			map[a][b]=atoi(p);
			b++;
			p=strtok(NULL,d);
		}
		b=0;
		a++;
	}

	/* output a matrix mapping being reading here*/
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			cout << map[i][j]<< "\t" ;
		}
		cout << endl;
	}
	/*create different matrix to store data in the future and also 2D matrix as for calibration and threshold cut*/
	Matrix energyA(n*10,Row(10));	
	MatrixPsd psdA(n*10,RowPsd(10));
	Matrix peakA(n*10,Row(10));
	Matrix timingA(n*10,Row(10));
	Matrix energyB(n*10,Row(10));	
	MatrixPsd psdB(n*10,RowPsd(10));
	Matrix peakB(n*10,Row(10));
	Matrix timingB(n*10,Row(10));	
	Matrix cube(n*10,Row(10));
        Matrix timingAB(n*10,Row(10));
	
	int lowthreshold[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int highthreshold[10][10]={	{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000}};	

	float calibration[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};	
	
	/*Start reading data from different file*/
	int i=0;	
	int j=0;
	while (fincondition.getline(str,lineMaximum))
	{
		
		p=strtok (str,d);
		if (i<11 && i>0)			
		{	while (p)
			{	
				lowthreshold[i-1][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		else if (i>11 && i < 22)			
		{	while (p)
			{	
				highthreshold[i-12][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		else if (i >22 && i<33)			
		{	while (p)
			{	
				calibration[i-23][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		j=0;
		i++;
		//cout << i << endl;
	}
	cout << "finish condition" << endl;

	i=0;
	while (finenergyA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{	
			energyA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpeakA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{			
			peakA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpsdA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			psdA[i][j]=atof(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fintimingA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
		i=0;
	while (finenergyB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{	
			energyB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpeakB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{			
			peakB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpsdB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			psdB[i][j]=atof(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fintimingB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fincube.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			cube[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
       	i=0;
	while (fintimingAB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingAB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
        //printMatrix (timingAB);
	
	if (option=="AB pick out" || option=="AB")
	{
		cout << "function AB pick" << endl;
		int xidA, yidA, zidA, xidB, yidB, zidB;
		string MapOption;
		int facepick;
		TString xAid, yAid, zAid, xBid, yBid, zBid,Aname, Bname;
		Aname="Event-A-";
		Bname="Event-B-";
		cout << "input xid for A event" << endl;
		cin >> xidA;
		cout << "input yid for A event" << endl;
		cin >> yidA;
		cout << "input zid for A event" << endl;
		cin >> zidA;
		cout << "input xid for B event" << endl;
		cin >> xidB;
		cout << "input yid for B event" << endl;
		cin >> yidB;
		cout << "input zid for B event" << endl;
		cin >> zidB;
		cout << "Do you need 2D events mapping? (y/n)" << endl;
		cin >> MapOption;
		xAid.Form("%d",xidA);
		xBid.Form("%d",xidB);
		yAid.Form("%d",yidA);
		yBid.Form("%d",yidB);
		zAid.Form("%d",zidA);
		zBid.Form("%d",zidB);
		TString rootid=Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid;
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
		ofstream foutevent;
		foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");
		ofstream foutdata;
		foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".csv");
		foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
		int eventcount=0;
		/*create category of histogram to hold spectrum*/
		TString energystr="energy ";
		TString peakstr="peak ";
		TString psdstr="psd ";
		TH1D peakhistA[128];
		TH1D integralhistA[128];
		TH1D peakhistB[128];
		TH1D integralhistB[128];
		TH2D psdhistA[128];
		TH2D psdhistB[128];
                TH2D psdhistAB[128];
		TString titleA=" A";
		TString titleB=" B";
		int xid=0;
		int yid=0;
		for (int i=0; i<128; i++)
		{
			TString xidname,yidname;
			xid=i%10;
			yid=i/10;
			xidname.Form ("%d",xid);
			yidname.Form ("%d",yid);
			TString titlename = yidname+" "+xidname;
			TString histname = yidname+" "+xidname;
			peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
			integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
			psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
			peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
			integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
			psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
                        psdhistAB[i]=TH2D(psdstr+histname+titleA+titleB,titlename,1500,0,350000,100,0,1);
		}
		int histTag=0;
		bool vetoA=true;
		bool vetoB=true;
		bool vetoApick=true;
		bool vetoBpick=true;
		for (int i=0; i<n; i++)
		{
			//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
			//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);
			if (true)
			{
				//cout << eventID[i] << "A 3 data"<< endl;
				vetoA=cubeveto(energyA, xidA,yidA,zidA,i);
				//cout << eventID[i] << "B 3 data"<< endl;
				vetoB=cubeveto(energyB, xidB,yidB,zidB,i);
				//cout << vetoA << "\t" << vetoB << endl;
				if (!vetoA && !vetoB)
				{
					foutevent << eventID[i] << endl;
					eventcount++;		
					foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << psdA[yidA+i*10][xidA] << "," <<  psdA[yidA+i*10][zidA+5] << "," << psdA[5+i*10+zidA][xidA] << "," << psdB[yidB+i*10][xidB] << "," <<  psdB[yidB+i*10][zidB+5] << "," << psdB[5+i*10+zidB][xidB]<<timingAB[yidB+i*10][xidB] << "," << timingAB[yidB+i*10][5+zidB] << "," << timingAB[5+zidB+i*10][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," << "\n";
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5||k<5)
							{
								histTag=j*10+k;								
								peakhistA[histTag].Fill(peakA[j+10*i][k]);
								integralhistA[histTag].Fill(energyA[j+10*i][k]);
								psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
								peakhistB[histTag].Fill(peakB[j+10*i][k]);
								integralhistB[histTag].Fill(energyB[j+10*i][k]);
								psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
                                                                psdhistAB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
                                                                psdhistAB[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
							
							}
						}
					}
					if (MapOption=="y")
					{			
						TString eventchar;						
						eventchar.Form ("%d",eventID[i]);	
						TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyA[j+i*10][k];l++)
									{
										temp2dhisA->Fill(k,9-j);						
									}
								}
							}
						}
						TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyB[j+i*10][k];l++)
									{
										temp2dhisB->Fill(k,9-j);						
									}
								}
							}
						}
						//cout << "wrtie the 2D map" << endl;
						temp2dhisA->Write();
						temp2dhisB->Write();
						delete temp2dhisA;
						delete temp2dhisB;
					}
				}
			}
		}
		cout << "total "<< eventcount << "\t events being found" << endl;
		foutdata.clear();
		foutdata.close();
		foutevent.clear();
		foutevent.close();
		f->Write();
		f->Close();
	}
        /*AB pick with a rewrite energy mapping by veto the energy event's that not time realted*/

        if (option=="AB pick out with timing check" || option=="ABT")
	{
		cout << "function AB pick" << endl;
		int xidA, yidA, zidA, xidB, yidB, zidB;
		string MapOption;
		int facepick;
		TString xAid, yAid, zAid, xBid, yBid, zBid,Aname, Bname;
		Aname="Event-A-";
		Bname="Event-B-";
		cout << "input xid for A event" << endl;
		cin >> xidA;
		cout << "input yid for A event" << endl;
		cin >> yidA;
		cout << "input zid for A event" << endl;
		cin >> zidA;
		cout << "input xid for B event" << endl;
		cin >> xidB;
		cout << "input yid for B event" << endl;
		cin >> yidB;
		cout << "input zid for B event" << endl;
		cin >> zidB;
		cout << "Do you need 2D events mapping? (y/n)" << endl;
		cin >> MapOption;
		xAid.Form("%d",xidA);
		xBid.Form("%d",xidB);
		yAid.Form("%d",yidA);
		yBid.Form("%d",yidB);
		zAid.Form("%d",zidA);
		zBid.Form("%d",zidB);
		TString rootid=Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid;
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
		ofstream foutevent;
		foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");
		ofstream foutdata;
		foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".csv");
		foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
		int eventcount=0;
		/*create category of histogram to hold spectrum*/
		TString energystr="energy ";
		TString peakstr="peak ";
		TString psdstr="psd ";
		TH1D peakhistA[128];
		TH1D integralhistA[128];
		TH1D peakhistB[128];
		TH1D integralhistB[128];
		TH2D psdhistA[128];
		TH2D psdhistB[128];
		TString titleA=" A";
		TString titleB=" B";
		int xid=0;
		int yid=0;                 
		int histTag=0;
		for (int i=0; i<128; i++)
		{
			TString xidname,yidname;
			xid=i%10;
			yid=i/10;
			xidname.Form ("%d",xid);
			yidname.Form ("%d",yid);
			TString titlename = yidname+" "+xidname;
			TString histname = yidname+" "+xidname;
			peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
			integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
			psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
			peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
			integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
			psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
		}
                /*Initial timing histogram */
                TH1D* timingA12hist=new TH1D("timing12A","timing12A",60,-30,30);
                TH1D* timingA13hist=new TH1D("timing13A","timing13A",60,-30,30);
                TH1D* timingA23hist=new TH1D("timing23A","timing23A",60,-30,30);
                TH1D* timingB12hist=new TH1D("timing12B","timing12B",60,-30,30);
                TH1D* timingB13hist=new TH1D("timing13B","timing13B",60,-30,30);
                TH1D* timingB23hist=new TH1D("timing23B","timing23B",60,-30,30);
                
                /*Rewrite the energy mapping by veto with the help of correlated timing*/
                
                Matrix energyAtm(n*10,Row(10));	
	        Matrix energyBtm(n*10,Row(10));
                int xtrigid, ytrigid;
                cout << "input xid for trigger event" << endl;
		cin >> xtrigid;
		cout << "input yid for trigger event" << endl;
		cin >> ytrigid;
                for (int i=0; i<n; i++)
                {
                        for (int j=0; j<10; j++)
                        {
                                for (int k=0; k<10; k++)
                                {
                                        if (j!=ytrigid || k!=xtrigid)
                                        {
                                                if ( abs(timingA[j+i*10][k]-timingA[ytrigid+i*10][xtrigid])>20)
                                                energyAtm[j+i*10][k]=0;
                                                else
                                                energyAtm[j+i*10][k]=energyA[j+i*10][k];
                                                if ( abs(timingB[j+i*10][k]-timingB[ytrigid+i*10][xtrigid])>20)
                                                energyBtm[j+i*10][k]=0;
                                                else
                                                energyBtm[j+i*10][k]=energyB[j+i*10][k];
                                        }
                                        else
                                        {
                                                energyAtm[j+i*10][k]=energyA[j+i*10][k];
                                                energyBtm[j+i*10][k]=energyB[j+i*10][k];
                                        }
                                }
                        }
                }            
        
                /*bool variable used for event pick*/
		bool vetoA=true;
		bool vetoB=true;
		bool vetoApick=true;
		bool vetoBpick=true;
		for (int i=0; i<n; i++)
		{
			//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
			//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);
			if (true)
			{
				//cout << eventID[i] << "A 3 data"<< endl;
				vetoA=cubeveto(energyAtm, xidA,yidA,zidA,i);
				//cout << eventID[i] << "B 3 data"<< endl;
				vetoB=cubeveto(energyBtm, xidB,yidB,zidB,i);
				//cout << vetoA << "\t" << vetoB << endl;
				if (!vetoA && !vetoB)
				{
					foutevent << eventID[i] << endl;
					eventcount++;		
                                        /*
foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
                                	*/
					foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," <<"\n";
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5||k<5)
							{
								histTag=j*10+k;								
								peakhistA[histTag].Fill(peakA[j+10*i][k]);
								integralhistA[histTag].Fill(energyA[j+10*i][k]);
								psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
								peakhistB[histTag].Fill(peakB[j+10*i][k]);
								integralhistB[histTag].Fill(energyB[j+10*i][k]);
								psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							
							}
						}
					}
                                        
                                        /*timing histogram filling*/
                                        timingA12hist->Fill(timingA[yidA+i*10][xidA]-timingA[yidA+i*10][5+zidA]); 
                                        timingA13hist->Fill(timingA[yidA+i*10][xidA]-timingA[5+zidA+i*10][xidA]); 
                                        timingA23hist->Fill(timingA[yidA+i*10][5+zidA]-timingA[5+zidA+i*10][xidA]); 
                                        timingB12hist->Fill(timingB[yidB+i*10][xidB]-timingB[yidB+i*10][5+zidB]); 
                                        timingB13hist->Fill(timingB[yidB+i*10][xidB]-timingB[5+zidB+i*10][xidB]); 
                                        timingB23hist->Fill(timingB[yidB+i*10][5+zidB]-timingB[5+zidB+i*10][xidB]);
                                        
        
					if (MapOption=="y")
					{			
						TString eventchar;						
						eventchar.Form ("%d",eventID[i]);	
						TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyAtm[j+i*10][k];l++)
									{
										temp2dhisA->Fill(k,9-j);						
									}
								}
							}
						}
						TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyBtm[j+i*10][k];l++)
									{
										temp2dhisB->Fill(k,9-j);						
									}
								}
							}
						}
						//cout << "wrtie the 2D map" << endl;
						temp2dhisA->Write();
						temp2dhisB->Write();
						delete temp2dhisA;
						delete temp2dhisB;
					}
				}
			}
		}
		cout << "total "<< eventcount << "\t events being found" << endl;


                timingA12hist->Write();
                timingA13hist->Write();
                timingA23hist->Write();
                timingB12hist->Write();
                timingB13hist->Write();
                timingB23hist->Write();

                delete timingA12hist;
                delete timingA13hist;
                delete timingA23hist;
                delete timingB12hist;
                delete timingB13hist;
                delete timingB23hist;
                
                foutdata.clear();
		foutdata.close();
		foutevent.clear();
		foutevent.close();
		f->Write();
		f->Close();
	}	
        

	/*mapping function*/
	if (option=="Event Energy mapping " || option=="m")
	{
		cout << "function mapping pick" << endl;		
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-energy-map"+rootap,"recreate");
		/*ask for event to be plot*/
		vector<int> eventIDplot;
		int eventid=0;
		string eventNew;
		do
		{
			cout<<"input the event you need to map"<<endl;
			cin >> eventid;
			eventIDplot.push_back(eventid);
			cout << "Do you need input another channel to be plot (y/n)" << endl;
			cin >> eventNew;
		} while (eventNew=='y');

		bool eventcheck=false;
		for (int i=0; i<n; i++)
		{	
			eventcheck=false;			
			for (int j=0; j<eventIDplot.size(); j++)
			{
				if (eventID[i]==eventIDplot[j])
					eventcheck=true;
			}
			if (eventcheck)
			{
				
				TString eventchar;						
				eventchar.Form ("%d",eventID[i]);	
				TH2D* temp2dhisA1=new TH2D("eventA1 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisB1=new TH2D("eventB1 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisA2=new TH2D("eventA2 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisB2=new TH2D("eventB2 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisA3=new TH2D("eventA3 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisB3=new TH2D("eventB3 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				for(int j=0;j<5;j++)
				{
					for(int k=0; k<5; k++)
					{
						{
							for (int l=0;l<energyA[j+i*10][k];l++)
							{
								temp2dhisA1->Fill(k,4-j);						
							}
							for (int l=0;l<energyB[j+i*10][k];l++)
							{
								temp2dhisB1->Fill(k,4-j);						
							}
						}
					}
				}
				for(int j=0;j<5;j++)
				{
					for(int k=5; k<10; k++)
					{
						{
							for (int l=0;l<energyA[j+i*10][k];l++)
							{
								temp2dhisA2->Fill(k-5,4-j);						
							}
							for (int l=0;l<energyB[j+i*10][k];l++)
							{
								temp2dhisB2->Fill(k-5,4-j);						
							}
						}
					}
				}
				for(int j=5;j<10;j++)
				{
					for(int k=0; k<5; k++)
					{
						{
							for (int l=0;l<energyA[j+i*10][k];l++)
							{
								temp2dhisA3->Fill(k,9-j);						
							}
							for (int l=0;l<energyB[j+i*10][k];l++)
							{
								temp2dhisB3->Fill(k,9-j);						
							}
						}
					}
				}				
			

				cout << "wrtie the 2D map" << endl;
				temp2dhisA1->Write();
				temp2dhisB1->Write();
				temp2dhisA2->Write();
				temp2dhisB2->Write();
				temp2dhisA3->Write();
				temp2dhisB3->Write();
				delete temp2dhisA1;
				delete temp2dhisB1;
				delete temp2dhisA2;
				delete temp2dhisB2;
				delete temp2dhisA3;
				delete temp2dhisB3;
			
			}
		}
		f->Write();
		f->Close();
	}
	

	if (option=="B event pick" || option=="B")
	{
		cout << "function AB pick" << endl;		
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-B-pickout"+rootap,"recreate");
		ofstream foutevent;
		foutevent.open ("B-Event.txt");
		int eventcount=0;
		/*create category of histogram to hold spectrum*/
		TString rowstr="row ";
		TString colstr="col ";
		TString chanlstr="channel ";
		TString energystr="energy ";
		TString peakstr="peak ";
		TString psdstr="psd ";
		TString eventstr="event";
		TH1D peakhistB[128];
		TH1D integralhistB[128];
		TH2D psdhistB[128];
		TString titleB=" B";
		int histcount=0;
		for (int i=1; i<129; i++)
		{
			TString rowname,colname,chanlname;
			int rowcount=(i-1)/32;
			int colcount=0;
			if (i<33)
			colcount=(i-1)/8;
			else if (i<65)
			colcount=(i-33)/8;
			else if (i<97)
			colcount=(i-65)/8;
			else if (i<129)
			colcount=(i-97)/8;
			int channelcount=(i-1)%8;
			rowname.Form ("%d",rowcount);
			colname.Form ("%d",colcount);
			chanlname.Form ("%d",channelcount);
			TString titlename = rowstr+" "+rowname+" "+colstr+" "+colname+" "+chanlstr+" "+chanlname;
			TString histname = rowname+" "+colname + " "+chanlname;
			peakhistB[i-1]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
			integralhistB[i-1]=TH1D(energystr+histname+titleB,titlename,1500,0,350000); // needs to be modified the size and maximum
			psdhistB[i-1]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
		}
		int xidA, yidA, zidA, xidB, yidB, zidB;
		string MapOption;
		cout << "input xid for B event" << endl;
		cin >> xidB;
		cout << "input yid for B event" << endl;
		cin >> yidB;
		cout << "input zid for B event" << endl;
		cin >> zidB;
		cout << "Do you need 2D events mapping? (y/n)" << endl;
		cin >> MapOption;
		bool vetoB=true;
		bool vetoBpick=true;
		for (int i=0; i<n; i++)
		{
			//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);
			if (true)
			{
				//cout << eventID[i] << "B 3 data"<< endl;
				vetoB=cubeveto(energyB, xidB,yidB,zidB,i);
				if (!vetoB)
				{
					foutevent << eventID[i] << endl;
					eventcount++;					
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5||k<5)
							{
								if (cube[j][k]>0)
								{								
									peakhistB[cube[j][k]-1].Fill(peakB[j+10*i][k]);
									integralhistB[cube[j][k]-1].Fill(energyB[j+10*i][k]);
									psdhistB[cube[j][k]-1].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
									//cout << "fill one channel" << endl;
								}
							}
						}
					}
					if (MapOption=="y")
					{			
						TString eventchar;						
						eventchar.Form ("%d",eventID[i]);	
						TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyB[j+i*10][k];l++)
									{
										temp2dhisB->Fill(k,9-j);						
									}
								}
							}
						}
						cout << "wrtie the 2D map" << endl;
						temp2dhisB->Write();
						delete temp2dhisB;
					}
				}
			}
		}
		cout << "total "<< eventcount << "\t events being found" << endl;
		foutevent.clear();
		foutevent.close();
		f->Write();
		f->Close();
	}

        if (option=="timing of AB" || option=="T")
	{
                cout << "function timing for AB picked" << endl;			
                TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-timing"+rootap,"recreate");
                /**/ 
                ifstream fevent;
                string eventcategory;
                cout << "type the event category needs to be read" << endl;
                cin >> eventcategory;
                fevent.open (eventcategory);
                int lineMaximum=1000;
	        char str[lineMaximum];
	        const char *d=" '\t'";
	        char* p;
                vector<int> eventIDpick;
	        while (fevent.getline(str,lineMaximum))
	        {
		        p=strtok (str,d);		
		        eventIDpick.push_back (atoi(p));
	        }
                for (int i=0; i<eventIDpick.size(); i++)
                {
                        cout << eventIDpick[i] << endl;
                }
                ofstream fout;
                fout.open ("timing.csv");
                fout << "eventID, xid, yid, zid, timingAB-1, timingAB-2, timingAB-3, timingA-1, timingA-2, timingA-3, timingB-1, timingB-2, timingB-3, \n";
                int xid, yid, zid;
                cout << "input xid for timing check" << endl;
		cin >> xid;
		cout << "input yid for timing check" << endl;
		cin >> yid;
                cout << "input yid for timing check" << endl;
		cin >> zid; 
	        bool eventcheck=false;
                TH1D* timingABhist=new TH1D("timing","timing",1000,0,10000);	   
                vector<int> eventnozerocheck;             
                for (int i=0; i<n; i++)
		{
                        eventcheck=false;
                        for (j=0; j<eventIDpick.size(); j++)
                        {
                                if (eventIDpick[j]==eventID[i])
                                eventcheck=true;
                        }			
                        if (eventcheck)
			{
                                cout << eventID[i] << endl;	
                                /*eventnozerocheck.push_back (timingA[yid+i*10][xid]);
                                eventnozerocheck.push_back (timingA[yid+i*10][5+zid]);
                               	eventnozerocheck.push_back (timingA[5+zid+i*10][xid]);
                               	eventnozerocheck.push_back (timingB[yid+i*10][xid]);
                               	eventnozerocheck.push_back (timingB[yid+i*10][5+zid]);
                                eventnozerocheck.push_back (timingB[5+zid+i*10][xid]);*/
                               
                                timingABhist->Fill(timingAB[yid+i*10][xid]);
                                fout <<eventID[i] << "," << xid << "," << yid << "," << zid << ","<< timingAB[yid+i*10][xid]  << "," << timingAB[yid+i*10][5+zid]  << "," << timingAB[5+zid+i*10][xid]  << "," << timingA[yid+i*10][xid] << "," << timingA[yid+i*10][5+zid] << "," << timingA[5+zid+i*10][xid] << "," << timingB[yid+i*10][xid] << "," << timingB[yid+i*10][5+zid] << "," << timingB[5+zid+i*10][xid] << "," <<"\n";	
			}
		}
                timingABhist->Write();
                delete timingABhist;
                fout.clear();
                fout.close();
                fevent.clear();
                fevent.close();
		f->Write();
		f->Close();
	}
	/*create a root file that recreate a root file hold data into different histogram*/
	cout << "good analysis " << n << " event being"<<endl;
	fincondition.clear();
	fincondition.close();
	finenergyA.clear();
	finenergyA.close();
	finpsdA.clear();
	finpsdA.close();
	finpeakA.clear();
	finpeakA.close();
	fintimingA.clear();
	fintimingA.close();
	finenergyB.clear();
	finenergyB.close();
	finpsdB.clear();
	finpsdB.close();
	finpeakB.clear();
	finpeakB.close();
	fintimingB.clear();
	fintimingB.close();
        fintimingAB.clear();
        fintimingAB.close();
	fincube.clear();
	fincube.close();
	finevent.clear();
	finevent.close();
	return 0;
}
