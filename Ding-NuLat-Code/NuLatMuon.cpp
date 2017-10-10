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

int main(int argc, char* argv[])
{
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",argv[argc-1]);
	TFile* f=new TFile (analysisfile+rootap,"recreate"); //summary root file being created, every time code being runned, file being recreated.
	cout << "root file being created here" << endl;
	/* The output txt file will be listed as follows
	 * muon event-> peakA.txt, energyA.txt, psdA.txt, timingA.txt
	 * summary-> event.txt, summary.txt, cubeID.txt
	*/
	string analysisname=argv[argc-1];	
	ofstream peakfileMuon;
	peakfileMuon.open (analysisname+" peakMuon.txt",fstream::app);
	ofstream energyfileMuon;
	energyfileMuon.open (analysisname+" energyMuon.txt",fstream::app);
	ofstream psdfileMuon;
	psdfileMuon.open (analysisname+" psdMuon.txt",fstream::app);
	ofstream timingfileMuon;
	timingfileMuon.open (analysisname+" timingMuon.txt",fstream::app);
	ofstream cubeIDfile;
	cubeIDfile.open (analysisname+" cubeID.txt",fstream::app);	
	ofstream eventfile;
	eventfile.open (analysisname+" event.txt",fstream::app);	
	ofstream anasummary;
	anasummary.open (analysisname+" summary.txt",fstream::app);
	
	/* variable to hold content for each line from data*/
	int lineMaximum=1000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*holding event information for the current line being analyzed*/
	int countline=0;
	int chan=0;
	int tempdatacount=0;
        int eventscrodnumber=0;
        int eventscrod=0;
	int eventnumber=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	int eventfirst=0; // variable for initialize the analysis for the first event	 
	int tempcondition[6]={0}; 	// {eventnumber rolnum colnum channelnum}
	int totalenergysumMuon=0;
	/* Tstring name and 64 histogram being created, right now the histogram is PSD,Energy by both peak and Integral, might get rid off peak method*/
	TString rowstr="row ";	
	TString colstr="col ";
	TString chanlstr="channel ";
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TString eventstr="event";
	TH1D peakhistMuon[128];
	TH1D integralhistMuon[128];
	TH2D psdhistMuon[128];
	TString title=" muon";
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
		peakhistMuon[i-1]=TH1D(peakstr+histname+title,titlename,2000,0,2000);	
		integralhistMuon[i-1]=TH1D(energystr+histname+title,titlename,1500,0,500000); // needs to be modified the size and maximum
		psdhistMuon[i-1]=TH2D(psdstr+histname+title,titlename,1500,0,500000,100,0,1);
		}
	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;	
	vector<double> pulseinfoMuon;
	vector<int> event;
	vector<int> energyspecMuon;
	vector<int> energypeakMuon;
	vector<double> psdanalysisMuon;
	vector<int> timingMuon;
	/*Geometry information being prcessed current version
				x=0	x=1	x=2	x=3	x=4	x=5	x=6	x=7	x=8	x=9			
		y=0		130	131	132	133	134	030	031	032	033	034
		y=1		135	137	124	125	126	035	037/106	020	021	022
		y=2		127	110	111	112	113	023	024	025	026	027			
		y=3		114	115	116	117	100	010	011	012	013	014
		y=4		101	102	103	104	105	015	016	017	106?	107?
		y=5		
		y=6		
		y=7										107	037	
		y=8		000	001	002	003	004	005	006	007	136	036	
		y=9		120	121	122	123	sum+				extra	eventID
		031 033 035 037 021 023 are the channel 1/10 attenuation 
		007 is the channel being trigged		
	*/
	
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
	vector<int> scrod;
	vector<int> row;
	vector<int> col;
	vector<int> channel;
	vector<int> cubeID;
	bool pulsecheck=false;
	int xID[100]={0};
	int yID[100]={0};
	cout << "Matrix initialize here" << endl;
	/*Matrix initialization 
		EnergyPeakMatrix	EnergyInteMatrix	CubeMatrix	PsdMatrix
	*/
	int TimingMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyPeakMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyInteMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};

	double PsdMatrixMuon[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}}; 
	int CubeMatrix[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};


	
	/* reading data file, the variable input when run the program has the struture : file1 , file2, file3,.... , filelast, output-file-name*/
	for (int i=1; i<argc-1; i++)
	{
		fin.open (argv[i]);
		countline=0; // initialize countline=0; data file is read line by line, use this number to give reading process proper tag.
		if (!fin.is_open())
		{
			cerr << "Could not open "<< argv[i]<< endl;
			fin.clear();
		}
		cout << argv[i] << "\t file being analysised here" << endl; // Output the file being analysised at current run
		while (fin.getline(str,lineMaximum))
		{	
			/* line by line processing */
			if (countline<3)
			cout << "note output"<<endl;	// the first three line just file information and structure
			else
			{
				p=strtok (str,d);
				while (p)
				{
					/* word by word frome one line processing, the first 8 number is the inforamtion belong to basic hardware information*/
					if (tempdatacount==1)
					{
						tempcondition[0]=atoi(p);
						if (countline==3)
						{
							eventnumber=tempcondition[0];
							eventfirst=tempcondition[0];
							cout << "event number being initialized here with value \t" << eventfirst << endl;
						}
						
					}
                                        else if (tempdatacount==2)
					{
						tempcondition[5]=atoi(p);
						if (countline==3)
						{
							eventscrodnumber=tempcondition[5];
							cout << "event row being initialized here with value \t" << eventrow << endl;
						}
					}
					else if (tempdatacount==3)
					{
						tempcondition[1]=atoi(p);
						if (countline==3)
						{
							eventrow=tempcondition[1];
							cout << "event row being initialized here with value \t" << eventrow << endl;
						}
					}
					else if (tempdatacount==4)
                                        {
						tempcondition[2]=atoi(p);
						if (countline==3)
						{
							eventcol=tempcondition[2];
							cout << "event col being initialized here with value \t" << eventcol << endl;
						}
                                        }
					else if (tempdatacount==5)
                                        {
						tempcondition[3]=atoi(p);
						if (countline==3)
						{
							eventchanl=tempcondition[3];
							cout << "event channel being initialized here with value \t" << eventchanl << endl;
						}
                                        }
					else if (tempdatacount==8)
					{
						/*After all the pre-process of the data, start deal with the issue from last event*/						
						/* if event number is not the same, means new event starts here, at this point, analysis all the information analysised before for the previous event.*/
											
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl  || tempcondition[5]!=eventscrod) 
						{
							chan=eventscrod*1000+eventrow*100+eventcol*10+eventchanl;
							pulsecheck=false;
						/*loop being calledl to check whether the channel is match with map, pick up all the signal which is matched the map and only doing analysis for them*/
							for (int i=0; i<10; i++)
							{
								for (int j=0; j<10; j++)
								{
									if (chan==map[i][j])
									{
										pulsecheck=true;
									}
								}
							}							
							if (pulsecheck)
							{							
								pulseinfoMuon=pulseProcess(pulse);
								histcount=eventchanl+eventcol*8+eventrow*32+1;
								cubeID.push_back(histcount);
								/*pulseinfo hase the structure that
								 * 	pulseinfo.[0](totalenergy);
									pulseinfo.[1](peakamp);
									pulseinfo.[2](psdratio);
									pulseinfo.[3](peakpos);	 // for the peakpos as timing, we need to redo this use CFD */
								energyspecMuon.push_back(pulseinfoMuon[0]);
								energypeakMuon.push_back(pulseinfoMuon[1]);
								psdanalysisMuon.push_back(pulseinfoMuon[2]);
								timingMuon.push_back(pulseinfoMuon[3]);
								event.push_back(eventnumber);
                                                                scrod.push_back(eventscrod);
								row.push_back(eventrow);
								col.push_back(eventcol);
								channel.push_back(eventchanl);	
							}
							if (tempcondition[0]!=eventnumber)
							{				
								//cout << "event=" << eventnumber << " being processed" << endl;
								for (int j=0; j<cubeID.size(); j++) // create a proper mapping for different channel to a 2D Matrix
								{								
									/*call the function to find the position certain channel should be located*/
									xID[j]=	mapXID (map , scrod[j], row[j], col[j], channel[j]);
									yID[j]= mapYID (map , scrod[j], row[j], col[j], channel[j]);	
									// write certain information into proper matrix
									EnergyPeakMatrixMuon[yID[j]][xID[j]]=energypeakMuon[j];
									EnergyInteMatrixMuon[yID[j]][xID[j]]=energyspecMuon[j];
									PsdMatrixMuon[yID[j]][xID[j]] = psdanalysisMuon[j];	
									TimingMatrixMuon[yID[j]][xID[j]]=timingMuon[j];
									CubeMatrix[yID[j]][xID[j]]=cubeID[j]; // might just get rid of this matrix and the txt file and just use the map
									if (xID[j]<5 && yID[j]<5)									
									{
										totalenergysumMuon+=energyspecMuon[j]; // totalenergy currently add whole blue face which is from 04-04
									}	
														
								}
								/*put sum into matrix [8][9], leave [9][9] for event number*/
								EnergyInteMatrixMuon[8][9]=totalenergysumMuon;
								CubeMatrix[8][9]=128; // currently overlap with 000, assign another number when fully functional board
								/*write result into a preanalysis root file so we can check briefly how the data looks like*/
								for(int j=0;j<10;j++)
								{
									for(int k=0; k<10; k++)								
									{
										if (CubeMatrix[j][k]>0 && CubeMatrix[j][k]<129)
										{											
											peakhistMuon[CubeMatrix[j][k]-1].Fill(EnergyPeakMatrixMuon[j][k]);
											integralhistMuon[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixMuon[j][k]);
											psdhistMuon[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixMuon[j][k],PsdMatrixMuon[j][k]);
										}
																				
									}
								}
								/*writing the pre-analysis result to txt file for future analysis*/
								for(int j=0;j<10;j++)
								{
									for(int k=0; k<10; k++)
									{			
										//cout << k ;					
										if ( j==9 && k==9)
										{
											/*event number being record every matrix at [9][9], currently it is easier for people to check*/
											peakfileMuon << event[0]  << "\t";
											energyfileMuon << event[0] << "\t";
											psdfileMuon << event[0] << "\t";	
											timingfileMuon << event[0] << "\t";
											cubeIDfile << event[0] << "\t" ;
											eventfile << event[0] << endl;	
										}
										else
										{
											//cout << j << k << endl;	
											peakfileMuon << EnergyPeakMatrixMuon[j][k] << "\t" ;
											energyfileMuon << 	EnergyInteMatrixMuon[j][k] << "\t";
											psdfileMuon << PsdMatrixMuon[j][k] << "\t";
											timingfileMuon << TimingMatrixMuon[j][k] << "\t";
											cubeIDfile << CubeMatrix[j][k] << "\t";
										}	
																				
									}
									energyfileMuon << endl;
									peakfileMuon << endl;
									psdfileMuon << endl;	
									timingfileMuon<<endl;
									cubeIDfile << endl;										
								}
								/*clear all the information from last event*/
								for (int j=0; j<10; j++)
								{
									for (int k=0; k<10; k++)
									{									
										EnergyPeakMatrixMuon[j][k]=0;
										EnergyInteMatrixMuon[j][k]=0;
										PsdMatrixMuon[j][k] = 0;
										TimingMatrixMuon[j][k]=0;
										CubeMatrix[j][k]=0;
									}
								}
								totalenergysumMuon=0;
								cubeID.clear();
								energyspecMuon.clear();
								energypeakMuon.clear();
								psdanalysisMuon.clear();
								timingMuon.clear();
								event.clear();
                                                                scrod.clear();
								row.clear();
								col.clear();
								channel.clear();
								//cout << "event process finished" << endl;
							} // end of last event processing
							pulse.clear();
							/*new event condition being initialized here*/
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];
                                                        eventscrod=tempcondition[5];		
							
						}
					}
					else if (tempdatacount>8)
					pulse.push_back(atoi(p));
					p=strtok(NULL,d);
					tempdatacount++;
				}
				tempdatacount=0;
			}
			countline++;
		}
		anasummary << "analysis file is=" << argv[i] << "\t first event is=" << eventfirst << "\t last event is=" << eventnumber << endl;
		fin.clear();
		fin.close();
	}
	anasummary.clear();	
	anasummary.close();
	cubeIDfile.clear();
	cubeIDfile.close();
	timingfileMuon.clear();
	timingfileMuon.close();
	peakfileMuon.clear();
	peakfileMuon.close();
	energyfileMuon.clear();
	energyfileMuon.close();
	psdfileMuon.clear();
	psdfileMuon.close();
	eventfile.clear();
	eventfile.close();
	f->Write();
	f->Close();
	return 0;
}
