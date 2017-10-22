#ifndef anafunction_H_
#define anafunction_H_
#include <stdlib.h>
#include <math.h>
#define WINDOW 23.5
#define SAMPLE 0.37
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;
/*function generall used when deal with the raw data */
vector<int> adjust (vector<int> &l);
vector<int> CutShock (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l,int number1, int number2, int number3, int number4);
int mapXID (Matrix &l , int row, int col,int chanl);
int mapYID (Matrix &l , int row, int col,int chanl);
vector<double> pulseProcess (vector<int>&l);
/*fuction being used directly for the AB trig event*/
vector<int> pulseAB (vector<int> &l, int a);
/*fuction being used for post anslysis*/
void recalibration (float a[10][10], int rl,int rh,int cl, int ch);
bool pick(int a[10][10],int b[10][10], Matrix &l , int n);
void print2Darray (int a[10][10]);
void printMatrix (Matrix &l);
void printfloat2Darray ( float a[10][10]);
int mapXID (Matrix &l , int channum);
int mapYID (Matrix &l , int channum);
bool pmtsideveto (Matrix &l , int a, int b, double c,int side, int n);
bool cubeveto (Matrix &l , int x, int y, int z, double c,int n);
vector<double> pulsePSDProcess (vector<int>&l,int push);
int timedif (int timebegin, int timeend);
bool twosideveto (Matrix &l , int x, int y, int z,double c, int n);
bool noiseveto (Matrix &l , int x, int y, int n);
int mapXIDtrig (Matrix &l , int ch);
int mapYIDtrig (Matrix &l , int ch);

/* function generall used for both NuLat-AB and NuLat-Muon and NuLat-PSD 
	adjust -> flip the negative pulse into positive so when do the integral, it directly becomes positive energy valuse
				It also readjust the base line by average the point before the pulse, treat that as baseline
	CutShock-> Function being called to cut all the shock signal which will mimic the pulse or confused peak finding mechanism
	maxfine ; minfind -> return the position where peak and valley of pulse located
	sum -> add point of pulse from certain range together, it used for different purpose like adjust baseline, find energy, and deal with psd	
	psd -> return the psd ratio value
	mapXID and mapYID -> return the xid and yid basid on the map of the experiment 
	pulseProcess-> deal with a pulse being cutshocked and adjusted, return a vector that contain meaningful informaiton about energy, peak-height, timing, psd
*/
vector<int> flip (vector<int> &l)
{
	vector<int> flip;
	for (int i=0; i<l.size(); i++)
	{
		flip.push_back(-l[i]); //adjust the baseline and flip the pulse
	}
	return flip;
}

vector<int> adjust (vector<int> &l)
{
	vector<int> adjust;
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	int minpos=minfind(l);
	int averageOfbegin=sum(l,0,30)/30; // treat the average of the first 30 point as the brief base point, another algorithum used before is find the maximum point, but that might come from accident shock and very sharp 
	int threshold=10;
	int leftzeropos=50; // set the inital value of leftzeropos as 50
	int adjustoff=0;
	for (int i=0; i<minpos;i++)
	{
		if (array[minpos-i]>(averageOfbegin-threshold))
		{
			leftzeropos=minpos-i;
			break;
		}		
	}
	int totalcount=leftzeropos;
	if (leftzeropos<5)
		leftzeropos=5;// to avoid condition that a purly noise pulse make the first point as the pulse peak
	adjustoff=sum(l,0,leftzeropos);
	int adjustoffset=int (adjustoff/totalcount);
	for (int i=0; i<l.size(); i++)
	{
		adjust.push_back(-(l[i]-adjustoffset)); //adjust the baseline and flip the pulse
	}
	delete []array;
	return adjust;
}

/*Funciton used to cut the shock of the signal, the shock will win the competition to the real pulse peak and cause the result confusing, the algorithum is search and tag all the sample higher then 50, then compare them to the neighbours, if the difference larger then 50, just rewrite them to their neighbours*/
vector<int> CutShock (vector<int> &l)
{
	vector<int> cutpulse;
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	int threshold=50; // creat a threshold , so all the data cross +50 and -50 will take into account.
	int arraypos=0;
	int lefttest=0;
	int righttest=0;
	vector<int> spik;
	for (int i=0; i<l.size(); i++)
	{
		if (abs(array[i])>1500) // deal with the condition that the begining of the pulse is really low , then a hugh step to real pulse
		array[i]=0;
		if (abs(array[i])>threshold)
		spik.push_back(i);
		
	} // store all the position that pulse higher the threhold
	for (int i=0; i<spik.size(); i++)
	{
		arraypos=spik[i];
		lefttest=abs (array[arraypos]-array[arraypos-1]);
		righttest=abs (array[arraypos]-array[arraypos+1]);	// search all the potential spike, identify them by compare to the neighbour, if both neighbour have 50 different, then identify them as spike
		if (lefttest>50&&righttest>50)     // currently make it cut at 50, so it will not cut the shock less then 50
		array[arraypos]=array[arraypos-1]; //rewrite the spike by make them equal to their neighbour.
	}
	for (int i=0; i<l.size(); i++)
	{
		cutpulse.push_back(array[i]); //rewrite the cutted spike pulse back to the pulse
	}
	delete []array;
	return cutpulse;
}
/*Function to implemented CFD method to find timing information*/
int timeCFD (vector<int>&l)
{
        int* arrayA=new int[l.size()];
        int* arrayB=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                if (i<30)
                arrayA[i]=0;
                else	        
                arrayA[i]=l[i-30];
        }
        for (int i=0; i<l.size(); i++)
        {
                arrayB[i]=-l[i]*0.4;
        }
        int* array=new int[l.size()];
        vector<int> modifiedpulse;
        for (int i=0; i<l.size(); i++)
	{
		modifiedpulse.push_back(arrayA[i]+arrayB[i]); //adjust the baseline and flip the pulse
	}
        for (int i=0; i<l.size(); i++)
        {
                array[i]=arrayA[i]+arrayB[i];
        } 
        //cout << "pulse being reorganized as array" << endl;
        int minpos=minfind(modifiedpulse);
        int maxpos=maxfind(modifiedpulse);
        int FirstNegative=minpos;
        int FirstPositive=maxpos;
        double zeropos=0;
        //cout << minpos << "\t" << maxpos << "\t" << modifiedpulse[minpos] << "\t" <<  modifiedpulse[maxpos] <<endl;
        if (minpos < maxpos && modifiedpulse[minpos]<0 && modifiedpulse[maxpos]>0)
        {        
                for (int i=0; i<abs(maxpos-minpos); i++)
                {
                        if (modifiedpulse[maxpos-i]<0)
                        {
                                FirstNegative=maxpos-i;
                                break;
                        }                
                }
                for (int i=0; i<abs(maxpos-minpos); i++)
                {
                        if (modifiedpulse[minpos+i]>0)
                        {
                                FirstPositive=minpos+i;
                                break;
                        }                
                }
                //zeropos=((modifiedpulse[FirstPositive]-modifiedpulse[FirstNegative])*FirstPositive-(FirstPositive-FirstNegative)*modifiedpulse[FirstPositive])/(modifiedpulse[FirstPositive]-modifiedpulse[FirstNegative]);
                zeropos=(FirstPositive+FirstNegative)/2;
        }
        return zeropos;
                     
}


/*Function used ot find the position where the sample is maximum*/
int maxfind (vector<int>&l)
{
	int* array=new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];
	}
	int tempvalue=0;
	int maxpos=0;
	for (int i=1; i<l.size(); i++)
	{
		if (tempvalue<array[i])
		{
			tempvalue=array[i];
			maxpos=i;
		}
		
	}
	delete []array;
	return maxpos;
}

/*Function used ot find the position where the sample is minmum*/
int minfind(vector<int> &l)
{
	int* array=new int[l.size()];
	int minpos=0;
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];
	}
	int tempdata=array[0];
	for (int i=0; i<l.size(); i++)
	{
		if (tempdata>array[i])
		{
			tempdata=array[i];
			minpos=i;
		}
	}
	delete []array;
	return minpos;
}

/*Function used to add sample together, majorly used as integral of the pulse to find the energy, variable will be the vector hold the pulse sample information, and the window for integration*/
int sum (vector<int>&l, int number1, int number2)
{
	int* array=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                array[i]=l[i];
        }
	int tempsum=0;
	if (number1==number2)
	{
		delete []array;		
		return 0;
	}
	else 
	{
		for (int i=0; i< number2-number1 ; i++)
		{
			tempsum=tempsum+array[number1+i];
		}
		delete []array;
		return tempsum;
	}
}

/*Fuction used to find psd, majorly call the function sum to get result with two different integration window and find the ratio of them*/
double psd (vector<int>&l,int number1,int number2, int number3, int number4)
{
	/*find the ratio between peakpush(number2), peakleftzero(number1), peakrightzero(number3 currently define as the end of pulse)*/	
	int* array=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                array[i]=l[i];
        }
	int peakleft=number1;
	int peakpush=number2;
	int peakright=number3;
        int peakend=number4;
	int totalenergy=sum(l,peakleft,peakright);
	int tailenergy=sum(l,peakpush,peakend);
	double psdana=0.;
	psdana=(double (tailenergy)/double (totalenergy));
	return psdana;	
}

/*Function to map the channel to get related X position*/
int mapXID (Matrix &l , int scrod, int row, int col,int chanl)
{
	int xid=0;
	int chan=scrod*1000+row*100+col*10+chanl;
        //cout << chan << endl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			{
				xid=j;
				break;
			}
		}
	}
	return xid;
}

/*Function to map the channel to get related Y position*/
int mapYID (Matrix &l ,int scrod, int row, int col,int chanl)
{
	int yid=0;
	int chan=scrod*1000+row*100+col*10+chanl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			{
				yid=i;
				break;
			}
		}
	}
	return yid;
}

/*Function to map the trig channel to get related X position, the difference between this and former function is this one read the channel number directly*/
int mapXIDtrig (Matrix &l , int ch)
{
	int xid=0;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (ch==l[i][j])
			{
				xid=j;
			}
		}
	}
	return xid;
}

/*Function to map the trig channel to get related Y position, the difference between this and former function is this one read the channel number directly*/
int mapYIDtrig (Matrix &l , int ch)
{
	int yid=0;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (ch==l[i][j])
			{
				yid=i;
			}
		}
	}
	return yid;
}

/*Function to be used to deal with the pulse, get useful informaiton for futher anlaysis
	First: Call the CutShock and adjust function, this two funciton basically modified the pulse shape, cut the shock, readjust the baseline and flip the pulse
	Second: Call the math function, majorly sum() to get integral, and also find peak information
	Third: store all the information into certain vector named pulseinfo;
		pulseinfo [0] -> Energy
		pulseinfo [1] -> Peak
		pulseinfo [2] -> PSD
		pulseinfo [3] -> peakpos ---realted to timing issue, needs more work to add CFD ability
*/
vector<double> pulseProcess (vector<int>&l)
{
        //cout<< "start " << endl;
	vector<double> pulseinfo;	
	vector<int> cuttedpulse=CutShock(l);					
	vector<int> adjustedpulse=flip(cuttedpulse);
	int *pulsey=new int[adjustedpulse.size()];
	int *xaxis=new int[adjustedpulse.size()];
	int totalenergy=0;		
	int peakamp=0;
	double psdratio=0.;	
	int peakpos=0;	
	int threshold=5;
	int leftzeropos=0;
	int rightzeropos=adjustedpulse.size();
        int pulsetiming=0;	
	for (int j=0; j<adjustedpulse.size(); j++)
	{
		pulsey[j]=adjustedpulse[j];
		xaxis[j]=j;
	}
	peakpos=maxfind(adjustedpulse);
	peakamp=pulsey[peakpos];
	for (int j=0; j<peakpos;j++)
	{
		if (pulsey[peakpos-j]<threshold)
		{
			leftzeropos=peakpos-j;
			break;
		}		
	}
	for (int j=peakpos; j<adjustedpulse.size();j++)
	{
		if (pulsey[j]<threshold)
		{
			rightzeropos=j;
			break;
		}		
	}
        int psdpos=peakpos+95;
	if (psdpos>adjustedpulse.size())
	psdpos=peakpos; 
	totalenergy=sum(adjustedpulse,leftzeropos,rightzeropos);
	psdratio = psd (adjustedpulse,leftzeropos,psdpos,rightzeropos,adjustedpulse.size());
        pulsetiming= timeCFD (adjustedpulse);
	/*event  information storage*/						
	pulseinfo.push_back(totalenergy);
	pulseinfo.push_back(peakamp);
	pulseinfo.push_back(psdratio);
	pulseinfo.push_back(pulsetiming);
	adjustedpulse.clear();
        //cout << "end" << endl;
	delete []pulsey;
	delete []xaxis;	
	return pulseinfo;						
}

/*Function to separate AB events, when variable a=0, A event being filled, otherwise fill B event*/
vector<int> pulseAB (vector<int> &l, int a)
{
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	
	int halfcount=l.size()/2;
	/*when a=0, fill the A event, when a=1, fill the B event, cut half and fill into both A and B*/
	vector<int> pulse;	
	if (a==0)
	{
		for (int i=0; i<halfcount; i++)
		{
			pulse.push_back(array[i]);
			
		}
	}
	else 
	{
		for (int i=halfcount; i<l.size(); i++)
		{
			pulse.push_back(array[i]);
		}	
	}
	delete []array;
	return pulse;
}

/*Function define the time difference between AB event time window*/
int timedif (int timebegin, int timeend)
{
	int timeinterval=0;
	if (timeend>timebegin)
		timeinterval=timeend-timebegin+1;
	else
		timeinterval=512-abs(timeend-timebegin)+1;
	return timeinterval;
}

/*function to be used as post analysis*/

void print2Darray (int a[10][10])
{
	for (int i=0; i< 10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << a[i][j] << endl;
			}
			else
			{
				cout << a[i][j] << "\t" ;
			}
		}
	}
}

void printfloat2Darray ( float a[10][10])
{
	for (int i=0; i< 10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << a[i][j] << endl;
			}
			else
			{
				cout << a[i][j] << "\t" ;
			}
		}
	}
}

void printMatrix (Matrix &l)
{
	for (int i=0; i< l.size(); i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << l[i][j] << endl;
			}
			else
			{
				cout << l[i][j] << "\t" ;
			}
		}
	}
}

void recalibration (float a[10][10], int rl,int rh,int cl, int ch)
{
	float multiply=100;
	for (int i=rl ; i< rh+1 ; i++)
	{
		for (int j=cl; j<ch+1; j++)
		{
			multiply=multiply*a[i][j]/10;
			cout << "multiply=" << multiply << endl;
		}
		
	}
	for (int i=rl ; i< rh+1 ; i++)
	{
		for (int j=cl; j<ch+1; j++)
		{
			a[i][j]=multiply/a[i][j]/1000000000;
		}
	}
}


bool pick(int a[10][10],int b[10][10], Matrix &l , int n)
{
	/*deal with the threshold*/
	bool eventveto=false;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (i<5|| j<5)
			{			
		
				if (l[i+10*n][j]>b[i][j] || l[i+10*n][j] < a[i][j] )
				{			
					//cout << "veto from pick \t" << l[i+10*n][j] << endl;			
					eventveto=true;
				}
			}
		}
	} 
	return eventveto;
}

bool noiseveto (Matrix &l , int x, int y, int n)
{
	bool eventveto=true;
	if (l[y+10*n][x]>1000)
	eventveto=false;
	return eventveto;
	
}

bool cubeveto (Matrix &l , int x, int y, int z, double c, int n)
{
	bool eventveto=false;
	bool sideone=pmtsideveto(l,x,y,c,0,n);
	bool sidetwo=pmtsideveto(l,x,z,c,1,n);
	bool sidethree=pmtsideveto(l,z,y,c,2,n);
	if (sideone || sidethree || sidetwo)	
		eventveto=true;
	return eventveto;
}

bool twosideveto (Matrix &l , int x, int y, int z,double c, int n)
{
	bool eventveto=false;
	bool sideone=pmtsideveto(l,x,y,c,0,n);
	bool sidethree=pmtsideveto(l,z,y,c,2,n);
	if (sideone || sidethree)	
		eventveto=true;
	return eventveto;
}


bool pmtsideveto (Matrix &l , int a, int b, double c, int side, int n)
{
	bool eventveto=false;

	int matrixpickx=0;
	int matrixpicky=0;
	int xid=0;
	int yid=0;
	/*
	side=0 -> x=0-4, y=0-4;
	side=1 -> x=0-4, y=5-9;
	side=2 -> x=5-9, y=0-4; 
	*/
	if (side==0)
	{
		matrixpickx=0;
		matrixpicky=0;
		xid=a;
		yid=b;
	}
	else if (side==1)
	{
		matrixpickx=5;
		matrixpicky=0;
		xid=a;
		yid=b+5;
	}
	else
	{
		matrixpickx=0;
		matrixpicky=5;
		xid=a+5;
		yid=b;
	}
	int trigchan, testchan;	
	//cout << l[yid+n*10][xid] << endl;
	for (int i=matrixpickx; i<matrixpickx+5; i++)
	{
		for (int j=matrixpicky; j<matrixpicky+5; j++)
		{
			
			if (i!=yid || j!=xid)
			{		
				if (l[i+10*n][j]>l[yid+n*10][xid]*c)
				{
					eventveto=true;
					//cout << "veto "<<testchan << "\t" << trigchan<<"\t" << l[yid+n*10][xid] <<endl;
				}
			}
				
		}
	} 
	return eventveto;
}

/*PSD*/

vector<double> pulsePSDProcess (vector<int>&l,int push)
{
	vector<double> pulseinfo;	
	vector<int> cuttedpulse=CutShock(l);					
	vector<int> adjustedpulse=flip(cuttedpulse);
	int *pulsey=new int[adjustedpulse.size()];
	int *xaxis=new int[adjustedpulse.size()];
	int totalenergy=0;		
	int peakamp=0;
	double psdratio=0.;	
	int peakpos=0;
	int threshold=0;
	int leftzeropos=0;
        int rightzeropos=adjustedpulse.size();
	//cout << "pulse being called" << endl;	
	for (int j=0; j<adjustedpulse.size(); j++)
	{
		pulsey[j]=adjustedpulse[j];
		//cout << pulsey[j] ;
		xaxis[j]=j;
	}
	peakpos=maxfind(adjustedpulse);
	int psdpos=peakpos+push;
	if (psdpos>adjustedpulse.size())
	psdpos=adjustedpulse.size()-2; 
	for (int j=0; j<peakpos;j++)
	{
		if (pulsey[peakpos-j]<threshold)
		{
			leftzeropos=peakpos-j;
			break;
		}		
	}
        for (int j=peakpos; j<adjustedpulse.size();j++)
	{
		if (pulsey[j]<threshold)
		{
			rightzeropos=j;
			break;
		}		
	}
	totalenergy=sum(adjustedpulse,leftzeropos,adjustedpulse.size());
	//cout << peakamp << endl;
	psdratio = psd (adjustedpulse,leftzeropos,psdpos, rightzeropos, adjustedpulse.size());
	/*event  information storage*/						
	pulseinfo.push_back(totalenergy);
	pulseinfo.push_back(psdratio);
	//cout << "totalenergy=" << totalenergy << "\t" << "psdratio=" << psdratio << endl;
	adjustedpulse.clear();
	delete []pulsey;
	delete []xaxis;	
	return pulseinfo;						
}

#endif
