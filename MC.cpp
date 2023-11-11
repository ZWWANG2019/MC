#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>
#include <iomanip>

using namespace std;

const int L = 100;                
const int S=L*L;
const int MCS=30000;           
const int MS=20000;             
const double J=1;     
const double J2=J;  
const double R=2;
double E=0;
double T=0.4;	


void InitNeigh(int *ne)
{
	for(int i=0;i<S;i++)
	{
		int l=i-1;
		int r=i+1;
		int u=i-L;
		int d=i+L;
		
		const int bc=i%L;
		if(bc==0)
			l+=L;
		else if(bc+1==L)
			r-=L;
		
		if(u<0)	u+=S;
		else if(d>=S)	d-=S;
		
		cout<<i<<"\t"<<l<<"\t"<<r<<"\t"<<u<<"\t"<<d<<endl;
		 
		*ne++=l;
		*ne++=r;
		*ne++=u;
		*ne++=d;
	}
}

void InitSpin(int *s)
{
	for(int i=0;i<S;i++)
	{
		//*s=0; 
		*s=rand()%3-1;
		if(i%L==0)	cout<<endl;		
		cout<<*s<<" ";
		s++;
	}
}

void InitRF(double *r)
{
	for(int i=0;i<S;i++)
	{
		*r=2*R*rand()/RAND_MAX-R;
		if(i%L==0)	cout<<endl;		
		cout<<setprecision(2)<<*r<<" ";
		r++;
	}
	cout<<"\n"<<endl;
}

//----------------------------------------
void Save(const int *spin)
{
	ofstream out("Spin.dat",ios::app);	
	for(int j=0;j<L;j++)
	{
		for(int i=0;i<L;i++)
			out<<*spin++<<' ';
			out<<"\n";
	}
	out<<"\n"<<endl; 
	out.close();
}

double Ex(const int si, const int sj)
{
	double e=0;
	const int s=si+sj;
	
	 if(s==1||s==-1)
		e=J;
	 else if(si*sj==-1)
	 	e=J2;
	
	 return e;
}


double Polarization(const int *s)
{
 	int p=0;
	for(int i=0;i<S;i++)
   		p+=*s++;
	return 1.0*p/S;
}

double Magneism(const int *s)
{
 	int m=0;
	for(int i=0;i<S;i++)
   		{
   			if(*s==0) m++;
   			s++;
		   }
		
	return 1.0*m/S;
}


void MonteCarlo(int *spin,const double *r,const int *ne)
{
 	int iaccept=0;
 	double sp=0;
 	double sm=0;
	  	
	 for(int mcs=MCS;mcs>0;mcs--)
	 {
  		for(int j=0;j<S;j++)
		{
			const int i=rand()%S;
			int *si=spin+i;
			const int *nei=ne+4*i;
			const int *sl=spin+*nei;
			const int *sr=spin+*(nei+1);			
			const int *su=spin+*(nei+2);			
			const int *sd=spin+*(nei+3);
			
			int snew=0;		
			if(*si==0)	snew=2*(rand()%2)-1;
			else if(*si==1)	snew=rand()%2-1;
			else	snew=rand()%2;
				
			double de=Ex(*si,*sl)+Ex(*si,*sr)+Ex(*si,*su)+Ex(*si,*sd)
					-Ex(snew,*sl)-Ex(snew,*sr)-Ex(snew,*su)-Ex(snew,*sd);
  			de+=*(r+i)*(abs(*si)-abs(snew));
			de-=E*(*si-snew);
			
			if(de>0||(RAND_MAX*exp(de/T)>rand()))
			{
				*si=snew;
             	iaccept++;
     		}
     	}
		if(mcs<MS)
	    {
		   	double p=Polarization(spin);
		   	sp+=p;
		   	double m=Magneism(spin);
		   	sm+=m;
		}	
	}
	
 	ofstream out("data.dat",ios::app); 
 	out <<T<<'\t'<<E<<'\t'<<sp/MS<<'\t'<<sm/MS<<endl;
 	out.close();
}

int main()
{
 	srand((unsigned)time(NULL));  
 	
 	int *s=new int[S];
 	int *ne=new int[4*S];
 	double *r=new double[S];
 	InitNeigh(ne); 	 	
 	InitSpin(s);
 	InitRF(r);
 	
 	for(E=0;E<=0.2;E=E+0.05)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
		Save(s);
	}
		
	for(E=0.25;E>=0.05;E-=0.05)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);
	}
	
		for(E=0;E>=-0.1;E-=0.02)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);
	}
	
		for(E=-0.15;E>=-0.3;E-=0.05)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);
	}
		
	for(E=-0.25;E<=0;E+=0.05)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);	
	}
	
	for(E=0.02;E<=0.1;E+=0.02)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);	
	}
	
	for(E=0.15;E<=0.25;E+=0.05)
	{
		cout<<"Temperature="<<T<<'\t'<<"E="<<E<<'\t'<<"R="<<R<<endl;
		MonteCarlo(s,r,ne);
				Save(s);	
	}
		
	delete []r;
	delete []ne;
	delete []s;
	return 0;
}


