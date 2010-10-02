#ifndef BAYESPSI3UTR_H_
#define BAYESPSI3UTR_H_

#include <cln/cln.h>
#include <cln/lfloat.h>
#include <iostream>
#include <iomanip>
using namespace cln;
using namespace std;

#ifndef  GETPRECISEFLOAT_H_
#define GETPRECISEFLOAT_H_

inline cl_R getPreciseFloat(float _f=0,int _prec=30)
{
	return cl_float(_f,float_format(_prec));
}

#endif

inline cl_I getBase(unsigned int Ne)
{
	return factorial(Ne);
}

inline cl_R getLTerm(float inside)
{
	return exp(getPreciseFloat(-1*inside));
}

inline cl_R getRTerm(float inside,cl_I Ne)
{	
	return expt(getPreciseFloat(inside),Ne);
}

inline float getNaivePsi_3UTR(unsigned int extension_length,unsigned int extension_count,unsigned int common_length,unsigned int common_count)
{
	float extension_density=float(extension_count)/float(extension_length);
	float common_density=float(common_count)/float(common_length);
	return (extension_density/(common_density));
}

inline float getBayesPsi_3UTR(unsigned int extension_length,unsigned int extension_count,unsigned int common_length,unsigned int common_count)
{
	//float extension_density=float(extension_count)/float(extension_length);
	//cerr<<"b1"<<endl;
	float common_density=float(common_count)/float(common_length);
	
	if(common_length==0 || extension_length==0)
	{
		return -10.0;
	}
	if(common_count==0)
	{
		return -15.0;
	}
	//now simulate psi= 0 to 1
	//cerr<<"b2"<<endl;
	
	cl_I base=getBase(extension_count);
	//cerr<<"base="<<base<<endl;
	//cerr<<"b3"<<endl;
	cl_R integration=getPreciseFloat();
	cl_R sump=getPreciseFloat();
	int times=100;
	float inc=float(1)/times;
	//cerr<<"b4"<<endl;
	for(float psi=0;psi<=1;psi+=inc)
	{
		//cerr<<"b5a"<<endl;
		float inside=float(extension_length)*float(psi)*common_density;
	//	cerr<<"b5b"<<endl;
		cl_R lTerm=(getLTerm(inside));
	//	cerr<<"b5c"<<endl;
		cl_R rTerm=(getRTerm(inside,extension_count));
	//	cerr<<"b5d"<<endl;
		//probability over the interval inc
		cl_R p=lTerm*(rTerm/base)*getPreciseFloat(inc);  //inc is there to integrate a portion of the graph covered by the increment
	//	cerr<<"b5d"<<endl;
		cl_R x=getPreciseFloat(psi)*p;
	//	cerr<<"b5e"<<endl;
		//cerr<<setprecision(3)<<psi<<"\t"<<p<<endl;
		
		sump+=p;
		integration+=x;
		
	}
	
	//cerr<<"b5f"<<endl;
	if(sump==0)
	{
		return -15.0;
	}
	
	return float_approx(integration/sump);
}


#endif /*BAYESPSI3UTR_H_*/
