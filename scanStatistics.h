#ifndef SCANSTATISTICS_H_
#define SCANSTATISTICS_H_

#include <stlib/src/numerical/random/poisson/PoissonPdf.h>
#include <stlib/src/numerical/random/poisson/PoissonCdf.h>
#include <math.h>
#include <cln/cln.h>
#include <cln/lfloat.h>
using namespace cln;


#ifndef  GETPRECISEFLOAT_H_
#define GETPRECISEFLOAT_H_

inline cl_R getPreciseFloat(float _f=0,int _prec=30)
{
	return cl_float(_f,float_format(_prec));
}

#endif

numerical::PoissonPdf<> Poisson_PDF;
numerical::PoissonCdf<> Poisson_CDF;

inline cl_R Poisson_Prob(unsigned int k,const cl_R& psi)
{
	//cl_R psi=getPreciseFloat(_psi);
	cl_R prob=getPreciseFloat(0.0);
	prob=Poisson_PDF(float_approx(psi),k);
	return prob;
}

inline cl_R __SS_p(unsigned int _k, const cl_R& _psi)
{	
	
	return Poisson_Prob(_k,_psi);
}

inline cl_R SS_Poisson_Prob_Approx3_3(unsigned int _k, double _psi, unsigned int _L)
{
	cl_R psi=getPreciseFloat(_psi);
	cl_R L=getPreciseFloat(_L);
	cl_R delta=(_k-1)*(L-1)*__SS_p(_k,psi);
	return 1-exp(-delta);
}

inline cl_R __SS_Fp( int k,const cl_R& psi)
{
	
	cl_R sum=getPreciseFloat(0);
	if(k<0)
	{
		return sum;
	}
	
	sum=Poisson_CDF(float_approx(psi),k);
	
	return sum;
}

inline cl_R __SS_Q2(unsigned int k, const cl_R& psi)
{
	return expt(__SS_Fp(k-1,psi),2)-(k-1)*__SS_p(k,psi)*__SS_p(k-2,psi)-(k-1-psi)*__SS_p(k,psi)*__SS_Fp(k-3,psi);
	
}

inline cl_R __SS_A1(unsigned int k, const cl_R& psi)
{
	////cerr<<"A1"<<endl;
	return 2*__SS_p(k,psi)*__SS_Fp(k-1,psi)*((k-1)*__SS_Fp(k-2,psi)-psi*__SS_Fp(k-3,psi));
}

inline cl_R __SS_A2(unsigned int k, const cl_R& psi)
{
	////cerr<<"A2"<<endl;
	return 0.5*expt(__SS_p(k,psi),2)*((k-1)*(k-2)*__SS_Fp(k-3,psi)-2*(k-2)*psi*__SS_Fp(k-4,psi)+expt(psi,2)*__SS_Fp(k-5,psi));
}

inline cl_R __SS_A3(unsigned int k, const cl_R& psi)
{
	////cerr<<"A3"<<endl;
	cl_R sum=getPreciseFloat(0);
	for (unsigned int r=1;r<=k-1;r++)
	{
		sum+=__SS_p(2*k-r,psi)*expt(__SS_Fp(r-1,psi),2);
	}
	////cerr<<"A3"<<endl;
	return sum;
}


inline cl_R __SS_A4(unsigned int k, const cl_R& psi)
{
	////cerr<<"A4"<<endl;
	cl_R sum=getPreciseFloat(0);
	for ( unsigned int r=2;r<=k-1;r++) ///Changed unsigned-signed problem
	{
		////cerr<<"A4"<<r<<(k-1)<<endl;
		sum+=__SS_p(2*k-r,psi)*__SS_p(r,psi)*((r-1)*__SS_Fp(r-2,psi)-psi*__SS_Fp(r-3,psi));
	}
	//////cerr<<"A4"<<endl;
	return sum;
}

inline cl_R __SS_Q3(unsigned int k, const cl_R& psi)
{
	return expt(__SS_Fp(k-1,psi),3)-__SS_A1(k,psi)+__SS_A2(k,psi)+__SS_A3(k,psi)-__SS_A4(k,psi);
}

inline cl_R SS_Poisson_Prob_Naus_Approx(unsigned int _k, double _psi, unsigned int _L)
{
	cl_R psi=getPreciseFloat(_psi);
	
	cl_R Q2=__SS_Q2(_k,psi);
	////cerr<<Q2<<endl;
	cl_R Q3=__SS_Q3(_k,psi);
	////cerr<<Q3<<endl;
	cl_R allQ=Q2*expt(Q3/Q2,_L-2);
	////cerr<<allQ<<endl;
	return 1-allQ;
}

inline cl_R SS_Poisson_Prob_Naus_Approx(unsigned int _k, cl_R _psi, cl_I _L)
{
	cl_R psi=_psi;
	
	cl_R Q2=__SS_Q2(_k,psi);
	////cerr<<Q2<<endl;
	cl_R Q3=__SS_Q3(_k,psi);
	////cerr<<Q3<<endl;
	cl_R allQ=Q2*expt(Q3/Q2,_L-2);
	////cerr<<allQ<<endl;
	return 1-allQ;
}


inline double ___SS_Poisson_Prob_Alm_Approx(unsigned int _k, cl_R _lambda, cl_R _w,cl_R _T)
{
	double psi=float_approx(_lambda*_w);
//	cerr<<"a"<<endl;
	cl_R a=_lambda*(_T-_w);
//	cerr<<"b"<<endl;
	cl_R b=-((_k-psi)/_k);
//	cerr<<"c "<<_k<<" "<<psi<<endl;
	cl_R c=(__SS_p(_k-1,psi));
//	cerr<<"a="<<a<<" b="<<b<<" c="<<c<<endl;
	cl_R d=float_approx(a)*float_approx(b)*float_approx(c);
//	cerr<<"d"<<endl;
//	cerr<<"after2"<<endl;
	double Fp=float_approx(__SS_Fp(_k-1,psi));
//	cerr<<"after3"<<endl;
	double p=float_approx(__SS_p(_k-1,psi));
//	cerr<<"after4 "<<float_approx(b)<<" "<<float_approx(a)<<" "<<p<<endl;
	double inside=float_approx(b)*float_approx(a)*p;
//	cerr<<"after5 "<<inside<<endl;
	double expp=exp(inside);
//	cerr<<"after6 "<<Fp<<" "<<expp<<endl;
	double QQ=Fp*expp;
//	cerr<<"after7 "<<QQ<<endl;
	if(std::isinf(QQ))
		return 1.0;
	
	return 1.0-QQ;
	//return 1-float
	//return 1-float_approx(__SS_Fp(_k-1,psi))*exp(float_approx(-(_k-psi)/_k)*float_approx(a)*float_approx(__SS_p(_k-1,psi)));
	
}

inline double SS_Poisson_Prob_Alm_Approx(unsigned int _k, cl_R _lambda, cl_R _w,cl_R _T)
{
	double psi=float_approx(_lambda*_w);
//	cerr<<"a"<<endl;
	cl_R a=_lambda*(_T-_w);
//	cerr<<"b"<<endl;
	cl_R b=-((_k-psi)/_k);
//	cerr<<"c "<<_k<<" "<<psi<<endl;
	cl_R c=(__SS_p(_k-1,psi));
	cerr<<"a="<<float_approx(a)<<" b="<<b<<" c="<<c<<" ";
	cl_R d=(a)*(b)*(c);
//	cerr<<"d"<<endl;
//	cerr<<"after2"<<endl;
	cl_R Fp=(__SS_Fp(_k-1,psi));
//	cerr<<"after3"<<endl;
	cl_R p=(__SS_p(_k-1,psi));
	cerr<<"after4 "<<float_approx(b)<<" "<<float_approx(a)<<" "<<p<<" ";
	cl_R inside=(b)*(a)*p;
	cerr<<"after5 "<<inside<<" ";
	cl_R expp=exp(inside);
	cerr<<"after6 "<<Fp<<" "<<expp<<" ";
	cl_R QQ=Fp*expp;
	cerr<<"after7 "<<QQ<<" ";
	//if(std::isinf(float_approx(QQ)))
	//	return 0.0;
	
	return float_approx(1.0-(QQ));
	//return 1-float
	//return 1-float_approx(__SS_Fp(_k-1,psi))*exp(float_approx(-(_k-psi)/_k)*float_approx(a)*float_approx(__SS_p(_k-1,psi)));
	
}

#endif /*SCANSTATISTICS_H_*/
