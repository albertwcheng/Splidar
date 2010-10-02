#ifndef SPLIDARGRAPH_H_
#define SPLIDARGRAPH_H_

#include "snipIncludes.h"
#include "Nucleic.h"

//#define INTRON5_LENGTH_TO_GET 250
//#define INTRON3_LENGTH_TO_GET 250
//#define EXON5_LENGTH_TO_GET 50
//#define EXON3_LENGTH_TO_GET 50

#define SPLIDARGRAPH_DEBUG_LEVEL 2

#define SG_DEBUG1(OSCHAIN) STDERR1(OSCHAIN,SPLIDARGRAPH_DEBUG_LEVEL)
#define SG_DEBUG2(OSCHAIN) STDERR2(OSCHAIN,SPLIDARGRAPH_DEBUG_LEVEL)
#define SG_DEBUG3(OSCHAIN) STDERR3(OSCHAIN,SPLIDARGRAPH_DEBUG_LEVEL)






class ExonPairWJnx
{
public:
	int fixedCoord;
	GffEntry::ExonPtr variantExon;
	GffEntry::ExonPtr commonExon;
	GffEntry::JnxPtr jnx;
	inline ExonPairWJnx(int _fixedCoord,GffEntry::ExonPtr _variantExon, GffEntry::ExonPtr _commonExon, GffEntry::JnxPtr _jnx):
		fixedCoord(_fixedCoord),variantExon(_variantExon),commonExon(_commonExon),jnx(_jnx){}
	
	inline bool operator < (const ExonPairWJnx& right) const
	{
		if(fixedCoord==right.fixedCoord)
		{
			if(variantExon==right.variantExon)
			{
				return commonExon<right.commonExon;
			}
			else
				return variantExon<right.variantExon;
		}
		else
			return fixedCoord<right.fixedCoord;
	}
	
	inline bool operator > (const ExonPairWJnx& right) const
	{
		if(fixedCoord==right.fixedCoord)
		{
			if(variantExon==right.variantExon)
			{
				return commonExon>right.commonExon;
			}
			else
				return variantExon>right.variantExon;
		}
		else
			return fixedCoord>right.fixedCoord;		
	}
	
	
	inline bool operator == ( const ExonPairWJnx& right) const
	{
		return (fixedCoord==right.fixedCoord) && (commonExon==right.commonExon) && (variantExon==right.variantExon);
	}
	
	
	inline bool operator !=(const ExonPairWJnx& right) const
	{
		return !(*this==right);
	}
	
	inline bool operator <=(const ExonPairWJnx& right) const
	{
		return !(*this>right);	
	}
	
	inline bool operator >=(const ExonPairWJnx& right) const
	{
		return !(*this<right);
	}
};

#define DEFAULT_JNX_READ_FLAG_VALUE 1
#define DEFAULT_JNX_POS_FLAG_VALUE 1

class FlankingSeq
{
public:
	static const string INTRON5;
	static const string INTRON3;
	static const string EXON5;
	static const string EXON3;
	string Type;
	string Sequence;
	inline FlankingSeq(string _Type,string _Sequence):Type(_Type),Sequence(_Sequence){}
	inline operator string() const{
		return Type+Sequence;
	}
};

class TrafficInfoAS
{
public:
	inline static double getDensity(const KeyPair<int,int>& densityVector)
	{
		return double(densityVector.k1)/double(densityVector.k2);
	}
	KeyPair<int,int> noFlankingInfo;
	KeyPair<int,int> withFlankingInfo;
	KeyPair<int,int> middleExonsFlow;
	KeyPair<int,int> flankingExonsFlow;
	KeyPair<int,int> jnxsFlow;
	set<KeyPair<int,int> > bounds;
	set<KeyPair<int,int> > boundsExonic;
	set<KeyPair<int,int> > specBounds;


	inline int getBoundsLength() const
	{
		int L=0;
		for(set<KeyPair<int,int> >::const_iterator i=bounds.begin();i!=bounds.end();i++)
		{
			L+=::len01(*i);
		}

		return L;
	}

	
	vector<KeyPair<int,int> > ibounds;
	

	
	int JRF;
	int JPF;
	string gsid;
	string jnxstring;
	string JRCheckString;
	string JPCheckString;
	vector<FlankingSeq> UFSeq;
	vector<FlankingSeq> MFSeq; //variable or middle
	vector<FlankingSeq> DFSeq; //if only one common, use this
	inline TrafficInfoAS():JRF(DEFAULT_JNX_READ_FLAG_VALUE),JPF(DEFAULT_JNX_POS_FLAG_VALUE),noFlankingInfo(0,0),withFlankingInfo(0,0), middleExonsFlow(0,0), flankingExonsFlow(0,0),jnxsFlow(0,0){}
	inline double getDensity(bool useNoFlanking=true) const
	{
		if(useNoFlanking)
			return TrafficInfoAS::getDensity(noFlankingInfo);
		else
			return TrafficInfoAS::getDensity(withFlankingInfo);
	}
	inline KeyPair<int,int> getFirstBound() const
	{
		return (*bounds.begin());
	}
	inline KeyPair<int,int> getLastBound() const
	{
		return (*bounds.rbegin());
	}
	inline string getCoordPathStringInner(const set<KeyPair<int,int> > &pathStore, KeyPair<int,int> *minMax=NULL,string sep=",",string sepIntra="-") const
	{
		string returnVal;

		if(minMax)
		{
			minMax->k1=INT_MAX;
			minMax->k2=INT_MIN;
		}

		if(pathStore.size()<1)
			return "";




		set<KeyPair<int,int> >::const_iterator i=pathStore.begin();
		returnVal=StringUtil::str(i->k1+1)+sepIntra+StringUtil::str(i->k2); //to 11

		if(minMax)
		{
			minMax->k1=MIN(minMax->k1,i->k1+1);
			minMax->k2=MAX(minMax->k2,i->k2);
		}
		i++;

		for(;i!=pathStore.end();i++)
		{
			if(minMax)
			{
				minMax->k1=MIN(minMax->k1,i->k1+1);
				minMax->k2=MAX(minMax->k2,i->k2);
			}

			returnVal+=sep+StringUtil::str(i->k1+1)+sepIntra+StringUtil::str(i->k2);
		}

		return returnVal;
	}
	inline string getCoordPathString(KeyPair<int,int> *minMax=NULL,string sep=",",string sepIntra="-") const
	{
		/*string returnVal;
		
		if(minMax)
		{	
			minMax->k1=INT_MAX;
			minMax->k2=INT_MIN;
		}
		
		if(bounds.size()<1)
			return "";
		
		
		
		
		set<KeyPair<int,int> >::const_iterator i=bounds.begin();
		returnVal=StringUtil::str(i->k1+1)+sepIntra+StringUtil::str(i->k2); //to 11
		
		if(minMax)
		{
			minMax->k1=MIN(minMax->k1,i->k1+1);
			minMax->k2=MAX(minMax->k2,i->k2);
		}
		i++;
		
		for(;i!=bounds.end();i++)
		{
			if(minMax)
			{
				minMax->k1=MIN(minMax->k1,i->k1+1);
				minMax->k2=MAX(minMax->k2,i->k2);
			}
			
			returnVal+=sep+StringUtil::str(i->k1+1)+sepIntra+StringUtil::str(i->k2);
		}
		
		return returnVal;*/
		return this->getCoordPathStringInner(bounds,minMax,sep,sepIntra);
	}

	inline string getSpecCoordPathString(KeyPair<int,int>*minMax=NULL,string sep=",",string sepIntra="-") const
	{
		return this->getCoordPathStringInner(this->specBounds,minMax,sep,sepIntra);
	}

	inline string getExonicCoordPathString(KeyPair<int,int>*minMax=NULL,string sep=",",string sepIntra="-") const
	{
		return this->getCoordPathStringInner(this->boundsExonic,minMax,sep,sepIntra);
	}
 };
	
class SplidarOutputExtraData
{
public:
	bool commonExonIsCobound;
	
};
class GenericSplidarOutputFormat
{
public:
	int readLength;
	GffEntry::Locus* locus;
	string locusName;
	ofstream* fout;
	ofstream* foutSeq;
	RandomAccessFile *raf;
	char strand;
	int eventID;
	int commonExonIsCobound;
	 string EXCELHyperLinkPrefix;
	 string EXCELHyperLinkSuffix;
	 int INTRON5_LENGTH_TO_GET;
	 int INTRON3_LENGTH_TO_GET;
	 int EXON5_LENGTH_TO_GET;
	 int EXON3_LENGTH_TO_GET;
	//static const string EXCELHyperLinkPrefix;
	//static const string EXCELHyperLinkSuffix;
	inline bool hasDataOut()
	{
		return fout;
	}
	inline bool hasSequenceOut()
	{
		return foutSeq;
	}
	inline void outData(
			const string&eventType,
			 TrafficInfoAS& inInfo,
			 TrafficInfoAS& exInfo,
			 SplidarOutputExtraData& extraDat



	)
	{
		KeyPair<int,int> minMax1,minMax2,minMax1Spec,minMax2Spec,minMax1Exonic,minMax2Exonic;
		string inCoordString=inInfo.getCoordPathString(&minMax1);
		string exCoordString=exInfo.getCoordPathString(&minMax2);
		string inSpecCoordString=inInfo.getSpecCoordPathString(&minMax1Spec);
		string exSpecCoordString=exInfo.getSpecCoordPathString(&minMax2Spec);
		string inCoordStringExonic=inInfo.getExonicCoordPathString(&minMax1Exonic);
		string exCoordStringExonic=exInfo.getExonicCoordPathString(&minMax2Exonic);
		
		int minBound=MIN(minMax1.k1,minMax2.k1);
		int maxBound=MAX(minMax1.k2,minMax2.k2);
		
		eventID++;
		
		if(fout)
		{
			//(*fout)<<eventID<<"\t";
			(*fout)<<eventType<<"\t";
			(*fout)<<locusName<<"\t";
			(*fout)<<inInfo.gsid<<"/"<<exInfo.gsid<<"\t";
			(*fout)<<inInfo.jnxstring<<"/"<<exInfo.jnxstring<<"\t";
			(*fout)<<locus->chr<<"\t";
			(*fout)<<locus->strand<<"\t";		
			(*fout)<<inCoordString<<"/"<<exCoordString<<"\t";
			(*fout)<<int(extraDat.commonExonIsCobound)<<"\t";
			(*fout)<<inCoordStringExonic<<"/"<<exCoordStringExonic<<"\t";

			(*fout)<<inSpecCoordString<<"\t"<<exSpecCoordString<<"\t";


			(*fout)<<EXCELHyperLinkPrefix<<locus->chr<<":"<<minBound<<"-"<<maxBound<<EXCELHyperLinkSuffix<<"\t";							
			
			(*fout)<<inInfo.flankingExonsFlow.k1<<"\t"<<inInfo.flankingExonsFlow.k2<<"\t";
			(*fout)<<inInfo.jnxsFlow.k1<<"\t"<<inInfo.jnxsFlow.k2<<"\t";			
			(*fout)<<inInfo.middleExonsFlow.k1<<"\t"<<inInfo.middleExonsFlow.k2<<"\t";		
			(*fout)<<inInfo.JRCheckString<<"\t"<<inInfo.JPCheckString<<"\t";
			(*fout)<<inInfo.JRF<<"\t"<<inInfo.JPF<<"\t";
			(*fout)<<exInfo.flankingExonsFlow.k1<<"\t"<<exInfo.flankingExonsFlow.k2<<"\t";
			(*fout)<<exInfo.jnxsFlow.k1<<"\t"<<exInfo.jnxsFlow.k2<<"\t";			
			(*fout)<<exInfo.middleExonsFlow.k1<<"\t"<<exInfo.middleExonsFlow.k2<<"\t";
			(*fout)<<exInfo.JRCheckString<<"\t"<<exInfo.JPCheckString<<"\t";
			(*fout)<<exInfo.JRF<<"\t"<<exInfo.JPF<<"\t";
			(*fout)<<inInfo.noFlankingInfo.k1<<"\t"<<inInfo.noFlankingInfo.k2<<"\t"; //inc
			(*fout)<<exInfo.noFlankingInfo.k1<<"\t"<<exInfo.noFlankingInfo.k2<<"\t"; //exc
			(*fout)<<(exInfo.noFlankingInfo.k1+inInfo.noFlankingInfo.k1)<<"\t"; //inc reads + exc reads
			
			double two_d=inInfo.getDensity();
			double one_d=exInfo.getDensity();
			
			(*fout)<<inInfo.noFlankingInfo.k1<<"\t"; //NI
			(*fout)<<exInfo.withFlankingInfo.k1<<"\t"; //NE+
			if(inInfo.noFlankingInfo.k1+exInfo.withFlankingInfo.k1==0)
				(*fout)<<"nan"<<"\t";
			else
				(*fout)<<(double(inInfo.noFlankingInfo.k1)/double(inInfo.noFlankingInfo.k1+exInfo.withFlankingInfo.k1))<<"\t"; //IR
									
			(*fout)<<two_d<<"\t"; //inc density
			(*fout)<<one_d<<"\t"; //exc density
	
			if(one_d+two_d==0)
				(*fout)<<"nan"<<endl;
			else
				(*fout)<<((two_d)/(one_d+two_d))<<endl; //EIL or Psi
		}
		
		if(!raf)
			cerr<<"raf not good"<<endl;
		
		if(!foutSeq)
			cerr<<"foutSeq not good"<<endl;
		
		if(raf && foutSeq)
		{
			//(*foutSeq)<<eventID<<"\t"
			(*foutSeq)<<eventType<<"\t";
			(*foutSeq)<<locusName<<"\t";
			(*foutSeq)<<inInfo.gsid<<"/"<<exInfo.gsid<<"\t";
			(*foutSeq)<<inInfo.jnxstring<<"/"<<exInfo.jnxstring<<"\t";
			(*foutSeq)<<locus->chr<<"\t";
			(*foutSeq)<<locus->strand<<"\t";
			(*foutSeq)<<inCoordString<<"/"<<exCoordString<<"\t"; 
			(*foutSeq)<<int(extraDat.commonExonIsCobound)<<"\t";
			(*foutSeq)<<inCoordStringExonic<<"/"<<exCoordStringExonic<<"\t";
			(*foutSeq)<<inSpecCoordString<<"\t"<<exSpecCoordString<<"\t";

			(*foutSeq)<<EXCELHyperLinkPrefix<<locus->chr<<":"<<minBound<<"-"<<maxBound<<EXCELHyperLinkSuffix<<"\t";							
			
			
			
			
			{
			//inc seq
				set<KeyPair<int,int> >::iterator begin=inInfo.bounds.begin(); ///####
				set<KeyPair<int,int> >::iterator meiyi=inInfo.bounds.end(); ///####
				set<KeyPair<int,int> >::iterator end=inInfo.bounds.end(); ///####
				meiyi--; ///####
				

						
				
				set<KeyPair<int,int> >::iterator i=begin;
				set<KeyPair<int,int> >::iterator j=i;
				j++;
				
				KeyPair<int,int> iBound;
				KeyPair<int,int> piBound;
				if(i!=end)
				{
					if(j==end)
					{
						raf->transfer(*foutSeq,i->k1,i->k2);
					}
					else
					{
						for(;j!=end;i++,j++)
						{
							iBound=KeyPair<int,int>(i->k2,j->k1);
							
							raf->transfer(*foutSeq,i->k1,i->k2);
							(*foutSeq)<<"|";
							
							if(i==begin)
							{
								if(this->strand==GffEntry::FORWARD)
								{
									
								
								inInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,raf->get(i->k2-MIN(EXON3_LENGTH_TO_GET,len01(*i)/2),i->k2)));
								inInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,raf->get(iBound.k1,iBound.k1+MIN(INTRON5_LENGTH_TO_GET,len01(iBound)/2))));
								}
								else
								{
									inInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,::reverse_complement(raf->get(i->k2-MIN(EXON5_LENGTH_TO_GET,len01(*i)/2),i->k2))));
									inInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,::reverse_complement(raf->get(iBound.k1,iBound.k1+MIN(INTRON3_LENGTH_TO_GET,len01(iBound)/2)))));
									
								}
								
						}else
						{
							if(this->strand==GffEntry::FORWARD)
							{
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,raf->get(piBound.k2-MIN(INTRON3_LENGTH_TO_GET,len01(piBound)/2),piBound.k2)));
								
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,raf->get(i->k1,i->k1+MIN(EXON5_LENGTH_TO_GET,len01(*i)/2))));
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,raf->get(i->k2-MIN(EXON3_LENGTH_TO_GET,len01(*i)/2),i->k2)));
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,raf->get(iBound.k1,iBound.k1+MIN(INTRON5_LENGTH_TO_GET,len01(iBound)/2))));
							}else
							{
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,::reverse_complement(raf->get(piBound.k2-MIN(INTRON5_LENGTH_TO_GET,len01(piBound)/2),piBound.k2))));
								
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,::reverse_complement(raf->get(i->k1,i->k1+MIN(EXON3_LENGTH_TO_GET,len01(*i)/2)))));
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,::reverse_complement(raf->get(i->k2-MIN(EXON5_LENGTH_TO_GET,len01(*i)/2),i->k2))));
								inInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,::reverse_complement(raf->get(iBound.k1,iBound.k1+MIN(INTRON3_LENGTH_TO_GET,len01(iBound)/2)))));
								
							}
							
						}
							
							piBound=iBound;
							
						}
						
						//now the last rigthmost exon!
						if(this->strand==GffEntry::FORWARD)
						{
						inInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,raf->get(piBound.k2-MIN(INTRON3_LENGTH_TO_GET,len01(piBound)/2),piBound.k2)));
						
						inInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,raf->get(i->k1,i->k1+MIN(EXON5_LENGTH_TO_GET,len01(*i)/2))));
						}
						else
						{
							inInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,::reverse_complement(raf->get(piBound.k2-MIN(INTRON5_LENGTH_TO_GET,len01(piBound)/2),piBound.k2))));
								
								inInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,::reverse_complement(raf->get(i->k1,i->k1+MIN(EXON3_LENGTH_TO_GET,len01(*i)/2)))));
									
						}
						
						raf->transfer(*foutSeq,i->k1,i->k2);
					
					}
				}
				
					
				(*foutSeq)<<"\t";
			
			}
			

			{
			//exc seq
				set<KeyPair<int,int> >::iterator begin=exInfo.bounds.begin(); ///####
				set<KeyPair<int,int> >::iterator meiyi=exInfo.bounds.end(); ///####
				set<KeyPair<int,int> >::iterator end=exInfo.bounds.end(); ///####
				meiyi--; ///####
				

						
				
				set<KeyPair<int,int> >::iterator i=begin;
				set<KeyPair<int,int> >::iterator j=i;
				j++;
				
				KeyPair<int,int> iBound;
				KeyPair<int,int> piBound;
				if(i!=end)
				{
					if(j==end)
					{
						raf->transfer(*foutSeq,i->k1,i->k2);
					}
					else
					{
						for(;j!=end;i++,j++)
						{
							iBound=KeyPair<int,int>(i->k2,j->k1);
							
							raf->transfer(*foutSeq,i->k1,i->k2);
							(*foutSeq)<<"|";
							
							if(i==begin)
							{
								
								if(this->strand==GffEntry::FORWARD)
								{
								exInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,raf->get(i->k2-MIN(EXON3_LENGTH_TO_GET,len01(*i)),i->k2)));
								exInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,raf->get(iBound.k1,iBound.k1+MIN(INTRON5_LENGTH_TO_GET,len01(iBound)))));
								}
								else  
								{
								exInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,::reverse_complement(raf->get(i->k2-MIN(EXON5_LENGTH_TO_GET,len01(*i)),i->k2))));
								exInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,::reverse_complement(raf->get(iBound.k1,iBound.k1+MIN(INTRON3_LENGTH_TO_GET,len01(iBound))))));
									
								}
								
							}else
							{
								if(this->strand==GffEntry::FORWARD)
								{
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,raf->get(piBound.k2-MIN(INTRON3_LENGTH_TO_GET,len01(piBound)),piBound.k2)));
								
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,raf->get(i->k1,i->k1+MIN(EXON5_LENGTH_TO_GET,len01(*i)))));
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,raf->get(i->k2-MIN(EXON3_LENGTH_TO_GET,len01(*i)),i->k2)));
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,raf->get(iBound.k1,iBound.k1+MIN(INTRON5_LENGTH_TO_GET,len01(iBound)))));
								}
								else
								{
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,::reverse_complement(raf->get(piBound.k2-MIN(INTRON5_LENGTH_TO_GET,len01(piBound)),piBound.k2))));
									
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,::reverse_complement(raf->get(i->k1,i->k1+MIN(EXON3_LENGTH_TO_GET,len01(*i))))));
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,::reverse_complement(raf->get(i->k2-MIN(EXON5_LENGTH_TO_GET,len01(*i)),i->k2))));
									exInfo.MFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,::reverse_complement(raf->get(iBound.k1,iBound.k1+MIN(INTRON3_LENGTH_TO_GET,len01(iBound))))));
									
								}
							}
							
							piBound=iBound;
							
						}
						
						//now the last rigthmost exon!
						if(this->strand==GffEntry::FORWARD)
						{
						exInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::INTRON3,raf->get(piBound.k2-MIN(INTRON3_LENGTH_TO_GET,len01(piBound)),piBound.k2)));
						
						exInfo.DFSeq.push_back(FlankingSeq(FlankingSeq::EXON5,raf->get(i->k1,i->k1+MIN(EXON5_LENGTH_TO_GET,len01(*i)))));
						}
						else
						{
							exInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::INTRON5,::reverse_complement(raf->get(piBound.k2-MIN(INTRON5_LENGTH_TO_GET,len01(piBound)),piBound.k2))));
								
							exInfo.UFSeq.push_back(FlankingSeq(FlankingSeq::EXON3,::reverse_complement(raf->get(i->k1,i->k1+MIN(EXON3_LENGTH_TO_GET,len01(*i))))));
							
						}
						
						raf->transfer(*foutSeq,i->k1,i->k2);
					
					}
				}
				
					
				(*foutSeq)<<"\t";
			
			}			
			
			
			
			//the flankings
			{
				for(vector<FlankingSeq>::const_iterator i=inInfo.UFSeq.begin();i!=inInfo.UFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				(*foutSeq)<<"\t";
				
				for(vector<FlankingSeq>::const_iterator i=inInfo.MFSeq.begin();i!=inInfo.MFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				(*foutSeq)<<"\t";				

				for(vector<FlankingSeq>::const_iterator i=inInfo.DFSeq.begin();i!=inInfo.DFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				(*foutSeq)<<"\t";
				for(vector<FlankingSeq>::const_iterator i=exInfo.UFSeq.begin();i!=exInfo.UFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				(*foutSeq)<<"\t";
				
				for(vector<FlankingSeq>::const_iterator i=exInfo.MFSeq.begin();i!=exInfo.MFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				(*foutSeq)<<"\t";				

				for(vector<FlankingSeq>::const_iterator i=exInfo.DFSeq.begin();i!=exInfo.DFSeq.end();i++)
				{
					(*foutSeq)<<string(*i)<<"|";
				}
				
				//(*foutSeq)<<"\t";				
				
				
			}
			
			
			
			(*foutSeq)<<endl;
		
		}
	}

	GenericSplidarOutputFormat(int _readLength,GffEntry::Locus* _locus, string _locusName, ofstream* _fout, ofstream* _foutSeq, RandomAccessFile* _raf,

			string _EXCELHyperLinkPrefix,
			 string _EXCELHyperLinkSuffix,
			 int _INTRON5_LENGTH_TO_GET,
			 int _INTRON3_LENGTH_TO_GET,
			 int _EXON5_LENGTH_TO_GET,
			 int _EXON3_LENGTH_TO_GET)
		:readLength(_readLength),locus(_locus),locusName(_locusName),fout(_fout),foutSeq(_foutSeq),eventID(0), raf(_raf), strand(locus->strand),

		EXCELHyperLinkPrefix(_EXCELHyperLinkPrefix),
		  EXCELHyperLinkSuffix(_EXCELHyperLinkSuffix),
		  INTRON5_LENGTH_TO_GET(_INTRON5_LENGTH_TO_GET),
		  INTRON3_LENGTH_TO_GET(_INTRON3_LENGTH_TO_GET),
		  EXON5_LENGTH_TO_GET(_EXON5_LENGTH_TO_GET),
		  EXON3_LENGTH_TO_GET(_EXON3_LENGTH_TO_GET)
	{
		
	}
};



class SplidarGraph: public SpliceTraversalGraph, public GenericSplidarOutputFormat
{
protected:
	Splidar_OpFlag op;
public:	

	//Splidar_OpFlag(bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	inline SplidarGraph(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NESuperGroupPtr _sgroup,GffEntry::Locus* _locus,string _locusName,bool _checkLocusName,int _lifeSpan,int _readLength
			):op(_op),SpliceTraversalGraph(_sgroup,_locusName,_checkLocusName,_lifeSpan),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{
		strand=locus->strand;
	}
	
	inline SplidarGraph(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::NExonGroupPtr _exongroup,GffEntry::Locus* _locus,string _locusName,bool _checkLocusName,int _lifeSpan,int _readLength
			 ):op(_op),SpliceTraversalGraph(_exongroup,_locusName,_checkLocusName,_lifeSpan),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf,   _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{
		strand=locus->strand;
	}
	
	inline SplidarGraph(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::IntExonPtrMMapDI _range,GffEntry::Locus* _locus,string _locusName,bool _checkLocusName,int _lifeSpan,int _readLength
			):op(_op),SpliceTraversalGraph(_range,_locusName,_checkLocusName,_lifeSpan),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf,   _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{
		strand=locus->strand;
	}	
	
	//old constructor
	/*inline SplidarGraph(ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::NExonGroupPtr _exongroup,GffEntry::Locus* _locus,string _locusName,bool _checkLocusName,int _lifeSpan,int _readLength
			 ):SpliceTraversalGraph(_exongroup,_locusName,_checkLocusName,_lifeSpan),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _EXCELHyperLinkPrefix,_EXCELHyperLinkSuffix,_INTRON5_LENGTH_TO_GET,_INTRON3_LENGTH_TO_GET,_EXON5_LENGTH_TO_GET)
	{
		strand=locus->strand;
	}
	
	inline SplidarGraph(ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::IntExonPtrMMapDI _range,GffEntry::Locus* _locus,string _locusName,bool _checkLocusName,int _lifeSpan,int _readLength
			string _EXCELHyperLinkPrefix,
			string _EXCELHyperLinkSuffix,
			int _INTRON5_LENGTH_TO_GET,
			int _INTRON3_LENGTH_TO_GET,
			int _EXON5_LENGTH_TO_GET,
			int _EXON3_LENGTH_TO_GET):SpliceTraversalGraph(_range,_locusName,_checkLocusName,_lifeSpan),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _EXCELHyperLinkPrefix,_EXCELHyperLinkSuffix,_INTRON5_LENGTH_TO_GET,_INTRON3_LENGTH_TO_GET,_EXON5_LENGTH_TO_GET)
	{
		strand=locus->strand;
	}	*/
	
	~SplidarGraph()
	{
		
	}
	
	 void startGraph()
	{
		cerr<<"starting graph";
	}
	 void discoverExon(const ExonPostItem& curItem)
	{
		
	}
	 void discoverExonGroup(const ExonGroupPostItem& curItem)
	{
		
	}
	 void arriveVisitedExon(const ExonPostItem&curItem, ExonPostRow* prevItems)
	{
		
	}
	 void arriveVisitedExonGroup(const ExonGroupPostItem&curItem, ExonGroupPostRow* prevItems)
	{
		
	}
	 

	

	

	
	TrafficInfoAS getTrafficAS(AESPathThread::backtrack_iterator beg,AESPathThread::backtrack_iterator end1,bool includeFlankExon=true)
	{
		AESPathThread::backtrack_iterator i=beg;

		TrafficInfoAS returnVal;
		//cerr<<"gt"<<"a"<<endl;
		for(i=beg;i!=end1;i++)
		{
			//if(excludeFlankingi==beg)
			//cerr<<"gt"<<"a0"<<endl;
			const InEdgeAndVertex& iev=*i;
			if(i==beg) //right flanking exon, and jnx before.
			{
				//cerr<<"gt"<<"a1"<<endl;
				returnVal.bounds.insert(iev.exon->getBound());
				returnVal.gsid=iev.exon->exonGroup->sid;
				
				if(this->op.getCount && this->hasDataOut())
				{
				KeyPair<int,int> jnxDen=iev.inJnx->getDensity(readLength,true);
				
				returnVal.JRCheckString=StringUtil::str(jnxDen.k1);
				returnVal.JPCheckString=StringUtil::str(jnxDen.k2);
				returnVal.JRF*=jnxDen.k1;
				returnVal.JPF*=jnxDen.k2;
				
				addDensityVectors(returnVal.jnxsFlow,jnxDen);
				addDensityVectors(returnVal.noFlankingInfo,jnxDen);				
				addDensityVectors(returnVal.withFlankingInfo,jnxDen);
				
				KeyPair<int,int> flankExonDen=iev.exon->getDensity(readLength,true);
				
				addDensityVectors(returnVal.flankingExonsFlow,flankExonDen);
				
				if(includeFlankExon)
					addDensityVectors(returnVal.withFlankingInfo,flankExonDen);
				}
				//cerr<<"gt"<<"<a1"<<endl;
			}
			else if(!i.hasAtLeastOneMore()) //left flanking exon (no jnx before)
			{
				//cerr<<"gt"<<"a3"<<endl;
				
				
				returnVal.bounds.insert(iev.exon->getBound());
				returnVal.gsid=iev.exon->exonGroup->sid+"_"+returnVal.gsid;
				
				if(this->op.getCount&&this->hasDataOut())
				{
				KeyPair<int,int> flankExonDen=iev.exon->getDensity(readLength,true);
				
				addDensityVectors(returnVal.flankingExonsFlow,flankExonDen);
				
				if(includeFlankExon)
					addDensityVectors(returnVal.withFlankingInfo,flankExonDen);				
				}
				//cerr<<"gt"<<"<a3"<<endl;
			}
			else //middle exon and jnx to the left
			{
				//cerr<<"gt"<<"a2"<<endl;
				returnVal.bounds.insert(iev.exon->getBound());
				returnVal.gsid=iev.exon->exonGroup->sid+"_"+returnVal.gsid;
				
				if(this->op.getCount&&this->hasDataOut())
				{
				
				KeyPair<int,int> jnxDen=iev.inJnx->getDensity(readLength,true);
				KeyPair<int,int> exonDen=iev.exon->getDensity(readLength,true);
				
				returnVal.JRCheckString=StringUtil::str(jnxDen.k1)+","+returnVal.JRCheckString;
				returnVal.JPCheckString=StringUtil::str(jnxDen.k2)+","+returnVal.JPCheckString;
				returnVal.JRF*=jnxDen.k1;
				returnVal.JPF*=jnxDen.k2;											
				
				addDensityVectors(returnVal.noFlankingInfo,jnxDen);	
				addDensityVectors(returnVal.noFlankingInfo,exonDen);				
				addDensityVectors(returnVal.withFlankingInfo,jnxDen);
				addDensityVectors(returnVal.withFlankingInfo,exonDen);	
				
				addDensityVectors(returnVal.middleExonsFlow,exonDen);
				addDensityVectors(returnVal.jnxsFlow,jnxDen);
				}
				//cerr<<"gt"<<"<a2"<<endl;
			}
		}
		
		return returnVal;
		
		
		
	}
	

	 void _endGraph()
	{
		cerr<<"end graph"<<endl;
		//now try to analyze the structure if ExonPostMatrix
		for(ExonPostMatrix::I i=this->exonPostMatrix.begin();i!=this->exonPostMatrix.end();i++)
		{
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;
		//	cerr<<"a:"<<endl;
			SrcExonPostStruct * srcS=epr->getSrcExonStructure();
	//		cerr<<"b"<<endl;
			//now here comes the easier but critical step!
	
			for(SrcExonPostStruct::I1 i1=srcS->begin();i1!=srcS->end();i1++)
			{
				//cerr<<"b.1"<<endl;
				GffEntry::ExonPtr srcExon=i1->first;
				SrcExonPostStruct::S2* s2=i1->second;
				//cerr<<"b.2"<<endl;
				for(SrcExonPostStruct::I2 i2=s2->begin();i2!=s2->end();i2++)
				{
					//cerr<<"b.3"<<endl;
					int time=i2->first;
					//cerr<<"b.3a"<<endl;
					SrcExonPostStruct::S3 *s3=i2->second;
					//cerr<<"b.3b"<<endl;
					//s3 holds the different routes
					//cerr<<"b.4:"<<srcExon<<","<<terminalExon<<endl;
					cerr<<"the different routes for exon transition from "<<srcExon->getBound();
					cerr<<" to "<<terminalExon->getBound()<<" with " <<time<<" steps:"<<endl;
					//cerr<<"b.5"<<endl;
					for(SrcExonPostStruct::I3 i3=s3->begin();i3!=s3->end();i3++)
					{
						ExonPostItem* epi=*i3;
						AESPathThread::backtrack_iterator bi=epi->exonThread->rbegin(epi->vertexPointer);
						AESPathThread::backtrack_iterator be=epi->exonThread->rend();
						for(;bi!=be;bi++)
						{
							const InEdgeAndVertex& iev=*bi;
							cerr<<iev.exon->getBound()<<"\t";
						}
						
						cerr<<endl;
					}
				}
			}
			
		//	cerr<<"c"<<endl;
			delete srcS;
		//	cerr<<"d"<<endl;
		}
		
	}
	 void endPath(AESPathThread* exonThread)
	{
		
	}
};

template<typename ContainerClass,typename KeyClass>
bool key_exists(const ContainerClass& container,const KeyClass& key )
{
	return container.find(key)!=container.end();
}

#include "AEP_Splidar.h"
#include "A53SS_Splidar.h"
#include "ATE_Splidar.h"
#include "RI_Splidar.h"
























#endif /*SPLIDARGRAPH_H_*/
