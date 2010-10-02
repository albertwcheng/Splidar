#ifndef _COLLAPSE_WIG
#define _COLLAPSE_WIG

#include <list>
#include "snipIncludes.h"

using namespace std;
#define LINE_BUFFERSIZE 1024


#define RESULT_ISLAND 0
#define RESULT_MERGELEFT 1
#define RESULT_MERGERIGHT 2
#define RESULT_MERGEBOTH 3
#define RESULT_IGNORED -1

class CollapseWig{
public:
	
	
	
	/*
		browser position chr1:-1--1
		browser hide all
		browser full acembly knownGene
		track type=wiggle_0 name="R2_chr1" description="Log2 Read Counts for R2_chr1" visibility=full color=0,0,255 priority=20
		chr1 3289859 3289888 1.0
		chr1 3410579 3410608 1.0
		chr1 3436649 3436678 1.0
		chr1 3627719 3627748 1.0
		chr1 3661229 3661258 1.0
		chr1 3661769 3661798 1.0
		chr1 3745229 3745258 1.0
		chr1 3753029 3753058 1.0
	 */
	class WigElement{
		public:
		string chrom;
		int start0;
		int end1;
		double score;
		inline string getChrom() const{
			return chrom;
		}
		inline double getScore() const{
			return score;
		}
		inline int getStart1() const{
			return start0+1;
		}
		inline int getStart0() const {
			return start0;
		}
		inline int getEnd0() const {
			return end1-1;
		}
		inline int getEnd1() const{
			return end1;
		}
		inline bool isValid() const{
			return start0<end1 && start0>=0 && end1>=1;
		}
		inline bool overlapWith(const WigElement& theOtherOne) const
		{
			return (this->getStart1()>=theOtherOne.getStart1() && this->getStart1()<= theOtherOne.getEnd1()) 
			|| (this->getEnd1()>=theOtherOne.getStart1() && this->getEnd1()<=theOtherOne.getEnd1()) 
			|| (theOtherOne.getStart1()>=this->getStart1() && theOtherOne.getStart1()<=this->getEnd1())
			|| (theOtherOne.getEnd1()>=this->getStart1() && theOtherOne.getEnd1()<=this->getEnd1());
	
		}
		
		inline bool isWithinOf(const WigElement& theOtherOne,int distance) const
		{
			return overlapWith(theOtherOne) || (this->getEnd1()<theOtherOne.getStart1() && this->getEnd1()+distance>=theOtherOne.getStart1()) || (this->getStart1()>theOtherOne.getEnd1() && this->getStart1()-distance<=theOtherOne.getEnd1());
		}
		
		inline void mergeWith(const WigElement& theOtherOne)
		{
			//for now just change the start and end
			if(!this->isValid() || !theOtherOne.isValid())
				return;
			
			this->start0=MIN(this->start0,theOtherOne.start0);
			this->end1=MAX(this->end1,theOtherOne.end1);
		}
		
		WigElement():start0(INT_MAX),end1(INT_MIN){
			
		}
		inline int getLength() const
		{
			return getEnd1()-getStart0();
		}
		WigElement(string _chrom,int _start0,int _end1):chrom(_chrom),start0(_start0),end1(_end1)
		{}
		WigElement(string _expr)
		{
			start0=INT_MAX;
			end1=INT_MIN;
			vector<string> fields;
			StringUtil::split(_expr," ",fields);
			//cerr<<_expr<<endl;
			if(_expr.substr(0,1)=="#")
			{
				chrom="comment";
			}
			if(_expr.substr(0,7)=="browser")
			{
				chrom="browser";
				return;
			}
			if(_expr.substr(0,5)=="track")
			{
				chrom="track";
			
				return;
			}
			if(fields.size()<4)
				return;
			chrom=fields[0];
			start0=StringUtil::atoi(fields[1]);
			end1=StringUtil::atoi(fields[2]);
			score=StringUtil::atof(fields[3]);
			
		}
	};
	
	
	
	
	//list<WigElement> wigs;
	map<string,list<WigElement> > chrwigs;
	
	list<WigElement>::iterator curI;
	list<WigElement>::iterator prevI;
	
	
	class ScoreRange{
	public:
		int distance;
		double score;
		double scoreDensity;
		ScoreRange(int _distance,double _score,double _scoreDensity):distance(_distance),score(_score),scoreDensity(_scoreDensity){}
	};
	
	vector<ScoreRange> scoreRanges;
	
	
	void addScoreRange(int _distance,double _score,double _scoreDensity)
	{
		scoreRanges.push_back(ScoreRange(_distance,_score,_scoreDensity));
	}
	
	double islandScore;
	double islandScoreDensity;
	//string prevChrom;
	
	void reset(list<WigElement>& wigs)
	{
		curI=wigs.begin();
		prevI=wigs.begin();
		//prevChrom="";
	}
	
	CollapseWig(double _islandScore,double _islandScoreDensity):islandScore(_islandScore),islandScoreDensity(_islandScoreDensity)
		
	{
		//reset();
		
		//scoreRanges.push_back(ScoreRange(_distanceToCollapseBlocks,_scoreToAcceptBlock,_scoreDensityToAcceptBlock));
	}
	
	
	void finalize(int minBlockSize)
	{
		for(map<string,list<WigElement> >::iterator i=chrwigs.begin();i!=chrwigs.end();i++)
		{
			finalize(i->second,minBlockSize);
		}
	}
	
	
	void finalize(list<WigElement>& wigs,int minBlockSize)
	{
		vector<list<WigElement>::iterator> toRemove;
		for(list<WigElement>::iterator i=wigs.begin();i!=wigs.end();i++)
		{
			if(i->getLength()<minBlockSize)
			{
				toRemove.push_back(i);
			}
		}
		
		for(vector<list<WigElement>::iterator >::iterator i=toRemove.begin();i!=toRemove.end();i++)
		{
			wigs.erase(*i);
		}
	}
	
	void writeBed(ostream& os,string namePrefix)
	{
		
		int indx=0;
		
		for(map<string,list<WigElement> >::iterator i=chrwigs.begin();i!=chrwigs.end();i++)
		{
			writeBed(os, namePrefix,i->second,indx);
		}		
	}
	
	void writeBed(ostream& os,string namePrefix,list<WigElement>& wigs,int& indx)
	{
		
		
		for(list<WigElement>::iterator i=wigs.begin();i!=wigs.end();i++)
		{
			indx++;
			os<<i->getChrom();
			os<<"\t";
			os<<i->getStart0();
			os<<"\t";
			os<<i->getEnd1();
			os<<"\t";
			os<<namePrefix;
			os<<indx;
			os<<"\t";
			os<<i->getScore();
			os<<endl;
		}
	}
	

	
	
	inline bool ifInMergeRange(const WigElement& biggerElement,const WigElement& targetElement)
	{
		double elementScore=targetElement.getScore();
		double elementScoreDensity=elementScore/targetElement.getLength();
		
		for(vector<ScoreRange>::iterator i=scoreRanges.begin();i!=scoreRanges.end();i++)
		{
			if(elementScore<i->score || elementScoreDensity<i->scoreDensity) //not match score requirement
				continue;
			
			if(biggerElement.isWithinOf(targetElement,i->distance))
				return true;
			
		}
		
		return false;
	}
	
	int addBlock(const WigElement& element,list<WigElement>& wigs)
	{
		//handle different cases:
		
		double elementScore=element.getScore();
		double elementScoreDensity=elementScore/element.getLength();
		
		if(curI==wigs.end())
		{
			//nothing in the list~!
			if(elementScore<islandScore || elementScoreDensity<islandScoreDensity)
				return RESULT_IGNORED;
			
			curI=wigs.insert(curI,element);
			return RESULT_ISLAND;
				
		}else
		{
			
			prevI=curI;
			
			if(curI->getChrom()!=element.getChrom())
			{
				while(curI!=wigs.end() && curI->getChrom()!=element.getChrom())
				{
					curI++;
				}
				
				
				prevI=curI;
				
			}
			
			while(curI!=wigs.end() && curI->getStart1()<element.getStart1())
			{
				prevI=curI;
				curI++;
				
				
			}
			
			//check if element overlaps with curI or prevI
			if (ifInMergeRange(*prevI,element))  //(prevI->isWithinOf(element,distanceToCollapseBlocks))
			{
				prevI->mergeWith(element);
				
				
				//extra merging
				if(prevI!=curI && curI!=wigs.end() && ifInMergeRange(*curI,element))//prevI->isWithinOf(*curI,distanceToCollapseBlocks))
				{
					prevI->mergeWith(*curI);
					wigs.erase(curI); //remove curI element update curI
					
					curI=prevI;
					return RESULT_MERGEBOTH;
				}
				else
				{
					curI=prevI;
					return RESULT_MERGELEFT;
				}
				
				
				
			}
			else
			{
				if(curI!=wigs.end() && ifInMergeRange(*curI,element))
				{
					curI->mergeWith(element);
					return RESULT_MERGERIGHT;
					
				}else
				{
					
					if(elementScore<islandScore || elementScoreDensity<islandScoreDensity)
						return RESULT_IGNORED;
					curI=wigs.insert(curI,element);
					return RESULT_ISLAND;
				}
				
			}
			
			
		}
		
		
	}
	
	void printWigs(list<WigElement>& wigs)
	{
		for(list<WigElement>::iterator i=wigs.begin();i!=wigs.end();i++)
		{
			cerr<<"["<<i->getChrom()<<":"<<i->getStart1()<<"-"<<i->getEnd1()<<"]";
		}
		cerr<<endl;
	}
	
	map<string, list<WigElement> >::iterator curWigs;
	
	void readAndCollapseWigs(string filename,bool resetOnTrack=true)
	{
		ifstream file(filename.c_str());
		
		
		char * lineBuffer=new char[LINE_BUFFERSIZE];
		int lino=0;
		

		string prevChrom;
		
		
		while(!file.eof())
		{
			
			if(file.eof())
				break;
			
			lino++;
			file.getline(lineBuffer,LINE_BUFFERSIZE); //getline
			
			if(lineBuffer[0]=='#') //comment ignore!!
				continue;
			
			WigElement newElement(lineBuffer); //feed line and construct a new wig element
			
			if(!newElement.isValid())
			{
				if(newElement.getChrom()=="track" && resetOnTrack)
				{
					cerr<<lineBuffer<<endl;
					//reset();
					prevChrom="";
				}
			}
			else //if(newElement.getScore()>=scoreToAcceptBlock && newElement.getScore()/newElement.getLength()>=scoreDensityToAcceptBlock) //check if the element is valid, to ignore headers and invalid lines
			{
				if(prevChrom!=newElement.getChrom()){
					cerr<<newElement.getChrom()<<endl;
					prevChrom=newElement.getChrom();
					
					curWigs=chrwigs.find(newElement.getChrom());
					if(curWigs==chrwigs.end())
					{
						chrwigs.insert(map<string,list<WigElement> >::value_type(newElement.getChrom(),list<WigElement>() ));
						curWigs=chrwigs.find(newElement.getChrom());
						
					}
					
					curI=curWigs->second.begin();
											prevI=curWigs->second.begin();
				}
				
				addBlock(newElement,curWigs->second);
			}
			
		}
		
		delete[] lineBuffer;
	}

	
};



#endif
