/**
 * \file GffEntry.h
 * \author W. Albert Cheng (awcheng@mit.edu).
 * \brief Parser, storer of Gff entries as exons and transcripts
 *
 */
#ifndef GFFENTRY_H_
#define GFFENTRY_H_

#include "snipIncludes.h"
#include <sstream>
#include <map>
//#include <multimap>
#include <vector>
#include "sam.h"


#define DEFAULT_ANNO_SOURCE "GffEntry"
#define pExonOf(i) ((*i).first)
#define exonOf(i) (*((*i).first))
#define gffsOf(i) (*((*i).second))
#define foreachGff(s,i) for(vector<RankedGffEntryPair>::iterator i=s.begin();i!=s.end();i++)

#define USAGE_UNKNOWN 0
#define USAGE_ASSIGNEXONGROUP 1
#define USAGE_COUNTBLOCKFREQ 2
#define MAQ_FORMAT 0
#define SAM_FORMAT 1
#define BAM_FORMAT 2


inline int getReadMapFormatCode(const string& format)
{
	string readMapFormat=StringUtil::toUpper(format);

	if(readMapFormat=="MAPVIEW")
	{
		return MAQ_FORMAT;
	}
	else if(readMapFormat=="SAM")
	{
		return SAM_FORMAT;
	}
	else
	{
		cerr<<"Unknown format "<<readMapFormat<<endl;
		dier("",0);

	}
}

#define COMPACTPOSITIONSUFFIX ".fold.s"
class GffEntry;

inline void addDensityVectors(KeyPair<int,int>& to, const KeyPair<int,int>& from)
{
	to.k1+=from.k1;
	to.k2+=from.k2;
}

inline bool isSmaller(const KeyPair<int,int>& left01, const KeyPair<int,int>& right01)
{
	return left01.k2<=right01.k1;
}

inline bool isAinB(const KeyPair<int,int>&a01, const KeyPair<int,int>&b01)
{
	return b01.k1<=a01.k1 && b01.k2>=a01.k2;
}


inline void overlapBoundMutFirst( KeyPair<int,int>& b1, const KeyPair<int,int> &b2)
{
	b1.k1=MAX(b1.k1,b2.k1);
	b1.k2=MIN(b1.k2,b2.k2);
}

inline void expandBoundMutFirst(KeyPair<int,int>&b1, const KeyPair<int,int>& b2)
{
	b1.k1=MIN(b1.k1,b2.k1);
	b1.k2=MAX(b1.k2,b2.k2);
}

inline bool isValid(const KeyPair<int,int>& bound)
{
	return (bound.k1<bound.k2 && bound.k1>=0);
}

inline bool areOverlapping(const KeyPair<int,int>& a01, const KeyPair<int,int>& b01)
{
	if(!isValid(a01) || !isValid(b01))
		return false;


	return ((b01.k1>=a01.k1 && b01.k1<a01.k2) ||
			(b01.k2>a01.k1 && b01.k2<=a01.k2) ||
			(a01.k1>=b01.k1 && a01.k1<b01.k2) ||
			(a01.k2>b01.k1 && a01.k2<=b01.k2)
	       );
}

inline void expandBound(KeyPair<int,int>& target, const KeyPair<int,int>& source)
{
	target.k1=MIN(target.k1,source.k1);
	target.k2=MAX(target.k2,source.k2);
}

inline KeyPair<int,int> overlapBound(const KeyPair<int,int>& a, const KeyPair<int,int>& b)
{
	return KeyPair<int,int>(MAX(a.k1,b.k1),MIN(a.k2,b.k2));
}


class BoundSet
{
public:
	typedef KeyPair<int,int> Bound;

	set<Bound> _set;

	typedef set<Bound>::iterator iterator;
	typedef set<Bound>::const_iterator const_iterator;


	inline void insertInternal(const Bound& _b)
	{
		if(!isValid(_b))
			return;

		_set.insert(_b);
	}

	inline BoundSet(const Bound& _b)
	{
		insertInternal(_b);
	}

	inline BoundSet()
	{

	}

	inline bool isEmpty()
	{
		return _set.size()==0;
	}




	inline void findOverlapIteratorPairs(const_iterator& overlapFirst, const_iterator& overlapLastInclusive, const_iterator& overlapLast,const Bound& _b,const_iterator starter,const_iterator ender) const
	{
		overlapFirst=_set.end();
		overlapLastInclusive=_set.end();
		for(const_iterator i=starter;i!=ender;i++)
		{
			if(::areOverlapping(*i,_b))
			{
				if(overlapFirst==_set.end())
				{
					overlapFirst=i;

				}

				overlapLastInclusive=i;


			}else if(overlapFirst!=_set.end())
			{
				//not overlapping and has been overlapped in previous. chain broken: quit
				break;
			}
		}

		if(overlapFirst!=_set.end()){
			overlapLast=overlapLastInclusive;
			overlapLast++;
		}
	}
	inline void findOverlapIteratorPairs(const_iterator& overlapFirst, const_iterator& overlapLastInclusive, const_iterator& overlapLast,const Bound& _b) const
	{
		return findOverlapIteratorPairs(overlapFirst, overlapLastInclusive, overlapLast, _b,begin(),end());
	}

	inline void addBound(const Bound& _b)
	{
		if(!isValid(_b))
			return;

		const_iterator overlapFirst;
		const_iterator overlapLast;
		const_iterator overlapLastInclusive;

		findOverlapIteratorPairs(overlapFirst,overlapLastInclusive,overlapLast,_b);

		if(overlapFirst==_set.end())
		{
			//no overlap, add directly
			insertInternal(_b);
		}
		else
		{
			KeyPair<int,int> firstOverlappedBound=*overlapFirst;
			KeyPair<int,int> lastOverlappedBound=*overlapLastInclusive;

			if(firstOverlappedBound==lastOverlappedBound && _b.k1>=firstOverlappedBound.k1 && _b.k2<=firstOverlappedBound.k2)
				return;

			_set.erase(overlapFirst,overlapLast);

			insertInternal(Bound(MIN(firstOverlappedBound.k1,_b.k1),MAX(lastOverlappedBound.k2,_b.k2)));



		}


	}

	inline void addBound(int k1,int k2)
	{
		addBound(Bound(k1,k2));
	}

	inline bool overlaps(const Bound& _b)
	{
		for(iterator i=_set.begin();i!=_set.end();i++)
		{
			if(areOverlapping(*i,_b))
				return true;
		}

		return false;
	}

	inline bool getOverlapLength(const Bound& _b)
	{
		int len=0;
		for(iterator i=_set.begin();i!=_set.end();i++)
		{
			len+=::overlaplen0101(*i,_b);
		}

		return len;
	}



	inline void overlapInPlace(const Bound& _b)
	{
		const_iterator overlapFirst;
		const_iterator overlapLast;
		const_iterator overlapLastInclusive;

		findOverlapIteratorPairs(overlapFirst,overlapLastInclusive,overlapLast,_b);

		//here!
		if(overlapFirst==_set.end())
		{
			//no overlap,  clear everything
			_set.clear();
		}
		else
		{
			KeyPair<int,int> firstOverlappedBound=*overlapFirst;
			KeyPair<int,int> lastOverlappedBound=*overlapLastInclusive;


			iterator overlapSecond=overlapFirst;
			overlapSecond++;
			_set.erase(_set.begin(),overlapSecond); //erase everything to the overlap first inclusive
			_set.erase(overlapLastInclusive,_set.end()); //erase everything from the last overlap to the end


			firstOverlappedBound.k1=MAX(firstOverlappedBound.k1,_b.k1);
			lastOverlappedBound.k2=MIN(lastOverlappedBound.k2,_b.k2);

			insertInternal(firstOverlappedBound);
			insertInternal(lastOverlappedBound);


		}

	}
	inline iterator begin()
	{
		return _set.begin();
	}
	inline iterator end()
	{
		return _set.end();

	}
	inline const_iterator begin() const
	{
		return _set.begin();
	}
	inline const_iterator end() const
	{
		return _set.end();

	}
	inline void subtractInPlaceBy( const BoundSet& subtractor)
	{
		set<Bound> new_set;

		const_iterator sub_i=subtractor._set.begin();
		if(sub_i==subtractor._set.end()) //subtractor is empty.
			return;

		const_iterator sub_e=subtractor.end();

		for(iterator i=_set.begin();i!=_set.end();i++)
		{
			//now find all overlap of *i in subtractor
			const_iterator overlapFirst;
			const_iterator overlapLast;
			const_iterator overlapLastInclusive;

			subtractor.findOverlapIteratorPairs(overlapFirst,overlapLastInclusive,overlapLast,*i,sub_i,sub_e);

			if(overlapFirst==subtractor._set.end())
			{
				//no overlap
				new_set.insert(*i);
			}
			else
			{
				int new_k1=i->k1;

				for(iterator j=overlapFirst;j!=overlapLast;j++)
				{
					if(j->k1<=i->k1)
					{
						if(j->k2>=i->k2)
						{
							// [i ]   subtractee
							//[ j  ]  subtractor
							//becomes null

							//invalidate new_k1
							new_k1=-1;
							break;
						}else
						{
							// [i  ]
							//[ j]
							new_k1=j->k2;
						}

					}else if(j->k1>i->k1)
					{
						//form new bound with previous
						// [ i  ]
						//  ( [j  ]
						//  ( [j]
						//  ^prev new_k1
						new_set.insert(Bound(new_k1,j->k1));
						if(j->k2<i->k2)
						{
							// [ i    ]
							//   [j](
							//      ^new_k1
							new_k1=j->k2;
						}
						else
						{
							//invalidate new_k1;
							new_k1=-1;
						}
					}


				}


				if(new_k1!=-1)
				{
					// [ i   ]
					//[ ] [](
					//      ^new_k1

					new_set.insert(Bound(new_k1,i->k2));
				}

				sub_i=overlapLastInclusive; //speed up next one;
			}


		}

		_set=new_set; //repalce set;


	}
};

/**
 * A class to give the information of the membership of an exon to a transcript and the corresponding rank of the exon in that transcript
 *
 */

class RankedGffEntryPair
{
public:
	int rank;
	GffEntry* gff;
	RankedGffEntryPair(int _rank,GffEntry *_gff):rank(_rank),gff(_gff){}

};

/**
 * A GffEntry is a data structure representing an entry (A line) in the gff file.
 * It also represent a transcript
 */

#define INDEXEDEXON_BUFFER_OVERHEAD 100

class NExonGroup;
class NESuperGroup;
//class ExonBoundGroup;

class GffEntry
{
private:
	GffEntry(vector<string>& splits,string annoSource=DEFAULT_ANNO_SOURCE);
	static int GffEntryCount;
	static int GffExonCount;

public:
	class Exon;
	class Loader;
	class TranscriptLoadingListener;
	class Locus;
	class GBlock;
	class Jnx;
	class JnxTag;
	class JnxTagBlock;
	class SelexaMatch;

	typedef SmartPtr<Jnx> JnxPtr;
	//typedef SmartPtr<JnxTagBlock> JnxTagBlockPtr;

	typedef SmartPtr<Exon> ExonPtr;
	typedef SmartPtr<GBlock> GBlockPtr;
	typedef SmartPtr<SelexaMatch> SelexaMatchPtr;

	static int gTotalExonReadsUsed;

	static bool isUMPLoaded;

	static inline bool isUniquelyMappablePositionLoaded()
	{
		return isUMPLoaded;
	}

	//class ExonPtr;

	/**
	 * A map that store the exonome. mapping each exon to n transcripts (several transcripts may share the same exon)
	 */
	static map<ExonPtr,vector<RankedGffEntryPair>* > exonome;
	static vector<ExonPtr> indexedExonome;
	static map<string, set<GBlockPtr> *> gblocks; //chr -> start-ordered Gblocks Ptr

	typedef map<string,set<GBlockPtr>* >::iterator ChrGBlockI;
	typedef set<GBlockPtr>::iterator GBlockI;

	//static map<string, vector<Jnx*>* > jnxs;
	//static map<string, vector<JnxTagBlock*>* > jnxtagblocks;
	static map<string, map<int,JnxTag*>* > jnxtags; //chr -> [ id -> jnxTag Ptr ]



	class JnxTagID
	{
	public:
		string chr;
		string setName;
		int id;
		inline operator string() const
		{
			return chr+"|"+setName+":"+StringUtil::str(id);
		}
		void set(const string& strID)
		{
			vector<string> splitons;
			StringUtil::split(strID,"|",splitons,true);
			if(splitons.size()<2)
				return;

			chr=splitons[0];
			StringUtil::split(splitons[1],":",splitons,true);
			if(splitons.size()<2)
				return;
			setName=splitons[0];
			id=StringUtil::atoi(splitons[1]);
		}
		JnxTagID(const string& strID)
		{
			set(strID);
		}
		JnxTagID(){}
		JnxTagID(const string& _chr,const string& _setName, int _id):chr(_chr),setName(_setName),id(_id){}
	};



	static inline void resetJnxTags(const map<string, map<int,GffEntry::JnxTag*>* > ::iterator& i)
	{
		cerr<<"reseting jnxs on "<<(*i).first<<endl;
		map<int,GffEntry::JnxTag*>* sset=(*i).second;
		for(map<int,GffEntry::JnxTag*>::iterator j=sset->begin();j!=sset->end();j++)
		{
			delete (*j).second;
		}
		delete sset;
	}

	static inline void resetJnxTags(string chr="")
	{
		if(chr=="")
		{
			for(map<string, map<int,GffEntry::JnxTag*>* >::iterator i=jnxtags.begin();i!=jnxtags.end();i++)
			{
				GffEntry::resetJnxTags(i);
			}
			//reset everything;
		}
		else
		{
			map<string, map<int,GffEntry::JnxTag*>* >::iterator i=jnxtags.find(chr);
			if(i==jnxtags.end())
				return;
			GffEntry::resetJnxTags(i);

			jnxtags.erase(i);

		}
	}

	/*static inline void resetJnxs(const map<string,vector<GffEntry::JnxPtr>* >::iterator& i)
	{
		cerr<<"reseting jnxs on "<<(*i).first<<endl;
		vector<GffEntry::JnxPtr>* sset=(*i).second;
		for(vector<GffEntry::JnxPtr>::iterator j=sset->begin();j!=sset->end();j++)
		{
			delete (*j);
		}
		delete sset;
	}*/

	/*static inline void resetJnxs(string chr="")
	{
		if(chr=="")
		{
			for(map<string, vector<GffEntry::JnxPtr>* >::iterator i=jnxs.begin();i!=jnxs.end();i++)
			{
				GffEntry::resetJnxs(i);
			}
			//reset everything;
		}
		else
		{
			map<string,vector<GffEntry::JnxPtr>* >::iterator i=jnxs.find(chr);
			if(i==jnxs.end())
				return;
			GffEntry::resetJnxs(i);
			
		}
	}*/




	/*static inline void resetJnxTagBlocks(const map<string,vector<GffEntry::JnxTagBlock* >* >::iterator& i)
	{
		cerr<<"reseting jnxs on "<<(*i).first<<endl;
		vector<JnxTagBlock*>* sset=(*i).second;
		for(vector<JnxTagBlock*>::iterator j=sset->begin();j!=sset->end();j++)
		{
			delete (*j);
		}
		delete sset;
	}

	static inline void resetJnxTagBlocks(string chr="")
	{
		if(chr=="")
		{
			for(map<string,vector<JnxTagBlock*> *>::iterator i=jnxtagblocks.begin();i!=jnxtagblocks.end();i++)
			{
				GffEntry::resetJnxTagBlocks(i);
			}
			//reset everything;
		}
		else
		{
			map<string,vector<JnxTagBlock*> *>::iterator i=jnxtagblocks.find(chr);
			if(i==jnxtagblocks.end())
				return;
			GffEntry::resetJnxTagBlocks(i);

		}
	}	*/
		
	static inline void resetAllJnxStuff()
	{
		//GffEntry::resetJnxTags();
		//GffEntry::resetJnxTagBlocks();
		GffEntry::resetJnxTags();
		//GffEntry::resetJnxs();
	}



	class Jnx ///modifying this consider modidying JxnTagSet as well
	{

	public:
		//string chr;
		char strand;
		//int gUp;
		//int gDown;
		GffEntry::ExonPtr leftExon;
		GffEntry::ExonPtr rightExon;
		vector<JnxTagBlock*> jblocks;
		Jnx(GffEntry::ExonPtr _left,GffEntry::ExonPtr _right,char _strand):leftExon(_left),rightExon(_right),strand(_strand){}
		
		inline bool getNaiveAmbiguity() const
		{
			return strand==GffEntry::BOTH_STRANDS;
		}
		
		
		inline int getgUpIntron0() const
		{
			return leftExon->getEnd1();
		}

		inline int getgDownIntron1() const
		{
			return rightExon->getStart0();
		}

		inline int getgUpExon0() const
		{
			return leftExon->getEnd0();
		}
		inline int getgDownExon1() const
		{
			return rightExon->getStart1();
		}

		inline KeyPair<int,int> getDensity(int readLength,bool useUniquelyMappablePosition)
		{
			KeyPair<int,int> result(0,0);

			for(vector<JnxTagBlock*>::iterator i=jblocks.begin();i!=jblocks.end();i++)
			{
				JnxTagBlock* jblock=(*i);
				result.k1+=jblock->getSelexaMatchFreq();
				
				if(useUniquelyMappablePosition)
					result.k2+=jblock->getUniquelyMappablePosition();
				else
					result.k2+=jblock->getMappablePositions();
			}

			return result;
		}
	};

	class JnxTagSet ////
	{

	public:
		//string chr;
		char strand;
		//int gUp;
		//int gDown;

		GffEntry::GBlockPtr leftBlock;
		GffEntry::GBlockPtr rightBlock;

		inline KeyPair<int,int> getBound11() const
		{
			return KeyPair<int,int>(leftBlock->getEnd1()+1,rightBlock->getStart1()-1);
		}

		inline bool operator > (const JnxTagSet& right) const
		{
			return this->getBound11()>right.getBound11();
		}

		inline bool operator < (const JnxTagSet& right) const
		{
			return this->getBound11()<right.getBound11();
		}

		inline bool operator == ( const JnxTagSet& right) const
		{
			return this->getBound11()==this->getBound11(); //is this sufficient??!
		}

		inline bool operator !=(const JnxTagSet& right) const
		{
			return !(*this==right);
		}

		inline bool operator >=(const JnxTagSet& right) const
		{
			return !(*this<right);

		}

		inline bool operator <=(const JnxTagSet& right) const
		{
			return !(*this>right);
		}



		static const int Union;
		static const int Intersection;

		set<JnxTagBlock*> jblocks;

		inline void unionWith(const Jnx& _jnx)
		{
			for(vector<JnxTagBlock*>::const_iterator i=_jnx.jblocks.begin();i!=_jnx.jblocks.end();i++)
			{
				this->jblocks.insert(*i);
			}
		}




		inline JnxTagSet(const Jnx& _jnx):strand(_jnx.strand),leftBlock(_jnx.leftExon->getGenomicLastBlock()),rightBlock(_jnx.rightExon->getGenomicFirstBlock())
		{

			this->unionWith(_jnx);
		}

		inline void insectWith(const Jnx& _jnx)
		{
			vector<JnxTagBlock*> toDelete;
			for(set<JnxTagBlock*>::iterator i=this->jblocks.begin();i!=this->jblocks.end();i++)
			{
				bool presentInForeignSet=false;
				for(vector<JnxTagBlock*>::const_iterator j=_jnx.jblocks.begin();j!=_jnx.jblocks.end();j++)
				{
					if ((*i)==(*j)){
							presentInForeignSet=true;
							break;
					}
				}

				if(!presentInForeignSet)
				{
					toDelete.push_back(*i);
				}

			}

			for(vector<JnxTagBlock*>::iterator i=toDelete.begin();i!=toDelete.end();i++)
			{
				this->jblocks.erase(*i);
			}
		}

		inline void addTags(const Jnx& _jnx,int mode)
		{
			if(mode==GffEntry::JnxTagSet::Union)
			{
				this->unionWith(_jnx);
			}
			else
			{
				this->insectWith(_jnx);
			}
		}


		inline bool getNaiveAmbiguity() const
		{
			return strand==GffEntry::BOTH_STRANDS;
		}


		inline KeyPair<int,int> getDensity(int readLength,bool useUniquelyMappablePosition)
		{
			KeyPair<int,int> result(0,0);

			for(set<JnxTagBlock*>::iterator i=jblocks.begin();i!=jblocks.end();i++)
			{
				JnxTagBlock* jblock=(*i);
				result.k1+=jblock->getSelexaMatchFreq();

				if(useUniquelyMappablePosition)
					result.k2+=jblock->getUniquelyMappablePosition();
				else
					result.k2+=jblock->getMappablePositions();
			}

			return result;
		}
	};

	typedef SmartPtr<JnxTagSet> JnxTagSetPtr;


	class JnxTagBlock
	{
	public:
		JnxTag* tag;
		int startOnTag1;
		int endOnTag1;
		int blockSpanStart;
		int blockSpanEndInc;
		int startIndexSelexaMatchOnTag;
		int endIndexSelexaMatchOnTagInc;
		int umpos;

		inline JnxTagBlock(JnxTag* _tag=NULL,int _startOnTag1=1,int _endOnTag1=1,int _blockSpanStart=0,int _blockSpanEndInc=0,int _startIndexSelexaMatchOnTag=0,int _endIndexSelexaMatchOnTagInc=-1):
			tag(_tag),startOnTag1(_startOnTag1),endOnTag1(_endOnTag1),blockSpanStart(_blockSpanStart),blockSpanEndInc(_blockSpanEndInc),startIndexSelexaMatchOnTag(_startIndexSelexaMatchOnTag),endIndexSelexaMatchOnTagInc(_endIndexSelexaMatchOnTagInc),umpos(0)
			{}

		inline int getUniquelyMappablePosition() const
		{
			return umpos;
		}

		inline int getMappablePositions() const
		{
			return endOnTag1-startOnTag1+1;
		}

		inline int getSelexaMatchFreq() const
		{
			if(endIndexSelexaMatchOnTagInc<startIndexSelexaMatchOnTag)
				return 0;

			return endIndexSelexaMatchOnTagInc-startIndexSelexaMatchOnTag+1;
		}


		inline KeyPair<int,int> getBound11() const
		{
			return KeyPair<int,int>(startOnTag1,endOnTag1);
		}

		inline int getEnd1() const
		{
			return endOnTag1;
		}

		inline int getStart1() const
		{
			return startOnTag1;
		}

		inline int getLength() const
		{
			return len11(getStart1(),getEnd1());
		}

	};
	
	class JnxTag
	{
	public:
		multiset<SelexaMatchPtr> selexaMatches;
		JnxTagID id;
		//blocks/bounds size
		char strand;
		int lengthTag;
		//skip 1 field
		//exonpath size
		int boundSize;
		int exonPathSize;
		KeyPair<int,int>* bounds;
		GffEntry::Exon*** exonPaths;
		vector<JnxTagBlock*> jblocks;

		inline ~JnxTag()
		{
			if(bounds)
				delete[] bounds;

			if(exonPaths)
			{
				for(int i=0;i<exonPathSize;i++)
				{
					delete[] exonPaths[i];
				}

				delete[] exonPaths;
			}

			for(vector<JnxTagBlock*>::iterator i=jblocks.begin();i!=jblocks.end();i++)
			{
				delete (*i);
			}

			for(multiset<SelexaMatchPtr>::iterator j=selexaMatches.begin();j!=selexaMatches.end();j++)
			{
				delete (*j);
			}
		}
		inline void allocBoundsAndExonPaths()
		{
			bounds=new KeyPair<int,int>[this->boundSize];
			exonPaths=new GffEntry::Exon**[this->exonPathSize];
			for(int i=0;i<exonPathSize;i++)
			{
				exonPaths[i]=new GffEntry::Exon*[this->boundSize];
			}
		}

		inline JnxTag(): bounds(NULL),exonPaths(NULL){}
		
		//new: Apr29, 2009: correct for start
		inline void updateJnxTagBlockWithSelexaData()
		{
				multiset<SelexaMatchPtr>::iterator i=selexaMatches.begin();

				if(i==selexaMatches.end())
					return;

				int start=0;
				int end=0;

				//care about end
				for(vector<JnxTagBlock*>::iterator j=jblocks.begin();j!=jblocks.end();j++)
				{
						JnxTagBlock* jblock=(*j);


						//new:
						jblock->endIndexSelexaMatchOnTagInc=0;

						if((*i)->getPos1()>jblock->endOnTag1)
							continue;

						while(i!=selexaMatches.end() && (*i)->getPos1()<=jblock->endOnTag1)
						{
							i++;
							end++;
						}


						//jblock->startIndexSelexaMatchOnTag=0//start;
						jblock->endIndexSelexaMatchOnTagInc=end-1;

						if(i==selexaMatches.end())
							break;
				}

				//care about start

				multiset<SelexaMatchPtr>::reverse_iterator ri=selexaMatches.rbegin();

				start=selexaMatches.size()-1;

				for(vector<JnxTagBlock*>::reverse_iterator j=jblocks.rbegin();j!=jblocks.rend();j++)
				{
					JnxTagBlock* jblock=(*j);

					jblock->startIndexSelexaMatchOnTag=10000;

					if((*ri)->getPos1()<jblock->startOnTag1)
						continue;

					while(ri!=selexaMatches.rend() && (*ri)->getPos1()>=jblock->startOnTag1)
					{
						ri++;
						start--;
					}

					jblock->startIndexSelexaMatchOnTag=start+1;

					if(ri==selexaMatches.rend())
						break;
				}
		}

		//before Apr29
		/*inline void updateJnxTagBlockWithSelexaData()
		{
			vector<SelexaMatchPtr>::iterator i=selexaMatches.begin();

			if(i==selexaMatches.end())
				return;

			int start=0;
			int end=0;

			for(vector<JnxTagBlock*>::iterator j=jblocks.begin();j!=jblocks.end();j++)
			{
				JnxTagBlock* jblock=(*j);

				if((*i)->getPos1()>jblock->endOnTag1)
					continue;

				while(i!=selexaMatches.end() && (*i)->getPos1()<=jblock->endOnTag1)
				{
					i++;
					end++;
				}


				jblock->startIndexSelexaMatchOnTag=start;
				jblock->endIndexSelexaMatchOnTagInc=end-1;

				if(i==selexaMatches.end())
					return;
			}
		}*/

	};






	static inline GffEntry::ExonPtr getExon(int exonID)
	{
		return indexedExonome[exonID];
	}
	static inline void registerIndexedExon(ExonPtr exon)
	{
		indexedExonome.push_back(exon);
		/*if(exon->exonID==313416)
		{	cerr<<"aaaaa"<<endl;
			cerr<<exon<<endl;
			cerr<<exon->chr<<endl;
		}*/
		//cerr<<0<<"::>"<<GffEntry::indexedExonome[0]<<endl;
	}
	/**
	 * A vector to store the collection of all transcripts
	 */
	static vector<GffEntry*> transcriptome;

	static inline GffEntry* getGffEntry(int transcriptID)
	{
		return transcriptome[transcriptID];
	}


	/**
	 * sibling vectors to store all overlapping exons
	 */

	static vector< vector<ExonPtr>* > siblingVectors;

	/**
	 * All loci in genome, ordered by first exon
	 *
	 */

	//static map<ExonPtr ,Locus* > loci;
	static map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * > loci;
	
	static multimap<string, GffEntry::Locus*> nameIndexedLoci;
	
	inline static void indexLociByName()
	{
		for(map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::iterator i=loci.begin();i!=loci.end();i++)
		{
			multimap<GffEntry::ExonPtr,GffEntry::Locus*> * m2=(*i).second;
			for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator j=m2->begin();j!=m2->end();j++)
			{
				GffEntry::Locus* locus=(*j).second;


				for(set<string>::iterator k=locus->names.begin();k!=locus->names.end();k++)
				{
					nameIndexedLoci.insert(multimap<string,GffEntry::Locus*>::value_type(*k,locus));

				//cerr<<"inserting loci name:"<<(*k)<<endl;
				}
			}
		}
	}



	inline static bool isJnxInfoLoaded(string chr)
	{
		map<string,map<int,JnxTag*>* >::iterator i=jnxtags.find(chr);
		return i!=jnxtags.end();
	}

	inline static bool isJnxInfoLoaded()
	{
		return GffEntry::jnxtags.size()>0;
	}

	inline static bool isExonomeLoaded()
	{
		return GffEntry::exonome.size()>0;
	}

	inline static bool isTranscriptomeLoaded()
	{
		return GffEntry::transcriptome.size()>0;
	}

	inline static bool isGffLoaded()
	{
		return GffEntry::isExonomeLoaded() && GffEntry::isTranscriptomeLoaded();
	}

	inline static bool isLocusLoaded()
	{
		return GffEntry::loci.size()>0;
	}

	inline static bool isGBlocksLoaded()
	{
		return GffEntry::gblocks.size()>0;
	}

	inline static bool isSelexaMatchesLoaded()
	{
		return GffEntry::selexaMatchLoaded;
	}

	inline static GffEntry::Locus* findLocusByName(string name)
	{
		multimap<string, GffEntry::Locus*>::iterator i=
		GffEntry::nameIndexedLoci.find(name);

		if(i==GffEntry::nameIndexedLoci.end())
		{
			return NULL;
		}

		return (*i).second;
	}

	static bool selexaMatchLoaded;


/*	class ExonPtr
	{
	private:
		GffEntry::Exon* pEx;
	public:
		inline ExonPtr(GffEntry::Exon* _pEx=NULL):pEx(_pEx){}
		inline void operator=(GffEntry::Exon *_pEx)
		{
			pEx=_pEx;
		}

		inline operator GffEntry::Exon*() const
		{
			return pEx;
		}

		inline GffEntry::Exon& operator *() const
		{
			return *pEx;
		}


		inline GffEntry::Exon* operator ->() const
		{
			return pEx;
		}

		inline bool operator<(const ExonPtr& right) const
		{
			return (**this)<(*right);
		}
		inline bool operator>(const ExonPtr& right) const
		{
			return (**this)>(*right);
		}
		inline bool operator==(const ExonPtr& right) const
		{
			return (**this)==(*right);
		}
		inline bool operator<=(const ExonPtr& right) const
		{
			return ((**this)==(*right)) || ((**this)<(*right));
		}

		inline bool operator>=(const ExonPtr& right) const
		{
			return ((**this)==(*right))  || ((**this)>(*right));
		}

		inline bool operator!=(const ExonPtr& right) const
		{
			return !((**this)==(*right));
		}

		inline bool ptrEquals(const ExonPtr& right) const
		{
			return pEx==right.pEx;
		}

	};*/



	int transcriptID; //global ID
	int localID;
	int bin;
	string name;
	string chrom;
	char strand;
	int txStart;
	int txEnd;
	int cdsStart;
	int cdsEnd;
	int exonCount;
	string annoSource;
	bool visited;

	typedef map<ExonPtr,vector<RankedGffEntryPair>* >::iterator ExonomeWalker;

	int getExonArrayIndex(int i) const;
	int getExonRankFromIndex(int index) const;
	ExonPtr getExonOfRank(int i) const;

	/**
	 *
	 *a static function to reset all gblocks
	 *
	 */
	static void resetGBlocks(map<string, set<GBlockPtr> * >::iterator i);

	static void resetGBlocks(string chr="");

	/**
	 * a static function to reset the all exonome to nothing
	 * @see GffEntry::exonome
	 */
	static void resetExonome();

	/**
	 * a static function to reset all transcriptome
	 * @see GffEntry::transcriptome
	 */
	static void resetTranscriptome();

	/**
	 * a static function to reset everything (transcriptome and exonome)
	 * @see GffEntry::transcriptome
	 * @see GffEntry::exonome
	 * @see GffEntry::resetExonome
	 * @see GffEntry::resetTranscriptome
	 */
	static void resetEverything();


	/**
	 * a static function to reset the sibling vectors
	 */

	static void resetSiblingVectors();

	/**
	 * a static function to reset the loci
	 */
	static void resetLoci();

	/**
	 * a static function to create and register a GffEntry (transcript) into the transcriptome
	 * @param splits an vector of string generated from the splitting of a line in the Gff file
	 * @return the pointer to the GffEntry (transcript) created
	 * @see GffEntry::transcriptome
	 */
	static GffEntry* createGffEntry(vector<string> &splits,string annoSource=DEFAULT_ANNO_SOURCE,int *transcriptID=NULL);

	/**
	 * a static function to create and register an exon into exonome, and associate it with an GffEntry (transcript)
	 * @param entry the pointer to the transcript to which the exon belongs to
	 * @param chr the chromosome of the transcript/exon
	 * @param rank the exon rank in the transcript
	 * @param start the genomic start (in the genomic strand) coordinate (starts from 1) of the exon
	 * @param end the genomic end (in the genomic strand) coordinate (starts from 1) of the exon inclusive
	 * @param frame the reading frame of the exon
	 * @return the pointer to the exon created
	 * @see GffEntry::exonome
	 * @see GffEntry::createGffEntry
	 */
	static ExonPtr getNewExon(GffEntry *entry,string chr,int rank,int start,int end,int frame,char strand);

	/**
	 * A class to describe a loci (transcript group)
	 */

	


	class Locus
	{
	public:
		set<string> names;
		vector<GffEntry*> transcripts;
		set<ExonPtr> exonSet;
		ExonPtr *exons;


		//NExonGroup* root;
		map<string,NExonGroup*> exongroups;
		map<string,NESuperGroup*> leftSuperGroups;
		map<string,NESuperGroup*> rightSuperGroups;

		int exonCount;
		

		string chr;
		




		string getFirstName()
		{
			if(names.size()==0)
				return "";

			return *names.begin();
		}
		
		

		Locus():chr(""),exons(NULL)
		{

		}
		~Locus();
		
		inline void finalize()
		{

			exonCount=exonSet.size();
			KeyPair<int,int> b=getBound();
			if(exonCount>1000)
			{
				cerr<<"exonCount> 200 ec="<<exonCount<<",name="<<(*names.begin())<<",coord="<<chr<<":"<<(b.k1+1)<<"-"<<b.k2<<endl;
				//exit(0);
			}
			exons=new ExonPtr[exonCount];
			int j=0;
			for(set<ExonPtr>::iterator i=exonSet.begin();i!=exonSet.end();i++)
			{
				exons[j++]=(*i);
			}
		}

		char strand;

		inline KeyPair<int,int> getBound()
		{
			ExonPtr head=getHeadExon();
			ExonPtr tail=getTailExon();
			return KeyPair<int,int>(head->start,tail->end);
		}

		inline ExonPtr getTailExon()
		{
			if(exons)
			{
				return exons[exonCount-1];
			}

			set<ExonPtr >::reverse_iterator i=exonSet.rbegin();
			return (*i);

		}
		inline ExonPtr getHeadExon()
		{
			if(exons)
			{
				return exons[0];
			}
			set<ExonPtr >::iterator i=exonSet.begin();
			return (*i);
		}

		inline void addTranscript(GffEntry * v)
		{
			names.insert(set<string>::value_type(v->name2));
			transcripts.push_back(v);
			int exonCount=v->exonCount;
			for(int i=0;i<exonCount;i++)
			{
				ExonPtr cp=v->exons[i];
				exonSet.insert(set<ExonPtr >::value_type(cp));
			}
		}

		//vector<ExonBoundGroup*> leftExonBoundGroups;
		//vector<ExonBoundGroup*> rightExonBoundGroups;


	    void assignExonBoundGroups_FlushState(char curState, set<
	    		GffEntry::ExonPtr>* Opens, map<int,set<GffEntry::ExonPtr>* >& Closes, map<string,set<GffEntry::ExonPtr>* > &leftBoundExons,map<string,set<GffEntry::ExonPtr>* > &rightBoundExons);

	    void assignExonBoundGroups_FlushTerminalState(char flusheeState, set<
	    		GffEntry::ExonPtr>* TerminalOpens, map<int,set<GffEntry::ExonPtr>* >& TerminalCloses, map<string,set<GffEntry::ExonPtr>* > &terminalLeftBoundExons,map<string,set<GffEntry::ExonPtr>* > &terminalRightBoundExons);
		void assignExonBoundGroups(int jobID,bool hasTerminal=false);


		void printExonBounds(ostream& os);

		friend ostream& operator<<(ostream& os,  GffEntry::Locus& loc);
		friend istream& operator>>(istream& os,  GffEntry::Locus& loc);
	};


	/**
	 *
	 * A class to represent a selexaMapped to a GBlock
	 *
	 */

	class SelexaMatch
	{
	public:
		




		inline static bool readMapviewMatch (istream&is, GffEntry::SelexaMatch& match)
		{
			string istr; //dummy string
			int iint;  //dummy int
			is>>istr; //1) readname (ignored)
			is>>match.chr; //2) chromosome;
			is>>match.gpos0; //3) this is gpos1
			match.gpos0--; //convert to gpos0;
			is>>match.strand; //4) strand
			is>>iint; //5) insert size from the outer coordinates of a pair (ignored)
			is>>istr; //6) paired flag (ignored)
			is>>match.quality; //7) quality
			is>>iint; //8) single mapping quality (ignored)
			is>>iint; //9) alternative mapping quality (ignored)
			is>>match.nMismatches; //10) number of mismatches of best hit, assumed to be the one
			is>>iint; //11) sum of qualities mismatched bases of the best hit (ignored)
			is>>iint; //12) sum of 0-mismatch hits of the first 24bp (ignored)
			is>>iint; //13) sum of 1-mismatch hits of the first 24bp (ignored)
			is>>iint; //14) read length (ignored) known
			is>>istr; //15) read sequence (ignored)
			is>>istr; //16) read quality (ignored)

			return true;
		}





		inline static int getSamStrand(int flag)
		{
			if((flag & 0x0010)==0)
				return GffEntry::FORWARD;
			else
				return GffEntry::REVERSE;
		}

		inline static void readOptFlag(map<string,pair<char,string> >&optflags,vector<string>&fields,int start)
		{
			vector<string> splits;
			for( unsigned int i=start;i<fields.size();i++)
			{
				StringUtil::split(fields[i],":",splits);
				if(splits.size()!=3)
				{
					cerr<<"strange error: opt flag is not properly formatted "<<fields[i]<<endl;
					continue;
				}

				optflags.insert(map<string,pair<char,string> >::value_type(fields[0],pair<char,string>(fields[1][0],fields[2])));

			}
		}

		inline static bool readSamMatch(istream& is,GffEntry::SelexaMatch& match,int  minMapQuality=0,int maxMismatches=INT_MAX, bool ifReadOptFlag=true) //return whther a match is read or not
		{

			string line=consumeStreamUntilNewLine(is);
			if (line[0]=='@')
			{
				//header directive is ignored for now
				match.chr="@";
				return false;
			}


			vector<string> fields;
			StringUtil::split(line,"\t",fields);

			if(fields.size()<14)
			{
				//cerr<<"Error: Not Enough fields "<<line<<endl;
				return false;
			}

			//0:QNAME/ readname is ignored.
			//1:FLAG bitwiseflag (where you get strand info)
			int flag=StringUtil::atoi(fields[1]);
			match.strand=getSamStrand(flag);

			//2:ref
			match.chr=fields[2]; ///ref
			//3:POS at 1-based
			match.gpos0=StringUtil::atoi(fields[3]); //now is 1--based
			match.gpos0--; //make 0-based
			//4:MAP QUALITY
			match.quality=StringUtil::atoi(fields[4]);
			if (match.quality<minMapQuality)
				return false;

			//5:CIGAR STRING,ignored for now
			//6:MATE reference,ignored for now
			//7:MATE POS, ignored for now
			//8:ISIZE insert size, ignored for now
			//9:SEQ: ingnored
			//10: query qulity,ignored
			//11:TAG ignored
			//12:VTYPE ignored
			//13:VALUE ignored
			match.nMismatches=0;
			map<string,pair<char,string> > optflags;
			readOptFlag(optflags,fields,11);
			map<string,pair<char,string> >::iterator nmi=optflags.find("NM");
			if(nmi!=optflags.end())
			{
				match.nMismatches=StringUtil::atoi(nmi->second.second);
			}

			if(match.nMismatches>maxMismatches)
				return false;

			return true;
		}
		
		
		inline static bool readMatchFromBam(samfile_t *fp,bam1_t* b,GffEntry::SelexaMatch& match,int  minMapQuality=0,int maxMismatches=INT_MAX, bool ifReadOptFlag=true)
		{
			
			
			
			uint32_t *cigar = bam1_cigar(b);
			const bam1_core_t *c = &b->core;
			int i, l;
			if (b->core.tid < 0) return 0;
			for (i = l = 0; i < c->n_cigar; ++i) {
				int op = cigar[i]&0xf;
				if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
					l += cigar[i]>>4;
			}
			
			//printf("%s\t%d\t%d\t%s\t%d\t%c\n", fp->header->target_name[c->tid],
			//	   c->pos, c->pos + l, bam1_qname(b), c->qual, (c->flag&BAM_FREVERSE)? '-' : '+');
			
			
			
			match.strand=((c->flag&BAM_FREVERSE)? '-' : '+');
			match.chr=fp->header->target_name[c->tid]; ///ref
			match.gpos0=c->pos; 
			match.quality=c->qual;
			
			if (match.quality<minMapQuality)
				return false;
			

			match.nMismatches=0;  //update this!!!
			//map<string,pair<char,string> > optflags;
			//readOptFlag(optflags,fields,11);
			//map<string,pair<char,string> >::iterator nmi=optflags.find("NM");
			//if(nmi!=optflags.end())
			//{
			//	match.nMismatches=StringUtil::atoi(nmi->second.second);
			//}
			//
			if(maxMismatches!=INT_MAX){
				cerr<<"warning: mismatch not implemented for BAM reading here"<<endl;
				
			}
			
			
			if(match.nMismatches>maxMismatches)
				return false;
			
			return true;
			
		}
		
		

		inline static bool readMatch(istream& is,GffEntry::SelexaMatch& match, int format)
		{
			switch(format){
			case MAQ_FORMAT:
				return readMapviewMatch(is,match);
			case SAM_FORMAT:
				return readSamMatch(is,match);

			}

			return false;
		}

		string chr;
		char strand;
		int gpos0; //genomic position0
		//int length; known,right?
		int *mispos; //mismatch positions (relative)
		char *misbases; //mismatchTo bases (relative)
		int nMismatches; //how many mismatches
		int quality; //quality of mapping
		inline int getPos0()
		{
			return gpos0;
		}
		inline int getPos1()
		{
			return gpos0+1;
		}
		inline SelexaMatch():mispos(NULL),misbases(NULL),nMismatches(0)
		{

		}
		inline ~SelexaMatch()
		{
			if(nMismatches>0 && mispos && misbases)
			{
				delete[] mispos;
				delete[] misbases;
			}
		}
		inline bool operator<(const SelexaMatch& right) const
		{
			if(chr==right.chr)
			{
				return gpos0<right.gpos0;
			}else
				return chr<right.chr;
		}
		inline bool operator>(const SelexaMatch& right) const
		{
			if(chr==right.chr)
			{
				return gpos0>right.gpos0;
			}else
				return chr>right.chr;
		}
		inline bool operator==(const SelexaMatch& right) const
		{
			return chr==right.chr && gpos0==right.gpos0;
		}
		inline bool operator!=(const SelexaMatch& right) const
		{
			return !(*this==right);
		}
		inline bool operator>=(const SelexaMatch& right) const
		{
			return !(*this<right);
		}
		inline bool operator<=(const SelexaMatch& right) const
		{
			return !(*this>right);
		}


	};


	typedef vector<GffEntry::GBlockPtr>::iterator BlockI;
	typedef pair<BlockI,BlockI> DBlockI;

	/**
	 *
	 * A class to remember the smallest blocks of the genome divided by elements (exons, introns, intergenics)
	 *
	 */

	class GBlock
	{
	public:
		string chr;
		int start0;
		int end1;
		char strand;
		int ID;
		string name;

		set<ExonPtr> exons;

		int umposL;
		int umposR;

		BoundSet uniquelyMappableBounds;

		int usageFrequency;

		multiset<GffEntry::SelexaMatchPtr> selexaMatches;

		typedef multiset<SelexaMatchPtr>::iterator SIterator;
		typedef pair<SIterator,SIterator> DIterator;


		KeyPair<int,int> getDensity(int readLength,bool useUniquelyMappablePosition,bool protrude=true,bool ambigCheck=true);

		inline KeyPair<int,int>  getBound(){ //01

			return KeyPair<int,int>(start0,end1);
		}
		
		template<class iteratorClass,class reverseIteratorClass>
		inline static KeyPair<int,int> getDensityOfContigBlocks(iteratorClass from,iteratorClass toEx,int readLength,bool letReadsToProtrudeLastBlock=false,bool useUniquelyMappablePosition=true,bool ambigCheck=true)
		{
			KeyPair<int,int> returnVal(0,0);

			if(toEx==from)
			{
				cerr<<"strange error: no blocks included"<<endl;
				return returnVal;
			}

			reverseIteratorClass rFrom(toEx);
		//	rFrom--;
			reverseIteratorClass rToEx(from);
		//	rToEx--;

			//cerr<<"a before getting density"<<endl;
			GBlockPtr curBlock=*rFrom;
			//cerr<<"b getting density  for "<<curBlock->getBound()<<endl;
			addDensityVectors(returnVal,curBlock->getDensity(readLength,useUniquelyMappablePosition,letReadsToProtrudeLastBlock,ambigCheck));
			//cerr<<"c getting density  for "<<curBlock->getBound()<<endl;
			rFrom++;

			for(;rFrom!=rToEx;rFrom++)
			{
			//	cerr<<"a before getting density"<<endl;
				curBlock=*rFrom;
			//	cerr<<"b getting density  for "<<curBlock->getBound()<<endl;
				addDensityVectors(returnVal,curBlock->getDensity(readLength,useUniquelyMappablePosition,true,ambigCheck));
			//	cerr<<"c getting density  for "<<curBlock->getBound()<<endl;
			}

			return returnVal;

		}

		template<class iteratorClass,class reverseIteratorClass>
		inline static KeyPair<int,int> getDensityOfContigBlocks(pair<iteratorClass,iteratorClass> blockRange, int readLength, bool letReadsToProtrudeLastBlock=false,bool useUniquelyMappablePosition=true,bool ambigCheck=true)
		{
			return getDensityOfContigBlocks<iteratorClass,reverseIteratorClass>(blockRange.first,blockRange.second,readLength,letReadsToProtrudeLastBlock,useUniquelyMappablePosition,ambigCheck);

		}

		inline int getSelexaMatchSize() const
		{
			return selexaMatches.size();
		}
		inline bool isExonic() const
		{
			return this->noOfExons()>0;
		}
		inline int noOfExons() const
		{
			return exons.size();
		}
		
		inline KeyPair<int,int> getBound11()
		{
			return KeyPair<int,int>(getStart1(),getEnd1());
		}

		inline KeyPair<int,int> getBound00OutsideRelPosForReadLength(int readLength)
		{
			int len=getLength();
			int ppp=len-readLength;
			return KeyPair<int,int>(MAX(0,ppp+1),len-1);
		}

		inline KeyPair<int,int> getBound00WithinRelPosForReadLength(int readLength)
		{
			int ppp=getLength()-readLength;
			return KeyPair<int,int>(0,ppp);
		}

		inline KeyPair<int,int> getBound00WithinGenomicPosForReadLength(int readLength)
		{
			int ppp=getLength()-readLength;
			if(ppp<=0)
			{
				return KeyPair<int,int>(-1,-2);
			}
			return KeyPair<int,int>(this->start0,this->start0+ppp);
		}

		inline KeyPair<int,int> getAbs11FromRel00(const KeyPair<int,int>& rel00)
		{
			if(!isBound00Valid(rel00))
				return KeyPair<int,int>(100,-100);
				
			return KeyPair<int,int>(start0+rel00.k1+1,start0+rel00.k2+1);
		}

		inline KeyPair<int,int> getBoundAbs11WithinPosForReadLength(int readLength)
		{
			return getAbs11FromRel00(getBound00WithinRelPosForReadLength(readLength));
		}

		inline KeyPair<int,int> getboundAbs11OutsidePosForReadLength(int readLength)
		{
			return getAbs11FromRel00(getBound00OutsideRelPosForReadLength(readLength));
		}
		
		
		inline DIterator getSelexaMatchesWithinGenomicBound00(int astart0,int aend0)
		{

			if(!::isBound00Valid(KeyPair<int,int>(astart0,aend0)))
			{
				//cerr<<"error: get SelexaMatch Within Genomic Bound 00 astart0 > aend0"<<endl;
				//die_exit("");
				return DIterator(selexaMatches.end(),selexaMatches.end());
			}

			SelexaMatch sm;
			sm.chr=chr;
			sm.gpos0=astart0;

			SIterator low=selexaMatches.lower_bound(&sm);

			sm.gpos0=aend0;
			SIterator high=selexaMatches.upper_bound(&sm);

			return DIterator(low,high);
		}

		inline DIterator getSelexaMatchesWithStartWithinRelPos(const KeyPair<int,int>& bound00)
		{
			return getSelexaMatchesWithStartWithinRelPos(bound00.k1,bound00.k2);
		}

		inline DIterator getSelexaMatchesWithStartWithinRelPos(int rstart0,int rend0)
		{
			if(rstart0>rend0)
			{
				return DIterator(selexaMatches.end(),selexaMatches.end());
			}

			int astart0=start0+rstart0;
			int aend0=start0+rend0;

			return getSelexaMatchesWithinGenomicBound00(astart0,aend0);
		}


		
		inline int getEnd1()
		{
			return end1;
		}
		
		inline int getStart1()
		{
			return start0+1;
		}


		inline int getStart0()
		{
			return start0;
		}
		
		inline int getEnd0()
		{
			return end1-1;
		}

		inline bool getNaiveAmbiguity()
		{
			return strand==GffEntry::BOTH_STRANDS;
		}
		
		
		inline int getLength()
		{
			return end1-start0;
		}
		
		inline int getUniquelyMappablePosL()
		{
			return this->umposL;
		}

		inline int getUniquelyMappablePosR()
		{
			return this->umposR;
		}

		inline int getNaiveMappablePos(int readLength=1)
		{
			return MAX(0,getLength()-readLength+1);
		}

		inline int getRelativePos(int gpos0)
		{
			return gpos0-start0;
		}

		inline int getRelativePos(SelexaMatchPtr pSm)
		{
			return getRelativePos(pSm->gpos0);
		}

		inline int getUniquelyMappablePositions(const KeyPair<int,int>& _bound)
		{
			//everything is still genomic coordinates
			return this->uniquelyMappableBounds.getOverlapLength(_bound);
		}


		inline ~GBlock()
		{
			if(selexaMatches.size()>0)
			{
				for(multiset<SelexaMatchPtr>::iterator i=selexaMatches.begin();i!=selexaMatches.end();i++)
				{
					delete *i;
				}
			}
		}

		inline GBlock(const string& _chr="",int _start0=0,int _end1=0)
			:chr(_chr),start0(_start0),end1(_end1),usageFrequency(0),umposL(0),umposR(0){}


		inline bool operator > (const GBlock& block) const
		{
			if(chr==block.chr)
			{
				if(start0==block.start0)
				{
				//	cerr<<"strange error: block starts are the same: "<<chr<<":"<<start0<<"-"<<end1<<endl;
				//	cerr<<"\t\t\t\t\t\t: "<<block.chr<<":"<<block.start0<<"-"<<block.end1<<endl;
					return end1>block.end1;
				}else
					return start0>block.start0;

			}
			else
				return chr>block.chr;
		}
		inline bool operator < (const GBlock& block) const
		{
			if(chr==block.chr)
			{
				if(start0==block.start0)
				{
				//	cerr<<"strange error: block starts are the same: "<<chr<<":"<<start0<<"-"<<end1<<endl;
				//	cerr<<"\t\t\t\t\t\t: "<<block.chr<<":"<<block.start0<<"-"<<block.end1<<endl;

					return end1<block.end1;
				}else
					return start0<block.start0;

			}
			else
				return chr<block.chr;
		}
		inline bool operator == (const GBlock& block) const
		{
			return (chr==block.chr && start0==block.start0 && end1==block.end1);
		}
		
		inline bool operator !=(const GBlock& block) const
		{
			return !(*this==block);
		}
		inline bool operator >=(const GBlock& block) const
		{
			return !(*this<block);
		}
		inline bool operator <=(const GBlock& block) const
		{
			return !(*this>block);
		}


	};

	//class ExonBoundGroup;



	/**
	 * A class to describe an exon recognized by chr, start, end
	 *
	 */

	class Exon
	{
	public:


		bool operator <(const Exon &right) const;
		bool operator >(const Exon &right) const;
		bool operator ==(const Exon& right) const;
		bool operator !=(const Exon& right) const;

		bool operator >=(const Exon& right) const;

		bool operator <=(const Exon& right) const;
		inline int getLength() const
		{
			return end-start; //end is 1based, start is 0based;
		}
		string chr;
		int rank;
		int start;
		int end;
		int frame;
		int exonID;
		bool visited;
		char strand;
		int siblingID;

		int maxCodingStart;
		int minCodingEnd;



		int inDegree;

		int usageFlag;
		/*
		 * reference to all transcripts with this exon
		 */
		vector<RankedGffEntryPair>* assTranscripts;
		vector<GffEntry::GBlockPtr>* blocks;
		map<ExonPtr,Jnx*>* outJnxs;

		string leftBoundID;
		string rightBoundID;


		//ExonBoundGroup* leftExonBoundGroup;
		//ExonBoundGroup* rightExonBoundGroup;


		typedef map<ExonPtr,Jnx*>::iterator OutJnxI;



		NExonGroup *exonGroup;
		string egsid;

		KeyPair<int,int> getDensity(int readLength,bool useUniquelyMappablePosition=true, bool ambigCheck=true);

		inline GffEntry::GBlockPtr getGenomicFirstBlock()
		{
			return *(blocks->begin());

		}

		inline GffEntry::GBlockPtr getGenomicLastBlock()
		{
			return *(blocks->rbegin());
		}

		template<class containerClass,class iteratorClass>
		inline static pair<iteratorClass,iteratorClass> getBlocksInRangeGeneric(/*const*/ containerClass& container,KeyPair<int,int> bound, bool completeIn=true)
		{
			pair<iteratorClass,iteratorClass> rangeI;

			rangeI.first=container.begin();
			rangeI.second=container.end();

			bool started=false;

			iteratorClass i=container.begin();

			for(;i!=container.end();i++)
			{
				if(completeIn)
				{
					if(::isAinB((*i)->getBound(),bound))
					{
						if(!started)
						{
							rangeI.first=i;
							started=true;
						//	cerr<<"set start at "<<(*i)->getBound()<<" for bound "<<bound<<endl;
						}
					}
					else
					{
						if(started)
						{
							rangeI.second=i;
						//	cerr<<"set end at "<<(*i)->getBound()<<" for bound "<<bound<<endl;
							return rangeI;
						}
					}
				}
				else
				{
					if(::areOverlapping((*i)->getBound(),bound))
					{
						if(!started)
						{
							rangeI.first=i;
							started=true;
						}
					}else
					{
						if(started)
						{
							rangeI.second=i;
							return rangeI;
						}
					}


				}
			}

		//	cerr<<"end get blocks"<<endl;
			rangeI.second=i;
			return rangeI;
		}




		inline DBlockI getBlocksInRange(KeyPair<int,int> bound,bool completeIn=true)
		{



			DBlockI rangeI;

			if(!this->blocks)
			{
				cerr<<"strange error die: exon with no blocks!?"<<endl;
				return rangeI;
			}

			rangeI.first=blocks->end();
			rangeI.second=blocks->end();

			bool started=false;

			vector<GffEntry::GBlockPtr>::iterator i=blocks->begin();

			for(;i!=blocks->end();i++)
			{
				if(completeIn)
				{
					if(::isAinB((*i)->getBound(),bound))
					{
						if(!started)
						{
							rangeI.first=i;
							started=true;
					//		cerr<<"set start at "<<(*i)->getBound()<<" for bound "<<bound<<endl;
						}
					}
					else
					{
						if(started)
						{
							rangeI.second=i;
					//		cerr<<"set end at "<<(*i)->getBound()<<" for bound "<<bound<<endl;
							return rangeI;
						}
					}
				}
				else
				{
					if(::areOverlapping((*i)->getBound(),bound))
					{
						if(!started)
						{
							rangeI.first=i;
							started=true;
						}
					}else
					{
						if(started)
						{
							rangeI.second=i;
							return rangeI;
						}
					}


				}
			}

		//	cerr<<"end get blocks"<<endl;
			rangeI.second=i;
			return rangeI;
		}





		/*
		 * reference to all overlapping exons
		 */
		vector<ExonPtr> *siblings;

		inline int getEnd0() const
		{
			return end-1;
		}

		inline int getStart0() const
		{
			return start;
		}
		inline int getStart1() const
		{
			return start+1;
		}
		inline int getEnd1() const
		{
			return end;
		}

		inline KeyPair<int,int> getBound() const
		{
			return KeyPair<int,int>(this->getStart0(),this->getEnd1());

		}

		Jnx* addNewJnxBlock(ExonPtr _rightExon,JnxTagBlock* _jblocks, char _strand);


		inline int getInDegree() const
		{
			return this->inDegree;
		}

		inline int getOutDegree() const
		{
			if(!this->outJnxs)
				return 0;

			return this->outJnxs->size();
		}

		inline bool hasOutJnxs() const
		{
			return getOutDegree()>0;//this->outJnxs;
		}

		inline bool hasInJnxs() const
		{
			return this->getInDegree()>0;
		}

		inline bool hasInOrOut() const
		{
			return this->hasInJnxs() || this->hasOutJnxs();//this->getInDegree()>0 || this->getOutDegree()>0;
		}

		inline bool hasBothInAndOut() const
		{
			return this->hasInJnxs() && this->hasOutJnxs();
		}

		~Exon();
		Exon();
		Exon(string _chr,int _rank,int _start,int _end,int _frame,int _exonID);
		void set(string _chr,int _rank,int _start,int _end,int _frame,int _exonID);
	};

	/*
	 * class Loader.
	 * A class to assist the loading of the GffFiles, single or in sequence
	 *
	 */

	class Loader
	{
	private:
		TranscriptLoadingListener *listener;
	public:
		void resetEverything();
		void loadGffFile(string filename,string annoSource=DEFAULT_ANNO_SOURCE,bool reset=true);
		void loadLoci(string filename);
		void loadGBlocks(string filename);
		void loadUniquelyMappablePosition(int readLength,string chrMapRefFile,string prefix,string suffix=COMPACTPOSITIONSUFFIX,bool toGBlock=true,bool toJnxBlock=true);
		void loadSelexaMatches(string filename,int format);
		void loadJnxInfo(string filename,int minSpan,int readLength);
		void loadJnxSelexaMatches(string filename,int format);
		void loadSiblings(string filename);
		vector<KeyPair<string,string> > loadGffFiles(string sourceFile,bool reset=true, string prefix="");
		Loader(TranscriptLoadingListener * _listener=NULL);
	};

	/*
	 * class TranscriptLoadingListener
	 * The listener (abstract) class for each line of transcript being loaded in GffFile
	 *
	 */

	class TranscriptLoadingListener
	{
	public:
		virtual void transcriptLoaded(GffEntry& entry,string annoSource)=0;
		virtual ~TranscriptLoadingListener(){}
	};

	ExonPtr *exons;
	string id;
	string name2;
	string cdsStartStat;
	string cdsEndStat;

	static const char FORWARD;
	static const char REVERSE;
	static const char UNKNOWN;
	static const char BOTH_STRANDS;

	
	static multiset<GffEntry::SelexaMatchPtr> blindSpotSelexaMatches;
	
	/*
	 * get first exon of this transcript on the genomic strand
	 */
	inline GffEntry::ExonPtr getGenomicFirstExon()
	{
		return this->exons[0];
	}
	/*
	 * get last exon of this transcript on the genomic strand
	 */
	inline GffEntry::ExonPtr getGenomicLastExon()
	{
		return this->exons[this->exonCount-1];
	}
	
	~GffEntry();
};




#define TRAVEL_DOWNTO_LEAVES INT_MAX
#define TRAVEL_DOWNONE_LEVEL 2
#define DONT_TRAVEL_DOWN 1
#define ONLY_LEAVEL_EXONS DONT_TRAVEL_DOWN

class NExonGroup;


#define ExonBoundGroup_LEFT  'o'
#define ExonBoundGroup_RIGHT  'c'
#define ExonBoundGroup_TLEFT  'O'
#define ExonBoundGroup_TRIGHT  'C'
/*class ExonBoundGroup: public set<GffEntry::ExonPtr>
{
public:
	KeyPair<int,int> coordRange; //11
	char mode;
	string name;

	inline string getID()
	{
		return name;
	}

	ExonBoundGroup(string _name,char _mode):mode(_mode),name(_name),coordRange(INT_MAX,INT_MIN)
	{}

	bool operator < (const ExonBoundGroup& right) const
	{
		if(coordRange.k1==right.coordRange.k1)
		{
			if(coordRange.k2==right.coordRange.k2)
			{
				if(mode==right.mode)
				{
					return name < right.name;
				}
				else
				{
					return mode<right.mode;
				}
			}
			else
			{
				return coordRange.k2<right.coordRange.k2;
			}
		}else
		{
			return coordRange.k1<right.coordRange.k1;
		}
	}

	bool operator > (const ExonBoundGroup& right) const
	{
		if(coordRange.k1==right.coordRange.k1)
		{
			if(coordRange.k2==right.coordRange.k2)
			{
				if(mode==right.mode)
				{
					return name > right.name;
				}
				else
				{
					return mode>right.mode;
				}
			}
			else
			{
				return coordRange.k2>right.coordRange.k2;
			}
		}else
		{
			return coordRange.k1>right.coordRange.k1;
		}
	}


};*/


class NExonGroup;
class NESuperGroup;


typedef SmartPtr<NExonGroup> NExonGroupPtr;
typedef SmartPtr<NESuperGroup> NESuperGroupPtr;

class NESuperGroup: public set<NExonGroupPtr>
{
public:
	typedef NExonGroup* NExonGroupP;
	//KeyPair<int,int> bound; //bound is the intersection of all member exon groups
	string id;

	 NESuperGroup(const string& _id);


	 void add(NExonGroupPtr exongroup);


	string getID() const;

	inline bool operator < (const NESuperGroup& sgroup) const
	{
		return this->getID()<sgroup.getID();
	}
	
	inline bool operator > (const NESuperGroup& sgroup) const
	{
		return this->getID()>sgroup.getID();
	}
	
	inline bool operator==(const NESuperGroup& sgroup) const
	{
		return this->getID()==sgroup.getID();
	}

	inline bool operator!=(const NESuperGroup& sgroup) const
	{
		return !(*this==sgroup);
	}
	
	inline bool operator<=(const NESuperGroup& sgroup) const
	{
		return !(*this>sgroup);
	}
	inline bool operator>=(const NESuperGroup& sgroup) const
	{
		return !(*this<sgroup);
	}
	

};


//typedef NESuperGroup* NESuperGroupPtr;

class NExonGroup
{
public:
	typedef SmartPtr<NExonGroup> NExonGroupPtr;
	set<GffEntry::ExonPtr> levelExons;
	//set<NExonGroupPtr>* children;
	string locusName;
	//NExonGroup::NExonGroupPtr relRoot;
	//string relRootString;
	//bool isExtremeRight;
	//bool isExtremeLeft;
	NESuperGroupPtr leftSuperGroup;
	NESuperGroupPtr rightSuperGroup;
	GffEntry::Locus *locus;

	//set<NExonGroupPtr>  *parents; //set when finalized.


	/*void becomeChildOf(NExonGroupPtr _parent)
	{
		if(!_parent->children)
		{
			_parent->children=new set<NExonGroupPtr>;

		}

		_parent->children->insert(this);

		if(!this->parents)
		{
			this->parents=new set<NExonGroupPtr>;
		}
	}*/


	string sid;

	//KeyPair<int,int> bound;

	KeyPair<int,int> levelExonBound;
	//int depth;


	typedef set<GffEntry::ExonPtr>::iterator ExonI;
	//typedef set<NExonGroupPtr>::iterator ChildI;
	typedef set<GffEntry::ExonPtr>::reverse_iterator RExonI;
	//typedef set<NExonGroupPtr>::reverse_iterator RChildI;

	typedef multimap<int,GffEntry::ExonPtr> IntExonPtrMMap;
	typedef IntExonPtrMMap::iterator IntExonPtrMMapI;
	typedef pair<IntExonPtrMMapI,IntExonPtrMMapI> IntExonPtrMMapDI;
	
	//copy constructor to handle the root condition
	/*NExonGroup(const NExonGroup& neg):  leftSuperGroup(NULL),rightSuperGroup(NULL), relRootString(neg.relRootString), children(neg.children) , locusName(neg.locusName), sid(neg.sid), bound(neg.bound), levelExons(neg.levelExons),levelExonBound(INT_MAX,INT_MIN),parent(NULL)
	{
		//now update membership of levelExons
		for(set<GffEntry::ExonPtr>::iterator i=levelExons.begin();i!=levelExons.end();i++)
			(*i)->exonGroup=this;
	}*/
	

	void getBoundariesMap(IntExonPtrMMap *uBoundMap,IntExonPtrMMap*dBoundMap);

	//void getSuperGroups(set<NESuperGroupPtr> *leftSGroups,set<NESuperGroupPtr> *rightSGroups);

	/*inline KeyPair<int,int> getBound() const
	{
		return bound;
	}	*/

	inline NExonGroup(): /*isExtremeRight(false), isExtreamLeft(false), */ leftSuperGroup(NULL),rightSuperGroup(NULL),levelExonBound(INT_MAX,INT_MIN) {}

	/*inline bool isLeave()
	{
		return !hasChildren();
	}

	inline NExonGroupPtr getLeftMostChild()
	{
		if(!this->children || this->children->size()<2)
		{
			return NULL;
		}

		return *this->children->begin();

	}
	inline NExonGroupPtr getRightMostChild()
	{
		if(!this->children || this->children->size()<2)
		{
			return NULL;
		}

		return *this->children->rbegin();
	}

	inline void expandBoundBy(GffEntry::ExonPtr exon)
	{
		expandBound(this->bound,exon->getBound());
	}
	inline void expandBoundBy(NExonGroup* exongroup)
	{
		expandBound(this->bound,exongroup->getBound());
	}*/

	inline ExonI insertExon(GffEntry::ExonPtr exon)
	{
		////cerr<<"insertExon "<<exon<<endl;
		pair<ExonI,bool> insertStat=levelExons.insert(exon);

		if(insertStat.second==false)
			return insertStat.first;

		exon->exonGroup=this;

		exon->egsid=this->sid+"."+StringUtil::str((unsigned int)levelExons.size());
		//this->expandBoundBy(exon);
		expandBound(this->levelExonBound,exon->getBound());

		return insertStat.first;
	}
	/*inline ChildI insertChild(NExonGroup* child)
	{
		if(!children)
			children=new set<NExonGroupPtr>;
		this->expandBoundBy(child);
		pair<ChildI,bool> insertStat=children->insert(child);
		return insertStat.first;
	}
	
	inline void insertChildren(ChildI first, ChildI last1)
	{
		for(ChildI i=first;i!=last1;i++)
		{
			this->insertChild(*i);
		}
	}
	
	inline void eraseChild(const NExonGroupPtr& child)
	{
		if(!hasChildren())
			return;
		
		this->children->erase(child);
	}
	
	inline void eraseChildren(ChildI first,ChildI last1)
	{
		if(!hasChildren())
			return;
		//just erase, nothing to do else? no need to shrink bound?
		this->children->erase(first,last1);
	}

	inline bool hasChildren()
	{
		return children;
	}
	
	void destroyTree()
	{
		if(children)
		{
			for(ChildI i=children->begin();i!=children->end();i++)
			{
				NExonGroupPtr eg=(*i);
				eg->destroyTree();
				delete eg;
			}
		}
	}*/

	inline ~NExonGroup()
	{
		/*delete children;

		if(leftSuperGroup && leftSuperGroup->at(0)==this)
		{
			delete leftSuperGroup;
		}

		if(rightSuperGroup && rightSuperGroup->at(0)==this)
		{
			delete rightSuperGroup;
		}*/

		//let locus destroy the left and right supergruops

	}
	
	/*inline NExonGroup(GffEntry::ExonPtr firstExon):  leftSuperGroup(NULL),rightSuperGroup(NULL),relRootString(""), children(NULL), bound(INT_MAX,INT_MIN) ,levelExonBound(INT_MAX,INT_MIN), parent(NULL)
	{
		insertExon(firstExon);
	}*/
	
	bool operator < (const NExonGroup& right) const
	{
		//NExonGroup* rootPtr=root;
		//NExonGroup* rightPtr=right.root;

		//if(rootPtr==rightPtr) //don't even need this because by sorting bounds and the fact that ExonGroup doesn't overlap, it is automatically sorted by exonGroup too.
			//return this->bound<right.bound; //isSmaller(this->bound,right.bound);
		//else
		//	return rootPtr<rightPtr; //just for the purpose of segregating entries by root exon group ptr, so ordering by address is ok

		if(this->levelExonBound==right.levelExonBound)
		{
			return this->sid<right.sid;
		}
		else
		{
			return this->levelExonBound<right.levelExonBound;
		}

	}
	
	bool operator == (const NExonGroup& right) const
	{
		return this->levelExonBound==right.levelExonBound && this->sid==right.sid; //they are bound to have the same root, because they are actually identical!
	}
	
	bool operator > (const NExonGroup& right) const
	{
		//NExonGroup* rootPtr=root;
		//NExonGRoup* rightPtr=right.root;
		
		//if(rootPtr==rightPtr)
			//return this->bound>right.bound;
		//else
		//	return rootPtr>rightPtr; //just for the purpose of segregating entries by root exon group ptr, so ordering by address is ok
		
		if(this->levelExonBound==right.levelExonBound)
		{
			return this->sid>right.sid;
		}
		else
		{
			return this->levelExonBound>right.levelExonBound;
		}

	}
	bool operator !=(const NExonGroup& right) const
	{
		//if(!root.ptrEquals(right))
		//	return false;

		return !(*this==right);
	}

	bool operator >=(const NExonGroup& right) const
	{
		return !(*this<right);
	}
	bool operator <=(const NExonGroup& right) const
	{
		return !(*this>right);
	}



};


typedef SmartPtr<NExonGroup> NExonGroupPtr;


istream& operator>>(istream& is, GffEntry::GBlock& block);

ostream& operator<<(ostream& os, const GffEntry::Locus& loc);

istream& operator >> (istream&is, GffEntry::SelexaMatch& match);

#define RJT_START 0
#define RJT_END 1
#define RJT_BOTHSTARTEND 2
#define __DLINE  //cerr<<__LINE__<<endl;

inline void _readJnxTag(istream&is, GffEntry::JnxTag& jnxtag, int minSpan,int readLength)
{
	string tstr;

	cerr<<endl;
	
	__DLINE
	is>>tstr;
	jnxtag.id.set(tstr);
	__DLINE
	is>>tstr; //chr, ignored: already gotten from id;
	
	//is>>boundSize;
	is>>jnxtag.strand;
	is>>jnxtag.lengthTag;
	is>>jnxtag.boundSize;
	is>>jnxtag.exonPathSize;
	jnxtag.allocBoundsAndExonPaths();
	__DLINE

	vector<int> boundsCoord;

	//now make JnxTagBlocks;
	//it's in 11 land!!
	//{

	KeyPair<int,int> *lifeSpansExons=new KeyPair<int,int>[jnxtag.boundSize];
	
	int lenInc=0;
	int lenExc=0;
	
	
	KeyPair<int,int> prevBo(INT_MIN,INT_MIN);
	
	__DLINE
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		is>>tstr;
		StringUtil::splitInt2(tstr,"-",boundsCoord,true);
		if(boundsCoord.size()!=2)
			continue;
		jnxtag.bounds[i].k1=boundsCoord[0]; //11
		jnxtag.bounds[i].k2=boundsCoord[1]; //11
		if(boundsCoord[0]<prevBo.k1 || boundsCoord[0]<prevBo.k2)
			cerr<<"WRONGB:"<<boundsCoord[0]<<"-"<<boundsCoord[1]<<endl;
		prevBo=jnxtag.bounds[i];
		
		
	}

	cerr<<"******* " <<jnxtag.id.chr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<jnxtag.bounds[i].k1<<"\t";
	}
	
	cerr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<jnxtag.bounds[i].k2<<"\t";
	}
	cerr<<endl;
	cerr<<"*******"<<endl;
	
	int maxPosition1GivenTagLength=jnxtag.lengthTag-readLength+1;
	__DLINE

	
	//now calculate the lifespan of each blocks in jnxtag;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		int thisLen=len11(jnxtag.bounds[i]);
		lenInc+=thisLen;
		int calLifeSpansLowBound=lenExc+minSpan-readLength+1;
		lifeSpansExons[i].k1=calLifeSpansLowBound-1; //1 before /***/
		lifeSpansExons[i].k2=lenInc-minSpan+1;				//1
		lenExc+=thisLen;
	}
	
	cerr<<"*******"<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<lifeSpansExons[i].k1<<"\t";
	}
	
	cerr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<lifeSpansExons[i].k2<<"\t";
		if(lifeSpansExons[i].k2<0)
			cerr<<"WRONGAAA"<<endl;
	}	
	cerr<<endl;
	cerr<<"*******"<<endl;
	
	__DLINE
	int iS=1;
	int iE=0;
	int nS=jnxtag.boundSize;
	int nE=jnxtag.boundSize-1;

	int prevB1=1;
	int flag;
	int coord1b4;
	
	int blockSpanStart;
	int blockSpanEndInc;
	__DLINE
	while(iE<nE && lifeSpansExons[iE].k2<2)
	{
		iE++;
	}
	__DLINE
	if(iE<nE)
	{
		iS=iE+1;
		blockSpanStart=iE;
		__DLINE
		while(iS<nS && lifeSpansExons[iS].k1<=0) /***/
		{
			iS++;
		}
		
		
		if(iS==nS)
		{
		prevB1=1;
		}
		else
		{
		prevB1=lifeSpansExons[iS].k1+1;
		iS++;
		}
		blockSpanEndInc=iS-1;
		__DLINE
		while(iS<nS || iE<nE)
		{
			if(iE==nE || (iS<nS && lifeSpansExons[iS].k1<lifeSpansExons[iE].k2) ) //START COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_START;
			}
			else if(iS==nS || (iE<nE && lifeSpansExons[iS].k1>lifeSpansExons[iE].k2)) //END COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iE].k2);
				flag=RJT_END;
			}
			else //start end equal => COME together
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_BOTHSTARTEND;
			}
			
		//	cerr<<flag<<" event:"<<prevB1<<"-"<<coord1b4<<" for blocks "<<blockSpanStart<<"-"<<blockSpanEndInc<<endl;
			
			//create a new JnxTagBlock!!
			if(blockSpanStart==blockSpanEndInc)
			{
				cerr<<"WRONG: blockSpanStart == blockSpanEndInc!"<<blockSpanStart<<","<<blockSpanEndInc<<endl;
			}
			else
			{
				cerr<<"blockSpanStart,blockSpanEndInc="<<blockSpanStart<<","<<blockSpanEndInc<<endl;
			GffEntry::JnxTagBlock* newJBlock=new GffEntry::JnxTagBlock(&jnxtag,prevB1,coord1b4,blockSpanStart,blockSpanEndInc);
			jnxtag.jblocks.push_back(newJBlock);
			}
			
			//update the blockSpanStart,end now!
			
			switch(flag)
			{
			case RJT_START:
				blockSpanEndInc=iS; //now iS from the right starts.
				iS++;
				break;
			case RJT_END:
				blockSpanStart=iE+1; //the iE is Ended
				iE++;
				break;
			case RJT_BOTHSTARTEND:
				blockSpanStart=iE+1;
				blockSpanEndInc=iS;
				iE++;
				iS++;
				break;
			}
			
			
			
			prevB1=coord1b4+1;
			
			if(coord1b4==maxPosition1GivenTagLength)
				break;
		}
	}
	
	__DLINE
	
	for(int i=0;i<jnxtag.exonPathSize;i++)
	{
		for(int j=0;j<jnxtag.boundSize;j++)
		{
			int tint;
			is>>tint;
		//	cerr<<"read "<<tint;
			jnxtag.exonPaths[i][j]=GffEntry::getExon(tint);
		//	cerr<<"; exon"<<jnxtag.exonPaths[i][j]<<endl;
			
			if(j>0)
			{//	cerr<<" exon j-1:"<<jnxtag.exonPaths[i][j-1]<<endl;
				__DLINE
				jnxtag.exonPaths[i][j-1]->addNewJnxBlock(jnxtag.exonPaths[i][j],NULL,jnxtag.strand);
				__DLINE
			}
		}
	}
	
	__DLINE
	is>>tstr; //read in the ";"
	if(tstr!=";"){
		//die("line not terminated with ;");
		//cerr<<"line "
	}
	
	__DLINE
	//now goto every of the JnxTagBlock and add to the Exons;
	
	
	for(vector<GffEntry::JnxTagBlock*>::iterator jtbi=jnxtag.jblocks.begin();jtbi!=jnxtag.jblocks.end();jtbi++)
	{
		
		GffEntry::JnxTagBlock* jblock=(*jtbi);
	//	cerr<<"jblock:"<<jblock->blockSpanStart<<"-"<<jblock->blockSpanEndInc<<"\ttotal"<<jnxtag.boundSize<<endl;
		
		for(int i=0;i<jnxtag.exonPathSize;i++)
		{
			for(int j=jblock->blockSpanStart;j<jblock->blockSpanEndInc;j++)
			{
				jnxtag.exonPaths[i][j]->addNewJnxBlock(jnxtag.exonPaths[i][j+1],jblock,jnxtag.strand);
			}
		}
	}
	__DLINE
	if(jnxtag.jblocks.size()>=2)
	{
		//GffEntry::JnxTagBlock * block;
		//block->startOnTag1==1?
		if(jnxtag.jblocks[0]->startOnTag1!=1)
		{
			cerr<<"-------------------ERROR: jblocks[0]->startOnTag1!=1"<<endl;
		}
		if(jnxtag.jblocks[jnxtag.jblocks.size()-1]->endOnTag1!=maxPosition1GivenTagLength)
		{
			cerr<<"-------------------ERROR: jblocks[N-1]->endOnTag1!=max"<<endl;
		}
	}
		
	__DLINE
	delete[] lifeSpansExons;
	//}
	
	
	//return is;
}

inline void _readJnxTagok(istream&is, GffEntry::JnxTag& jnxtag, int minSpan,int readLength)
{
	string tstr;

	cerr<<endl;

	__DLINE
	is>>tstr;
	jnxtag.id.set(tstr);
	__DLINE
	is>>tstr; //chr, ignored: already gotten from id;
	
	//is>>boundSize;
	is>>jnxtag.strand;
	is>>jnxtag.lengthTag;
	is>>jnxtag.boundSize; 
	is>>jnxtag.exonPathSize;
	jnxtag.allocBoundsAndExonPaths();
	__DLINE
	
	vector<int> boundsCoord;
	
	//now make JnxTagBlocks;
	//it's in 11 land!!
	//{
	
	KeyPair<int,int> *lifeSpansExons=new KeyPair<int,int>[jnxtag.boundSize];
	
	int lenInc=0;
	int lenExc=0;
	
	
	KeyPair<int,int> prevBo(INT_MIN,INT_MIN);
	
	__DLINE
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		is>>tstr;
		StringUtil::splitInt2(tstr,"-",boundsCoord,true);
		if(boundsCoord.size()!=2)
			continue;
		jnxtag.bounds[i].k1=boundsCoord[0]; //11
		jnxtag.bounds[i].k2=boundsCoord[1]; //11
		if(boundsCoord[0]<prevBo.k1 || boundsCoord[0]<prevBo.k2)
			cerr<<"WRONGB:"<<boundsCoord[0]<<"-"<<boundsCoord[1]<<endl;
		prevBo=jnxtag.bounds[i];
		
		
	}
	
	cerr<<"******* " <<jnxtag.id.chr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<jnxtag.bounds[i].k1<<"\t";
	}

	cerr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<jnxtag.bounds[i].k2<<"\t";
	}
	cerr<<endl;
	cerr<<"*******"<<endl;

	int maxPosition1GivenTagLength=jnxtag.lengthTag-readLength+1;
	__DLINE


	//now calculate the lifespan of each blocks in jnxtag;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		int thisLen=len11(jnxtag.bounds[i]);
		lenInc+=thisLen;
		int calLifeSpansLowBound=lenExc+minSpan-readLength+1;
		lifeSpansExons[i].k1=calLifeSpansLowBound-1; //1 before /***/
		lifeSpansExons[i].k2=lenInc-minSpan+1;				//1 
		lenExc+=thisLen;
	}
	
	cerr<<"*******"<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<lifeSpansExons[i].k1<<"\t";
	}

	cerr<<endl;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		cerr<<lifeSpansExons[i].k2<<"\t";
		if(lifeSpansExons[i].k2<0)
			cerr<<"WRONGAAA"<<endl;
	}
	cerr<<endl;
	cerr<<"*******"<<endl;

	__DLINE
	int iS=1;
	int iE=0;
	int nS=jnxtag.boundSize;
	int nE=jnxtag.boundSize-1;
	
	int prevB1=1;
	int flag;
	int coord1b4;
	
	int blockSpanStart;
	int blockSpanEndInc;
	__DLINE
	while(iE<nE && lifeSpansExons[iE].k2<2)
	{
		iE++;
	}
	__DLINE
	if(iE<nE)
	{
		iS=iE+1;
		blockSpanStart=iE;
		__DLINE
		
		if(lifeSpansExons[iS].k1<=0)
		{
			prevB1=1;
			
		
			while(iS<nS && lifeSpansExons[iS].k1<=0) /***/
			{
				iS++;
			}
		
		}
		else
		{
			prevB1=lifeSpansExons[iS].k1+1;
			iS++;
		}
		
		
		blockSpanEndInc=iS-1;
		__DLINE
		
		while(iS<nS || iE<nE)
		{
			if(iE==nE || (iS<nS && lifeSpansExons[iS].k1<lifeSpansExons[iE].k2) ) //START COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_START;
			}
			else if(iS==nS || (iE<nE && lifeSpansExons[iS].k1>lifeSpansExons[iE].k2)) //END COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iE].k2);
				flag=RJT_END;
			}
			else //start end equal => COME together /*iS<nS && iE<nE && lifeSpansExons[iS].k1==lifeSpansExons[iE].k2*/
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_BOTHSTARTEND;
			}
			
		//	cerr<<flag<<" event:"<<prevB1<<"-"<<coord1b4<<" for blocks "<<blockSpanStart<<"-"<<blockSpanEndInc<<endl;
			
			//create a new JnxTagBlock!!
			if(blockSpanStart==blockSpanEndInc)
			{
				cerr<<"WRONG: blockSpanStart == blockSpanEndInc!"<<blockSpanStart<<","<<blockSpanEndInc<<endl;
			}
			else
			{
				cerr<<"blockSpanStart,blockSpanEndInc="<<blockSpanStart<<","<<blockSpanEndInc<<";pos=["<<prevB1<<","<<coord1b4<<"]\t";
				
				bool containExonSmallerThanSpan=false;
				for(int e=blockSpanStart;e<=blockSpanEndInc;e++)
				{
					if(len11(jnxtag.bounds[e])<minSpan)
					{
						containExonSmallerThanSpan=true;
					}
				}
				
				if(containExonSmallerThanSpan)
				{
					cerr<<":ignored: contain Exon Smaller than Span"<<endl;
				}
				else
				{
					cerr<<":created"<<endl;
					GffEntry::JnxTagBlock* newJBlock=new GffEntry::JnxTagBlock(&jnxtag,prevB1,coord1b4,blockSpanStart,blockSpanEndInc);
					jnxtag.jblocks.push_back(newJBlock);
				}
			}
			
			//update the blockSpanStart,end now!
			
			switch(flag)
			{
			case RJT_START:
				blockSpanEndInc=iS; //now iS from the right starts.
				iS++;
				break;
			case RJT_END:
				blockSpanStart=iE+1; //the iE is Ended
				iE++;
				break;
			case RJT_BOTHSTARTEND:
				blockSpanStart=iE+1;
				blockSpanEndInc=iS;
				iE++;
				iS++;
				break;
			}
			
			
			
			prevB1=coord1b4+1;
			
			if(coord1b4==maxPosition1GivenTagLength)
				break;
		}
	}
	
	__DLINE
	
	for(int i=0;i<jnxtag.exonPathSize;i++)
	{
		for(int j=0;j<jnxtag.boundSize;j++)
		{
			int tint;
			is>>tint;
		//	cerr<<"read "<<tint;
			jnxtag.exonPaths[i][j]=GffEntry::getExon(tint);
		//	cerr<<"; exon"<<jnxtag.exonPaths[i][j]<<endl;
			
			if(j>0)
			{//	cerr<<" exon j-1:"<<jnxtag.exonPaths[i][j-1]<<endl;
				__DLINE
				jnxtag.exonPaths[i][j-1]->addNewJnxBlock(jnxtag.exonPaths[i][j],NULL,jnxtag.strand);
				__DLINE
			}
		}
	}
	
	__DLINE
	is>>tstr; //read in the ";"
	if(tstr!=";"){
		//die("line not terminated with ;");
	}
	
	__DLINE
	//now goto every of the JnxTagBlock and add to the Exons;
	
	
	for(vector<GffEntry::JnxTagBlock*>::iterator jtbi=jnxtag.jblocks.begin();jtbi!=jnxtag.jblocks.end();jtbi++)
	{
		
		GffEntry::JnxTagBlock* jblock=(*jtbi);
	//	cerr<<"jblock:"<<jblock->blockSpanStart<<"-"<<jblock->blockSpanEndInc<<"\ttotal"<<jnxtag.boundSize<<endl;
		
		for(int i=0;i<jnxtag.exonPathSize;i++)
		{
			for(int j=jblock->blockSpanStart;j<jblock->blockSpanEndInc;j++)
			{
				jnxtag.exonPaths[i][j]->addNewJnxBlock(jnxtag.exonPaths[i][j+1],jblock,jnxtag.strand);
			}
		}
	}
	__DLINE
	if(jnxtag.jblocks.size()>=2)
	{
		//GffEntry::JnxTagBlock * block;
		//block->startOnTag1==1?
		if(jnxtag.jblocks[0]->startOnTag1!=1)
		{
			//cerr<<"-------------------ERROR: jblocks[0]->startOnTag1!=1"<<endl;
		}
		if(jnxtag.jblocks[jnxtag.jblocks.size()-1]->endOnTag1!=maxPosition1GivenTagLength)
		{
			//cerr<<"-------------------ERROR: jblocks[N-1]->endOnTag1!=max"<<endl;
		}
	}
		
	__DLINE
	delete[] lifeSpansExons;
	//}
	

	//return is;
}

inline void readJnxTag(istream&is, GffEntry::JnxTag& jnxtag, int minSpan,int readLength)
{
	string tstr;


	is>>tstr;
	jnxtag.id.set(tstr);

	is>>tstr; //chr, ignored: already gotten from id;

	is>>jnxtag.strand;
	is>>jnxtag.lengthTag;
	is>>jnxtag.boundSize;
	is>>jnxtag.exonPathSize;
	jnxtag.allocBoundsAndExonPaths();


	vector<int> boundsCoord;

	//now make JnxTagBlocks;
	//it's in 11 land!!
	//{

	KeyPair<int,int> *lifeSpansExons=new KeyPair<int,int>[jnxtag.boundSize];

	int lenInc=0;
	int lenExc=0;


	KeyPair<int,int> prevBo(INT_MIN,INT_MIN);


	for(int i=0;i<jnxtag.boundSize;i++)
	{
		is>>tstr;
		StringUtil::splitInt2(tstr,"-",boundsCoord,true);
		if(boundsCoord.size()!=2)
			continue;
		jnxtag.bounds[i].k1=boundsCoord[0]; //11
		jnxtag.bounds[i].k2=boundsCoord[1]; //11
		if(boundsCoord[0]<prevBo.k1 || boundsCoord[0]<prevBo.k2)
			cerr<<"WRONGB:"<<boundsCoord[0]<<"-"<<boundsCoord[1]<<endl;
		prevBo=jnxtag.bounds[i];


	}



	int maxPosition1GivenTagLength=jnxtag.lengthTag-readLength+1;



	//now calculate the lifespan of each blocks in jnxtag;
	for(int i=0;i<jnxtag.boundSize;i++)
	{
		int thisLen=len11(jnxtag.bounds[i]);
		lenInc+=thisLen;
		int calLifeSpansLowBound=lenExc+minSpan-readLength+1;
		lifeSpansExons[i].k1=calLifeSpansLowBound-1; //1 before /***/
		lifeSpansExons[i].k2=lenInc-minSpan+1;				//1
		lenExc+=thisLen;
	}


	int iS=1;
	int iE=0;
	int nS=jnxtag.boundSize;
	int nE=jnxtag.boundSize-1;

	int prevB1=1;
	int flag;
	int coord1b4;

	int blockSpanStart;
	int blockSpanEndInc;

	while(iE<nE && lifeSpansExons[iE].k2<2)
	{
		iE++;
	}

	if(iE<nE)
	{
		iS=iE+1;
		blockSpanStart=iE;


		if(lifeSpansExons[iS].k1<=0)
		{
			prevB1=1;


			while(iS<nS && lifeSpansExons[iS].k1<=0) /***/
			{
				iS++;
			}

		}
		else
		{
			prevB1=lifeSpansExons[iS].k1+1;
			iS++;
		}


		blockSpanEndInc=iS-1;


		while(iS<nS || iE<nE)
		{
			if(iE==nE || (iS<nS && lifeSpansExons[iS].k1<lifeSpansExons[iE].k2) ) //START COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_START;
			}
			else if(iS==nS || (iE<nE && lifeSpansExons[iS].k1>lifeSpansExons[iE].k2)) //END COMES first
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iE].k2);
				flag=RJT_END;
			}
			else //start end equal => COME together /*iS<nS && iE<nE && lifeSpansExons[iS].k1==lifeSpansExons[iE].k2*/
			{
				coord1b4=MIN(maxPosition1GivenTagLength,lifeSpansExons[iS].k1);
				flag=RJT_BOTHSTARTEND;
			}


			//create a new JnxTagBlock!!
			if(blockSpanStart==blockSpanEndInc)
			{
				//cerr<<": ingnored blockSpanStart == blockSpanEndInc!"<<blockSpanStart<<","<<blockSpanEndInc<<endl;
			}
			else
			{
				//cerr<<"blockSpanStart,blockSpanEndInc="<<blockSpanStart<<","<<blockSpanEndInc<<";pos=["<<prevB1<<","<<coord1b4<<"]\t";

				bool containExonSmallerThanSpan=false;
				for(int e=blockSpanStart;e<=blockSpanEndInc;e++)
				{
					if(len11(jnxtag.bounds[e])<minSpan)
					{
						containExonSmallerThanSpan=true;
					}
				}

				if(containExonSmallerThanSpan)
				{
					//cerr<<":ignored: contain Exon Smaller than Span"<<endl;
				}
				else
				{
					//cerr<<":created"<<endl;
					GffEntry::JnxTagBlock* newJBlock=new GffEntry::JnxTagBlock(&jnxtag,prevB1,coord1b4,blockSpanStart,blockSpanEndInc);
					jnxtag.jblocks.push_back(newJBlock);
				}
			}

			//update the blockSpanStart,end now!

			switch(flag)
			{
			case RJT_START:
				blockSpanEndInc=iS; //now iS from the right starts.
				iS++;
				break;
			case RJT_END:
				blockSpanStart=iE+1; //the iE is Ended
				iE++;
				break;
			case RJT_BOTHSTARTEND:
				blockSpanStart=iE+1;
				blockSpanEndInc=iS;
				iE++;
				iS++;
				break;
			}



			prevB1=coord1b4+1;

			if(coord1b4==maxPosition1GivenTagLength)
				break;
		}
	}


	for(int i=0;i<jnxtag.exonPathSize;i++)
	{
		for(int j=0;j<jnxtag.boundSize;j++)
		{
			int tint;
			is>>tint;
			jnxtag.exonPaths[i][j]=GffEntry::getExon(tint);

			if(j>0)
			{
				jnxtag.exonPaths[i][j-1]->addNewJnxBlock(jnxtag.exonPaths[i][j],NULL,jnxtag.strand);

			}
		}
	}


	is>>tstr; //read in the ";"
	if(tstr!=";"){
		die("line not terminated with ;");
	}


	//now goto every of the JnxTagBlock and add to the Exons;


	for(vector<GffEntry::JnxTagBlock*>::iterator jtbi=jnxtag.jblocks.begin();jtbi!=jnxtag.jblocks.end();jtbi++)
	{

		GffEntry::JnxTagBlock* jblock=(*jtbi);
	//	cerr<<"jblock:"<<jblock->blockSpanStart<<"-"<<jblock->blockSpanEndInc<<"\ttotal"<<jnxtag.boundSize<<endl;

		for(int i=0;i<jnxtag.exonPathSize;i++)
		{
			for(int j=jblock->blockSpanStart;j<jblock->blockSpanEndInc;j++)
			{
				jnxtag.exonPaths[i][j]->addNewJnxBlock(jnxtag.exonPaths[i][j+1],jblock,jnxtag.strand);
			}
		}
	}




	delete[] lifeSpansExons;
	//}


}





#endif /*GFFENTRY_H_*/
