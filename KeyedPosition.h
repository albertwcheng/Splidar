#ifndef KEYEDPOSITION_H_
#define KEYEDPOSITION_H_

#include "snipIncludes.h"
#include "Nucleic.h"
#include "NucleicBit.h"


#define uint unsigned int

class KeyedPosition;

class Position
{
public:
	uint chrID;
	int pos;
	Position(int _chrID=-1,int _pos=INT_MIN) :chrID(_chrID), pos(_pos)
	{
		
	}
	inline int getPos() const
	{
		if(isForward())
			return pos;
		else
			return -1*pos;
	}
	inline bool isValid() const
	{
		return pos!=INT_MIN;
	}
	inline bool isForward() const
	{
		return pos>0;
	}
	inline bool operator > (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() > right.getPos();
		}
		else
			return chrID > right.chrID;
	}
	inline bool operator < (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() < right.getPos();
		}
		else
			return chrID < right.chrID;
	}
	
	inline void readFromKeyedPosition(istream& is);
	inline void printAsText(ostream& os) const
	{
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
	}
	
	
};

inline ostream& operator << (ostream& os, const Position& pos)
{
	os.write((char*)&pos.chrID,sizeof(uint));
	os.write((char*)&pos.pos,sizeof(int));
	return os;
}

inline istream& operator >> (istream& is, Position& pos)
{
	is.read((char*)&pos.chrID,sizeof(uint));
	
	if(is.eof())
	{
		pos.pos=INT_MIN;
	}

	is.read((char*)&pos.pos,sizeof(int));
	
	return is;
}

class CompactPosition: public Position
{
public:
	
	int end0;
	inline CompactPosition(uint _chrID=-1,int _start0=INT_MIN,int _end0=INT_MIN): Position(_chrID,_start0),end0(_end0)
	{
		
	}
	inline bool operator > (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() > right.getPos();
		}
		else
			return chrID > right.chrID;
	}
	inline bool operator < (const Position& right) const
	{
		if(chrID==right.chrID)
		{
			return getPos() < right.getPos();
		}
		else
			return chrID < right.chrID;
	}
	
	inline int getStart0() const
	{
		return this->getPos();
	}
	
	inline int getEnd0() const
	{
		if(isForward())
			return end0;
		else
			return -1*end0;
	}
	
	inline int getStart1() const
	{
		return getStart0()+1;
	}
	
	inline int getEnd1() const
	{
		return getEnd0()+1;
	}
	
	inline KeyPair<int,int> getBound11() const
	{
		return KeyPair<int,int>(getStart1(),getEnd1());
	}
	
	inline KeyPair<int,int> getBound01() const
	{
		return KeyPair<int,int>(getStart0(),getEnd1());
	}

	inline int getLength() const
	{
		return getEnd0()-getStart0()+1;
	}
	
	
	inline void printAsText(ostream& os) const
	{
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
		os<<"\t";
		if(end0>0)
		{
		os<<"F:";
		os<<end0;
		}
		else
		{
		os<<"R:";
		os<<(-1*end0);
		}		
		
	}
	
};

inline ostream& operator << (ostream& os, const CompactPosition& pos)
{
	os.write((char*)&pos.chrID,sizeof(uint));
	os.write((char*)&pos.pos,sizeof(int));
	os.write((char*)&pos.end0,sizeof(int));	
	return os;
}

inline istream& operator >> (istream& is, CompactPosition& pos)
{
	is.read((char*)&pos.chrID,sizeof(uint));
	
	if(is.eof())
	{
		pos.pos=INT_MIN;
	}

	is.read((char*)&pos.pos,sizeof(int));
	is.read((char*)&pos.end0,sizeof(int));	
	
	return is;
}


class KeyedPosition: public Position
{
public:
	
	NucKey3b a;
	NucKey3b b;
	
	KeyedPosition()
	{
		a=NULL_KEY;
		b=NULL_KEY;
	}


	inline void writeAsPosition(ostream& os) const
	{
		os<<(const Position&)(*this);
	}
	
	inline bool operator > (const KeyedPosition& right) const
	{
		if(b==right.b)
		{
			return a<right.a;
		}
		else
			return b<right.b;
	}


	inline bool absEquals(const KeyedPosition& right) const
	{
		return a==right.a && b==right.b && chrID==right.chrID && pos==right.pos;
	}
	inline bool operator < (const KeyedPosition& right) const
	{
		if(b==right.b)
		{
			return a>right.a;
		}
		else
			return b>right.b;
	}
	
	inline bool operator == (const KeyedPosition& right) const
	{
		return a==right.a && b==right.b;
	}
	
	inline bool operator !=(const KeyedPosition& right) const
	{
		return !(*this==right);
	}
	
	inline bool operator >=(const KeyedPosition& right) const
	{
		return !(*this<right);
		
	}
	inline bool operator <=(const KeyedPosition& right) const
	{
		return !(*this>right);
	}

	inline KeyedPosition(const char* seq,int length,uint _chrID,int _pos,bool forward):Position(_chrID,_pos*(forward?1:-1))
	{
		a=NULL_KEY;
		b=NULL_KEY;
		
		if(length>KEY_CAPACITY*2)
		{
			cerr<<"error try to encode nuc > "<<(KEY_CAPACITY*2)<<endl;
			return;
		}
		
		if(length>KEY_CAPACITY)
		{
			b=Nuc2Key3b(seq,KEY_CAPACITY);
			a=Nuc2Key3b(seq+KEY_CAPACITY,length-KEY_CAPACITY);
		}
		else
		{
			b=Nuc2Key3b(seq,length);
		}
	}
	inline void printAsText(ostream& os) const
	{
		bin(a,os);
		bin(b,os);
		os<<"\t";
		os<<chrID;
		os<<"\t";
		
		if(pos>0)
		{
		os<<"F:";
		os<<pos;
		}
		else
		{
		os<<"R:";
		os<<(-1*pos);
		}
	}
};

inline ostream& operator << (ostream& os, const KeyedPosition& kpos)
{
	os.write((char*)&kpos.a,sizeof(NucKey3b));
	os.write((char*)&kpos.b,sizeof(NucKey3b));
	os.write((char*)&kpos.chrID,sizeof(uint));
	os.write((char*)&kpos.pos,sizeof(int));
	return os;
}

inline istream& operator >> (istream& is,  KeyedPosition& kpos)

{
	is.read((char*)&kpos.a,sizeof(NucKey3b));
	if(is.eof())
	{
		kpos.pos=INT_MIN;
	}
	is.read((char*)&kpos.b,sizeof(NucKey3b));
	is.read((char*)&kpos.chrID,sizeof(uint));
	is.read((char*)&kpos.pos,sizeof(int));
	return is;
}


#define UNREAD 0
#define KEYEDPOSITION 1
#define POSITION 2
#define COMPACTPOSITION 3



class GenomeNmersDecoder
{
public:
	ifstream fin;
	
	int typeRead;
	
	KeyedPosition kpos;
	Position pos;
	CompactPosition cpos;
	
	inline GenomeNmersDecoder(string filename):fin(filename.c_str(),ios::in|ios::binary),typeRead(UNREAD)
	{
		if(!fin.good())
		{
			cerr<<"file "<<filename<<" cannot be opened"<<endl;
			return;
		}
	}
	
	inline int numEntriesPending(int typeToRead=KEYEDPOSITION)
	{
		std::streampos origpos=fin.tellg();
		fin.seekg(0,ios::end);
		std::streampos endpos=fin.tellg();
		std::streampos bytesPending=endpos-origpos;
		

		int sizeItem;
		switch(typeToRead)
		{
	
		case KEYEDPOSITION:
			sizeItem=sizeof(KeyedPosition);
			break;
		case COMPACTPOSITION:
			sizeItem=sizeof(CompactPosition);
			break;
		default:
			sizeItem=sizeof(Position);
			
		}

		
		fin.seekg(origpos,ios::beg);
		
		if(bytesPending%sizeItem!=0)
		{
			cerr<<"Error: the size doesn't match Format"<<endl;
		}
		
		return bytesPending/sizeItem;
	}
	
	
	inline bool readEntry(int typeToRead=KEYEDPOSITION)
	{
		if(fin.eof())
			return false;
		
		//cerr<<"format="<<typeToRead<<endl;
		
		switch(typeToRead)
		{
		case KEYEDPOSITION:
			fin>>kpos;
			typeRead=KEYEDPOSITION;
			return kpos.isValid();
		case COMPACTPOSITION:
			fin>>cpos;
			typeRead=COMPACTPOSITION;
			return cpos.isValid();
		default:
			//cerr<<"readPosition"<<endl;
			fin>>pos;	
			//pos.printAsText(cerr);
			typeRead=POSITION;
			return pos.isValid();
		}
	}
	
	inline const Position& getPosition()
	{
		switch(typeRead)
		{	
		case KEYEDPOSITION:
			return (const Position&)kpos;
		case COMPACTPOSITION:
			return (const Position&)cpos;
		default:
			return pos;
		}
	}
	
	inline const CompactPosition& getCompactPosition()
	{
		return cpos;
	}
	
	inline const KeyedPosition& getKeyedPosition()
	{
		return kpos;
	}
	
	inline ~GenomeNmersDecoder()
	{
		fin.close();
	}
};

inline void Position::readFromKeyedPosition(istream& is)
{
	KeyedPosition kpos;
	is>>kpos;
	(*this)=(Position)kpos;
}

class PositionSorter
{
public:
	inline static void sort(string filebinPos,string sortedOut)
	{
		GenomeNmersDecoder dec(filebinPos);
		
		int numEntries=dec.numEntriesPending(POSITION);
		
		Position *poss=new Position[numEntries];
		
		int ind=0;
		while(dec.readEntry(POSITION))
		{
			poss[ind++]=dec.getPosition();
		}
		
		std::sort(poss,poss+numEntries);
		
		ofstream fout(sortedOut.c_str(),ios::binary|ios::out);
		
		for(int i=0;i<numEntries;i++)
		{
			fout<<poss[i];
		}
		
		cout<<numEntries<<" entries sorted"<<endl;
		
		delete[] poss;
		fout.close();
		
	}

	inline static void sortCompact(string filebinPos,string sortedOut)
	{
		GenomeNmersDecoder dec(filebinPos);
		
		int numEntries=dec.numEntriesPending(POSITION);
		
		Position *poss=new Position[numEntries];
		
		int ind=0;
		while(dec.readEntry(POSITION))
		{
			poss[ind++]=dec.getPosition();
		}
		
		std::sort(poss,poss+numEntries);
		
		ofstream fout(sortedOut.c_str(),ios::binary|ios::out);
		
		CompactPosition cpp;
		
		int nCompact=0;
		
		for(int i=0;i<numEntries;i++)
		{
			Position& curPos=poss[i];
			
			//cerr<<"compare"<<endl;
			//cerr<<"\t";
			//cpp.printAsText(cerr);
			//cerr<<endl;
			//cerr<<"\t";
			//curPos.printAsText(cerr);
			//cerr<<endl;			
			
			if(cpp.chrID!=curPos.chrID || curPos.getPos()>cpp.getEnd0()+1)
			{
				if(cpp.isValid())
				{	
					fout<<cpp;
					nCompact++;
				}
				
				cpp.chrID=curPos.chrID;
				cpp.pos=curPos.pos;
				cpp.end0=curPos.pos;
			}else
				cpp.end0=curPos.pos;
			
		}
		
		if(cpp.isValid())
		{
			
			fout<<cpp;
			nCompact++;
		}
		
		cout<<numEntries<<" entries sorted to "<<nCompact<<" contiguous positions"<<endl;
		delete[] poss;
		fout.close();
		
	}
	
};


class SmartChrMap : public ChrMap
{
public:
	inline SmartChrMap(string filename): ChrMap(filename)
	{
		
	}
	
	typedef KeyPair<string,int> ChrMapInfo;
	ChrMapInfo getChrMapInfoFromPChrID(int id)
	{
		string chrName=this->_map[id];
		vector<string> chrPart;
		StringUtil::split(chrName,"|",chrPart);
		if(chrPart.size()>1)
		{
			vector<string> chrPart1;
			StringUtil::split(chrPart[1],":",chrPart1);
			return KeyPair<string,int>(chrPart[0],StringUtil::atoi(chrPart1[1]));
		}else
			return KeyPair<string,int>(chrPart[0],-1);
		
		
		//return KeyPair<string,int>(chrPart[0],(chrPart.size()>1?StringUtil::atoi(chrPart[1]):-1));
	}
};


class PositionChrPartitioner
{
public:
	
	ChrMap indexedChrMap;
	
	string outputPrefix;
	string outputSuffix;
	bool ignoreReverse;
	typedef buffered_ofstream<Position,vector<Position>,vector<Position>::iterator> BufOut;
	map<string,BufOut*> bfOuts; //binned by chr
	typedef map<string,BufOut*>::iterator I;
	typedef map<string,BufOut*>::value_type V;
	
	inline PositionChrPartitioner(string _indexedChrFile,string _outputPrefix,string _outputSuffix,bool _ignoreReverse=true): indexedChrMap(_indexedChrFile),outputPrefix(_outputPrefix),outputSuffix(_outputSuffix),ignoreReverse(_ignoreReverse)
	{
		/*for(map<int,string>::iterator i=indexedChrMap._map.begin();i!=indexedChrMap._map.end();i++)
		{
			cerr<<i->first<<"=>"<<i->second<<endl;
		}*/
	}
	
	inline void partition(string filebin,int formatBin=KEYEDPOSITION)
	{
		GenomeNmersDecoder dec(filebin);
		
		//cerr<<"filebin="<<filebin<<endl;
		//cerr<<"format="<<formatBin<<endl;
		
		int nPart=0;
		int nPass=0;
		
		while(dec.readEntry(formatBin))
		{
			nPass++;
			//cerr<<"a"<<nPass<<endl;
			const Position& pos=dec.getPosition();
			
			//if(nPass%100)
			//cerr<<"folding "<<nPass<<endl;
			
			if(ignoreReverse && !pos.isForward())
				continue;
			//cerr<<"b"<<nPass<<endl;
			nPart++;
			uint chr=pos.chrID;
			vector<string> chrPart;
			//cerr<<"c"<<nPass<<"chr="<<chr<<endl;
			string chrName=indexedChrMap._map[chr];
			//cerr<<"d"<<nPass<<"got "<<chrName<<endl;
			StringUtil::split(chrName,"|",chrPart);
			//cerr<<"e"<<nPass<<endl;
			string realName=chrPart[0];
			//cerr<<"f"<<nPass<<endl;
			BufOut* bfout=NULL;
			I i=bfOuts.find(realName);
			//cerr<<"g"<<nPass<<endl;
			if(i==bfOuts.end())
			{
				bfout=new BufOut(outputPrefix+realName+outputSuffix,false);
				bfOuts.insert(V(realName,bfout));
			}else
				bfout=i->second;
			//cerr<<"h"<<nPass<<endl;
			bfout->push(pos,ios::binary|ios::out|ios::app,false);
			//cerr<<"i"<<nPass<<endl;
			
			
		}
		
		cerr<<"finish loop"<<endl;
		
		for(I i=bfOuts.begin();i!=bfOuts.end();i++)
		{
			BufOut* bfout=i->second;
			bfout->flush(ios::binary|ios::out|ios::app,false);
			delete bfout;
		}
		
		cout<<nPass<<" entries passed among which "<<nPart<<" were partitioned"<<endl;
		
	}
	
	inline ~PositionChrPartitioner()
	{
		
	}
	
	
};

class GenomeNmersEncoder
{
public:

	ofstream fmap;
	int nmersize;
	int prefixLength;
	string outputPrefix;
	string outputSuffix;
	
	
	typedef buffered_ofstream<KeyedPosition,vector<KeyedPosition>,vector<KeyedPosition>::iterator> BufOut;
	map<string,BufOut*> bfouts;
	
	inline ~GenomeNmersEncoder()
	{	
		for(map<string,BufOut*>::iterator i=bfouts.begin();i!=bfouts.end();i++)
		{
			BufOut* fout=(*i).second;
			fout->flush(ios::out|ios::app|ios::binary,false);
			delete fout;
		}
		fmap.close();
	}
	inline GenomeNmersEncoder(string _outputPrefix,string _outputSuffix,int _prefixLength,string filenameMap,int _nmersize):outputPrefix(_outputPrefix),outputSuffix(_outputSuffix),prefixLength(_prefixLength),fmap(filenameMap.c_str()),nmersize(_nmersize)
	{
		
	}
	inline const char* getStrPtr(int i,const string& str)
	{
		return str.c_str()+i;
	}
	inline string getPrefix(int i,const string& str)
	{
		return str.substr(i,prefixLength);
	}
	inline void transferFromFastaFile(string fastaFileName,string prefixConstraint="")
	{
		int curChrID=0;
		FastaFile ffile(fastaFileName);
		int nReads=0;
		
		int constraintLength=prefixConstraint.length();
		
		while(ffile.readEntry())
		{
			int nCurReads=0;
			curChrID++;
			
			fmap<<curChrID<<"\t"<<ffile.seqName<<endl;
			cout<<"Encoding "<<curChrID<<":"<<ffile.seqName<<" of length "<<ffile.seq.length();
			
			int seqLength=ffile.seq.length();
			if(seqLength<nmersize)
			{
				cerr<<"Ignored: sequence of "<<ffile.seqName<<" has length smaller than nmersize"<<endl;
				cout<<": ignored: sequence has length smaller then nmersize "<<endl;
				continue;
			}
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{
				BufOut* bout=NULL;
				string prefix=getPrefix(i,ffile.seq);
				
				if(prefixConstraint!="" && prefix.substr(0,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				map<string,BufOut*>::iterator bi=bfouts.find(prefix);
				if(bi==bfouts.end())
				{
					bout=new BufOut(outputPrefix+prefix+outputSuffix);
					bfouts.insert(map<string,BufOut*>::value_type(prefix,bout));
				}else
					bout=(*bi).second;
				
				KeyedPosition kpos(getStrPtr(i,ffile.seq),nmersize,curChrID,i+1,true);
				//kpos.printAsText(cerr);
				//cerr<<endl;
				bout->push(kpos,ios::out|ios::binary|ios::app,false);
				nCurReads++;
				
				
				
				
				
			}
			
			
			
			string rseq=reverse_complement(ffile.seq);
			
			
			for(int i=0;i<=seqLength-nmersize;i++)
			{
				BufOut* bout=NULL;
				string prefix=getPrefix(i,rseq);
				
				if(prefixConstraint!="" && prefix.substr(0,constraintLength)!=prefixConstraint)
				{
					continue;
				}
				
				map<string,BufOut*>::iterator bi=bfouts.find(prefix);
				if(bi==bfouts.end())
				{
					bout=new BufOut(outputPrefix+prefix+outputSuffix);
					bfouts.insert(map<string,BufOut*>::value_type(prefix,bout));
				}else
					bout=(*bi).second;
				KeyedPosition kpos(getStrPtr(i,rseq),nmersize,curChrID,seqLength-nmersize-i+1,false);
				//kpos.printAsText(cerr);
				//cerr<<endl;
				bout->push(kpos,ios::out|ios::binary|ios::app,false);
				nCurReads++;
				
				
				
				
			}			
			
			cout<<" outputing "<<nCurReads<<" reads with prefixConstraint "<<prefixConstraint<<endl;
			nReads+=nCurReads;
			
		}
		
		cerr<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl;
		cout<<curChrID<<" sequence(s) processed."<<nReads<<" simulated read(s) outputed "<<endl; 
		
	}
	
	
};

#endif /*KEYEDPOSITION_H_*/
