//#include "GffEntry.h"
#include "snipIncludes.h"

using namespace std;



GffEntry::GffEntry(vector<string>& splits,string _annoSource) :
	exons(NULL) {

	visited=false;
//	cerr<<"num splits="<<splits.size()<<endl;
	bin=StringUtil::atoi(splits[0]);
	name=splits[1];
	chrom=splits[2];
	strand=splits[3][0];
	txStart=StringUtil::atoi(splits[4]);
	txEnd=StringUtil::atoi(splits[5]);
	cdsStart=StringUtil::atoi(splits[6]);
	cdsEnd=StringUtil::atoi(splits[7]);
	exonCount=StringUtil::atoi(splits[8]);
	vector<int> starts, ends, frames;
	starts=StringUtil::splitInt(splits[9], ",");
	ends=StringUtil::splitInt(splits[10], ",");
	id=splits[11];
	name2=splits[12];
	if(name2=="")
	{
		cerr<<_annoSource<<":"<<name<<":"<<id<<":"<<name2<<" does not have a name2"<<endl;
	}
	cdsStartStat=splits[13];
	cdsEndStat=splits[14];
	frames=StringUtil::splitInt(splits[15], ",");
	annoSource=_annoSource;

	exons=new ExonPtr [exonCount];

	GffEntry::GffEntryCount++;

	for (int i=0; i<exonCount; i++) {
		exons[i]=GffEntry::getNewExon(this, chrom, (strand==FORWARD)?i:(exonCount-i-1), starts[i], ends[i],
				frames[i],strand);
		exons[i]->maxCodingStart=MAX(exons[i]->maxCodingStart,cdsStart);
		exons[i]->minCodingEnd=MIN(exons[i]->minCodingEnd,cdsEnd);
	}

}

ostream& operator<<(ostream& os, GffEntry::Locus& loc)
{
	GffEntry::ExonPtr head=loc.getHeadExon();
	string chr=head->chr;
	KeyPair<int,int> bound=loc.getBound();
	os<<chr<<"\t"<<bound.k1<<"\t"<<bound.k2<<"\t";
	os<<loc.transcripts[0]->strand<<"\t";
	os<<chr<<":"<<(bound.k1+1)<<"-"<<bound.k2<<"\t";
	os<<loc.names.size()<<"\t";
	for(set<string>::iterator i=loc.names.begin();i!=loc.names.end();i++)
	{
		os<<(*i)<<"\t";
	}
	os<<loc.transcripts.size()<<"\t";

	for(vector<GffEntry*>::iterator i=loc.transcripts.begin();i!=loc.transcripts.end();i++)
	{
		os<<(*i)->localID<<"\t";
		os<<(*i)->annoSource<<":"<<(*i)->transcriptID<<"\t";
	}
	os<<loc.exonSet.size()<<"\t";
	for(set<GffEntry::ExonPtr >::iterator i=loc.exonSet.begin();i!=loc.exonSet.end();i++)
	{
		os<<(*i)->exonID<<"\t";
	}

	os<<endl;
	return os;
}


istream& operator >>(istream &is, GffEntry::Locus& loc)
{
	string istr;
	int iint,nint;
	string chr;
	int start;
	int end;

	chr="";

	is>>chr;
	loc.chr=chr;

	if(chr=="")
		return is;

	is>>start;

	is>>end;

	is>>loc.strand;
	is>>istr;

	is>>nint;
	for(int i=0;i<nint;i++)
	{
	//	cerr<<"getting name "<<(i+1)<<" of "<<nint<<endl;
		is>>istr;
		loc.names.insert(set<string>::value_type(istr));
	}

	is>>nint;
	for(int i=0;i<nint;i++)
	{
	//	cerr<<"getting transcript "<<(i+1)<<" of "<<nint<<endl;
		is>>iint; //localID
		is>>istr; //annoSource:transcriptID discarded
		//cerr<<"\tget gffentry "<<iint<<endl;
		loc.transcripts.push_back(GffEntry::getGffEntry(iint));
	}

	is>>nint;
	for(int i=0;i<nint;i++)
	{
	//	cerr<<"getting exon "<<(i+1)<<" of "<<nint<<endl;
		is>>iint;
		//cerr<<"\tget exon "<<iint<<endl;
		GffEntry::ExonPtr pE=GffEntry::getExon(iint);
		//cerr<<"\tgot exon "<<pE<<endl;
		loc.exonSet.insert(set<GffEntry::ExonPtr >::value_type(pE));
		//cerr<<"\tinsert exon done"<<endl;
	}

	//cerr<<"finish reading"<<endl;

	GffEntry::ExonPtr pHead=loc.getHeadExon();

	KeyPair<int,int> bound=loc.getBound();

	KeyPair<int,int> readBound(start,end);

	if(pHead->chr!=chr || bound!=readBound)
	{
		cerr<<"unmatched coord!"<<endl;
	}
	//cerr<<"finish compare"<<endl;

	return is;
}

istream& operator >>(istream &is, GffEntry::GBlock& block)
{
	string istr;
	int iint,nint;
	string chr;
//	int start;
//	int end;

	chr="";

	is>>istr;//the id is ignored
	
	vector<string> vstr;
	StringUtil::split(istr,":",vstr,true);
	
	if(vstr.size()==2)
	{
		//right size now get the ID
		block.name=vstr[0];
		block.ID=StringUtil::atoi(vstr[1]);
	}
	
	
	is>>chr;
	block.chr=chr;

	if(chr=="")
		return is;
     
	is>>block.start0;
	is>>block.end1;
	is>>block.strand;


	is>>nint;
	for(int i=0;i<nint;i++)
	{
		is>>iint;
		GffEntry::ExonPtr pE=GffEntry::getExon(iint);
		block.exons.insert(pE);
	}
	
	is>>chr; //to absorb the stupid ;
	
	return is;
}


int GffEntry::getExonRankFromIndex(int i) const
{
	return (strand==FORWARD)?i+1:exonCount-i;
}

 int GffEntry::getExonArrayIndex(int i) const
{
	return (strand==FORWARD)?i:exonCount-i-1;
}

 GffEntry::ExonPtr GffEntry::getExonOfRank(int i) const
{
	return exons[getExonArrayIndex(i)];
}

GffEntry::ExonPtr GffEntry::getNewExon(GffEntry *entry, string chr, int rank,
		int start, int end, int frame,char strand) {

	int exonID=indexedExonome.size();

	Exon * pe=new Exon(chr, rank, start, end, frame, exonID);
	ExonPtr e(pe);

//	if(frame<-1 || frame>3)
//		cerr<<"invalid frame "<<entry->name<<" "<<entry->name2<<" "<<chr<<":"<<start<<"-"<<end<<" with frame "<<frame<<endl;
	vector<RankedGffEntryPair> * v;


	map<ExonPtr,vector<RankedGffEntryPair>* >::iterator i=exonome.find(e);
	if(i==exonome.end()){
		v=new vector<RankedGffEntryPair>;
		e->assTranscripts=v;
		registerIndexedExon(e);
		exonome.insert(map<ExonPtr,vector<RankedGffEntryPair>* >::value_type(e, v));
	}else
	{

		delete pe;
		e=(*i).first;
		v=(*i).second;
	}
	
	if(e->strand==GffEntry::UNKNOWN)
	{
		e->strand=strand;
	}else if(e->strand!=strand)
	{
		e->strand=GffEntry::BOTH_STRANDS;
	}



	v->push_back(RankedGffEntryPair(rank,entry));




	return e;

}

bool GffEntry::Exon::operator <(const Exon &right) const
{
	if(this->chr==right.chr)
	{
		if(this->start==right.start)
		{
			return this->end<right.end;
		}
		else
		{
			return this->start<right.start;
		}
	} else
	{
		return this->chr<right.chr;
	}
}
bool GffEntry::Exon::operator>(const Exon &right) const
{
	if(this->chr==right.chr)
	{
		if(this->start==right.start)
		{
			return this->end>right.end;
		}
		else
		{
			return this->start>right.start;
		}
	} else
	{
		return this->chr>right.chr;
	}
}
bool GffEntry::Exon::operator ==(const Exon& right) const
{
	return this->chr==right.chr && this->start == right.start && this->end == right.end;
}

bool GffEntry::Exon::operator !=(const Exon& right) const
{
	return !((*this)==right);
}
bool GffEntry::Exon::operator >=(const Exon& right) const
{
	return (*this)==right || (*this)>right;
}

bool GffEntry::Exon::operator <=(const Exon& right) const
{
	return (*this)==right || (*this)<right;
}

GffEntry::Exon::Exon(): blocks(NULL), outJnxs(NULL),siblings(NULL), visited(false), strand(GffEntry::UNKNOWN), siblingID(0), usageFlag(USAGE_UNKNOWN), exonGroup(NULL), maxCodingStart(INT_MIN),minCodingEnd(INT_MAX) {
	GffEntry::GffExonCount++;
	this->inDegree=0;
	
}
GffEntry::Exon::~Exon(){
	GffEntry::GffExonCount--;
	if(blocks)
		delete blocks;
	
	if(outJnxs)
	{
		for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator i=outJnxs->begin();i!=outJnxs->end();i++)
		{
			GffEntry::Jnx* jnx=(*i).second;
			delete jnx;
		}
		delete outJnxs;
	}
}
GffEntry::Exon::Exon(string _chr, int _rank, int _start, int _end, int _frame,int _exonID) :
	blocks(NULL),chr(_chr), rank(_rank), start(_start), end(_end), frame(_frame), exonID(_exonID),siblings(NULL), siblingID(0), outJnxs(NULL) , strand(GffEntry::UNKNOWN),visited(false), usageFlag(USAGE_UNKNOWN), exonGroup(NULL), maxCodingStart(INT_MIN),minCodingEnd(INT_MAX){
	GffEntry::GffExonCount++;
	this->inDegree=0;

}
void GffEntry::Exon::set(string _chr, int _rank, int _start, int _end,
		int _frame,int _exonID) {
	chr=_chr;
	rank=_rank;
	start=_start;
	end=_end;
	frame=_frame;
	exonID=_exonID;

	//if(frame<-1 || frame>3)
	//	cerr<<"invalid frame "<<exonID<<" "<<chr<<":"<<start<<"-"<<end<<" with frame "<<frame<<endl;
}

GffEntry::Jnx* GffEntry::Exon::addNewJnxBlock(GffEntry::ExonPtr _rightExon, GffEntry::JnxTagBlock* _jblock,char _strand)
{
	//cerr<<"a"<<endl;
	if(!this->outJnxs)
		this->outJnxs=new map<ExonPtr,GffEntry::Jnx*>;
	
	//	cerr<<"b:outjnx *:"<<this->outJnxs<<endl;	
	//	cerr<<"b:outjnx.size : "<<this->outJnxs->size()<<endl;
	//	cerr<<"b:rightExon:"<<_rightExon<<endl;
	//	cerr<<"b:retrieveRightExon:"<<_rightExon->exonID<<endl;
		
	Jnx* jnx;
	map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator i=this->outJnxs->find(_rightExon);
//	cerr<<"c"<<endl;
	if(i==this->outJnxs->end())
	{
		jnx=new Jnx(this,_rightExon,_strand);
		this->outJnxs->insert(map<GffEntry::ExonPtr,GffEntry::Jnx*>::value_type(_rightExon,jnx));
		_rightExon->inDegree++;
	}
	else
	{
		jnx=(*i).second;
		if(jnx->strand!=_strand) //may have problems, but who fxxking care now?
		{
			jnx->strand=GffEntry::BOTH_STRANDS;
		}
	}
//	cerr<<"d"<<endl;
	
	if(_jblock)
		jnx->jblocks.push_back(_jblock);
//	cerr<<"e"<<endl;
	return jnx;
}

GffEntry::Locus::~Locus()
{
	delete[] exons;
	//destry root, the exon trees
	
	/*if(root)
	{
		root->destroyTree();
		delete root;
	}*/

	//now destroy the exon bound groups

	for(map<string,NExonGroup*>::iterator i=this->exongroups.begin();i!=this->exongroups.end();i++)
	{
		delete i->second;
	}

	for(map<string,NESuperGroup*>::iterator i=this->leftSuperGroups.begin();i!=this->leftSuperGroups.end();i++)
	{
		delete i->second;
	}

	for(map<string,NESuperGroup*>::iterator i=this->rightSuperGroups.begin();i!=this->rightSuperGroups.end();i++)
	{
		delete i->second;
	}


	/*for(vector<ExonBoundGroup*>::iterator i=this->leftExonBoundGroups.begin();i!=leftExonBoundGroups.end();i++)
	{
		delete *i;
	}

	for(vector<ExonBoundGroup*>::iterator i=this->rightExonBoundGroups.begin();i!=rightExonBoundGroups.end();i++)
	{
		delete *i;
	}*/
}


const char GffEntry::FORWARD='+';
const char GffEntry::REVERSE='-';
const char GffEntry::UNKNOWN='.';
const char GffEntry::BOTH_STRANDS='=';
map<GffEntry::ExonPtr,vector<RankedGffEntryPair>* > GffEntry::exonome;
map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*>* > GffEntry::loci;
//map<GffEntry::ExonPtr ,GffEntry::Locus* > GffEntry::loci;
vector< vector<GffEntry::ExonPtr> *> GffEntry::siblingVectors;
vector<GffEntry::ExonPtr> GffEntry::indexedExonome;
vector<GffEntry*> GffEntry::transcriptome;
map<string, set<GffEntry::GBlockPtr> * > GffEntry::gblocks;
multiset<GffEntry::SelexaMatchPtr> GffEntry::blindSpotSelexaMatches;
int GffEntry::GffEntryCount=0;
int GffEntry::GffExonCount=0;
bool GffEntry::selexaMatchLoaded=false;
multimap<string, GffEntry::Locus*> GffEntry::nameIndexedLoci;
//map<string, set<GffEntry::JnxPtr>* > GffEntry::jnxs;
//map<string, vector<GffEntry::JnxTagBlock*>* > GffEntry::jnxtagblocks;	
map<string, map<int,GffEntry::JnxTag*>* > GffEntry::jnxtags;
int GffEntry::gTotalExonReadsUsed=0;

bool GffEntry::isUMPLoaded=false;

GffEntry::~GffEntry() {
	GffEntry::GffEntryCount--;

	//the exon arrays containing the exon Pointers
	if (exons) {
		delete[] exons;
	}
}

void GffEntry::resetEverything()
{
	GffEntry::resetTranscriptome();
	GffEntry::resetExonome();
	GffEntry::resetLoci();
	GffEntry::resetSiblingVectors();
	GffEntry::resetGBlocks();
	GffEntry::resetAllJnxStuff();
}
void GffEntry::resetTranscriptome()
{

	//cerr<<"pre Transcript Count="<<GffEntry::GffEntryCount<<endl;
	for(vector<GffEntry*>::iterator i=GffEntry::transcriptome.begin();i!=GffEntry::transcriptome.end();i++)
		delete (*i);
	//cerr<<"Check Transcript Count="<<GffEntry::GffEntryCount<<endl;
	GffEntry::transcriptome.clear();
}

void GffEntry::resetLoci()
{
	//for(map<ExonPtr,Locus* >::iterator i=GffEntry::loci.begin();i!=GffEntry::loci.end();i++)
	//	delete (*i).second;
	
	for(map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator i=GffEntry::loci.begin();i!=GffEntry::loci.end();i++)
	{
		multimap<GffEntry::ExonPtr,GffEntry::Locus*> * localMap=(*i).second;
		for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator j=localMap->begin();j!=localMap->end();j++)
		{
			delete (*j).second;
		}
		delete localMap;
	}	
	
	GffEntry::loci.clear();
}

void GffEntry::resetSiblingVectors()
{
	for(vector< vector<ExonPtr>* >::iterator i=GffEntry::siblingVectors.begin();i!=GffEntry::siblingVectors.end();i++)
		delete (*i);

	GffEntry::siblingVectors.clear();
}


GffEntry* GffEntry::createGffEntry(vector<string>& splits,string annoSource,int *transcriptID)
{
//	cerr<<"pre2"<<endl;
	GffEntry* entry=new GffEntry(splits,annoSource);
//	cerr<<"post2"<<endl;
	entry->localID=GffEntry::transcriptome.size();
	if(!transcriptID)
	{
		cerr<<"using size "<<endl;
		entry->transcriptID=entry->localID;
	}
	else
		entry->transcriptID=((*transcriptID)++); //GffEntry::transcriptome.size();


	GffEntry::transcriptome.push_back(entry);
	return entry;
}

//exons are not needed to be deleted because mantained by the map Exonome;

void GffEntry::resetGBlocks(map<string, set<GBlockPtr> * >::iterator i)
{
	cerr<<"reseting gblocks for chr "<<(*i).first<<endl;
		set<GBlockPtr> * sset=(*i).second;
		for(set<GBlockPtr>::iterator j=sset->begin();j!=sset->end();j++)
		{
			delete *j;
		}
		delete sset;

}


void GffEntry::resetGBlocks(string chr)
{
	if(chr=="")
	{
		for(map<string, set<GBlockPtr> * >::iterator i=GffEntry::gblocks.begin();i!=gblocks.end();i++)
		{
			GffEntry::resetGBlocks(i);
		}

		GffEntry::gblocks.clear();
	}
	else
	{
		map<string, set<GBlockPtr> * >::iterator i=GffEntry::gblocks.find(chr);
		if(i==GffEntry::gblocks.end())
			return;

		GffEntry::resetGBlocks(i);
		GffEntry::gblocks.erase(i);

	}
	
}
void GffEntry::resetExonome()
{
	for(vector<GffEntry::ExonPtr>::iterator i=GffEntry::indexedExonome.begin();
		i!=GffEntry::indexedExonome.end();
		i++)
	{
		delete (*i);
	}

	GffEntry::indexedExonome.clear();

	GffEntry::exonome.clear();
	/*cerr<<"pre Exon Count="<<GffEntry::GffExonCount<<endl;
	for(ExonomeWalker w=exonome.begin();w!=exonome.end();w++)
	{
		 set<GffEntry*>* entries=&gffsOf(w);

		 delete entries;
	}
	GffEntry::exonome.clear();
	cerr<<"Check Exon Count="<<GffEntry::GffExonCount<<endl;*/
}

GffEntry::Loader::Loader(TranscriptLoadingListener * _listener): listener(_listener)
{

}

void GffEntry::Loader::resetEverything()
{
	GffEntry::resetEverything();
}

vector<KeyPair<string,string> > GffEntry::Loader::loadGffFiles(string sourceFile, bool reset, string prefix)
{
	if(reset)
		resetEverything();

	vector<KeyPair<string,string> > annoMapping;

	if(StringUtil::isSuffix(".config",sourceFile.c_str()))
	{
		SuperConfigFile conf(sourceFile);
		vector<string> annotations=conf.readArray("annotations");
		for(vector<string>::iterator i=annotations.begin();i!=annotations.end();i++)
		{
			const string& annoSource=*i;
			string filename=conf.readSuperString(annoSource);

			annoMapping.push_back(KeyPair<string,string>(annoSource,filename));
			cerr<<"call load file "<<filename<<endl;
			loadGffFile(filename,annoSource,false);
		}
		/*set<string> keys=conf.keys();
		for(set<string>::iterator i=keys.begin();i!=keys.end();i++)
		{
			const string& name=*i;
			if(name[0]=='.')
			{
				string annoSource=name.substr(1);
				string filename=conf.readSuperString(name);
				annoMapping.push_back(KeyPair<string,string>(annoSource,filename));
				cerr<<"call load file "<<filename<<endl;
				loadGffFile(filename,annoSource,false);
			}
		}*/

	}
	else
	{
		ifstream sourceF(sourceFile.c_str());
		if(!sourceF.good())
		{
			cerr<<"Cannot read source file"<<endl;

			sourceF.close();
			return vector<KeyPair<string,string> >();
		}

		string annoSource;
		string filename;

		while(!sourceF.eof())
		{
			sourceF>>annoSource;
			sourceF>>filename;
			cerr<<annoSource<<" \t "<<filename<<endl;
			if(annoSource.length()<1 || annoSource.length()<1)
			{
				continue;
			}

			if(annoSource=="#prefix")
			{
				prefix=filename;
				continue;
			}

			annoMapping.push_back(KeyPair<string,string>(annoSource,filename));
			cerr<<"call load file "<<filename<<endl;
			loadGffFile(prefix+filename,annoSource,false);

		}

		sourceF.close();
	}
//	cerr<<0<<":"<<GffEntry::indexedExonome[0]<<endl;
//	cerr<<313416<<":"<<GffEntry::indexedExonome[313416]<<endl;
	return annoMapping;
}

void GffEntry::Loader::loadGffFile(string filename,string annoSource, bool reset)
{
	if(reset)
		resetEverything();

	ifstream gff(filename.c_str());

	if (!gff.good()) {

		cerr<<"Cannot open Gff"<<endl;
		return;

	}
	cerr<<"loading Gff File "<<filename<<" as "<<annoSource<<endl;
	char line[GFF_BUFFER_SIZE]; //the buffer for reading lines from gff file
	int id=0;

	while (!gff.eof()) //read all entries from gff file
	{


		gff.getline(line, GFF_BUFFER_SIZE);

		if(!strcmp(line,""))
			break;

		vector<string> spliton; //vector to carry the split tokens from each line of gff file
		StringUtil::split(string(line), string("\t"), spliton);

		if (spliton.size()<1)
			continue;
		if (spliton[0]=="#bin") //ignore commented lines from gff files
			continue;


		GffEntry &entry = *(GffEntry::createGffEntry(spliton,annoSource,&id)); //construct a GffEntry from the tokens

		if(listener)
			listener->transcriptLoaded(entry,annoSource);
	}
	gff.close();
//	cerr<<0<<":"<<GffEntry::indexedExonome[0]<<endl;
	//cerr<<313416<<":"<<GffEntry::indexedExonome[313416]<<endl;
}

void GffEntry::Loader::loadSiblings(string filename)
{
	ifstream fil(filename.c_str());
	vector<ExonPtr> * vep;
	
	if(!fil.good())
	{
		cerr<<"cannot open "<<filename<<endl;
	}
	
	int sibID=0;
	
	while(!fil.eof())
	{
		string eIDs("");
		fil>>eIDs;
		if(eIDs=="")
			break;
		
		if(eIDs=="[")
		{
			vep=new vector<ExonPtr>;
			continue;
		}
		
		if(eIDs=="]")
		{
			GffEntry::siblingVectors.push_back(vep);
			sibID++;
		}
		
		int exonID=StringUtil::atoi(eIDs);
		ExonPtr pE=GffEntry::getExon(exonID);
		pE->siblings=vep;
		pE->siblingID=sibID;
		vep->push_back(pE);
		
	}
	
	
	fil.close();
}


void GffEntry::Loader::loadLoci(string filename)
{
	ifstream fil(filename.c_str());

/*	cerr<<0<<":"<<GffEntry::indexedExonome[0]<<endl;
	cerr<<313416<<":"<<GffEntry::indexedExonome[313416]<<endl;
	cerr<<"chr="<<GffEntry::indexedExonome[313416]->chr<<endl;
*/
	
	if(!fil.good())
	{
		cerr<<"Cannot open "<<filename<<endl;
		return;
	}


	int c=0;
	
	
	string prevChr="";
	
	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * curMapLocus=NULL;
	
	while(!fil.eof())
	{
		GffEntry::Locus * loc=new GffEntry::Locus;
		c++;
	//	if(c%100==1)
	//		cerr<<"loading Locus #"<<c<<endl;

		fil>>(*loc);

		if(loc->chr=="")
		{
			delete loc;
			break;
		}
		
		if(prevChr!=loc->chr)
		{
			prevChr=loc->chr;
			map<string, multimap<GffEntry::ExonPtr, GffEntry::Locus*>* >::iterator i=GffEntry::loci.find(prevChr);
			if(i==GffEntry::loci.end())
			{
				curMapLocus=new multimap<GffEntry::ExonPtr,GffEntry::Locus*>;
				GffEntry::loci.insert(map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*> *>::value_type(prevChr,curMapLocus));
			}
			else
				curMapLocus=(*i).second;
		}
		
		loc->finalize();

		GffEntry::ExonPtr pHead=loc->getHeadExon();
		//GffEntry::loci.
		curMapLocus->insert(multimap<GffEntry::ExonPtr,GffEntry::Locus* >::value_type(pHead,loc));

	}

	fil.close();

}

void GffEntry::Loader::loadGBlocks(string filename)
{
	ifstream fil(filename.c_str());
	
	if(!fil.good())
	{
		cerr<<"Cannot open "<<filename<<endl;
		return;
	}


	int c=0;

	set<GffEntry::GBlockPtr> *chrBlock=NULL;
	
	string curChr="";
	
	while(!fil.eof())
	{
		
		GffEntry::GBlockPtr block=new GffEntry::GBlock;
		
		c++;
		
		//if(c%100==1)
		//	cerr<<"loading GBlock #"<<c<<endl;

		fil>>(*block);

		if(block->chr=="")
		{
			delete block;
			break;
		}
		if(block->chr!=curChr)
		{
			curChr=block->chr;
			map<string,set<GffEntry::GBlockPtr>* >::iterator i=GffEntry::gblocks.find(curChr);
			if(i==GffEntry::gblocks.end())
			{
				chrBlock=new set<GffEntry::GBlockPtr>;
				GffEntry::gblocks.insert(map<string,set<GffEntry::GBlockPtr> * >::value_type(curChr,chrBlock));
			}
			else
				chrBlock=(*i).second;
		}
		
		chrBlock->insert(block);
		for(set<GffEntry::ExonPtr>::iterator be=block->exons.begin();be!=block->exons.end();be++)
		{
			GffEntry::ExonPtr pE=(*be);
			if(!pE->blocks)
			{
				pE->blocks=new vector<GffEntry::GBlockPtr>;
			}
			pE->blocks->push_back(block);
			//cerr<<"adding gblock "<<block->chr<<":"<<block->start0<<"-"<<block->end1<<" to exon "<<pE->chr<<":"<<pE->start<<"-"<<pE->end<<endl;
			if(block->start0<pE->start || block->end1>pE->end)
			{
				cerr<<"horrible_error:"<<"adding gblock "<<block->chr<<":"<<block->start0<<"-"<<block->end1<<" to exon "<<pE->chr<<":"<<pE->start<<"-"<<pE->end<<endl;
			}
		}

	}

	fil.close();

}

void GffEntry::Loader::loadJnxSelexaMatches(string filename,int format)
{
	ifstream fil(filename.c_str());
	
	if(!fil.good())
	{
		cerr<<"Cannot open "<<filename<<endl;
		return;
	}

	
	///////BAM FORMAT SUPPORT: BLOCK 1
	samfile_t *samfp=NULL;
	bam1_t *bamt=bam_init1();
	
	
	if(format==BAM_FORMAT)
	{
		//close file and open as bam instead
		fil.close();
		if((samfp=samopen(filename.c_str(),"rb",0))==0){
			cerr<<"Cannot open bam file "<<filename<<endl;
			return;
		}
	}
	//////////////////
	
	
	map<int,GffEntry::JnxTag* > * chrJnxTags=NULL;
	
	int c=0;

	string curChr="";
	
	
	int prevPos=-1;
	int prevJID=-1;
	GffEntry::JnxTag* jnxtag=NULL;
	
	
	////BLOCK 2
	bool readRead=true;
	
	
	while(true)  //(!fil.eof())
	{
		
		if(format==BAM_FORMAT && samread(samfp,bamt)<0) //if no more to read.
			break;
		
		if(format!=BAM_FORMAT && fil.eof()) //if not, end of file reached, break
			break;
		////
		
		
		GffEntry::SelexaMatchPtr sm=new GffEntry::SelexaMatch;
		
		c++;
		
		if(c%10000==1)
			cerr<<"loading JnxSelexaMatch line #"<<c<<endl;


		////BAM format support : BLOCK 3
		
		if(format==BAM_FORMAT)
			readRead=GffEntry::SelexaMatch::readMatchFromBam(samfp,bamt,*sm);
		else
			readRead=GffEntry::SelexaMatch::readMatch(fil,*sm,format);
		///////
		

		if(sm->chr=="")
		{
			cerr<<"ending at "<<c<<endl;
			delete sm;
			break;
		}
		
		if(!readRead)
		{
			delete sm;
			continue;
		}


		GffEntry::JnxTagID jnxtagID(sm->chr);
		
		
		if(jnxtagID.chr!=curChr)
		{
			prevPos=-1;
			prevJID=-1;
			jnxtag=NULL;
			curChr=jnxtagID.chr;
			map<string,map<int,GffEntry::JnxTag*>* >::iterator i=GffEntry::jnxtags.find(curChr);
			if(i==GffEntry::jnxtags.end())
			{
				fil.close();
				cerr<<string("jnxtags for chr ")<<curChr<<" not found";
				die("");
			}
			else
			{
				chrJnxTags=(*i).second;
			}	
		
		}
		
		if(prevJID!=jnxtagID.id)
		{
			prevPos=-1;
			map<int,GffEntry::JnxTag*>::iterator i=chrJnxTags->find(jnxtagID.id);
			
			prevJID=jnxtagID.id;
			
			if(jnxtag)
				jnxtag->updateJnxTagBlockWithSelexaData();
			
			jnxtag=NULL;
			
			
			if(i==chrJnxTags->end())
			{
				//fil.close();
				cerr<<"Tags "<<jnxtagID.id<<"not added"<<endl;
				delete sm;
				
				//continue;
				//die("tags not added");
			}else
			{
				jnxtag=(*i).second;
			}
			
			if(jnxtag)
			{
				if(jnxtag->selexaMatches.size()>0)
				{
					//fil.close();
					//delete sm;
					//die("strange: jnxtag has already been added with selexa matches: Selexa File Not sorted by chr/id, require a sorted file");
				}
			
				
			}
			
			
		}
		
		if(sm)
			if(sm->gpos0<prevPos)
			{
				fil.close();
				die("Selexa File Not Sorted, require a sorted file");
			}
		
		if(jnxtag)
			jnxtag->selexaMatches.insert(sm);
		
		
	}
	
	cerr<<c<<" lines of jnx read loaded"<<endl;

	if(jnxtag)
	{
		jnxtag->updateJnxTagBlockWithSelexaData();
	}
	
	

	//GffEntry::selexaMatchLoaded=true;
	
	////BAM Support BLOCK4
	if(format==BAM_FORMAT){
		bam_destroy1(bamt);
		samclose(samfp);
	}
	else{
		fil.close();
	}
	//////
	

}


void GffEntry::Loader::loadJnxInfo(string filename,int minSpan,int readLength) //to JnxTags
{
	
	if(GffEntry::isJnxInfoLoaded())
		return;
	
	ifstream fil(filename.c_str());
	if(!fil.good())
	{
		cerr<<"Cannot open "<<filename<<endl;
		return;
	}
	
	map<int,JnxTag*> * chrTags=NULL;
		
	int c=0;

	string curChr="";
		

		
	while(!fil.eof())
	{
			
			GffEntry::JnxTag* newTag=new GffEntry::JnxTag;
			
			c++;
			
			if(c%1000==1)
				cerr<<"loading JnxTag #"<<c<<"\t";

			readJnxTag(fil,*newTag,minSpan,readLength);
			//cerr<<endl;
			if(newTag->id.chr=="")
			{
				delete newTag;
				break;
			}
			
			if(newTag->id.chr!=curChr)
			{

				curChr=newTag->id.chr;
				map<string,map<int,GffEntry::JnxTag*>* >::iterator i=GffEntry::jnxtags.find(curChr);
				if(i==GffEntry::jnxtags.end())
				{
					chrTags=new map<int,GffEntry::JnxTag*>;
					GffEntry::jnxtags.insert(map<string,map<int,GffEntry::JnxTag*>* >::value_type(curChr,chrTags));
				}
				else
				{
					chrTags=(*i).second;
				}	
				
			}
			
			chrTags->insert(map<int,GffEntry::JnxTag*>::value_type(newTag->id.id,newTag));
	}
	
	//cerr<<chrTags->size()<<"junction tags loaded"<<endl;
	
	
	fil.close();
}



void GffEntry::Loader::loadSelexaMatches(string filename,int format)
{
	ifstream fil(filename.c_str());
	
	if(!fil.good())
	{
		cerr<<"Cannot open "<<filename<<endl;
		return;
	}
	
	
	///////BAM FORMAT SUPPORT: BLOCK 1
	samfile_t *samfp=NULL;
	bam1_t *bamt=bam_init1();
	
	
	if(format==BAM_FORMAT)
	{
		//close file and open as bam instead
		fil.close();
		if((samfp=samopen(filename.c_str(),"rb",0))==0){
			cerr<<"Cannot open bam file "<<filename<<endl;
			return;
		}
	}
	//////////////////

	set<GffEntry::GBlockPtr> * chrBlock=NULL;
	
	int c=0;

	string curChr="";
	
	set<GffEntry::GBlockPtr>::iterator ib;
	
	int prevPos=-1;
	
	
	////BLOCK 2
	bool readRead=true;

	
	while(true)  //(!fil.eof())
	{
		
		if(format==BAM_FORMAT && samread(samfp,bamt)<0) //if no more to read.
			break;
		
		if(format!=BAM_FORMAT && fil.eof()) //if not, end of file reached, break
			break;
		////
		
		
		GffEntry::SelexaMatchPtr sm=new GffEntry::SelexaMatch;
		
		c++;
		
		//if(c%100==1)
		//	cerr<<"loading SelexaMatch #"<<c<<endl;
		
		
		////BAM format support : BLOCK 3
		
		if(format==BAM_FORMAT)
			readRead=GffEntry::SelexaMatch::readMatchFromBam(samfp,bamt,*sm);
		else
		   readRead=GffEntry::SelexaMatch::readMatch(fil,*sm,format);
		///////

		if(sm->chr=="")
		{
			delete sm;
			break;
		}
		

		if(!readRead)
		{
			delete sm;
			continue;
		}

		if(sm->chr!=curChr)
		{
			prevPos=-1;
			curChr=sm->chr;
			map<string,set<GffEntry::GBlockPtr>* >::iterator i=GffEntry::gblocks.find(curChr);
			if(i==GffEntry::gblocks.end())
			{
				fil.close();
				cerr<<string("chr ")<<curChr<<" not found"<<endl;
				die("");
			}
			else
			{
				chrBlock=(*i).second;
			}	
			
			ib=chrBlock->begin();
		}
		
		if(sm->gpos0<prevPos)
		{
			die("Selexa File Not Sorted, require a sorted file");
		}
		
		while(ib!=chrBlock->end() && (*ib)->getEnd0()<sm->gpos0)
		{

			//  * ib  ]  [ ]  [   ] (read)
			// ib go right until [last block] *end
			// or (read) *ib ]

			ib++;
		}
		
		//here either has reached the end of all blocks in that chromosome [  last block ] *end  (read), or [ (read) *ib ] or (read)  [ *ib ]


		//added 1/10/2010
		prevPos=sm->gpos0;
		///
		



		if(ib==chrBlock->end() || (*ib)->getStart0()>sm->gpos0)
		{

			// (read) [ *ib = first block]
			//or [  ]   (read)   [* ib ]
			//or [  last block ] *end  (read)
			GffEntry::blindSpotSelexaMatches.insert(sm);
			continue;
		}
		
		
		
		GffEntry::GBlockPtr block=(*ib);
		
		//added 4/21/2009: For case where GBlock is not contiguous
		//deleted 1/7/2010 again because it has already been handled above
		/*if(block->getStart0()>sm->gpos0)
		{
			GffEntry::blindSpotSelexaMatches.insert(sm);
		}	
		else*/

		{
			block->selexaMatches.insert(sm);
		//in this case (*ib).start0<=sm->gpos0<=(*ib).end0
		}
	}

	cerr<<c<<" lines of selexa matches loaded"<<endl;
	
	
	
	GffEntry::selexaMatchLoaded=true;
	
	
	
	//cerr<<"blindSportSelexa="<<GffEntry::
	
	////BAM Support BLOCK4
	if(format==BAM_FORMAT){
		bam_destroy1(bamt);
		samclose(samfp);
	}
	else{
		fil.close();
	}
	//////
}


void GffEntry::Loader::loadUniquelyMappablePosition(int readLength,string chrMapRefFile,string prefix,string suffix,bool toGBlock,bool toJnxBlock)
{
	if(GffEntry::isUniquelyMappablePositionLoaded())
	{
		return;
	}
	
	GffEntry::isUMPLoaded=true;
	
	SmartChrMap chrmap(chrMapRefFile);
	
	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlock not loaded: cannot load uniquely mappable position onto");
	}
	
	if(toJnxBlock && !GffEntry::isJnxInfoLoaded())
	{
		die("Jnx info not loaded: cannot load uniquely mappable position onto")
	}
	
	
	
	for(map<string,set<GffEntry::GBlockPtr>* >::iterator chrI=GffEntry::gblocks.begin();chrI!=GffEntry::gblocks.end();chrI++)
	{
		string chr=chrI->first;
		set<GffEntry::GBlockPtr>* gblockset=chrI->second;
		cerr<<"loading uniquely mappable position for chr "<<chr<<endl;
		::GenomeNmersDecoder dec(prefix+chr+suffix);
		map<int, GffEntry::JnxTag*> * jnxtags=NULL;
		GffEntry::JnxTag* jnxtag=NULL;
		vector<GffEntry::JnxTagBlock*>::iterator jblockI;
		
		if(toJnxBlock)
		{
			map<string, map<int,GffEntry::JnxTag*>* >::iterator jti=GffEntry::jnxtags.find(chr);
			
			if(jti!=GffEntry::jnxtags.end())
				jnxtags=jti->second;
			
			
			if(!jnxtags)
			{
					cerr<<"Error: jnxtags not available for current chromosome"<<endl;

			}
							
			
		}
		
		KeyPair<string,int> chrInfo("",-1);
		int umChrID=-1;
		
		
		set<GffEntry::GBlockPtr>::iterator pGBlock=gblockset->begin();
		
		while(dec.readEntry(COMPACTPOSITION))
		{
			const CompactPosition& cpos=dec.getCompactPosition();
			
			if(umChrID!=int(cpos.chrID)) ///Changed unsigned-signed problem
			{
				umChrID=cpos.chrID;
				
				chrInfo=chrmap.getChrMapInfoFromPChrID(cpos.chrID);
				
				if(chrInfo.k1!=chr)
				{
					cerr<<"Error, mispartitioned"<<endl;
					continue;
				}
				
				
				jnxtag=NULL;
				
				if(chrInfo.k2>=0)
				{
					if(!jnxtags)
						continue;
					
					map<int,GffEntry::JnxTag*>::iterator iJnxTag =jnxtags->find(chrInfo.k2);
					if(iJnxTag==jnxtags->end())
					{
						cerr<<"Error: jnxtag ID not found:"<< chrInfo.k2 <<endl;
					}else
					{
						jnxtag=iJnxTag->second;
						jblockI=jnxtag->jblocks.begin();
					}
				}
			}
			
			KeyPair<int,int> bUm11=cpos.getBound11();
			
			if(chrInfo.k2<0 && toGBlock) //incorporate uniquely mapping position to gblocks
			{
				if(pGBlock==gblockset->end())
					continue;
				
				while((*pGBlock)->getEnd1()<=cpos.getEnd1())
				{
					
					if((*pGBlock)->getEnd1()>=cpos.getStart1()) //can they overlap at all?
					{
						//get overlap, add unique positions
						KeyPair<int,int> bL11=(*pGBlock)->getBoundAbs11WithinPosForReadLength(readLength);
						KeyPair<int,int> bR11=(*pGBlock)->getboundAbs11OutsidePosForReadLength(readLength);
						
						
						//cerr<<(*pGBlock)->chr<<":"<<(*pGBlock)->getStart1()<<"-"<<(*pGBlock)->getEnd1()<<"("<<bL11.k1<<"-"<<bL11.k2<<","<<bR11.k1<<"-"<<bR11.k2<<")"<<" vs "<<bUm11.k1<<"-"<<bUm11.k2<<"\t";
						
						int lenL=overlaplen1111(bL11,bUm11);
						int lenR=overlaplen1111(bR11,bUm11);
						
						(*pGBlock)->umposL+=lenL;
						(*pGBlock)->umposR+=lenR;
						
						(*pGBlock)->uniquelyMappableBounds.addBound(overlapBound((*pGBlock)->getBound(),cpos.getBound01()));

						//cerr<<"addL="<<lenL<<";"<<"addR="<<lenR<<endl;
						
					}
					
					//advance pGBlock;
					
					//{
					//cerr<<"**** OUTPUT:";
					//GffEntry::GBlock* block=(*pGBlock);
					//int posU=(block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR());
					//int posL=block->getLength();
					//if(posU>posL) 
					//		cerr<<"errrr:";
								
					//cerr<<block->chr<<":"<<block->getStart1()<<"-"<<block->getEnd1()<<"\t"<<posU<<"\t"<<posL<<endl;
					//}
					
					pGBlock++;
					
					if(pGBlock==gblockset->end())
						break;
				}
				
				if(pGBlock==gblockset->end())
					continue;
				
				//here gblock->end > um.end
				
				if((*pGBlock)->getStart1()<=cpos.getEnd1()) //can they overlap at all?
				{
					//find overlap and then advance um (COPY)
					KeyPair<int,int> bL11=(*pGBlock)->getBoundAbs11WithinPosForReadLength(readLength);
					KeyPair<int,int> bR11=(*pGBlock)->getboundAbs11OutsidePosForReadLength(readLength);
					
					
					//cerr<<(*pGBlock)->chr<<":"<<(*pGBlock)->getStart1()<<"-"<<(*pGBlock)->getEnd1()<<"("<<bL11.k1<<"-"<<bL11.k2<<","<<bR11.k1<<"-"<<bR11.k2<<")"<<" vs "<<bUm11.k1<<"-"<<bUm11.k2<<"\t";
					
					int lenL=overlaplen1111(bL11,bUm11);
					int lenR=overlaplen1111(bR11,bUm11);
					
					(*pGBlock)->umposL+=lenL;
					(*pGBlock)->umposR+=lenR;
					
					(*pGBlock)->uniquelyMappableBounds.addBound(overlapBound((*pGBlock)->getBound(),cpos.getBound01()));

					//cerr<<"addL="<<lenL<<";"<<"addR="<<lenR<<endl;
				}
				
			}
			else if(chrInfo.k2>=0 && toJnxBlock) //also incorporate uniquely mapping position to jnxblocks
			{
				//similar to GBlock but now the JBlock is traversed and that jBlock is now pos to pos w.r.t. readStart, not bound
				if(jblockI==jnxtag->jblocks.end())
					continue;
				
				while((*jblockI)->getEnd1()<=cpos.getEnd1()) 
				{
					if((*jblockI)->getEnd1()>=cpos.getStart1()) //can they overlap at all?
					{
						KeyPair<int,int> b11=(*jblockI)->getBound11();
						(*jblockI)->umpos+=overlaplen1111(b11,bUm11);
					}
					
					jblockI++;
					
					if(jblockI==jnxtag->jblocks.end())
						break;
				}
				
				if(jblockI==jnxtag->jblocks.end())
					continue;
				
				//here jblock->end > um.end
				
				if((*jblockI)->getStart1()<=cpos.getEnd1()) //can they overlap at all?
				{
					KeyPair<int,int> b11=(*jblockI)->getBound11();
					(*jblockI)->umpos+=overlaplen1111(b11,bUm11);
				}
				
			}
			
			
		}
		
	}
}



KeyPair<int,int> GffEntry::GBlock::getDensity(int readLength,bool useUniquelyMappblePosition,bool protrude,bool ambigCheck)
{
	KeyPair<int,int> returnVal(0,0);
	
	if(ambigCheck && this->getNaiveAmbiguity())
		return returnVal;
	
	if(protrude)
		returnVal.k1=this->selexaMatches.size();
	else
		returnVal.k1=calculateConstituteExonFreq_inner_numOfReads(this,this->getSelexaMatchesWithStartWithinRelPos(this->getBound00WithinRelPosForReadLength(readLength)));

	//cerr<<"defcall"<<endl;
	
	if(useUniquelyMappblePosition) //wtf! this was wrong?!
	{
//		cerr<<"using uniq pos"<<endl;
		if(protrude)
			returnVal.k2=this->getUniquelyMappablePosL()+this->getUniquelyMappablePosR();
		else
			returnVal.k2=this->getUniquelyMappablePosL();
		
	}
	else
	{

		if(protrude)
			returnVal.k2=this->getLength();
		else
			returnVal.k2=MAX(0,this->getLength()-readLength+1); //changed: 4/21/2009
		
		
	}
	

	
	return returnVal;
}

KeyPair<int,int> GffEntry::Exon::getDensity(int readLength,bool useUniquelyMappablePosition,bool ambigCheck)
{
	
	KeyPair<int,int> result(0,0);
	
	if(!blocks)
		return result;
	
	return GffEntry::GBlock::getDensityOfContigBlocks<vector<GffEntry::GBlockPtr>::iterator,vector<GffEntry::GBlockPtr>::reverse_iterator >(blocks->begin(),blocks->end(),readLength,false,useUniquelyMappablePosition,ambigCheck);
	
}
//change to above code: 4/21/2009
/*KeyPair<int,int> GffEntry::Exon::getDensity(int readLength,bool useUniquelyMappablePosition,bool ambigCheck)
{
	
	
	KeyPair<int,int> result(0,0);
	
	if(!blocks)
		return result;
	
	vector<GffEntry::GBlockPtr>::reverse_iterator i=blocks->rbegin();
	
//	cerr<<"exon getnd: exon coordinate:"<<this->getStart1()<<"-"<<this->getEnd1()<<"\t";
	GffEntry::GBlock* g;
	g=(*i);
//	cerr<<(g->start0+1)<<"-"<<(g->end1)<<"\t";
	
	
	
	if(ambigCheck && (*i)->getNaiveAmbiguity())
	{
		
	}
	else
	{
		result.k1+=calculateConstituteExonFreq_inner_numOfReads((*i)->getSelexaMatchesWithStartWithinRelPos((*i)->getBound00WithinRelPosForReadLength(readLength)));
	
	
		if(useUniquelyMappablePosition)
			result.k2+=(*i)->getUniquelyMappablePosL();
		else
			result.k2+=(*i)->getNaiveMappablePos(readLength);
	}
	
	i++;
	
	for(;i!=blocks->rend();i++)
	{
		g=(*i);
//		cerr<<(g->start0+1)<<"-"<<(g->end1)<<"\t";
		
		
		if(ambigCheck && (*i)->getNaiveAmbiguity())
		{
			
		}else
		{
		
		result.k1+=(*i)->selexaMatches.size();
		
		if(useUniquelyMappablePosition)
			result.k2+=(*i)->getUniquelyMappablePosL()+(*i)->getUniquelyMappablePosR();
		else
			result.k2+=(*i)->getNaiveMappablePos();
		}
	}
	
	cerr<<endl;
	return result;
}*/

/*void NExonGroup::getBoundariesMap(IntExonPtrMMap* uBoundMap,IntExonPtrMMap* dBoundMap)
{


	NExonGroup::ExonIterator i=this->getAllLevelsExonIterator();
	GffEntry::ExonPtr exonNext;
	while((exonNext=i.nextItem()))
	{
		if(uBoundMap)
			uBoundMap->insert(IntExonPtrMMap::value_type(exonNext->start,exonNext));
		
		if(dBoundMap)
			dBoundMap->insert(IntExonPtrMMap::value_type(exonNext->end,exonNext));
	}



}*/

/*void NExonGroup::getSuperGroups(set<NESuperGroupPtr> *leftSGroups,set<NESuperGroupPtr> *rightSGroups)
{
	NExonGroup::NExonGroupIterator i=this->getExonGroupIterator();
	pair<int,NExonGroup::NExonGroupPtr> exongNext;
	while(NExonGroup::NExonGroupIterator::isValidItem(exongNext=i.nextItem()))
	{
		cerr<<this->locusName<<" arriving at exon group "<<exongNext.second->sid<<endl;
		
		if(leftSGroups && exongNext.second->leftSuperGroup){
			leftSGroups->insert(exongNext.second->leftSuperGroup);

		}

		if(rightSGroups && exongNext.second->rightSuperGroup)
		{
			cerr<<this->locusName<<" inserting right super group"<<exongNext.second->rightSuperGroup->getID()<<endl;
			rightSGroups->insert(exongNext.second->rightSuperGroup);
		}
	}
}*/

string NESuperGroup::getID() const
{
	return id;
}

NESuperGroup::NESuperGroup(const string& _id): id(_id)
{
	this->id=_id;

	//bound=firstExonGroup->levelExonBound;

}

void NESuperGroup::add(NExonGroupPtr exongroup)
{
	this->insert(exongroup);

	//::overlapBoundMutFirst(this->bound,exongroup->levelExonBound);
}

const int GffEntry::JnxTagSet::Union=1;
const int GffEntry::JnxTagSet::Intersection=2;


void GffEntry::Locus::assignExonBoundGroups_FlushState(char curState, set<
		GffEntry::ExonPtr>* Opens, map<int,set<GffEntry::ExonPtr>* >& Closes, map<string,set<GffEntry::ExonPtr>* > &leftBoundExons,map<string,set<GffEntry::ExonPtr>* > &rightBoundExons) {


	if (curState == ExonBoundGroup_LEFT && Opens &&  Opens->size()>0) {
		//flush the open state
		string bid=string("L")+StringUtil::str((unsigned int)leftBoundExons.size()+1);



		for (set<GffEntry::ExonPtr>::iterator i = Opens->begin(); i
				!= Opens->end(); i++) {


			GffEntry::ExonPtr exon=*i;
			//cerr<<"flushing exon state open "<<exon->getBound()<<" within leftbound group "<<bid<<endl;
			exon->leftBoundID=bid;
		}

		//Opens.clear(); don't clear because it will be transferred.
		leftBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(bid,Opens));

	} else if (curState == ExonBoundGroup_RIGHT && Closes.size()>0) {



		string bid=string("R")+StringUtil::str((unsigned int)rightBoundExons.size()+1);

		set<GffEntry::ExonPtr>* thisCloses=new set<GffEntry::ExonPtr>;






		for(map<int,set<GffEntry::ExonPtr>*>::reverse_iterator i=Closes.rbegin();i!=Closes.rend();i++)
		{
			NExonGroupPtr neg=new NExonGroup;
			neg->locusName=this->getFirstName();

			//remember this in locus
			int iOpenNum=i->first;
			set<GffEntry::ExonPtr>* thisOpenNumCloses=i->second;
			//iOpenNum is the same as the L number when inserted in LEFT flush action

			if(iOpenNum>0)
			{
				neg->sid=string("L")+StringUtil::str(iOpenNum)+bid;
			}
			else
			{
				neg->sid=string("S")+StringUtil::str(-iOpenNum)+bid;
			}

			this->exongroups.insert(map<string,NExonGroupPtr>::value_type(neg->sid,neg));



			for(set<GffEntry::ExonPtr>::iterator j=thisOpenNumCloses->begin();j!=thisOpenNumCloses->end();j++)
			{
				//add this exon to exon group and rightbound exons
				GffEntry::ExonPtr exon=*j;
				exon->rightBoundID=bid;

				//cerr<<"flushing exon state close "<<exon->getBound()<<" within rightbound group "<<bid<<endl;

				thisCloses->insert(exon);
				neg->insertExon(exon);

			}

			//has to delete this because its contents are channeled but not the whole set transferred.
			delete thisOpenNumCloses;




		}


		rightBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(bid,thisCloses));

		//this needs to clear because it is not transferred.
		Closes.clear();
	}
}


void GffEntry::Locus::assignExonBoundGroups_FlushTerminalState(char flusheeState, set<
		GffEntry::ExonPtr>* TerminalOpens, map<int,set<GffEntry::ExonPtr>* >& TerminalCloses, map<string,set<GffEntry::ExonPtr>* > &terminalLeftBoundExons,map<string,set<GffEntry::ExonPtr>* > &terminalRightBoundExons) {


	if (flusheeState == ExonBoundGroup_TLEFT && TerminalOpens && TerminalOpens->size()>0) {
		//flush the open state
		string bid=string("S")+StringUtil::str((unsigned int)terminalLeftBoundExons.size()+1);



		for (set<GffEntry::ExonPtr>::iterator i = TerminalOpens->begin(); i
				!= TerminalOpens->end(); i++) {


			GffEntry::ExonPtr exon=*i;
			//cerr<<"flushing exon state open "<<exon->getBound()<<" within leftbound group "<<bid<<endl;
			exon->leftBoundID=bid;
		}

		//Opens.clear(); don't clear because it will be transferred.
		terminalLeftBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(bid,TerminalOpens));

	} else if (flusheeState == ExonBoundGroup_TRIGHT && TerminalCloses.size()>0) {



		string bid=string("E")+StringUtil::str((unsigned int)terminalRightBoundExons.size()+1);

		set<GffEntry::ExonPtr>* thisCloses=new set<GffEntry::ExonPtr>;






		for(map<int,set<GffEntry::ExonPtr>*>::reverse_iterator i=TerminalCloses.rbegin();i!=TerminalCloses.rend();i++)
		{
			NExonGroupPtr neg=new NExonGroup;
			neg->locusName=this->getFirstName();

			//remember this in locus
			int iOpenNum=i->first;
			set<GffEntry::ExonPtr>* thisOpenNumCloses=i->second;
			//iOpenNum is the same as the L number when inserted in LEFT flush action

			if(iOpenNum>0)
			{
				neg->sid=string("L")+StringUtil::str(iOpenNum)+bid;
			}
			else
			{
				neg->sid=string("S")+StringUtil::str(-iOpenNum)+bid;
			}


			this->exongroups.insert(map<string,NExonGroupPtr>::value_type(neg->sid,neg));



			for(set<GffEntry::ExonPtr>::iterator j=thisOpenNumCloses->begin();j!=thisOpenNumCloses->end();j++)
			{
				//add this exon to exon group and rightbound exons
				GffEntry::ExonPtr exon=*j;
				exon->rightBoundID=bid;

				//cerr<<"flushing exon state close "<<exon->getBound()<<" within rightbound group "<<bid<<endl;

				thisCloses->insert(exon);
				neg->insertExon(exon);

			}

			//has to delete this because its contents are channeled but not the whole set transferred.
			delete thisOpenNumCloses;




		}


		terminalRightBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(bid,thisCloses));

		//this needs to clear because it is not transferred.
		TerminalCloses.clear();
	}
}

void GffEntry::Locus::assignExonBoundGroups(int jobID,bool haveTerminal) {

	//hash all exons into this hash
	multimap<int, pair<char, GffEntry::ExonPtr> > boundHash;

	//for each exon
	for (set<GffEntry::ExonPtr>::iterator i = this->exonSet.begin(); i
			!= this->exonSet.end(); i++) {
		GffEntry::ExonPtr exon = *i;

		exon->usageFlag=jobID;
		//cerr<<"ebg assignment: hashing exon"<<exon->egsid<<endl;


		char leftClass=ExonBoundGroup_LEFT;

		if(haveTerminal && exon->getInDegree()==0)
		{
			leftClass=ExonBoundGroup_TLEFT;
		}

		char rightClass=ExonBoundGroup_RIGHT;

		if(haveTerminal && exon->getOutDegree()==0)
		{
			rightClass=ExonBoundGroup_TRIGHT;
		}

		boundHash.insert(
				multimap<int, pair<char, GffEntry::ExonPtr> >::value_type(
						exon->getBound().k1 + 1,
						pair<char, GffEntry::ExonPtr> (leftClass,
								exon))); //hash left bound


		boundHash.insert(
				multimap<int, pair<char, GffEntry::ExonPtr> >::value_type(
						exon->getBound().k2, pair<char, GffEntry::ExonPtr> (
								rightClass, exon)));  //hash right bound


	}

	//now preceed from left to right
	char curState = ExonBoundGroup_LEFT; //has to start with open/left //useless statement

	set<GffEntry::ExonPtr>* Opens=NULL;


	set<GffEntry::ExonPtr>* TerminalOpens=NULL;

	map<int, set<GffEntry::ExonPtr>* > TerminalCloses;

	map<GffEntry::ExonPtr,int> opennums;

	map<int,set<GffEntry::ExonPtr>*> Closes;

	map<string,set<GffEntry::ExonPtr>* > leftBoundExons; //L's
	map<string,set<GffEntry::ExonPtr>* > rightBoundExons; //R's
	map<string,set<GffEntry::ExonPtr>* > terminalLeftBoundExons; //S's
	map<string,set<GffEntry::ExonPtr>* > terminalRightBoundExons; //E's


	multimap<int, pair<char, GffEntry::ExonPtr> >::iterator hi =
			boundHash.begin();
	multimap<int, pair<char, GffEntry::ExonPtr> >::iterator hi2;


	//int curOpenNum=1;



	while (hi != boundHash.end()) {



		int coord = hi->first;
		hi2 = hi;

		while (hi2 != boundHash.end() && hi2->first == coord) {
			hi2++;
		}


		//now range [hi,hi2) contains the exons with the same coord-key
		//add these to Opens or Closes according to their mode
		while (hi != hi2) {
			pair<char, GffEntry::ExonPtr> mode_exon = hi->second;
			//cerr<<"receiving "<<mode_exon.first<<"@"<<mode_exon.second->getBound()<<endl;
			map<GffEntry::ExonPtr,int>::iterator open_infoI;
			int openNum;
			set<GffEntry::ExonPtr> *thisOpenNumCloses;
			map<int,set<GffEntry::ExonPtr>* >::iterator thisOpenNumClosesI;
			switch(mode_exon.first)
			{
				case ExonBoundGroup_LEFT:
					if(!Opens)
					{
						Opens=new set<GffEntry::ExonPtr>;
					}
					//cerr<<"insert opennum for "<<mode_exon.second->getBound()<<" = "<<(leftBoundExons.size()+1)<<endl;
					opennums.insert(map<GffEntry::ExonPtr,int>::value_type(mode_exon.second,leftBoundExons.size()+1));
					Opens->insert(mode_exon.second);
				break;

				case ExonBoundGroup_RIGHT:
					open_infoI=opennums.find(mode_exon.second);
					if(open_infoI==opennums.end())
					{
						//this shouldn't happen!
						die_exit("open info not found!!!");
					}

					openNum=open_infoI->second;


					thisOpenNumClosesI=Closes.find(openNum);
					if(thisOpenNumClosesI==Closes.end())
					{
						thisOpenNumCloses=new set<GffEntry::ExonPtr>;
						Closes.insert(map<int,set<GffEntry::ExonPtr>* >::value_type(openNum,thisOpenNumCloses));
					}
					else
					{
						thisOpenNumCloses=thisOpenNumClosesI->second;
					}

					thisOpenNumCloses->insert(mode_exon.second);
				break;

				case ExonBoundGroup_TLEFT:

					//flush TerminalCloses

					assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TRIGHT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);

					if(!TerminalOpens)
					{
						TerminalOpens=new set<GffEntry::ExonPtr>;
					}

					opennums.insert(map<GffEntry::ExonPtr,int>::value_type(mode_exon.second,-1*(terminalLeftBoundExons.size()+1)));
					TerminalOpens->insert(mode_exon.second);

					break;

				case ExonBoundGroup_TRIGHT:

					//flush TerminalOpens
					assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TLEFT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);
					TerminalOpens=NULL;

					open_infoI=opennums.find(mode_exon.second);
					if(open_infoI==opennums.end())
					{
						//this shouldn't happen!
						die_exit("open info not found!!!");
					}

					openNum=open_infoI->second;


					thisOpenNumClosesI=TerminalCloses.find(openNum);
					if(thisOpenNumClosesI==TerminalCloses.end())
					{
						thisOpenNumCloses=new set<GffEntry::ExonPtr>;
						TerminalCloses.insert(map<int,set<GffEntry::ExonPtr>* >::value_type(openNum,thisOpenNumCloses));
					}
					else
					{
						thisOpenNumCloses=thisOpenNumClosesI->second;
					}

					thisOpenNumCloses->insert(mode_exon.second);

					break;

			}

			hi++;
		}

		//now if state==open, but close has something, then state transition!
		//on the other hand, if state==close, but open has something then state transition!
		if (curState == ExonBoundGroup_LEFT && Closes.size() > 0) {

			//an end displaces a terminal start
			assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TLEFT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);
			TerminalOpens=NULL;

			assignExonBoundGroups_FlushState(curState,Opens,Closes,leftBoundExons,rightBoundExons);

			Opens=NULL;
			//complete state transition
			curState = ExonBoundGroup_RIGHT;
		} else if (curState == ExonBoundGroup_RIGHT && Opens && Opens->size() > 0) {

			//a start displaces a terminal end
			assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TRIGHT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);

			assignExonBoundGroups_FlushState(curState,Opens,Closes,leftBoundExons,rightBoundExons);

			//complete state transition
			//curOpenNum++;
			curState = ExonBoundGroup_LEFT;
		}


	}

	//at the end finalize, has to be the close state! because nothing is open-ended!!
	//ignore this check: for single-exon case!~!
	/*if(curState!=ExonBoundGroup_RIGHT)
	{
		die_exit("ending state is not ExonBoundGroup_RIGHT");
	}*/

	//assignExonBoundGroups_FlushState(ExonBoundGroup_LEFT,Opens,Closes);

	//flush everything
	assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TLEFT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);
	assignExonBoundGroups_FlushTerminalState(ExonBoundGroup_TRIGHT,TerminalOpens,TerminalCloses,terminalLeftBoundExons,terminalRightBoundExons);
	assignExonBoundGroups_FlushState(ExonBoundGroup_LEFT,Opens,Closes,leftBoundExons,rightBoundExons);
	assignExonBoundGroups_FlushState(ExonBoundGroup_RIGHT,Opens,Closes,leftBoundExons,rightBoundExons);

	//first check all exon groups that every exons are consistent in their leftBoundID and rightBoundID
	for(map<string,NExonGroup*>::iterator i=this->exongroups.begin();i!=this->exongroups.end();i++)
	{
		NExonGroupPtr exongroup=i->second;

		string leftBID;
		string rightBID;
		for(NExonGroup::ExonI exonI=exongroup->levelExons.begin();exonI!=exongroup->levelExons.end();exonI++)
		{
			GffEntry::ExonPtr exon=*exonI;
			string thisLeftBID=exon->leftBoundID;
			string thisRightBID=exon->rightBoundID;
			//cerr<<"checking exon BID consistency for exon with Bound ID ["<<thisLeftBID<<","<<thisRightBID<<"]"<<endl;
			if(thisLeftBID.length()==0 || thisRightBID.length()==0)
			{
				cerr<<exon->getBound()<<" ";
				die_exit("fatal error, BID not set");
			}




			if(leftBID.length()==0)
			{
				leftBID=thisLeftBID;
				rightBID=thisRightBID;
			}
			else
			{
				if(leftBID!=thisLeftBID || rightBID!=thisRightBID)
				{
					cerr<<exon->getBound()<<" ";
					die_exit("bound ID not consistent");
				}
			}

		}
	}


	//transfer~!:
	for(map<string,set<GffEntry::ExonPtr>* >::iterator i=terminalLeftBoundExons.begin();i!=terminalLeftBoundExons.end();i++)
	{
		leftBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(i->first,i->second));
	}
	//now convert left bound exons and right bound exons into super groups and free
	for(map<string,set<GffEntry::ExonPtr>* >::iterator i=leftBoundExons.begin();i!=leftBoundExons.end();i++)
	{
		string bid=i->first;
		set<GffEntry::ExonPtr>* _exons=i->second;
		//set<NExonGroupPtr> _exongroups;
		NESuperGroup* sgroup=new NESuperGroup(bid);

		for(set<GffEntry::ExonPtr>::iterator j=_exons->begin();j!=_exons->end();j++)
		{
			GffEntry::ExonPtr exon=*j;
			if(!exon->exonGroup)
			{
				die_exit("exon group not set for exon");
			}
			exon->exonGroup->leftSuperGroup=sgroup;
			sgroup->add(exon->exonGroup);
		}


		this->leftSuperGroups.insert(map<string,NESuperGroup*>::value_type(bid,sgroup));



		delete _exons;
	}

	//transfer~!:
	for(map<string,set<GffEntry::ExonPtr>* >::iterator i=terminalRightBoundExons.begin();i!=terminalRightBoundExons.end();i++)
	{
		rightBoundExons.insert(map<string,set<GffEntry::ExonPtr>* >::value_type(i->first,i->second));
	}

	for(map<string,set<GffEntry::ExonPtr>* >::iterator i=rightBoundExons.begin();i!=rightBoundExons.end();i++)
	{
		string bid=i->first;
		set<GffEntry::ExonPtr>* _exons=i->second;
		//set<NExonGroupPtr> _exongroups;
		NESuperGroup* sgroup=new NESuperGroup(bid);

		for(set<GffEntry::ExonPtr>::iterator j=_exons->begin();j!=_exons->end();j++)
		{
			GffEntry::ExonPtr exon=*j;
			if(!exon->exonGroup)
			{
				die_exit("exon group not set for exon");
			}
			exon->exonGroup->rightSuperGroup=sgroup;
			sgroup->add(exon->exonGroup);
		}

		this->rightSuperGroups.insert(map<string,NESuperGroup*>::value_type(bid,sgroup));

		delete _exons;

	}


}

void GffEntry::Locus::printExonBounds(ostream& os)
{
	for(set<GffEntry::ExonPtr>::iterator i=this->exonSet.begin();i!=this->exonSet.end();i++)
	{
		GffEntry::ExonPtr exon=*i;
		os<<this->chr<<"\t"<<exon->getBound().k1<<"\t"<<exon->getBound().k2<<"\t"<<"EXON_"<<this->getFirstName()+"_"+exon->egsid<<endl;
	}
}

