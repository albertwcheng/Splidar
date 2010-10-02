#include "snipIncludes.h"
#include <math.h>
#include "SplidarGraph.h"
#include <iomanip>

inline bool within000(int pos0,const KeyPair<int,int>& bound00)
{
	return pos0>=bound00.k1 && pos0<=bound00.k2;
}



class BlockCountStruct{

	public:
		int count;
		KeyPair<int,int> Bound;
		BlockCountStruct(int _count,const KeyPair<int,int>& _Bound):count(_count),Bound(_Bound){}
		void operator++(int)
		{
			count++;
		}

		int getCount() const
		{
			return count;
		}

		void intersect(const KeyPair<int,int>& _bound)
		{
			Bound.k1=MAX(Bound.k1,_bound.k1);
			Bound.k2=MIN(Bound.k2,_bound.k2);
		}

		int getStart0() const{
			return Bound.k1;
		}
		int getStart1() const{
			return Bound.k1+1;
		}
		int getEnd1() const{
			return Bound.k2;
		}
		int getLength() const{
			return Bound.k2-Bound.k1;
		}
		int getEnd0() const{
			return Bound.k2-1;
		}
	};

class rj_state
{
public:
	bool used;
	int curExon;
	vector<GffEntry::Exon*> path;
	vector<KeyPair<int,int> > blocks;
	int remflank;
	string extra;
	int gBound;
	char strand;
	static int count;

	rj_state(const rj_state& src)
	{
		(*this)=src;
		count++;
	}

	~rj_state()
	{
		count--;
	}

	/*inline rj_state( const rj_state& src): gDown(gUp)
	{
		cerr<<"start"<<endl;
		used=src.used;
		curExon=src.curExon;

		cerr<<"1"<<endl;

		if(src.path.size()>0)
			copy(src.path.begin(),src.path.end(),path.begin());

		cerr<<"2"<<endl;

		if(src.blocks.size()>0)
			copy(src.blocks.begin(),src.blocks.end(),blocks.begin());
		cerr<<"3"<<endl;
		remflank=src.remflank;
		extra=src.extra;
		gUp=src.gUp;
		strand=src.strand;
		cerr<<"end"<<endl;
	}*/
	inline rj_state( char _strand, int _gBound, int _curExon, int _remflank, string _extra=""):
		strand(_strand), gBound(_gBound), curExon(_curExon),remflank(_remflank),extra(_extra), used(false){ count++;}
};

int rj_state::count=0;

class RankExonPair
{
public:
	int rank;
	GffEntry::ExonPtr pEx;
	RankExonPair(int _rank, GffEntry::ExonPtr _pEx):rank(_rank),pEx(_pEx){}
};

class gUNE_state;

class gUNE_states
{
public:

	gUNE_state *R; //the remembered state
	gUNE_state *Q; //the editing state
	set<GffEntry*> sameGffRQ;
	gUNE_states();

	inline void addExonToQ(GffEntry::ExonomeWalker wex);

	inline int curEnd();

	inline void setCurEnd(int x);

	inline void setCurStart(int x);

	inline void outputRQpairs(ostream& os,string chr,int&gstart,int&gend,RandomAccessFile&raf,ConservedBlockQueue&cbq);

	inline void shift();

	~gUNE_states();

};

class gUNE_state
{
public:

	map<GffEntry*,RankExonPair> mGE;
	int curEnd;
	int curStart;


	gUNE_state():curEnd(-1),curStart(-1)
	{

	}
	inline bool isEmpty()
	{
		return mGE.size()==0;
	}
	inline void reset()
	{
		mGE.clear();
		curEnd=-1;
	}

	void addExonToState(GffEntry::ExonomeWalker wex/*,gUNE_state *R=NULL,gUNE_states* state=NULL*/)
	{
		//replace any exons with the new one if coming from the same transcript, else add it to Queue
		//go through each Gff in wex

		vector<RankedGffEntryPair>& hGffs=gffsOf(wex);
		 GffEntry::ExonPtr pEx=&exonOf(wex);
		//see if a previous exon has been associated with the GffEntry
		//cerr<<hGffs.size()<<endl;
		foreachGff(hGffs,g)
		{
			//cerr<<"he"<<endl;
			map<GffEntry*,RankExonPair>::iterator iGE=mGE.find((*g).gff);
			if(iGE!=mGE.end())
			{
				mGE.erase(iGE);
			}

			mGE.insert(map<GffEntry*,RankExonPair>::value_type((*g).gff,RankExonPair((*g).rank,pEx)));
			/*if(R&&state)
			{
				if(R->mGE.find((*g).gff)!=R->mGE.end())
					state->sameGffRQ.insert(set<GffEntry*>::value_type((*g).gff));
			}*/

		}


	}
};

gUNE_states::~gUNE_states() {
	delete R;
	delete Q;
}
gUNE_states::gUNE_states() {
	R=new gUNE_state;
	Q=new gUNE_state;
}

inline void gUNE_states::addExonToQ(GffEntry::ExonomeWalker wex) {
	Q->addExonToState(wex);

}

inline int gUNE_states::curEnd() {
	return Q->curEnd;
}

inline void gUNE_states::setCurEnd(int x) {
	Q->curEnd=x;
}

inline void gUNE_states::setCurStart(int x) {
	Q->curStart=x;
}



#define COMMONEINFO2 entry.name2+";"+entry.name+";"+entry.id+";"+entry.strand+";"+StringUtil::str(entry.txStart)+";"+StringUtil::str(entry.txEnd)+";"
#define COMMONEINFO entry->name2+";"+entry->name+";"+entry->id+";"+entry->strand+";"+StringUtil::str(entry->txStart)+";"+StringUtil::str(entry->txEnd)+";"
const string TYPE_PROMOTER ("Promoter");
const string TYPE_3UTR ("3'UTR");
const string TYPE_INTRON ("Intron");
const string TYPE_TXSTART("TxStart");

#define Filldot(x) if(x=="") x=".";

inline void gUNE_states::outputRQpairs(ostream& os, string chr, int &gleft,
		int &gright, RandomAccessFile&raf, ConservedBlockQueue&cbq) {
	bool mtQ = Q->isEmpty();
	bool mtR = R->isEmpty();
	if (mtQ && mtR) //both Q and R are empty
		return;

	if (mtQ)
		Q->curStart = -1;
	if (mtR)
		R->curEnd = -1;

	gleft = R->curEnd + 1;
	gright = Q->curStart - 1; //bound is inclusive so -1.

	os << ">\t" << chr << "\t" << gleft << "\t" << gright << "\t";

	string _53elements;
	string _5elements;
	string _3elements;

	for (map<GffEntry*, RankExonPair>::iterator i = R->mGE.begin(); i
			!= R->mGE.end(); i++) {

		GffEntry *entry = (*i).first;

		int rank = (*i).second.rank;
		map<GffEntry*, RankExonPair>::iterator j;

		if (entry->strand == GffEntry::FORWARD && rank == entry->exonCount - 1) {
			_5elements += COMMONEINFO + TYPE_3UTR + ";" + StringUtil::str(
					(*i).second.pEx->end) + "|";
		} else if (entry->strand == GffEntry::REVERSE && rank == 0) {
			_5elements += COMMONEINFO + TYPE_PROMOTER + ";" + StringUtil::str(
					(*i).second.pEx->end) + "|";

		} else {

			_53elements += COMMONEINFO + TYPE_INTRON;
			int gnextrank;

			if (entry->strand == GffEntry::FORWARD) {
				_53elements += StringUtil::str(rank) + "," + StringUtil::str(
						rank + 1);
				gnextrank = rank + 1;
			} else {
				_53elements += StringUtil::str(rank - 1) + ","
						+ StringUtil::str(rank);
				gnextrank = rank - 1;
			}

			_53elements += ";" + StringUtil::str(
					entry->getExonOfRank(rank)->end) + "," + StringUtil::str(
					entry->getExonOfRank(gnextrank)->start) + "|";

		}

	}

	for (map<GffEntry*, RankExonPair>::iterator i = Q->mGE.begin(); i
			!= Q->mGE.end(); i++) {

		GffEntry *entry = (*i).first;

		int rank = (*i).second.rank;

		if (entry->strand == GffEntry::REVERSE && rank == entry->exonCount - 1) {
			_3elements += COMMONEINFO + TYPE_3UTR + ";" + StringUtil::str(
					(*i).second.pEx->start) + "|";
		} else if (entry->strand == GffEntry::FORWARD && rank == 0) {
			_3elements += COMMONEINFO + TYPE_PROMOTER + ";" + StringUtil::str(
					(*i).second.pEx->start) + "|";
		} else if (R->mGE.find(entry) != R->mGE.end()) {
			//skip, already outputed.
			continue;
		} else {

			_53elements += COMMONEINFO + TYPE_INTRON;
			int gprevrank;
			if (entry->strand == GffEntry::FORWARD) {
				_53elements += StringUtil::str(rank - 1) + ","
						+ StringUtil::str(rank);
				gprevrank = rank - 1;
			} else {
				_53elements += StringUtil::str(rank) + "," + StringUtil::str(
						rank + 1);
				gprevrank = rank + 1;

			}

			_53elements += ";" + StringUtil::str(
					entry->getExonOfRank(gprevrank)->end) + ","
					+ StringUtil::str(entry->getExonOfRank(rank)->start) + "|";

		}

	}

	//fill dots to empty elements for easier visualization in txt format.
	Filldot(_53elements);
	Filldot(_5elements);
	Filldot(_3elements);

	os << _53elements << "\t" << _5elements << "\t" << _3elements << endl;

	raf.transfer(os, gleft - 1, gright); //gright +1 becoz inc->exclusive; becoz annotation start from 0, -1 to both.

	os << endl;

	cbq.output(gleft, gright);

}

inline void gUNE_states::shift() {
	//shift the states, forget the R state, swap pointed targets of Q and R,
	R->reset(); //first reset the R state
	sameGffRQ.clear(); //reset the sameGff set
	//here goes the swapping
	gUNE_state *temp=R;
	R=Q;
	Q=temp;
}

/*void getUltimateGBlockAnnotation(string gfffile,string snipfile)
{
	ifstream gff(gfffile.c_str()); //input file stream for gff file
	ofstream snip(snipfile.c_str()); //output file stream for snip file
	char line[GFF_BUFFER_SIZE]; //the buffer for reading lines from gff file

	GffEntry::resetEverything(); //reset the exonome

	int gstart,gend;


	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file

	if(!gff.good())
	{
		cerr<<"cannot open "<<gfffile<<endl;
		gff.close();
		return;
	}

	while(!gff.eof()) //read all entries from gff file, and add (exon,Gff) to a sorted map exonome by (exon.start,exon.end)
	{

		gff.getline(line,GFF_BUFFER_SIZE);



		vector<string> spliton;  //vector to carry the split tokens from each line of gff file
		StringUtil::split(string(line),string("\t"),spliton);

		if(spliton.size()<1)
			continue;
		if(spliton[0]=="#bin") //ignore commented lines from gff files
			continue;

		GffEntry &entry = *(GffEntry::createGffEntry(spliton)); //construct a GffEntry from the tokens


	}

	//do the real things here!!
	int count=GffEntry::exonome.size();
	cerr<<"Exonome Size for "<<gfffile<<"="<<count<<endl;

	//for each exon in the exonome

	gUNE_states state;

	int c=0;

	for(GffEntry::ExonomeWalker ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{
			//new chr, do something?

			curChr=exonOf(ei).chr; //and then open a new seq

		}

		c++;
		if(c%500==0 || c>=count-10)
			cerr<<"handling:"<<c<<endl;

		//get start and end of current exon
		int estart=exonOf(ei).start;
		int eend=exonOf(ei).end;

		//divide into three cases:

		if(estart>state.curEnd())
		{

			//case D: need to output and start new curExons vector
			state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);
			state.shift();
			//and do similar stuff.
			state.addExonToQ(ei);
			state.setCurStart(estart);
			state.setCurEnd(eend);

		}
		else
		{
			//case B and C
			state.addExonToQ(ei);
			if(eend>state.curEnd()){
				//case B
				state.setCurEnd(eend);
			}

		}

	}


	state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);
	//the end; output the end;
	state.shift();
	state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);

	GffEntry::resetEverything(); //free memory
	gff.close();
	snip.close();
}*/

void getUltimateNonExonic(string gfffile,string chrDir,string snipfile,string conservedblockfile,string chrLabel)
{
	ifstream gff(gfffile.c_str()); //input file stream for gff file
	ofstream snip(snipfile.c_str()); //output file stream for snip file
	char line[GFF_BUFFER_SIZE]; //the buffer for reading lines from gff file

	ifstream cbfin(S(conservedblockfile)); //the conserved block file stream (input)

	ConservedBlockQueue cbq(cbfin,snip); //the conserved block file reader (input from the conserved file stream, output to the snip file)

	GffEntry::resetEverything(); //reset the exonome

	int gstart,gend;


	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file

	if(!gff.good())
	{
		curChr=chrLabel;
		cerr<<"no_info for "<<chrLabel<<endl;
		raf=new RandomAccessFile(chrDir+curChr+".seq");
		snip<<">\t"<<CHROMOSOME_NOINFO<<"\t"<<chrLabel;
		snip<<endl;
		raf->transfer(snip,0,-1);
		snip<<endl;
		cbq.output(0,INT_MAX);
		snip<<endl;
		snip.close();
		raf->close();
		delete raf;
		gff.close();
		return;
	}
	while(!gff.eof()) //read all entries from gff file, and add (exon,Gff) to a sorted map exonome by (exon.start,exon.end)
	{

		gff.getline(line,GFF_BUFFER_SIZE);



		vector<string> spliton;  //vector to carry the split tokens from each line of gff file
		StringUtil::split(string(line),string("\t"),spliton);

		if(spliton.size()<1)
			continue;
		if(spliton[0]=="#bin") //ignore commented lines from gff files
			continue;

		/*GffEntry &entry = *(*/GffEntry::createGffEntry(spliton)/*)*/; //construct a GffEntry from the tokens


	}

	//do the real things here!!
	int count=GffEntry::exonome.size();
	cerr<<"Exonome Size for "<<chrLabel<<"="<<count<<endl;

	//for each exon in the exonome


	gUNE_states state;

	int c=0;

	for(GffEntry::ExonomeWalker ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			if(raf) //remember to close and delete the previous RandomAccessFile Object
			{
				raf->close();
				delete raf;
			}


			curChr=exonOf(ei).chr; //and then open a new seq
			raf=new RandomAccessFile(chrDir+curChr+".seq");

		}

		c++;
		if(c%500==0 || c>=count-10)
			cerr<<"handling:"<<c<<endl;

		//get start and end of current exon
		int estart=exonOf(ei).start;
		int eend=exonOf(ei).end;

		//divide into three cases:

		if(estart>state.curEnd())
		{

			//case D: need to output and start new curExons vector
			state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);

			state.shift();

			//and do similar stuff.
			state.addExonToQ(ei);

			state.setCurStart(estart);

			state.setCurEnd(eend);

		}
		else
		{

			//case B and C
			state.addExonToQ(ei);
			if(eend>state.curEnd()){
				//case B
				state.setCurEnd(eend);
			}

		}


		//cerr<<"loopend"<<endl;





	}


	state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);
	//the end; output the end;
	state.shift();
	state.outputRQpairs(snip,curChr,gstart,gend,*raf,cbq);


	if(raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
	{
		raf->close();
		delete raf;
	}


	GffEntry::resetEverything(); //free memory
	gff.close();
	snip.close();
}

#define MAX(a,b) (((a)>(b))?(a):(b))

void writeExonXRef(string dstFile,string gfffile, bool gffSrcFile, bool resetEverything)
{
	ofstream dst(dstFile.c_str());

	if(gfffile!="")
	{
		GffEntry::resetEverything();
			GffEntry::Loader loader(NULL);

				if(gffSrcFile)
				{
					cerr<<"try to load "<<endl;
					loader.loadGffFiles(gfffile);
				}else
				{
					loader.loadGffFile(gfffile,"",true);
				}
	}else
	{
		cerr<<"using loaded Gffs"<<endl;
	}

	cerr<<"exonome size="<<GffEntry::exonome.size()<<endl;
	GffEntry::ExonomeWalker ei;

	for(ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{
		const GffEntry::Exon& exon=exonOf(ei);
		dst<<exon.exonID<<"\t"<<exon.chr<<":"<<exon.start<<"-"<<exon.end
		<<"\t"<<exon.chr<<"\t"<<exon.start<<"\t"<<exon.end<<endl;

	}
	if(resetEverything)
		GffEntry::resetEverything();

	dst.close();
}

void writeTranscriptXRef(string dstFile,string gfffile, bool gffSrcFile, bool resetEverything)
{
	ofstream dst(dstFile.c_str());

	if(gfffile!="")
	{
		GffEntry::resetEverything();
			GffEntry::Loader loader(NULL);

				if(gffSrcFile)
				{
					cerr<<"try to load "<<endl;
					loader.loadGffFiles(gfffile);
				}else
				{
					loader.loadGffFile(gfffile,"",true);
				}
	}else
	{
		cerr<<"using loaded Gffs"<<endl;
	}

	cerr<<"transcriptome size="<<GffEntry::transcriptome.size()<<endl;




	for(vector<GffEntry*>::iterator i=GffEntry::transcriptome.begin();
		i!=GffEntry::transcriptome.end();i++)
	{
		const GffEntry& transcript=*(*i);
		dst<<transcript.localID<<"\t"<<transcript.annoSource<<":"<<transcript.transcriptID<<"\t"<<transcript.exons[0]->chr<<endl;

		//<<":"<<exon.start<<"-"<<exon.end
		//<<"\t"<<exon.chr<<"\t"<<exon.start<<"\t"<<exon.end<<endl;

	}
	if(resetEverything)
		GffEntry::resetEverything();

	dst.close();
}

class IsoNameGroup
{
public:
	vector<GffEntry*> members;
	bool visited;
	IsoNameGroup(bool _visited=false): visited(_visited)
	{
		
	}
};

void getLoci(string gfffile,string snipfile, string siblingVectorFile ,bool gffSrcFile,string gffLabel,bool linkViaName2, bool linkViaExon, bool linkViaExonGroup ,bool writeExonXref,bool writeTranscriptXref,bool resetEverything)
{

	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream sib(siblingVectorFile.c_str());


//	int gstart,gend;

	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file
	int numExonsSib=0;
	int singleton=0;
	int maxGroupSize=0;
	int numVector=0;



	GffEntry::Loader loader(NULL);

	if (gfffile != "") {
		GffEntry::resetEverything(); //reset the exonome
		if (gffSrcFile) {

			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, gffLabel, true);
		}
	}


	//do the real things here!!
	int count=GffEntry::exonome.size();
	cerr<<"Exonome Size "<<"="<<count<<endl;

	//for each exon in the exonome

	int c=0;

	vector<GffEntry::ExonPtr>* sibling=NULL;

	int curEnd=-1;
	int curStart=-1;
	char curStrand=GffEntry::UNKNOWN;
	int estart,eend,estrand;

	int pbs=0;

	int kadd=0;
	GffEntry::ExonomeWalker ei;

	int strandThisTime;

	strandThisTime=GffEntry::FORWARD;
	
	map<string,IsoNameGroup* > isoNameTranscripts;
	
	for(vector<GffEntry*>::iterator i=GffEntry::transcriptome.begin();i!=GffEntry::transcriptome.end();i++)
	{
		GffEntry* entry=*i;
		map<string,IsoNameGroup* >::iterator isoNameGroupI=isoNameTranscripts.find(entry->name2);
		IsoNameGroup* isoNameGroup;
		if(isoNameGroupI==isoNameTranscripts.end())
		{
			isoNameGroup=new IsoNameGroup;
			isoNameTranscripts.insert(map<string,IsoNameGroup* >::value_type(entry->name2,isoNameGroup));
		}else
			isoNameGroup=(*isoNameGroupI).second;
		
		if(isoNameGroup->members.size()>0)
			if(entry->chrom!=isoNameGroup->members[0]->chrom)
			{
				cerr<<"Error :diff chr for gene existence " <<entry->name2<<" chr[0]= "<< isoNameGroup->members[0]->chrom <<" cur chr= "<<entry->chrom<<endl;
			}
		
		isoNameGroup->members.push_back(entry);
	}
	
	
	for(ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(exonOf(ei).siblings)
			continue; //connected already

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			curChr=exonOf(ei).chr; //and then open a new seq
			cerr<<"handling chr "<<curChr<<endl;
			curEnd=-1;
			curStart=-1;
//			char curStrand=GffEntry::UNKNOWN;
		}

		 estart=exonOf(ei).start;
		 eend=exonOf(ei).end;
		 
		 estrand= exonOf(ei).strand;//((*((exonOf(ei)).assTranscripts))[0]).gff->strand;
				
		 if(estrand==GffEntry::BOTH_STRANDS)
		 {
			 cerr<<"Exon of both strands\t"<<curChr<<":"<<(estart+1)<<"-"<<eend<<endl;
		 }
		 if(estrand!=GffEntry::BOTH_STRANDS && estrand!=strandThisTime)
		 {
			 //cerr<<"Exon of both strands\t"<<curChr<<":"<<(estart+1)<<"-"<<eend<<endl;
			 continue;
		 }
		//divide into three cases:

		if(estart>curEnd)
		{

			//case A:

			if(sibling)
			{
				kadd+=sibling->size();


				if(sibling->size()==1)
				{
					//numExonsSib--;
					singleton++;
					//cerr<<"singleton"<<endl;
				}
				else
				{
					numVector++;
					int ss=sibling->size();
					//if(ss>100)
					//{
					//snip<<ss<<"\t"<<curChr<<":"<<curStart<<"-"<<curEnd<<endl;
					//}
					maxGroupSize=MAX(ss,maxGroupSize);
					numExonsSib+=ss;
					for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
					{
						if((*i)->siblings)
							cerr<<"error: sibling exists already";

						(*i)->siblings=sibling;

					}
				}

			}else
			{
				cerr<<"NULL"<<endl;
			}





			sibling=new vector<GffEntry::ExonPtr>;

			GffEntry::siblingVectors.push_back(sibling);


			sibling->push_back((&exonOf(ei)));
			pbs++;
			curEnd=eend;
			curStart=estart;
			curStrand=estrand;
		}
		else
		{

			//case B and C
			sibling->push_back((&exonOf(ei)));
			pbs++;

			if(eend>curEnd)
				curEnd=eend;

		}

		c++;

		//cerr<<"loopend"<<endl;





	}


	if(sibling)
	{


		kadd+=sibling->size();
		numVector++;
		int ss=sibling->size();
		maxGroupSize=MAX(ss,maxGroupSize);
		numExonsSib+=ss;

		for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
		{
			(*i)->siblings=sibling;
		}


	}



	strandThisTime=GffEntry::REVERSE;
	curEnd=-1;
	curStart=-1;
	sibling=NULL;

	for(ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(exonOf(ei).siblings)
			continue; //connected already

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			curChr=exonOf(ei).chr; //and then open a new seq
			cerr<<"handling chr "<<curChr<<endl;
			curEnd=-1;
			curStart=-1;
		//	char curStrand=GffEntry::UNKNOWN;
		}

		 estart=exonOf(ei).start;
		 eend=exonOf(ei).end;
		 
		 
		 estrand= exonOf(ei).strand;//((*((exonOf(ei)).assTranscripts))[0]).gff->strand;
		 
		 if(estrand!=GffEntry::BOTH_STRANDS && estrand!=strandThisTime)
			 continue;

		//divide into three cases:

		if(estart>curEnd)
		{

			//case A:

			if(sibling)
			{
				kadd+=sibling->size();


				if(sibling->size()==1)
				{
					//numExonsSib--;
					singleton++;
					//cerr<<"singleton"<<endl;
				}
				else
				{
					numVector++;
					int ss=sibling->size();
					//snip<<ss<<"\t"<<curChr<<":"<<curStart<<"-"<<curEnd<<endl;
					maxGroupSize=MAX(ss,maxGroupSize);
					numExonsSib+=ss;
					for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
					{
						if((*i)->siblings)
							cerr<<"error: sibling exists already";

						(*i)->siblings=sibling;

					}
				}

			}else
			{
				cerr<<"NULL"<<endl;
			}





			sibling=new vector<GffEntry::ExonPtr>;

			GffEntry::siblingVectors.push_back(sibling);


			sibling->push_back((&exonOf(ei)));
			pbs++;
			curEnd=eend;
			curStart=estart;
			curStrand=estrand;
		}
		else
		{

			//case B and C
			sibling->push_back((&exonOf(ei)));
			pbs++;

			if(eend>curEnd)
				curEnd=eend;

		}

		c++;

		//cerr<<"loopend"<<endl;





	}


	if(sibling)
	{


		kadd+=sibling->size();
		numVector++;
		int ss=sibling->size();
		maxGroupSize=MAX(ss,maxGroupSize);
		numExonsSib+=ss;

		for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
		{
			(*i)->siblings=sibling;
		}


	}
	//cerr<<" #sibling vectors="<<GffEntry::siblingVectors.size()<<endl;
	cerr<<" exonome size="<<GffEntry::exonome.size()<<","<<count<<";handled="<<c<<endl;
	cerr<<" push-backs="<<pbs<<","<<kadd<<endl;
	cerr<<" #exons involved="<<numExonsSib<<endl;
	cerr<<" avg group size=" << (double(numExonsSib)/numVector)<<endl;
	cerr<<" max group size="<< maxGroupSize<<endl;
	cerr<<" #singleton ="<<singleton<<endl;
	
	int sibID=0;
	//output sibling vectors
	for(vector< vector<GffEntry::ExonPtr>* >::iterator i=GffEntry::siblingVectors.begin();
		i!=GffEntry::siblingVectors.end();
		i++)
	{
		sib<<"["<<"\t";
		
		vector<GffEntry::ExonPtr>& vep=**i;
		for(vector<GffEntry::ExonPtr>::iterator j=vep.begin();j!=vep.end();j++)
		{
			sib<<(*j)->exonID<<"\t";
			(*j)->siblingID=sibID;
		}
		
		sib<<"]"<<endl;
		sibID++;
	}
	
	
	
	//now do graph algorithms. BFS to all GffEntry;
	int procno=0;
	int linesoutput=0;

	list <GffEntry*> tasks(GffEntry::transcriptome.begin(),GffEntry::transcriptome.end());
	
	string prevChr="";
	
	multimap<GffEntry::ExonPtr, GffEntry::Locus*> * curMapLocus=NULL;
	
	while(!tasks.empty())
	{
		GffEntry *src=tasks.front();
		//tasks.pop_front();

		GffEntry::Locus *locus=new GffEntry::Locus;
		locus->chr=src->exons[0]->chr;
		
		if(locus->chr!=prevChr)
		{
			prevChr=locus->chr;
			map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::iterator i=GffEntry::loci.find(prevChr);
			if(i==GffEntry::loci.end())
			{
				curMapLocus=new multimap<GffEntry::ExonPtr,GffEntry::Locus*>;
				GffEntry::loci.insert(map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::value_type(prevChr,curMapLocus));
			}
			else
				curMapLocus=(*i).second;
		}
		

		queue<GffEntry*> Q;
		Q.push(src);
		char strand=src->strand;

		while(!Q.empty())
		{
			GffEntry* u=Q.front();
			u->visited=true;
			Q.pop();
			procno++;


			tasks.remove(u);
			if(procno%1000==1)
			{
				cerr<<"processing transcript "<<procno<<" remaining: "<<tasks.size()<<	endl;
			}
			locus->addTranscript(u);

			//now traverse transcript with same names;
			
			if (linkViaName2) {
				map<string,IsoNameGroup*>::iterator isoNameGroupI=
						isoNameTranscripts.find(u->name2);

				if (isoNameGroupI!=isoNameTranscripts.end()) {
					IsoNameGroup *isoNameGroup=(*isoNameGroupI).second;
					if (!isoNameGroup->visited) {
						for (vector<GffEntry*>::iterator it=
								isoNameGroup->members.begin(); it
								!=isoNameGroup->members.end(); it++) {
							GffEntry* tr=*it;
							if (tr->visited)
								continue;

							tr->visited=true;
							Q.push(tr);

						}

						isoNameGroup->visited=true;
					}
				}
			}
			
			//traverse exons
			

				int exonCount=u->exonCount;
				for (int i=0; i<exonCount; i++) {

				GffEntry::ExonPtr exon=u->exons[i];

				if (exon->visited)
					continue;

				exon->visited=true;
				if (linkViaExon) {
					vector<RankedGffEntryPair> *assTranscripts=
							exon->assTranscripts;
					for (vector<RankedGffEntryPair>::iterator j=
							assTranscripts->begin(); j!=assTranscripts->end(); j++) {
						GffEntry* v=(*j).gff;
						if (v->strand!=strand) {
							cerr<<"1) not same strand: ignore: "<<src->name2
									<<" vs "<<v->name2<<endl;
							continue;
						}
						if (!v->visited) {
							v->visited=true;
							Q.push(v);
						}

					}
				}

				//traverse the overlapped Exons
				if (linkViaExonGroup && exon->siblings) {
					for (vector<GffEntry::ExonPtr>::iterator j =
							exon->siblings->begin(); j != exon->siblings->end(); j++) {
						GffEntry::Exon *nei = (*j);
						if (nei->visited)
							continue;
						nei->visited = true;
						vector<RankedGffEntryPair> *assTranscripts =
								nei->assTranscripts;
						for (vector<RankedGffEntryPair>::iterator k =
								assTranscripts->begin(); k
								!= assTranscripts->end(); k++) {
							GffEntry* v = (*k).gff;
							if (v->strand != strand) {
								cerr<<"2) not same strand: ignore: "
										<<src->name2<<" vs "<<v->name2 <<endl;
								continue;
							}
							if (!v->visited) {
								v->visited = true;
								Q.push(v);
							}
						}
					}
				} //fi LinkViaExonGroup

			}
			





		}


		locus->finalize();

		//GffEntry::loci
		curMapLocus->insert(multimap<GffEntry::ExonPtr,GffEntry::Locus* >::value_type(locus->getHeadExon(),locus));
		
		snip<<(*locus);
		linesoutput++;
		if(linesoutput%100==1)
			cerr<<"outputing locus "<<linesoutput<<endl;



	}


	cerr<<"total loci="<<GffEntry::loci.size()<<endl;

	if(writeExonXref)
	{
		cerr<<"Writing Exon Xref to "<<snipfile<<".exref"<<endl;
		writeExonXRef(snipfile+".exref");
	}


	if(writeTranscriptXref)
	{
		cerr<<"Writing Transcript Xref to "<<snipfile<<".txref"<<endl;
		writeTranscriptXRef(snipfile+".txref");
	}

	if(resetEverything)
		GffEntry::resetEverything(); //free memory
	
	for(map<string,IsoNameGroup*>::iterator i=isoNameTranscripts.begin();i!=isoNameTranscripts.end();i++)
	{
		delete (*i).second;
	}
	
	
	//gff.close();
	sib.close();
	snip.close();
}


void _getLoci(string gfffile,string snipfile, string siblingVectorFile ,bool gffSrcFile,string gffLabel,bool writeExonXref,bool writeTranscriptXref,bool resetEverything)
{

	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream sib(siblingVectorFile.c_str());


//	int gstart,gend;

	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file
	int numExonsSib=0;
	int singleton=0;
	int maxGroupSize=0;
	int numVector=0;



	GffEntry::Loader loader(NULL);

	if (gfffile != "") {
		GffEntry::resetEverything(); //reset the exonome
		if (gffSrcFile) {

			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, gffLabel, true);
		}
	}


	//do the real things here!!
	int count=GffEntry::exonome.size();
	cerr<<"Exonome Size "<<"="<<count<<endl;

	//for each exon in the exonome

	int c=0;

	vector<GffEntry::ExonPtr>* sibling=NULL;

	int curEnd=-1;
	int curStart=-1;
	char curStrand=GffEntry::UNKNOWN;
	int estart,eend,estrand;

	int pbs=0;

	int kadd=0;
	GffEntry::ExonomeWalker ei;

	int strandThisTime;

	strandThisTime=GffEntry::FORWARD;

	for(ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(exonOf(ei).siblings)
			continue; //connected already

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			curChr=exonOf(ei).chr; //and then open a new seq
			cerr<<"handling chr "<<curChr<<endl;
			curEnd=-1;
			curStart=-1;
//			char curStrand=GffEntry::UNKNOWN;
		}

		 estart=exonOf(ei).start;
		 eend=exonOf(ei).end;
		 
		 estrand= exonOf(ei).strand;//((*((exonOf(ei)).assTranscripts))[0]).gff->strand;
				
		 if(estrand==GffEntry::BOTH_STRANDS)
		 {
			 cerr<<"Exon of both strands\t"<<curChr<<":"<<(estart+1)<<"-"<<eend<<endl;
		 }
		 if(estrand!=GffEntry::BOTH_STRANDS && estrand!=strandThisTime)
		 {
			 //cerr<<"Exon of both strands\t"<<curChr<<":"<<(estart+1)<<"-"<<eend<<endl;
			 continue;
		 }
		//divide into three cases:

		if(estart>curEnd)
		{

			//case A:

			if(sibling)
			{
				kadd+=sibling->size();


				if(sibling->size()==1)
				{
					//numExonsSib--;
					singleton++;
					//cerr<<"singleton"<<endl;
				}
				else
				{
					numVector++;
					int ss=sibling->size();
					//if(ss>100)
					//{
					//snip<<ss<<"\t"<<curChr<<":"<<curStart<<"-"<<curEnd<<endl;
					//}
					maxGroupSize=MAX(ss,maxGroupSize);
					numExonsSib+=ss;
					for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
					{
						if((*i)->siblings)
							cerr<<"error: sibling exists already";

						(*i)->siblings=sibling;

					}
				}

			}else
			{
				cerr<<"NULL"<<endl;
			}





			sibling=new vector<GffEntry::ExonPtr>;

			GffEntry::siblingVectors.push_back(sibling);


			sibling->push_back((&exonOf(ei)));
			pbs++;
			curEnd=eend;
			curStart=estart;
			curStrand=estrand;
		}
		else
		{

			//case B and C
			sibling->push_back((&exonOf(ei)));
			pbs++;

			if(eend>curEnd)
				curEnd=eend;

		}

		c++;

		//cerr<<"loopend"<<endl;





	}


	if(sibling)
	{


		kadd+=sibling->size();
		numVector++;
		int ss=sibling->size();
		maxGroupSize=MAX(ss,maxGroupSize);
		numExonsSib+=ss;

		for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
		{
			(*i)->siblings=sibling;
		}


	}



	strandThisTime=GffEntry::REVERSE;
	curEnd=-1;
	curStart=-1;
	sibling=NULL;

	for(ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{

		if(exonOf(ei).siblings)
			continue; //connected already

		if(curChr!=exonOf(ei).chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			curChr=exonOf(ei).chr; //and then open a new seq
			cerr<<"handling chr "<<curChr<<endl;
			curEnd=-1;
			curStart=-1;
		//	char curStrand=GffEntry::UNKNOWN;
		}

		 estart=exonOf(ei).start;
		 eend=exonOf(ei).end;
		 
		 
		 estrand= exonOf(ei).strand;//((*((exonOf(ei)).assTranscripts))[0]).gff->strand;
		 
		 if(estrand!=GffEntry::BOTH_STRANDS && estrand!=strandThisTime)
			 continue;

		//divide into three cases:

		if(estart>curEnd)
		{

			//case A:

			if(sibling)
			{
				kadd+=sibling->size();


				if(sibling->size()==1)
				{
					//numExonsSib--;
					singleton++;
					//cerr<<"singleton"<<endl;
				}
				else
				{
					numVector++;
					int ss=sibling->size();
					//snip<<ss<<"\t"<<curChr<<":"<<curStart<<"-"<<curEnd<<endl;
					maxGroupSize=MAX(ss,maxGroupSize);
					numExonsSib+=ss;
					for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
					{
						if((*i)->siblings)
							cerr<<"error: sibling exists already";

						(*i)->siblings=sibling;

					}
				}

			}else
			{
				cerr<<"NULL"<<endl;
			}





			sibling=new vector<GffEntry::ExonPtr>;

			GffEntry::siblingVectors.push_back(sibling);


			sibling->push_back((&exonOf(ei)));
			pbs++;
			curEnd=eend;
			curStart=estart;
			curStrand=estrand;
		}
		else
		{

			//case B and C
			sibling->push_back((&exonOf(ei)));
			pbs++;

			if(eend>curEnd)
				curEnd=eend;

		}

		c++;

		//cerr<<"loopend"<<endl;





	}


	if(sibling)
	{


		kadd+=sibling->size();
		numVector++;
		int ss=sibling->size();
		maxGroupSize=MAX(ss,maxGroupSize);
		numExonsSib+=ss;

		for(vector<GffEntry::ExonPtr>::iterator i=sibling->begin();i!=sibling->end();i++)
		{
			(*i)->siblings=sibling;
		}


	}
	//cerr<<" #sibling vectors="<<GffEntry::siblingVectors.size()<<endl;
	cerr<<" exonome size="<<GffEntry::exonome.size()<<","<<count<<";handled="<<c<<endl;
	cerr<<" push-backs="<<pbs<<","<<kadd<<endl;
	cerr<<" #exons involved="<<numExonsSib<<endl;
	cerr<<" avg group size=" << (double(numExonsSib)/numVector)<<endl;
	cerr<<" max group size="<< maxGroupSize<<endl;
	cerr<<" #singleton ="<<singleton<<endl;
	
	int sibID=0;
	//output sibling vectors
	for(vector< vector<GffEntry::ExonPtr>* >::iterator i=GffEntry::siblingVectors.begin();
		i!=GffEntry::siblingVectors.end();
		i++)
	{
		sib<<"["<<"\t";
		
		vector<GffEntry::ExonPtr>& vep=**i;
		for(vector<GffEntry::ExonPtr>::iterator j=vep.begin();j!=vep.end();j++)
		{
			sib<<(*j)->exonID<<"\t";
			(*j)->siblingID=sibID;
		}
		
		sib<<"]"<<endl;
		sibID++;
	}
	
	
	
	//now do graph algorithms. BFS to all GffEntry;
	int procno=0;
	int linesoutput=0;

	list <GffEntry*> tasks(GffEntry::transcriptome.begin(),GffEntry::transcriptome.end());
	
	string prevChr="";
	
	multimap<GffEntry::ExonPtr, GffEntry::Locus*> * curMapLocus=NULL;
	
	while(!tasks.empty())
	{
		GffEntry *src=tasks.front();
		//tasks.pop_front();

		GffEntry::Locus *locus=new GffEntry::Locus;
		locus->chr=src->exons[0]->chr;
		
		if(locus->chr!=prevChr)
		{
			prevChr=locus->chr;
			map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::iterator i=GffEntry::loci.find(prevChr);
			if(i==GffEntry::loci.end())
			{
				curMapLocus=new multimap<GffEntry::ExonPtr,GffEntry::Locus*>;
				GffEntry::loci.insert(map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::value_type(prevChr,curMapLocus));
			}
			else
				curMapLocus=(*i).second;
		}
		

		queue<GffEntry*> Q;
		Q.push(src);
		char strand=src->strand;

		while(!Q.empty())
		{
			GffEntry* u=Q.front();
			u->visited=true;
			Q.pop();
			procno++;


			tasks.remove(u);
			if(procno%1000==1)
			{
				cerr<<"processing transcript "<<procno<<" remaining: "<<tasks.size()<<	endl;
			}
			locus->addTranscript(u);



			int exonCount=u->exonCount;
			for(int i=0;i<exonCount;i++)
			{

				 GffEntry::ExonPtr exon=u->exons[i];

				if(exon->visited)
					continue;

				exon->visited=true;

				vector<RankedGffEntryPair> *assTranscripts=exon->assTranscripts;
				for(vector<RankedGffEntryPair>::iterator j=assTranscripts->begin();
					j!=assTranscripts->end();
					j++)
				{
					GffEntry* v=(*j).gff;
					if(v->strand!=strand)
					{
						cerr<<"1) not same strand: ignore: "<<src->name2<<" vs "<<v->name2<<endl;
						continue;
					}
					if(!v->visited)
					{
						v->visited=true;
						Q.push(v);
					}
				}


				//traverse the overlapped Exons
				if (exon->siblings) {
					for (vector<GffEntry::ExonPtr>::iterator j =
							exon->siblings->begin(); j != exon->siblings->end(); j++) {
						GffEntry::Exon *nei = (*j);
						if (nei->visited)
							continue;
						nei->visited = true;
						vector<RankedGffEntryPair> *assTranscripts =
								nei->assTranscripts;
						for (vector<RankedGffEntryPair>::iterator k =
								assTranscripts->begin(); k
								!= assTranscripts->end(); k++) {
							GffEntry* v = (*k).gff;
							if (v->strand != strand) {
								cerr<<"2) not same strand: ignore: "<<src->name2<<" vs "<<v->name2<<endl;
								continue;
							}
							if (!v->visited) {
								v->visited = true;
								Q.push(v);
							}
						}
					}
				}

			}





		}


		locus->finalize();

		//GffEntry::loci
		curMapLocus->insert(multimap<GffEntry::ExonPtr,GffEntry::Locus* >::value_type(locus->getHeadExon(),locus));
		
		snip<<(*locus);
		linesoutput++;
		if(linesoutput%100==1)
			cerr<<"outputing locus "<<linesoutput<<endl;



	}


	cerr<<"total loci="<<GffEntry::loci.size()<<endl;

	if(writeExonXref)
	{
		cerr<<"Writing Exon Xref to "<<snipfile<<".exref"<<endl;
		writeExonXRef(snipfile+".exref");
	}


	if(writeTranscriptXref)
	{
		cerr<<"Writing Transcript Xref to "<<snipfile<<".txref"<<endl;
		writeTranscriptXRef(snipfile+".txref");
	}

	if(resetEverything)
		GffEntry::resetEverything(); //free memory
	//gff.close();
	sib.close();
	snip.close();
}

class Keylet
{
public:
	int k1;
	int k2;
	string extra;
	int JnxID;

	inline Keylet(int _k1,int _k2): k1(_k1), k2(_k2) {
	}

	inline Keylet(int _k1, int _k2, const string& _extra): k1(_k1),k2(_k2),extra(_extra)
	{

	}

	inline bool operator <(const Keylet& obj) const
	{
		if(this->k1==obj.k1)
		{
			if(this->k2==obj.k2)
			{
				return (this->extra<obj.extra);
			}

			return (this->k2<obj.k2);
		}

		return this->k1<obj.k1;

	}

	inline bool operator !=(const Keylet& obj) const
	{
		return !(*this==obj);
	}

	inline bool operator >=(const Keylet& obj) const
	{
		return (*this)>obj || (*this)==obj;
	}

	inline bool operator <=(const Keylet& obj) const
	{
		return (*this)<obj || (*this)==obj;
	}

	inline bool operator ==(const Keylet& obj) const
	{
		return this->k1==obj.k1 && this->k2==obj.k2 && this->extra==obj.extra;
	}

	bool operator >(const Keylet& obj) const
	{
		if(this->k1==obj.k1)
		{
			if(this->k2==obj.k2)
			{
				return (this->extra>obj.extra);
			}
			return (this->k2>obj.k2);
		}

		return this->k1>obj.k1;

	}
};

class RankGffEntryExonPair: public RankedGffEntryPair
{
public:
	 GffEntry::ExonPtr pExon;
	 GffEntry::ExonPtr pExon2;
	inline RankGffEntryExonPair(int _rank, GffEntry::ExonPtr _pExon,  GffEntry::ExonPtr _pExon2, GffEntry* _entry):RankedGffEntryPair(_rank,_entry),pExon(_pExon),pExon2(_pExon2){}


};

/*class BiheadExoPair: public RankGffEntryExonPair
{
public:

	//vector<GffEntry*> transcripts1;
	//vector<GffEntry*> transcripts2;


	const GffEntry::Exon *pExon2;
	inline BiheadExoPair( const GffEntry::Exon*_pExon1, const GffEntry::Exon*_pExon2):
		RankGffEntryExonPair(0,_pExon1,NULL),pExon2(_pExon2){}

};*/

class novoInfo
{
public:
	KeyPair< GffEntry::ExonPtr, GffEntry::ExonPtr> exons;
	bool inFrame;
};

class Mapper
{
public:
	virtual string getValue(string ID)=0;
	virtual void freeMem(){

	}
	virtual ~Mapper(){}
};

class acemblyMapper: Mapper
{
public:
	string getValue(string ID)
	{
		vector<string> splitons;
		StringUtil::split(ID,".",splitons);
		return splitons[0];
	}
};

#define RETURN_CAPITAL 	0
#define RETURN_SMALL	1
#define RETURN_ORI_CASE	2

class OneMapper: Mapper
{
public:
	map<string,string> ma;
	bool caseSensitiveKey;
	bool returnCase;
	string empty;
	inline OneMapper(string filename,bool _caseSensitiveKey=false, bool _returnCase=RETURN_ORI_CASE, string _empty=""):
		caseSensitiveKey(_caseSensitiveKey),returnCase(_returnCase),empty(_empty)
	{
		ifstream f(filename.c_str());
		if(!f.good())
		{
			cerr<<"cannot open file:"<<filename<<endl;
			f.close();

		}

		string key;
		string value;

		while(!f.eof())
		{
			key="";
			value="";
			f>>key;
			f>>value;
			if(key.length()<1)
				continue;

			if(!caseSensitiveKey) //not case sensitive, change all key to CAPITAL;
				key=StringUtil::toUpper(key);

			switch(returnCase)
			{
			case RETURN_CAPITAL:
				value=StringUtil::toUpper(value);
				break;
			case RETURN_SMALL:
				value=StringUtil::toLower(value);
				break;
				//case RETURN_ORI_CASE: do nothing;
			}

			ma.insert(map<string,string>::value_type(key,value));
		}

		f.close();
	}

	string getValue(string ID)
	{
		if(!caseSensitiveKey)
			ID=StringUtil::toUpper(ID);

		map<string,string>::iterator i=ma.find(ID);
		if(i==ma.end())
			return empty;

		return (*i).second;

	}

	void freeMem(){
		delete (OneMapper*)this;
	}
};

class JnxSys
{
	static map<string,Mapper*> mappers;

};

class Keylet2
{
public:
	int jnxID;
	KeyPair<int,int> leftBlock;
	KeyPair<int,int>* rightBlocks;
	int nRightBlocks;
	int gUp;
	int gDown;
	Keylet2(int _gUp,int _gDown,const KeyPair<int,int>&_left, const vector<KeyPair<int,int> >&_right)
		:gUp(_gUp),gDown(_gDown),leftBlock(_left)
	{
		nRightBlocks=_right.size();
		rightBlocks=new KeyPair<int,int>[nRightBlocks];
		for(int i=0;i<nRightBlocks;i++)
			rightBlocks[i]=_right[i];
		
	}

	
	Keylet2(const Keylet2& src)
		:jnxID(src.jnxID),leftBlock(src.leftBlock),nRightBlocks(src.nRightBlocks),gUp(src.gUp),gDown(src.gDown)
	{
		rightBlocks=new KeyPair<int,int>[src.nRightBlocks];
		for(int i=0;i<src.nRightBlocks;i++)
			rightBlocks[i]=src.rightBlocks[i];
		
	}
	
	

	~Keylet2()
	{
		delete[] rightBlocks;
	}
	inline bool operator < (const Keylet2& right) const
	{
		if(gUp==right.gUp)
		{
			if(gDown==right.gDown)
			{
				if(leftBlock==right.leftBlock)
				{
					KeyPair<int,int>* x=rightBlocks;
					KeyPair<int,int>* y=right.rightBlocks;
					
					KeyPair<int,int> *xend=rightBlocks+nRightBlocks;
					KeyPair<int,int> *yend=right.rightBlocks+right.nRightBlocks;
					
					for(;x!=xend;x++,y++)
					{
						if(y==yend) //y exhausted before x
							return false;
						
						if(*x==*y)
						{
							
						}else    //if x,y unequal, compare
							return *x<*y;
					}
					
					return!(y==yend); //if x,y exhaused together, then false. if x exhausted before y, then true
					
				}else
					return leftBlock<right.leftBlock;
				
				
			}else
				return gDown<right.gDown;
			
		}else
			return gUp<right.gUp;
	}
	inline bool operator > (const Keylet2& right) const
	{
		if(gUp==right.gUp)
		{
			if(gDown==right.gDown)
			{
				if(leftBlock==right.leftBlock)
				{
					KeyPair<int,int>* x=rightBlocks;
					KeyPair<int,int>* y=right.rightBlocks;
					
					KeyPair<int,int> *xend=rightBlocks+nRightBlocks;
					KeyPair<int,int> *yend=right.rightBlocks+right.nRightBlocks;
					
					for(;x!=xend;x++,y++)
					{
						if(y==yend) //y exhausted before x
							return true;
						
						if(*x==*y)
						{
							
						}else    //if x,y unequal, compare
							return *x>*y;
					}
					
					return false; //if x,y exhaused together, then false. if x exhausted before y, then false. Anycase false!
					
				}else
					return leftBlock>right.leftBlock;
				
				
			}else
				return gDown>right.gDown;
			
		}else
			return gUp>right.gUp;
	}
	inline bool operator ==(const Keylet2& right) const {
		if (gUp==right.gUp && gDown==right.gDown && leftBlock==right.leftBlock
				&& nRightBlocks==right.nRightBlocks) {

			KeyPair<int,int>* x=rightBlocks;
			KeyPair<int,int>* y=right.rightBlocks;

			KeyPair<int,int> *xend=rightBlocks+nRightBlocks;

			for (; x!=xend; x++, y++) {
				if (*x==*y) {

				} else
					//if x,y unequal, compare
					return false;
			}

			return true;

		} else
			return false;

	}
	
	inline bool operator!=(const Keylet2& right) const{
		return !(*this==right);
	}
	
	inline bool operator>=(const Keylet2&right) const{
		return !(*this<right);
	}
	
	inline bool operator<=(const Keylet2& right) const{
		return !(*this>right);
	}
};

class gjk_path
{
public:
	
	GffEntry::Exon* leftExon;
	vector<GffEntry::Exon*> rightExons;
	inline bool operator <(const gjk_path& right) const
	{
		if(leftExon->exonID==right.leftExon->exonID)
		{	
			
			if(rightExons.size()==right.rightExons.size())
			{
				vector<GffEntry::Exon*>::const_iterator i=rightExons.begin();
				vector<GffEntry::Exon*>::const_iterator j=right.rightExons.begin();
				for(;i!=rightExons.end();i++,j++)
				{
					GffEntry::Exon* ie=(*i);
					GffEntry::Exon* je=(*j);
					if(ie->exonID!=je->exonID)
						return ie->exonID<je->exonID;
				}
				
				return false;
			}
			else
				return rightExons.size()<right.rightExons.size();
		}else
			return leftExon->exonID<right.leftExon->exonID;
	}
	inline bool operator >(const gjk_path& right) const
	{
		if(leftExon->exonID==right.leftExon->exonID)
		{	
			
			if(rightExons.size()==right.rightExons.size())
			{
				vector<GffEntry::Exon*>::const_iterator i=rightExons.begin();
				vector<GffEntry::Exon*>::const_iterator j=right.rightExons.begin();
				for(;i!=rightExons.end();i++,j++)
				{
					GffEntry::Exon* ie=(*i);
					GffEntry::Exon* je=(*j);
					if(ie->exonID!=je->exonID)
						return ie->exonID>je->exonID;
				}
				
				return false;
			}
			else
				return rightExons.size()>right.rightExons.size();
		}else
			return leftExon->exonID>right.leftExon->exonID;		
	}
	inline bool operator==(const gjk_path& right) const
	{
		if(leftExon->exonID==right.leftExon->exonID)
		{	
			
			if(rightExons.size()==right.rightExons.size())
			{
				vector<GffEntry::Exon*>::const_iterator i=rightExons.begin();
				vector<GffEntry::Exon*>::const_iterator j=right.rightExons.begin();
				for(;i!=rightExons.end();i++,j++)
				{
					GffEntry::Exon* ie=(*i);
					GffEntry::Exon* je=(*j);
					if(ie->exonID!=je->exonID)
						return false;
				}
				
				return true;
			}
			else
				return false;
		}else
			return false;	
	} 
	inline bool operator >=(const gjk_path& right) const
	{
		return !((*this)<right);
	}
	
	inline bool operator <=( const gjk_path& right) const
	{
		return !((*this)>right);
	}
	
	inline bool operator !=(const gjk_path& right) const
	{
		return !((*this)==right);
	}
	char strand;
	inline gjk_path(char _strand):strand(_strand){} 
};

typedef SmartPtr<gjk_path> gjk_pathPtr;




inline KeyPair<int,int>  divideAbyB_sorted(set<int>::iterator a1, set<int>::iterator a2,set<int>::iterator b1,set<int>::iterator b2,set<int>& AandB,set<int>& AsubB, bool clear=true)
{
	if(clear)
	{
		AandB.clear();
		AsubB.clear();
		
	}
	
	for(;a1!=a2;a1++)
	{

		
		while(b1!=b2 && *b1<*a1)
		{
			//BsubA
			b1++;					
		}		
		
		if(b1==b2) //end of b
		{
			AsubB.insert(*a1);
			continue;
		}		
		
		if(*a1<*b1)
		{
			//AsubB
			AsubB.insert(*a1);
		}
		else //*a1==*b1
		{
			//AandB
			AandB.insert(*a1);
			b1++;
		}
			
	}
	
	//BsubA by looping b1 until b2; ignored in this case
	
	return KeyPair<int,int>(AandB.size(),AsubB.size());
}





class OverlapState
{
public:
	int start;
	int length;
	int end0;
	int endExon;
	int startOffLength;
	set<int> jnxRowIDs;
	bool active;
	int prev;
	inline void setEnd0(int _endExon=-2,int _end0=-1)
	{
		//cerr<<"In setEnd0:prev="<<prev<<",end0="<<
		if(_end0==-1)
			return;
		
		end0=_end0;
		endExon=_endExon;
	}
	
	inline void activate(int curPos0)	
	{
		if(!active)
		{
			prev=curPos0;
		}
		active=true;
	}
	inline void inactivate(int curExon,int curPos0){
		active=false;
		setEnd0(curExon,curPos0);
		
	}
	
	inline bool walk(int curPos0)
	{
		if(!active || prev==curPos0)
			return false;
		
		length+=curPos0-prev+1;
		prev=curPos0+1;
		return true;
	}
	
	inline bool hasMovedSinceSplit()
	{
		return (length>startOffLength);
	}
	
	
	inline OverlapState(int _start, int _length=0, int _endExon=1, int _end0=0 , int _prev=-1,  bool _active=true)
		:start(_start),length(_length),active(_active),end0(_end0),startOffLength(_length),endExon(_endExon)
	{
		prev=MAX(_start,_prev);
	}
	inline bool operator<(const OverlapState& right) const
	{
		set<int>::iterator j=right.jnxRowIDs.begin();
		
		for(set<int>::iterator i=jnxRowIDs.begin();i!=jnxRowIDs.end();i++,j++)
		{
			if(j==right.jnxRowIDs.end()) //right exhausted before left
				return false;
			
			if(*i==*j) //same, do nothing
			{
				
			}else
				return *i<*j;

		}
		
		return !(j==right.jnxRowIDs.end());//right is still available=> left is exhausted before right.=>true;left right exhausted together=>false
	}
	inline bool operator>(const OverlapState& right) const
	{
		set<int>::iterator j=right.jnxRowIDs.begin();
		
		for(set<int>::iterator i=jnxRowIDs.begin();i!=jnxRowIDs.end();i++,j++)
		{
			if(j==right.jnxRowIDs.end()) //right exhausted before left
				return true;
			
			if(*i==*j) //same, do nothing
			{
				
			}else
				return *i>*j;
			
		}
		
		return false;//return false anyway: right is still available=> left is exhausted before right.=>false;left right exhausted together=>false
	}	
	inline bool operator==(const OverlapState& right)
	{
		if(jnxRowIDs.size()!=right.jnxRowIDs.size())
			return false;
		set<int>::iterator i;
		set<int>::iterator j;
		for(i=jnxRowIDs.begin(),j=right.jnxRowIDs.begin();
			i!=jnxRowIDs.end();i++,j++)
		{
			if(*i!=*j)
				return false;
		}
		
		return true;
	}
	inline bool operator!=(const OverlapState& right)
	{
		return !(*this==right);
	}
	inline bool operator>=(const OverlapState& right)
	{
		return !(*this<right);
	}
	inline bool operator<=(const OverlapState& right)
	{
		return !(*this>right);
	}


};



#define CT_START 0
#define CT_END 1

#define BASE0   0
#define BASE1	1
#define ORI_BASE	2

class JnxRow
{
public:
	const Keylet2 *keylet;
	const vector<gjk_path*>* vgp;
	int cursor;
	int cursorType;
	bool addedOnce;
	KeyPair<int,int> outBound;
	int outBExon;
	bool pendProgress;
	bool reachedEnd;
	
	class outlet
	{
	public:
		
		Keylet2* bounds;
		vector<gjk_path*>* path;
		outlet():bounds(NULL),path(NULL){}
		void free()
		{
			if(bounds)
				delete bounds;
			if(path)
				delete path;
		}
	};
	

	
	inline list<KeyPair<int,int> > searchOutBoundFromLeft(OverlapState* os,int readLength,KeyPair<int,int> &outBoundUpdater, int &outBExonUpdater)
	{
		list<KeyPair<int,int> > ov;
		
		outBoundUpdater.k1=INT_MAX;
		outBoundUpdater.k2=INT_MIN;
		
		int seeda0=os->start;
		int seedb0=outBound.k1-1;
		
//		cerr<<"start left search from "<<seeda0<<" "<<seedb0<<endl;		
		
		int rightBound0=os->end0;
		
		
		if(seeda0>seedb0)
		{	
//			cerr<<" seeda0 > seedb0 "<<endl;
			return ov;
		}	
		
		if(seeda0<=keylet->leftBlock.k1)
		{
			cerr<<"finished: ignored"<<endl;
		}
		
		if(!within000(seeda0,b01tob00(keylet->leftBlock)) || ! within000(seedb0,b01tob00(keylet->leftBlock)))
		{
//			cerr<<"left seeds not in first block"<<endl;
			return ov;
		}
		
		int remainingLength=readLength;
		int firstlen=len01(seedb0,keylet->leftBlock.k2);
//		cerr<<"firstlen="<<firstlen<<endl;
		remainingLength-=firstlen;
		
		ov.push_back(KeyPair<int,int>(seeda0,keylet->leftBlock.k2)); //in 01
		
		int rightI=0;
		
		int farEnd0=0; /*cw*/
		while(remainingLength>0)
		{
			KeyPair<int,int>& tBound=keylet->rightBlocks[rightI];
			int curLen=len01(tBound);
			if(rightI==os->endExon)
			{
//				cerr<<"remainingLength="<<remainingLength<<endl;
				farEnd0=tBound.k1+remainingLength-1;
				if(farEnd0>rightBound0)
				{
					cerr<<"farEnd0>end0"<<endl;
					farEnd0=rightBound0;
				}
				remainingLength=0;
//				cerr<<"1) farEnd0="<<farEnd0<<endl;
				ov.push_back(KeyPair<int,int>(tBound.k1,farEnd0+1)); //01
				break;
			}
			if(remainingLength<=curLen)
			{
				farEnd0=tBound.k1+remainingLength-1;
				if(farEnd0>rightBound0)
				{
//					cerr<<"farEnd0>end0"<<endl;
				}
				remainingLength=0;
//				cerr<<"2) farEnd0="<<farEnd0<<endl;
				ov.push_back(KeyPair<int,int>(tBound.k1,farEnd0+1)); //01
				break;
			}
			
		
			remainingLength-=curLen;
			ov.push_back(tBound);
			rightI++;
		}
		
		/*outBound.k1=os->start;
		
		if(farEnd0>outBound.k2) //protect the written rightEnd;
		{
			this->outBExon=rightI;
//			cerr<<"setting k2=farEnd0:"<<farEnd0<<endl;
			this->outBound.k2=farEnd0;
		}*/
		
		//update later, to request protection of left and right ends
		outBoundUpdater.k1=os->start;
		outBoundUpdater.k2=farEnd0;
		outBExonUpdater=rightI;
		
		
		return ov;
		
	}
	
	inline list<KeyPair<int,int> > searchOutBoundFromRight(OverlapState* os,int readLength,KeyPair<int,int> &outBoundUpdater, int &outBExonUpdater)
	{
		list<KeyPair<int,int> > ov;
		
		outBoundUpdater.k1=INT_MAX;
		outBoundUpdater.k2=INT_MIN;		
		
		int seeda0=outBound.k2+1;
		int seedaExon=outBExon;
		
		int seedb0=os->end0;
		int seedbExon=os->endExon;
		
		if(seeda0+1>=keylet->rightBlocks[seedaExon].k2) //0-1;
		{
			//end
			if(seedaExon==keylet->nRightBlocks-1)
			{
				cerr<<"ended: ignore"<<endl;
				return ov;
			}
			
			seedaExon++;
			seeda0=keylet->rightBlocks[seedaExon].k1; //0=0;
		}
		
		
//		cerr<<"start right search from "<<seeda0<<" "<<seedb0<<endl;
		
	//	int leftBound0=os->start;
		
		
		if(seeda0>seedb0)
		{	
//			cerr<<" seeda0 > seedb0 "<<endl;
			return ov;
		}	
		
		if(seedb0<=outBound.k2)
		{
//			cerr<<"seeb0 <= outBound.k2"<<endl;
			return ov;
		}
//		cerr<<"k"<<endl;
		
		if(seedbExon>seedaExon)
		{
			//at seedbExon
			cerr<<"seedbExon="<<seedbExon<<",seedaExon="<<seedaExon<<endl;
			ov.push_front(KeyPair<int,int>(keylet->rightBlocks[seedbExon].k1,seedb0+1)); //01
			
			
			for(int i=seedbExon-1;i>seedaExon;i--)
			{	cerr<<"doing "<<i<<endl;
				ov.push_front(keylet->rightBlocks[i]);
			}
		
		}
//		cerr<<"l"<<endl;
		int startExon=seedaExon;
		
		
		cerr<<"outBExon="<<outBExon<<endl;
		int remainingLength=readLength;
		cerr<<"m: startExon.k1="<<keylet->rightBlocks[startExon].k1<<",seeda0="<<seeda0<<endl;
		int LEN00=len00(keylet->rightBlocks[startExon].k1,seeda0);
		remainingLength-=LEN00; //00 => len
		if(LEN00<0)
			cerr<<"URGGGG: LEN00<0"<<endl;
		cerr<<"len00"<<LEN00<<endl;
		cerr<<"n: remainingLEngth:"<<remainingLength<<endl;
		if(remainingLength<=0)
		{
			cerr<<"right pass error remaining Length"<<endl;
		}
	//	cerr<<"o"<<endl;
		ov.push_front(KeyPair<int,int>(keylet->rightBlocks[startExon].k1,MIN(seedb0+1,keylet->rightBlocks[startExon].k2))); //01
	//	cerr<<"p"<<endl;
		for(int j=startExon-1;j>=0;j--)
		{
			int len=len01(keylet->rightBlocks[j]);
			
			remainingLength-=len;
			cerr<<"len01="<<len<<",remainingLength="<<remainingLength<<endl;
			
			ov.push_front(keylet->rightBlocks[j]);
			
			if(remainingLength<=0)
			{
				cerr<<"right pass error remaining Length smaller than 0 at "<<j<<endl;
			}			
		}
	//	cerr<<"q"<<endl;
		//now exon1
		int offset0=keylet->leftBlock.k2-remainingLength; //1->0
	//	cerr<<"r"<<endl;
		ov.push_front(KeyPair<int,int>(offset0,keylet->leftBlock.k2));
	//	cerr<<"s"<<endl;
		
		
		//must pass through exon1
		/*
		if(offset0<outBound.k1) //protect the written rightEnd;
		{
			
			outBound.k1=offset0;
		}
		
		outBound.k2=seedb0;
		outBExon=seedbExon;
		*.
		*/
		
		//update later, to request protection of left and right ends
		outBoundUpdater.k1=offset0;
		outBoundUpdater.k2=seedb0;
		outBExonUpdater=seedbExon;
		
		
		
		return ov;
		
	}	
	
	inline JnxRow(const Keylet2*_keylet=NULL,const vector< gjk_path*>* _vgp=NULL,int _minExonSpan=4)
		:outBound(INT_MAX,INT_MIN),keylet(_keylet),vgp(_vgp),cursor(-1),cursorType(CT_START),addedOnce(false),outBExon(-2), pendProgress(false), reachedEnd(false)
	{
		if(keylet && vgp)
			set(keylet,vgp,_minExonSpan);
	}
	
		inline void pendForProgress()
		{
			pendProgress=true;
		}
	
	inline void setOutBound(KeyPair<int,int> _outBound,int _outBExon)
	{
		cerr<<"request "<<_outBound<<";BExon@"<<_outBExon<<" while having "<<outBound<<"; BExon@"<<outBExon<<"\t";
		
		if(_outBound.k1<outBound.k1)
		{
			cerr<<"update k1:"<<_outBound.k1<<"\t";
			outBound.k1=_outBound.k1;
		}
		if(_outBound.k2>outBound.k2) //outBound.k2 and outBExon are updated as a pair
		{
			cerr<<"update k2,BE:"<<_outBound.k2<<","<<_outBExon;
			outBound.k2=_outBound.k2;
			outBExon=_outBExon;
		}
		
		cerr<<endl;
	}
	
	inline bool hasReachedEnd()
	{
		return reachedEnd;
	}
	
	inline /*outlet*/void extract(int a0,int b0,int bExon,KeyPair<int,int>& outBoundNew,int &outBExonNew) //base 0
	{
		//a0 is in the first exon;
		//outlet o;
		KeyPair<int,int> leftBlock(a0,outBound.k1);
		vector<KeyPair<int,int> > rightBlocks;
		
		
		
		int b1=b0+1;
		
		rightBlocks.push_back(KeyPair<int,int>(outBound.k2,(bExon==outBExon)?b1:(keylet->rightBlocks[outBExon].k2)));
		
		for(int t=outBExon+1;t<bExon;t++)
		{
			rightBlocks.push_back(KeyPair<int,int>(keylet->rightBlocks[t].k1,keylet->rightBlocks[t].k2));
		}
		
		if(bExon>outBExon)
		{
			rightBlocks.push_back(KeyPair<int,int>(keylet->rightBlocks[bExon].k1,b1));
		}
		
		
		
		//outBExon=bExon;
		
		//t==bExon
		
		
		
		outBound.k1=a0+1; //base 1 maintained
		outBound.k2=b0; //base 0 maintained
		outBExon=bExon;
		outBoundNew=outBound;
		outBExonNew=outBExon;
		
	}
	

	
	inline void set(const Keylet2*_keylet,const vector<gjk_path*>* _vgp,int minExonSpan)
	{
		keylet=_keylet;
		vgp=_vgp;
		
		outBExon=1;
		outBound.k1=keylet->gUp-minExonSpan+1; //base1 -> base0 ex
		outBound.k2=keylet->gDown+minExonSpan-2;; //base0 ex
		
	}
	inline bool progress()
	{
		if(!pendProgress)
			return true;
		
		pendProgress=false;
		
		switch(cursorType)
		{
		case CT_START:
			cursorType=CT_END;
			return true;
		case CT_END:
			if(cursor==keylet->nRightBlocks-1)
			{
				reachedEnd=true;
				return false;
			}
			cursor++;
			cursorType=CT_START;
			return true;
			//break;
		}
		
		/*dummy*/ return false;
		  
	}
	inline int getEventType()
	{
		return cursorType;
	}
	
	
	
	inline int getEventCoord(int base=BASE0)
	{
		if(reachedEnd)
			return INT_MAX;
		const KeyPair<int,int>*p;
		if(cursor==-1)
		{
			p=&keylet->leftBlock;
		}
		else
		{
			p=&keylet->rightBlocks[cursor];
		}
		
		switch(cursorType)
		{
		case CT_START:	
			//orginally base 0
			return (base==BASE1)?(p->k1+1):p->k1; 
		case CT_END:
			//originally base 1
			return (base==BASE0)?(p->k2-1):p->k2;
		}
		
		/*dummy*/return 0;
	}
};








 bool compareSet(set<int>& set1,set<int>& set2)
{
	if(set1.size()!=set2.size())
		return false;
	
	set<int>::iterator i;
	set<int>::iterator j;
	
	for(i=set1.begin(),j=set2.begin();i!=set1.end();i++,j++)
	{
		if(*i!=*j)
			return false;
	}
	
	return true;
}


inline bool insertIfLonger(set<SmartPtr<OverlapState> > &setToAdd, OverlapState* stateToAdd, bool deleteOldState=true, bool deleteUnusedNewState=true)
{
	set<SmartPtr<OverlapState> >::iterator fit=setToAdd.find(stateToAdd);
	
	if(fit==setToAdd.end()) //prev entry not found, add
	{
		setToAdd.insert(stateToAdd);
		return true;
		
		
	}else if((*fit)->length<stateToAdd->length){ //prev entry found, new longer, remove old
		
		if(deleteOldState)
			delete *fit;
		
		setToAdd.erase(fit);
		setToAdd.insert(stateToAdd);
		return true;
	}else { //new is not longer than the old, ignored.
		
		if(deleteUnusedNewState)
			delete stateToAdd;
		
		return false;
	}
}

class RangeResult
{
public:
	KeyPair<int,int> bound00;
	bool within;
	inline bool isValid()
	{
		return bound00.k1>=0 && bound00.k2>=0;
	}
	inline bool completeExonReturned()
	{
		return isValid() && !within;
	}
	RangeResult(int k1,int k2,bool _within)
		:bound00(k1,k2),within(_within){}
};





inline RangeResult rightRangeOf(int pos0,const KeyPair<int,int>& bound00)
{
	//RangeResult rr;
	if(pos0>bound00.k2)
	{
		return RangeResult(-1,-1,false);
	}else if(pos0<bound00.k1)
	{
		return RangeResult(bound00.k1,bound00.k2,false);
	}else // bound00.k2>=pos0>=bound00.k1
	{
		return RangeResult(pos0,bound00.k2,true);
	}
	
	
}

inline RangeResult leftRangeOf(int pos0,const KeyPair<int,int>& bound00)
{
	//RangeResult rr;
	if(pos0>bound00.k2)
	{
		return RangeResult(bound00.k1,bound00.k2,false);
	}else if(pos0<bound00.k1)
	{
		return RangeResult(-1,-1,false);
	}else // bound00.k2>=pos0>=bound00.k1
	{
		return RangeResult(pos0,bound00.k2,true);
	}
	
	
}
ostream &operator<<(ostream&os, set<int>& inset)
{
	for(set<int>::iterator i=inset.begin();i!=inset.end();i++)
		os<<"\t"<<(*i);
	return os;
}

ostream& operator << (ostream&os,  OverlapState& state)
{
	os<<"state.start,endExon,end,length"<<state.start<<","<<state.endExon<<","<<state.end0<<","<<state.length<<",active="<<(state.active?"true":"false")<<",prev="<<state.prev
		<<"\tstate.jxnrowIDs=";
	os<<state.jnxRowIDs;
	return os;
	//cerr<<endl;
}

ostream& operator <<(ostream &os, set<SmartPtr<OverlapState> >& setOfStates)
{
	cerr<<"set{"<<endl;
	for(set<SmartPtr<OverlapState> >::iterator i=setOfStates.begin();i!=setOfStates.end();i++)
		os<<**i<<endl;
	
	os<<"}"<<endl;
	return os;
}


template<class T>
inline const T& firstOf(set<T>& sset)
{
	return *sset.begin();
}


inline vector<KeyPair<int,int> > reverse_vector(vector<KeyPair<int,int> >& fv)
{
	return vector< KeyPair<int,int> >(fv.rbegin(),fv.rend());
}


inline void gkj_outJnx(ofstream &snipfile, ofstream &jnxinfo, RandomAccessFile* raf, const string&setName, int jID, const string& chr, char strand , int gUp, int gDown, list<KeyPair<int,int> > bounds, vector<gjk_path*> &exonPaths)
{
	if(bounds.size()<1)
		return;
	
//	cerr<<"outjnx"<<endl;
	
	//output the sequence file in fasta format
	snipfile<<">"<<chr<<"|"<<setName<<":"<<jID<<endl;
	
	list<KeyPair<int,int> >::iterator i=bounds.begin();
	
//	cerr<<"1"<<endl;
	
	raf->transfer(snipfile,(*i).k1,(*i).k2);
	
//	cerr<<"2"<<endl;
	snipfile<<endl;
	i++;
//	cerr<<"3"<<endl;	
	for(;i!=bounds.end();i++)
		raf->transfer(snipfile,(*i).k1,(*i).k2);
	snipfile<<endl;
//	cerr<<"4"<<endl;	
	//output the jnxinfo file
	
	jnxinfo<<chr<<"|"<<setName<<":"<<jID<<"\t";
	jnxinfo<<chr<<"\t";
	//jnxinfo<<bounds.size()<<"\t";
	jnxinfo<<strand<<"\t";
	
	int length=0;
//	cerr<<"5"<<endl;	
	for(i=bounds.begin();i!=bounds.end();i++) //stupid, but just make it work. who cares?
		length+=len01(*i);
	
	jnxinfo<<length<<"\t";
	if(length>76 || length<8 )
		cerr<<"length > 76 || < 8 : "<<length<<endl;
	jnxinfo<<bounds.size()<<"\t";
	
	//cerr<<"bbb"<<endl;
	set<gjk_pathPtr> exonPathSetNonRedundant;
	for(vector<gjk_path*>::iterator i=exonPaths.begin();i!=exonPaths.end();i++)
	{
		exonPathSetNonRedundant.insert(*i);
	}
	
	//cerr<<"pbbb"<<endl;
	
	jnxinfo<<exonPathSetNonRedundant.size()<<"\t";
	//jnxinfo<<exonPaths.size()<<"\t";
//	cerr<<"6"<<endl;	
	for(i=bounds.begin();i!=bounds.end();i++)
		jnxinfo<<((*i).k1+1)<<"-"<<(*i).k2<<"\t";   //11
//	cerr<<"7"<<endl;	
	for(/*vector<gjk_path*>*/set<gjk_pathPtr>::iterator j=/*exonPaths*/exonPathSetNonRedundant.begin();j!=/*exonPaths*/exonPathSetNonRedundant.end();j++)
	{
	//	cerr<<"7.1"<<endl;
		gjk_path * path=(*j);
	//	cerr<<"7.2"<<endl;
		if(!path)
		{
			cerr<<"path is invalid"<<endl;
			continue;
		}
		jnxinfo<<path->leftExon->exonID<<"\t";
	//	cerr<<"7.3"<<endl;
		if(path->rightExons.size()<bounds.size()-1)
		{
			cerr<<"wrong size"<<endl;
			break;
		}
	// cerr<<"7.4"<<endl;
		for(unsigned int k=0;k<bounds.size()-1;k++) /////Changed unsigned-signed problem
			jnxinfo<<path->rightExons[k]->exonID<<"\t";
	//	cerr<<"7.5"<<endl;
	}
//	cerr<<"8"<<endl;
	jnxinfo<<";"<<endl;
	
//	cerr<<"out jnx "<<jID<<" complete"<<endl;
	
	
}

inline void copyPathVector(vector<gjk_path*> &dst, const vector<gjk_path*> &src)
{
	for(vector<gjk_path*>::const_iterator i=src.begin();i!=src.end();i++)
		dst.push_back(*i);
}


inline ostream& operator<<(ostream& os, list<KeyPair<int,int> > &ov)
{
	for(list<KeyPair<int,int> >::iterator i=ov.begin();i!=ov.end();i++)
		os<<"\t"<<(*i).k1<<","<<(*i).k2;
	
	return os;
}

inline void gkj_inner_outJnx(ofstream &snipfile, ofstream &jnxinfo, RandomAccessFile* raf, const string& setName, int& jID, int readLength, int minExonSpan, string& curChr, char strand, int gUp, int gDown, int&cumF,KeyPair<int,int>&lastPrimJnx,map<Keylet2,vector<gjk_path*>* >::iterator &sit,int &ambig, int &ambigGroup, int &maxAmbig)
{
	JnxRow* jnxrows=NULL;
	set<SmartPtr<OverlapState> > overlapStates;
	set<SmartPtr<OverlapState> > outQueue;
	
//	cerr<<"gi1"<<endl;
	if(cumF>=2)
	{
		
//		cerr<<"gia1"<<endl;
		jnxrows=new JnxRow[cumF];
//		cerr<<"gia2"<<endl;		
		
		
		ambig+=cumF;
		maxAmbig=MAX(cumF,maxAmbig);
		ambigGroup++;
	
		cerr<<" cumF: "<<cumF<<":"<<curChr<<":"<<(lastPrimJnx.k1)<<"-"<<lastPrimJnx.k2<<endl;
		
		map<Keylet2,vector<gjk_path*>* >::iterator it2=sit;
		
		
		for(int j=0;j<cumF;j++)
		{
			jnxrows[j].set(&(*it2).first,(*it2).second,minExonSpan);
			
			const Keylet2& kk=(*it2).first;
			
		cerr<<"\t"<<(kk.leftBlock.k1+1)<<"-"<<kk.leftBlock.k2<<" ";

		cerr<<curChr<<":"<<(kk.gUp)<<"-"<<kk.gDown;
			
	//		cerr<<"  "<<kk.nRightBlocks;
			for(int k=0;k<kk.nRightBlocks;k++)
			{
				cerr<<" "<<(kk.rightBlocks[k].k1+1)<<"-"<<kk.rightBlocks[k].k2;
			}
		//	cerr<<endl;
			it2++;
		}
		
		//now solve the conflicts:
		
		set<int> in;
		int ended=0;
		
		
		ended=0;
		
		while(ended<cumF)
		{
			
			bool needToAddNewState=false;
			set<int> starts;
			set<int> ends;
			
			
			int curPos0=INT_MAX;
			
			
			for(int j=0;j<cumF;j++) //first pass, find min of coord
			{
				
				int rowCurCoord0=jnxrows[j].getEventCoord(BASE0);
				curPos0=MIN(curPos0,rowCurCoord0);
		

				
			}
			
			cerr<<"**********"<<endl;
			cerr<<"curPos0="<<curPos0;
			
			//second pass, construct signals
			for(int j=0;j<cumF;j++)
			{
				int rowCurCoord0=jnxrows[j].getEventCoord(BASE0);
				if(curPos0==rowCurCoord0)
				{
					jnxrows[j].pendForProgress();
					switch(jnxrows[j].cursorType)
					{
					case CT_START:
						if(!jnxrows[j].addedOnce)
						{
							needToAddNewState=true;
							jnxrows[j].addedOnce=true;
						}
						starts.insert(j);
						in.insert(j);
						break;
					case CT_END:
						ends.insert(j);
						in.erase(j);
						break;
					}
				}
			}
			
			
			//Process start signals
			cerr<<"\tin.size()="<<in.size()<<",overlapStates.size()="<<overlapStates.size()<<endl;
			if(in.size()>0)
			{
				cerr<<"IN:";
				cerr<<in<<endl;
			}
			if(starts.size()>0)
			{
			cerr<<"START:";
			cerr<<starts<<endl;
			}
			if(ends.size()>0)
			{
			cerr<<"ENDS:";
			cerr<<ends<<endl;
			}
			
			if(starts.size()>0 && overlapStates.size()>0)
			{
			cerr<<"process STARTS"<<endl;
			
			for(set<SmartPtr<OverlapState> >::iterator i=overlapStates.begin();
				i!=overlapStates.end();)
			{
				OverlapState* os=(*i);
				
				
				//try splitting!!
				OverlapState* AandB=new OverlapState(os->start,os->length,os->endExon,os->end0,os->prev,os->active);
				OverlapState* AsubB=new OverlapState(os->start,os->length,os->endExon,os->end0,os->prev,os->active);
				
				KeyPair<int,int> sizes=divideAbyB_sorted(os->jnxRowIDs.begin(),os->jnxRowIDs.end(),starts.begin(),starts.end(),AandB->jnxRowIDs,AsubB->jnxRowIDs);
				unsigned int AandBSize=sizes.k1; /////Changed unsigned-signed problem
				unsigned int AsubBSize=sizes.k2; /////Changed unsigned-signed problem
				
				//AandB are those reactivated
				//AsubB are left inactive
				cerr<<"encounter state "<<*os<<endl;
				cerr<<"AandB "<<*AandB<<endl;
				cerr<<"AsubB "<<*AsubB<<endl;
				
				bool eraseI=false;
				
				

				if(AandBSize==os->jnxRowIDs.size()) //complete overlap. all jnxs (re)start
				{
					os->activate(curPos0);
					delete AandB;
					delete AsubB;
					
				}
				else if(AandBSize>0 && AsubBSize>0) //A and B overlaps but not completely, need a split
				{
					//remove the original state first
					
					os->setEnd0(); //everyone should be in STOP state, so nevermind
					
					if(os->hasMovedSinceSplit())
						insertIfLonger(outQueue,os);
					else
						delete os;
					
					cerr<<"OUTQUEUE ";
					cerr<<outQueue;
					
					eraseI=true;
					//overlapStates.erase(i);
					
					if(AandBSize==1)
					{
						delete AandB;
						
					}else
					{
						AandB->activate(curPos0);
						insertIfLonger(overlapStates,AandB);
						
					}
					
				//irrelavent AsubB
					if(AsubBSize==1)
					{
						delete AsubB;
						
					}else
					{
						
						insertIfLonger(overlapStates,AsubB);
					}
					
					
				}else
				{
				//////
					delete AandB;
					delete AsubB;
				//////
				}
				
				if(eraseI)
					overlapStates.erase(i++);
				else
					i++;
				
		
			}
			}
			
			if(overlapStates.size()>0)
			{	
			cerr<<"add"<<endl;
			//here add the length;
			for(set<SmartPtr<OverlapState> >::iterator i=overlapStates.begin();
				i!=overlapStates.end();
				i++)
			{
				OverlapState* os=(*i);
				
				os->walk(curPos0);
				
				cerr<<*os<<endl;
			}
			}
			
			if(ends.size()>0 && overlapStates.size()>0)
			{
			cerr<<"process ENDS"<<endl;
		//Process end signals			
			for(set<SmartPtr<OverlapState> >::iterator i=overlapStates.begin();
				i!=overlapStates.end();)
			{
				bool eraseI=false;
				
				OverlapState* os=(*i);
				
				//try splitting!!
				OverlapState* AandB=new OverlapState(os->start,os->length,os->endExon,os->end0,os->prev,os->active);
				OverlapState* AsubB=new OverlapState(os->start,os->length,os->endExon,os->end0,os->prev,os->active);
				
				KeyPair<int,int> sizes=divideAbyB_sorted(os->jnxRowIDs.begin(),os->jnxRowIDs.end(),ends.begin(),ends.end(),AandB->jnxRowIDs,AsubB->jnxRowIDs);
				unsigned int AandBSize=sizes.k1; //changed signed-unsigned
				unsigned int AsubBSize=sizes.k2; //changed signed-unsigned
				
				cerr<<"encounter state "<<*os<<endl;
				cerr<<"AandB "<<*AandB<<endl;
				cerr<<"AsubB "<<*AsubB<<endl;				
				if(AandBSize==os->jnxRowIDs.size()) //complete overlap. all jnxs end
				{
					int fid=firstOf(os->jnxRowIDs);
					os->inactivate(jnxrows[fid].cursor,curPos0);
					cerr<<"\tinactivate "<<curPos0<<endl;
					delete AandB;
					delete AsubB;
				}
				else if(AandBSize>0 && AsubBSize>0) //A and B overlaps but not completely, need a split
				{
					//remove the original state first
					int fid=firstOf(os->jnxRowIDs);
					os->setEnd0(jnxrows[fid].cursor,curPos0);
					
					if(os->hasMovedSinceSplit())
						insertIfLonger(outQueue,os);
					else
						delete os;
					
				cerr<<"OUTQUEUE ";
					cerr<<outQueue;
						
					eraseI=true;
					//overlapStates.erase(i);
					
					if(AandBSize==1)
					{
						//insertIfLonger(outQueue,AandB);
						delete AandB;
					}else
					{
						fid=firstOf(os->jnxRowIDs);
						AandB->inactivate(jnxrows[fid].cursor,curPos0);
						insertIfLonger(overlapStates,AandB);
					}
					
				//irrelavent AsubB
					if(AsubBSize==1)
					{
						//just let it go not
						//insertIfLonger(outQueue,AsubB);
						delete AsubB;
						
					}else
					{
						insertIfLonger(overlapStates,AsubB);
					}
					
					
				}else
				{
				//////
					delete AandB;
					delete AsubB;
				//////
				}
				
				if(eraseI)
						overlapStates.erase(i++);
				else
						i++;
				
		
			}			
			
			}
			
			
			
			
			if(in.size()>1)
			{
		

				if(needToAddNewState)
				{
					OverlapState* newState=new OverlapState(curPos0);
					newState->jnxRowIDs.insert(in.begin(),in.end());
					newState->activate(curPos0);
					overlapStates.insert(newState);
				}
				
			} 
			
			
			for(int j=0;j<cumF;j++)
				if(!jnxrows[j].progress()) //reached the end this jnxrow
					ended++;	
		
		}//while (ended<cumF)
		
		for(set<SmartPtr<OverlapState> >::iterator i=overlapStates.begin();i!=overlapStates.end();i++)
		{
			//should have been STOP signaled
			OverlapState* os=*i;
			//os->setEnd0(curPos0);
			insertIfLonger(outQueue,os);
			
			cerr<<"OUTQUEUE ";
			cerr<<outQueue;
			//overlapStates.erase(i);
		}
		
		
		multimap<int,SmartPtr<OverlapState> > outQueue2; //int is the number of jnx involved. i.e., jnxRowIDs.length()
		
		
		for(set<SmartPtr<OverlapState> >::iterator i=outQueue.begin();i!=outQueue.end();i++)
		{
			
			SmartPtr<OverlapState> os=*i;
			
			if(os->length>=readLength)	
			{
				outQueue2.insert(multimap<int,SmartPtr<OverlapState> >::value_type(os->jnxRowIDs.size(),os));
			}
			
			//delete os;
		}
		
		KeyPair<int,int> outBoundUpdater;
		int outBExonUpdater;
		
		//now do the real output. and delete the overlapStates 
		cerr<<"the real output"<<endl;
		for(multimap<int,SmartPtr<OverlapState> >::reverse_iterator i=outQueue2.rbegin();i!=outQueue2.rend();i++)
		{
			OverlapState* os=(*i).second;
			//everything would have same coordinate, so get any one;
			int firstRowID=firstOf(os->jnxRowIDs);
			JnxRow& fjr=jnxrows[firstRowID];
			//first one has to be in left exon;
			
			
			const Keylet2* keylet=fjr.keylet;
			if(!within000(os->start,keylet->leftBlock))
			{
				cerr<<"assumption first encounter in left exon broken"<<endl;
				delete os;
				continue;
			}
			
			cerr<<"out overlapStates{"<<endl;
			cerr<<*os<<endl;
			cerr<<"}"<<endl;
			cerr<<"out OverlapJnx{"<<endl;
			
			cerr<<os->start<<" "<<(os->end0+1)<<endl;
			
			cerr<<"search bounds from first jnx row: "<<firstRowID<<endl;
			
			
			

			//os->start and os->end0 are now in base 0
			list<KeyPair<int,int> > ovl=fjr.searchOutBoundFromLeft(os,readLength,outBoundUpdater,outBExonUpdater);
			
			vector<gjk_path*> exonPaths;
			
			for(set<int>::iterator j=os->jnxRowIDs.begin();j!=os->jnxRowIDs.end();j++)
			{
				JnxRow& cjr=jnxrows[*j];
				cerr<<(*j)<<"\t";
				//cerr<<"copying bounds to jnx row: "<<(*j)<<" as "<<fjr.outBExon<<"@"<<fjr.outBound<<endl;
				cjr.setOutBound(outBoundUpdater,outBExonUpdater);
				copyPathVector(exonPaths,*cjr.vgp);
			}
			
			list<KeyPair<int,int> > ovr=fjr.searchOutBoundFromRight(os,readLength,outBoundUpdater,outBExonUpdater);

			for(set<int>::iterator j=os->jnxRowIDs.begin();j!=os->jnxRowIDs.end();j++)
			{
				JnxRow& cjr=jnxrows[*j];
				cerr<<(*j)<<"\t";
				//cerr<<"copying bounds to jnx row: "<<(*j)<<" as "<<fjr.outBExon<<"@"<<fjr.outBound<<endl;
				cjr.setOutBound(outBoundUpdater,outBExonUpdater);
				//copyPathVector(exonPaths,*cjr.vgp); //copied already
			}
			
			if(ovl.size()>0)
			{
			cerr<<"Left:";
			cerr<<ovl<<endl;
			
			jID++;
			gkj_outJnx(snipfile, jnxinfo, raf, setName, jID, curChr, strand , gUp, gDown, ovl,  exonPaths);

			}
			
			if(ovr.size()>0)
			{
			cerr<<"Right:";
			cerr<<ovr<<endl;
			
			jID++;
			gkj_outJnx(snipfile, jnxinfo, raf, setName, jID, curChr, strand , gUp, gDown, ovr,  exonPaths);

			}
			cerr<<"}"<<endl;
			//form new 
			
			
			
			delete os;
		}
		
		cerr<<"out rem Jnx{"<<endl;
		for(int i=0;i<cumF;i++)
		{
			cerr<<"****remaining Jnx for "<<i<<endl;
			JnxRow& jr=jnxrows[i];
		//	cerr<<"a"<<endl;
						
			OverlapState *os =new OverlapState(jr.keylet->leftBlock.k1 /*0*/, 10000 , jr.keylet->nRightBlocks-1, jr.keylet->rightBlocks[jr.keylet->nRightBlocks-1].k2-1 /*0*/ , 1000 ,false);
		//	cerr<<"b"<<endl;
						
			cerr<<"search from left"<<endl;
			list<KeyPair<int,int> > ovl=jr.searchOutBoundFromLeft(os,readLength,outBoundUpdater,outBExonUpdater); 
			jr.setOutBound(outBoundUpdater,outBExonUpdater);
		//	cerr<<"c"<<endl;
			cerr<<"search from right"<<endl;
			list<KeyPair<int,int> > ovr=jr.searchOutBoundFromRight(os,readLength,outBoundUpdater,outBExonUpdater);
			jr.setOutBound(outBoundUpdater,outBExonUpdater);
			//	cerr<<"d"<<endl;
			
			vector<gjk_path*> exonPaths;
			copyPathVector(exonPaths,*jr.vgp);
		//	cerr<<"e"<<endl;
			if(ovl.size()>0)
			{
			cerr<<"Left:";
			cerr<<ovl<<endl;
			
			jID++;
			gkj_outJnx(snipfile, jnxinfo, raf, setName, jID, curChr, strand , gUp, gDown, ovl,  exonPaths);

			}
			
			if(ovr.size()>0)
			{
			cerr<<"Right:";
			cerr<<ovr<<endl;
			
			jID++;
			gkj_outJnx(snipfile, jnxinfo, raf, setName, jID, curChr, strand , gUp, gDown, ovr,  exonPaths);

			}
			

			
			delete os;
		}
		
	//	cerr<<"}"<<endl;
		
	}else /* cumF==1 */
	{
	//	cerr<<"gib1"<<endl;
		//cerr<<"cumF = 1"<<endl;
		
	//	//map<Keylet2,vector<gjk_path*>* >::iterator &sit
		const Keylet2&keylet=(*sit).first;
	//	cerr<<"a"<<endl;
		
		list<KeyPair<int,int> > bounds;
	//	cerr<<"a1"<<endl;
		
		//cerr<<keylet.leftBlock.k1<<endl;
		
		bounds.push_back(keylet.leftBlock);
	//	cerr<<"a2"<<endl;
	//	cerr<<keylet.nRightBlocks<<endl;
		for(int ci=0;ci<keylet.nRightBlocks;ci++)
			bounds.push_back(keylet.rightBlocks[ci]);
		//cerr<<"b"<<endl;
	//	cerr<<"gib2"<<endl;		
		vector<gjk_path*> &exonPaths=*((*sit).second);
		
	//	cerr<<"exonPaths.size()="<<exonPaths.size()<<endl;
		
	/*	for(vector<gjk_path*>::iterator abc=exonPaths.begin();abc!=exonPaths.end();abc++)
		{
			//cerr<<(*abc)<<endl;
		}*/
		//cerr<<"c"<<endl;
		
		jID++;
		
		gkj_outJnx(snipfile, jnxinfo, raf, setName, jID, curChr, strand , gUp, gDown, bounds,  exonPaths);
		//cerr<<"d"<<endl;
				
	}
	
	if(jnxrows)
		delete[] jnxrows;
}

/*void clipRegion0(KeyPair<int,int> clipper, vector< KeyPair<int,int> > clippee, vector< KeyPair<int,int> >&orRegion, vector<KeyPair<int,int> >&newRegion, bool clear=true)
{
	if(clear)
	{
		orRegion.clear();
		newRegion.clear();
		
	}
	
	bool on=false;
	bool prevOn=false;
	
	int prevEnd=-10000;
	
	int orStart;
	
	for(vector<KeyPair<int,int> >::iterator i=clippee.begin();i!=clippee.end();i++)
	{
		
	
		
		if(within000(clipper.k1,*i))
		{
			
			if(!on)
			{//starting
				orStart=(*i).k1;
			}
			
			on=true;
			
		}
		
		if(!on)
		{
			orRegion.push_back(*i);
		}
		
		if(within000(clipper.k2,*i))
		{
			if(on)
			{
				orRegion.push_back(KeyPair<int,int>(orStart,(*i).k2));
			}
			
			on=false;
			
		}
		
		if(prevOn)
		{
			if(prevEnd>=0)
			{
				newRegion.push_back(KeyPair<int,int>(prevEnd+1,(*i).k1-1));
			}
		}	
		
		prevEnd=(*i).k2;
		prevOn=on;
	}
	
	
	
}
*/



class getGBlockState
{
public:
	typedef map<int,GffEntry::ExonPtr>::iterator SIterator;
	typedef pair<SIterator,SIterator> DIterator;
	multimap<int,GffEntry::ExonPtr> end0ExonMap;
	string curChr;
	string setName;
	int last0;
	int bID;
	inline getGBlockState(string _setName):setName(_setName),last0(0),bID(0)
	{
		
	}
	
	class EarlyExons
	{
	public:
		DIterator dIterator;
		int end0;
		
		inline EarlyExons(int _end0,const getGBlockState::DIterator& di)
			:end0(_end0),dIterator(di)
		{
			
		}
		
		inline getGBlockState::SIterator begin()
		{
			return dIterator.first;
		}
		
		inline getGBlockState::SIterator end()
		{
			return dIterator.second;
		}
		
		inline bool isValid()
		{
			return end0>=0;
		}
	};

	inline getGBlockState::SIterator begin()
	{
		return end0ExonMap.begin();
	}
	
	inline getGBlockState::SIterator end()
	{
		return end0ExonMap.end();
	}
	
	
	
	inline void removeExons( EarlyExons &eexon)
	{
		
		end0ExonMap.erase(eexon.begin(),eexon.end());
	}
	
	inline void addExon(GffEntry::ExonPtr exon)
	{
		end0ExonMap.insert(multimap<int,GffEntry::ExonPtr>::value_type(exon->getEnd0(),exon));
	}
	
	inline void reset()
	{
		last0=0;
		end0ExonMap.clear();
	}
	
	inline EarlyExons getEarlyExon()
	{
		if(end0ExonMap.size()<1)
		{
			return EarlyExons(-1,DIterator(end0ExonMap.end(),end0ExonMap.end()));
		}
		SIterator be=end0ExonMap.begin();
		int earlyExonEnd0=(*be).first;
		return EarlyExons(earlyExonEnd0,end0ExonMap.equal_range(earlyExonEnd0));
	}
};





inline void getGBlock_inner_addBlock(ofstream& snip,set<GffEntry::GBlockPtr> *sset,getGBlockState& state,int end1=INT_MAX)
{
	
	if(state.last0>=end1)
	{
		cerr<<"wtf: chr="<<state.curChr<<",last0>=end1:"<<state.last0<<">="<<(end1-1)<<endl;
		return;
	}
	
	
	
	GffEntry::GBlockPtr block=new GffEntry::GBlock(state.curChr,state.last0,end1);
	state.bID++;

	block->strand=GffEntry::UNKNOWN;
	
	int nExons=0;
	
	
	string exonsIDStrings="";
	
	for(getGBlockState::SIterator i=state.begin();i!=state.end();i++)
	{
		GffEntry::ExonPtr ep=(*i).second;
		if(block->strand==GffEntry::UNKNOWN)
			block->strand=ep->strand;
		else if(block->strand!=ep->strand)
			block->strand=GffEntry::BOTH_STRANDS;
		
		if(!ep->blocks)
		{
			ep->blocks=new vector<GffEntry::GBlockPtr>;
		}
		
		ep->blocks->push_back(block);
		
		block->exons.insert(ep);
		nExons++;
		exonsIDStrings+=StringUtil::str(ep->exonID)+"\t";
	}
	
	
	snip<<state.setName<<":"<<state.bID;
	snip<<"\t";
	snip<<block->chr<<"\t";
	snip<<block->start0<<"\t"<<block->end1<<"\t";
	snip<<block->strand<<"\t";
	snip<<nExons<<"\t";
	snip<<exonsIDStrings;
	snip<<";"<<endl;
	
	state.last0=end1;
	cerr<<"last0="<<state.last0<<endl;
	sset->insert(block);
	
	
}

inline void treatCurExons(ofstream &snip,set<GffEntry::GBlockPtr> *sset,getGBlockState& state,int start0)
{
	getGBlockState::EarlyExons eexon=state.getEarlyExon();
	while(eexon.isValid() && eexon.end0<start0) //something in the current set of exons has a smaller end0 then the new exon's start
	{
		//end the cur Block,write stuff
		cerr<<"ENDADD"<<endl;
		getGBlock_inner_addBlock(snip,sset,state,eexon.end0+1);
		//cerr<<"1b"<<endl;
		//remove the earlyExons
		state.removeExons(eexon);
		//prepare for next round
		eexon=state.getEarlyExon();
		
	}
}

void getGBlock(string snipfile, string gfffile, string setName,
		bool gffSrcFile, bool resetEverything) {

	
	ofstream snip(snipfile.c_str()); //output file stream for snip file



	if (gfffile != "") {
		GffEntry::resetEverything();
		GffEntry::Loader loader(NULL);

		if (gffSrcFile) {
			cerr << "try to load " << endl;
			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, "", true);
		}
	} else {
		cerr << "using loaded Gffs" << endl;
	}

	getGBlockState state(setName);

	int lc=0;

	set<GffEntry::GBlockPtr> * chrBlock=NULL;
	//walk through exonome
	
	int prevStart0=-1;
	
	int count=GffEntry::exonome.size();
	cerr<<"Exonome Size for ="<<count<<endl;

	//for each exon in the exonome
	
	set<GffEntry::ExonPtr> curExons;
	
	GffEntry::GBlockPtr curBlock=NULL;

	for(GffEntry::ExonomeWalker ei=GffEntry::exonome.begin();ei!=GffEntry::exonome.end();ei++)
	{
		GffEntry::ExonPtr curExon=pExonOf(ei);
		int start0=curExon->getStart0();
		
		lc++;
		
		if(lc%1000==1)
		{
			cerr<<"loading line "<<lc<<endl;
		}
		
		
		//if(curExon->getEnd0()-start0<1 )
		//	cerr<<"wtf2: end0 near start0 chr="<<state.curChr<<"start0="<<start0<<",end0="<<curExon->getEnd0()<<endl;
		
		if(state.curChr!=curExon->chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{
			//end the curBlock, write stuff
			treatCurExons(snip,chrBlock,state,INT_MAX);
			
			//reset state
			state.reset();
			state.curChr=curExon->chr;			
			
			map<string, set<GffEntry::GBlockPtr> * >::iterator i=GffEntry::gblocks.find(state.curChr);
			if(i==GffEntry::gblocks.end())
				
			{
				chrBlock=new set<GffEntry::GBlockPtr>;
				GffEntry::gblocks.insert(map<string, set<GffEntry::GBlockPtr>* >::value_type(state.curChr,chrBlock));
			}
			else
				chrBlock=(*i).second;
			
			
			prevStart0=-1;
			
		}
		else if(start0==prevStart0)
		{
			//same Start, just add
			state.addExon(curExon);
			continue;
		}
		
		//end cur Block, write stuff
		treatCurExons(snip,chrBlock,state,start0);
		
		cerr<<"STARTADD"<<endl;
		getGBlock_inner_addBlock(snip,chrBlock,state,start0);
		//add new Exon		
		state.addExon(curExon);
		prevStart0=start0;

		
	}
	
	//end the curBlock, write stuff
	treatCurExons(snip,chrBlock,state,INT_MAX);
	//reset state no need
	
	
	cerr<<"total "<<lc<<" lines processed"<<endl;



	if(resetEverything)
		GffEntry::resetEverything(); //free memory

	snip.close();
}

void getKnownJnxome(int readLength, int minExonSpan, string snipfile,
		string jnxinfofile, string chrDir, string gfffile, string setName,
		bool gffSrcFile, bool resetEverything) {

	int flank=readLength-minExonSpan;  //get the flanking length
	
	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream jnxinfo(jnxinfofile.c_str()); //output file stream for jnx info file

	int bufsize=flank*2+1; //buffer for the transfer of sequence. +1 to include \0

	if (gfffile != "") { //if gfffile is given, load it
		GffEntry::resetEverything();
		GffEntry::Loader loader(NULL);

		if (gffSrcFile) {
			cerr << "try to load " << endl;
			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, "", true);
		}
	} else { //else use loaded Gff
		cerr << "using loaded Gffs" << endl;
	}

	char *transferBuffer=new char[bufsize]; //init transfer buffer

	map<string,int> names;
	int jnxID=0;

	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file

	int tjnx=0;
	int cjnx=0;

	int ambig=0;
	int ambigGroup=0;
	int maxAmbig=0;
	
	int lc=0;
	int jID=0;
	char strands[]= { GffEntry::FORWARD, GffEntry::REVERSE };

	char curStrand;

	for (int stri=0; stri<2; stri++) {
		curStrand=strands[stri];
		map<string, map<Keylet2, vector< gjk_path*>*>*> alreadyin; //to remember which regions are already included (because multiple transcripts (of same gene) can have same promoter
		map<Keylet2,vector<gjk_path*> *> *collisionCheck=NULL;

		for (vector<GffEntry*>::iterator gfe=GffEntry::transcriptome.begin(); gfe
				!=GffEntry::transcriptome.end(); gfe++) //read all entries from gff file
		{
			GffEntry& entry=**gfe;

			if (entry.strand!=curStrand)
				continue;

			lc++;
			if (lc%1000==1)
				cerr<<"processing entry "<<lc<<endl;

			map<string,int>::iterator i;

			if (curChr!=entry.chrom) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
			{

				curChr=entry.chrom; //and then open a new seq
				map<string, map<Keylet2, vector<gjk_path*>*>*>::iterator alr;

				if ((alr=alreadyin.find(curChr))==alreadyin.end()) {
					//no record of curChr
					collisionCheck =new map<Keylet2, vector <gjk_path*>*>;
					alreadyin.insert(map<string,
							map< Keylet2, vector<gjk_path*>*>*>::value_type(
							curChr, collisionCheck));

				} else {
					collisionCheck=(*alr).second;
				}

			}
			
//			cerr<<"b"<<endl;
			
			 for (int i=0; i<entry.exonCount-1; i++) 
			 {
//				cerr<<"c0"<<endl;
				 
				
				
				 
				 
				int gUp=entry.exons[i]->end; //1 index incl
				int leftLength=entry.exons[i]->getLength();
				int gDown=entry.exons[i+1]->start; //0 index incl
				int rightLength=entry.exons[i+1]->getLength(); //length of the right exon
				int lFlank=flank;
				int rFlank=flank;
				
				
				//cerr<<"Jnxmaking Exons "<<entry.exons[i]->exonID<<"\t"<<entry.exons[i]->exonID<<endl;
				
//				cerr<<"c1"<<endl;

				if (leftLength<minExonSpan || rightLength<minExonSpan) {
					cerr<<"leftLength="<<leftLength<<",rigthLength="
							<<rightLength<<":either length < minExonSpan="
							<<minExonSpan<<",thus ignored"<<endl;
					continue;
				}
//				cerr<<"c2"<<endl;
				
				int curI=i;
				KeyPair<int,int> leftBlock(0, 0);
				vector<KeyPair<int,int> > rightBlocks;
//				cerr<<"c3"<<endl;
				gjk_path *path=new gjk_path(curStrand);
//				cerr<<"c4"<<endl;
				int lenThisTime=MIN(lFlank,leftLength);
				leftBlock.k1=gUp-lenThisTime;
				leftBlock.k2=gUp;

				//cerr<<"leftBock="<<leftBlock<<endl;
//				cerr<<"c5"<<endl;
				path->leftExon=entry.exons[curI];

				//only go right

				int rightUsed=0;
				//cerr<<"leftBock=";
				curI=i+1; //start off from the exon to the right of the jnx
//				cerr<<"c6"<<endl;
				while (rFlank>0) {
//					cerr<<"c7.1"<<endl;	
					if (curI<entry.exonCount-1)
						lenThisTime=MIN(rFlank,entry.exons[curI]->getLength());
					else
						lenThisTime=rFlank; //protrude into outside of TxEnd!!
					
//					cerr<<"c7.2"<<endl;	
					path->rightExons.push_back(entry.exons[curI]);
//					cerr<<"c7.3"<<endl;	
					KeyPair<int, int> rightBlk(
							entry.exons[curI]->start, //0-based end in exon; no need to +1
							entry.exons[curI]->start+lenThisTime); //0-based to 1-based
					
					rightBlocks.push_back(rightBlk); 

				//	cerr<<rightBlk<<"\t";
//					cerr<<"c7.4"<<endl;	
					rFlank-=lenThisTime;
					curI++;
					rightUsed++;

				}
				
				//cerr<<endl;
				
//				cerr<<"c8"<<endl;	

				Keylet2 limits(gUp, gDown, leftBlock, rightBlocks);

				map< Keylet2, vector<gjk_path*>*>::iterator j;
				vector<gjk_path*> *vgep;
//				cerr<<"c9"<<endl;
				if ((j=collisionCheck->find(limits))==collisionCheck->end()) //if the region hasn't been outputed yet, do
				{
				//	cerr<<"region not out yet"<<endl;
					vgep=new vector<gjk_path*>;

					limits.jnxID=(jnxID++);
//					cerr<<"preinsert"<<endl;
					collisionCheck->insert(map<Keylet2, vector<gjk_path*>*>::value_type(
							limits, vgep));
//					cerr<<"postinsert"<<endl;
					
				} else {
				//	cerr<<"region in, append exon"<<endl;
					cjnx++;
					vgep=(*j).second;
				}
//				cerr<<"c10"<<endl;
				vgep->push_back(path);
//				cerr<<"c11"<<endl;
			}

		}
		
		cerr<<"finish input"<<endl;

		
		
		

		
		for(map<string, map<Keylet2, vector< gjk_path*>*>*>::iterator chrRec=alreadyin.begin();
			chrRec!=alreadyin.end();
			chrRec++)
		{
			
			string curChr=(*chrRec).first;
			map<Keylet2,vector<gjk_path*>* >* curRecMap=(*chrRec).second;
			

			if(raf) //remember to close and delete the previous RandomAccessFile Object
			{
				raf->close();
				delete raf;
			}

			raf=new RandomAccessFile(chrDir+curChr+".seq");
			
			KeyPair<int,int> lastPrimJnx(-1,-1);
			int cumF=1;			
			
			map<Keylet2,vector<gjk_path*>* >::iterator sit=curRecMap->end();
			
			for(map<Keylet2,vector<gjk_path*>* >::iterator it=curRecMap->begin();
				it!=curRecMap->end();
				it++)
			{
				
				
				
				tjnx++;
				
				
				const Keylet2& keylet=(*it).first;
				//int gUp=keylet.gUp;
				//int gDown=keylet.gDown;
				
				
				
				if(lastPrimJnx.k1==keylet.gUp && lastPrimJnx.k2==keylet.gDown)
				{
					cumF++;
					
				}
				else
				{
					if(sit!=curRecMap->end())	
						gkj_inner_outJnx(snip,jnxinfo, raf, setName,jID,readLength,minExonSpan,curChr , curStrand, lastPrimJnx.k1, lastPrimJnx.k2,cumF,lastPrimJnx,sit,ambig,ambigGroup,maxAmbig);
					
					cumF=1;
					sit=it;
					lastPrimJnx.k1=keylet.gUp;
					lastPrimJnx.k2=keylet.gDown;
				}
				
				
				
			
			}
			
			if(sit!=curRecMap->end())
				gkj_inner_outJnx(snip,jnxinfo,raf, setName,jID,readLength,minExonSpan,curChr, curStrand, lastPrimJnx.k1, lastPrimJnx.k2, cumF,lastPrimJnx,sit,ambig,ambigGroup,maxAmbig);
								
			
			//delete curRecMap;
		}
		
		
		cerr<<"finish out, clearing up memory"<<endl;
		
		
		
		
		//clean up memory
		for(map<string, map<Keylet2, vector< gjk_path*>*>*>::iterator chrRec=alreadyin.begin();
			chrRec!=alreadyin.end();
			chrRec++)
		{
			map<Keylet2,vector<gjk_path*>* >* curRecMap=(*chrRec).second;
	
			
			for(map<Keylet2,vector<gjk_path*>* >::iterator it=curRecMap->begin();
				it!=curRecMap->end();
				it++)
			{

				
				
				vector<gjk_path*>* vgp=(*it).second;
				for(vector<gjk_path*>::iterator vi=vgp->begin();
					vi!=vgp->end();
					vi++)
				{
					gjk_path* pPath=(*vi);
					delete pPath;
				}
				
				delete vgp;
			
			}
			
			
			delete curRecMap;
		}
		
		
		cerr<<"finish clearing memory"<<endl;
		
		
		
		if (raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
		{
			raf->close();
			delete raf;
			raf=NULL;
		}
		
	}

	cerr<<"total "<<lc<<" lines processed"<<endl;



	cerr<<"Total jnx: "<<tjnx<<endl;
	cerr<<"collision jnx:" <<cjnx<<endl;
	cerr<<"Ambig="<<ambig<<",AmbigGroup="<<ambigGroup<<",AmbigMax="<<maxAmbig<<",AmbigAve="<<((double)ambig/ambigGroup)<<endl;
	

	delete[] transferBuffer;
	GffEntry::resetEverything(); //free memory

	jnxinfo.close();
	snip.close();
}

void _getKnownJnxome(int flank, string snipfile, string jnxinfofile , string chrDir,  string gfffile,  string setName,bool gffSrcFile,bool resetEverything) {

	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream jnxinfo(jnxinfofile.c_str()); //output file stream for jnx info file

	int bufsize=flank*2+1;

	if (gfffile != "") {
		GffEntry::resetEverything();
		GffEntry::Loader loader(NULL);

		if (gffSrcFile) {
			cerr << "try to load " << endl;
			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, "", true);
		}
	} else {
		cerr << "using loaded Gffs" << endl;
	}

	
	char *transferBuffer=new char[bufsize];

		map<string,int> names;
	int JnxID=0;

//	int boundOutCount=0;
	

	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file

	int tjnx=0;
	//GffEntry::resetEverything(); //reset the exonome

	//int tjnx=0;

	int lc=0;
	
	char strands[]={GffEntry::FORWARD, GffEntry::REVERSE};
	
	char curStrand;
	
	for (int stri=0; stri<2; stri++) {
		curStrand=strands[stri];

		map<string, map<Keylet, vector<KeyPair<rj_state*,rj_state*> >* >* > alreadyin; //to remember which regions are already included (because multiple transcripts (of same gene) can have same promoter


		map<Keylet,vector<KeyPair<rj_state*,rj_state*> >*  > *collisionCheck=NULL;
		
		for (vector<GffEntry*>::iterator gfe=GffEntry::transcriptome.begin(); gfe
				!=GffEntry::transcriptome.end(); gfe++) //read all entries from gff file
		{
			GffEntry& entry=**gfe;
			
			if(entry.strand!=curStrand)
				continue;
			
			lc++;
			if(lc%1000==1)
				cerr<<"processing entry "<<lc<<endl;

			map<string,int>::iterator i;
			map< Keylet, vector<KeyPair<rj_state*,rj_state*> >*>::iterator j;

			if (curChr!=entry.chrom) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
			{
				if (raf) //remember to close and delete the previous RandomAccessFile Object
				{
					raf->close();
					delete raf;
				}

				curChr=entry.chrom; //and then open a new seq
				raf=new RandomAccessFile(chrDir+curChr+".seq",bufsize,transferBuffer);
				map<string, map<Keylet, vector<KeyPair<rj_state*,rj_state*> >*>*>::iterator
						alr;

				if ((alr=alreadyin.find(curChr))==alreadyin.end()) {
					//no record of curChr
					collisionCheck
							=new map<Keylet, vector <KeyPair<rj_state*,rj_state*> >*>;
					alreadyin.insert(map<string,
							map< Keylet, vector<KeyPair<rj_state*,rj_state*> >*>* >::value_type(
							curChr, collisionCheck));

				} else {
					collisionCheck=(*alr).second;
				}

			}

			for (int i=0; i<entry.exonCount-1; i++) {
				int gUp=entry.exons[i]->end; //1 index incl
				int gDown=entry.exons[i+1]->start; //0 index incl
				int lFlank=flank;
				int rFlank=flank;

				string extra;

				string lExtra="";
				string rExtra="";
				//care about upper
				int curI=i;

				KeyPair<int,int> coord(0,0);
				rj_state* leftState=new rj_state(curStrand,gUp,0,0);
				rj_state* rightState=new rj_state(curStrand,gDown,0,0);
				
				int leftUsed=0;
				int rightUsed=0;
				
				while (lFlank>0) {
					int lenThisTime;

					if (curI>0)
						lenThisTime=MIN(lFlank,entry.exons[curI]->getLength());
					else
						lenThisTime=lFlank; //protrude into outside of TxStart!
						
					
					coord.k1=entry.exons[curI]->end-lenThisTime; //1-based end in exon => 0 based coord; no need to +1
					coord.k2=entry.exons[curI]->end; //1-based
					
					leftState->blocks.push_back(coord);
					leftState->path.push_back(entry.exons[curI]);
					
					lFlank-=lenThisTime;
					curI--;
					leftUsed++;
					leftState->extra=StringUtil::str(coord.k1)+"/"
							+StringUtil::str(coord.k2)+((leftUsed==1) ? ""
							: "/")+leftState->extra;
				}

				curI=i+1; //start off from the exon to the right of the jnx

				while (rFlank>0) {
					int lenThisTime;
					
					if (curI<entry.exonCount-1)
						lenThisTime=MIN(rFlank,entry.exons[curI]->getLength());
					else
						lenThisTime=rFlank; //protrude into outside of TxEnd!!

					coord.k1=entry.exons[curI]->start; //0-based end in exon; no need to +1
					coord.k2=entry.exons[curI]->start+lenThisTime; //0-based to 1-based
					
					rightState->blocks.push_back(coord);
					rightState->path.push_back(entry.exons[curI]);
					
					rFlank-=lenThisTime;
					curI++;
					rightUsed++;
					rightState->extra+=((rightUsed==1) ? "" : "/")
							+StringUtil::str(coord.k1)+"/"
							+StringUtil::str(coord.k2);

				}

				extra=leftState->extra+":"+rightState->extra;

				//strand-incensitive
				Keylet limits(gUp, gDown, extra);
				vector<KeyPair<rj_state*,rj_state*> > *vgep;

				if ((j=collisionCheck->find(limits))==collisionCheck->end()) //if the region hasn't been outputed yet, do
				{
					
					tjnx++;
					//int id=collisionCheck->size();

					vgep=new vector<KeyPair<rj_state*,rj_state*> >;

					limits.JnxID=(JnxID++);
					collisionCheck->insert(map<Keylet,
							vector<KeyPair<rj_state*,rj_state*> >* >::value_type(limits,
							vgep));

					snip << ">" << curChr << "|" << setName << ":"
							<< limits.JnxID << endl;

					for (int k=leftUsed-1; k>=0; k--)
						raf->transfer(snip, leftState->blocks[k].k1, leftState->blocks[k].k2);
					
					//snip<<endl;
					
					for (int k=0; k<rightUsed; k++)
						raf->transfer(snip, rightState->blocks[k].k1, rightState->blocks[k].k2);

					snip<<endl; //write the sequence to snip file

				} else {

					vgep=(*j).second;
				}

				vgep->push_back(KeyPair<rj_state*,rj_state*>(leftState,rightState));

			}

		}

		//here the jnxinfo>>>>
		//clear the collision thing.
		for (map<string, map<Keylet,vector<KeyPair<rj_state*,rj_state*> >* >*>::iterator i=alreadyin.begin();
				i!=alreadyin.end(); i++)
			{
				int numJnx=(*i).second->size();


				cerr<<numJnx<<"  "<<curStrand<<"ve strand jnxs in "<<(*i).first<<endl;
				//tjnx+=numJnx;

				
				map<Keylet, vector<KeyPair<rj_state*,rj_state*> >* > *pJnxEntry=(*i).second;

				string curChr=(*i).first;

				for(map<Keylet, vector<KeyPair<rj_state*,rj_state*> >* >::iterator j=pJnxEntry->begin();
					j!=pJnxEntry->end();
					j++)
				{

					//if(StringUtil::split((*j).first.extra,"/",tmp).size()>3)
					//	boundOutCount++;
					
					vector<KeyPair<rj_state*,rj_state*> >* vkp=(*j).second;
					
					
					
					jnxinfo<<setName<<":"<<(*j).first.JnxID<<"\t"<<curChr<<"\t"<<(*j).first.extra<<"\t";
					jnxinfo<<curStrand<<"\t";
					jnxinfo<<(*vkp)[0].k1->blocks.size()<<"\t";
					jnxinfo<<(*vkp)[0].k2->blocks.size()<<"\t";
					jnxinfo<<vkp->size()<<"\t";
					
					for(vector<KeyPair<rj_state*,rj_state*> >::iterator it=vkp->begin();
						it!=vkp->end();
						it++)
					{
						KeyPair<rj_state*,rj_state*> & kps=(*it);
						jnxinfo<<kps.k1->path[0]->exonID<<"\t";
						jnxinfo<<kps.k2->path[0]->exonID<<"\t";
						
						jnxinfo<<kps.k1->path.size()<<"\t";
						
						for(vector<GffEntry::Exon*>::iterator iter2=kps.k1->path.begin();
							iter2!=kps.k1->path.end();
							iter2++)
						{
							jnxinfo<<(*iter2)->exonID<<"\t";
						}

						jnxinfo<<kps.k2->path.size()<<"\t";
						
						for(vector<GffEntry::Exon*>::iterator iter2=kps.k2->path.begin();
							iter2!=kps.k2->path.end();
							iter2++)
						{
							jnxinfo<<(*iter2)->exonID<<"\t";
						}
						
						delete kps.k1;
						delete kps.k2;
						
						
					}
					
					jnxinfo<<"ok";
					jnxinfo<<endl;
				

					delete vkp;
				}

				delete pJnxEntry;
				
			}

	}

	cerr<<"total "<<lc<<" lines processed"<<endl;

	if (raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
	{
		raf->close();
		delete raf;
	}

	//int tjnx=0;
	//vector<string> tmp;

/*	for (map<string, map<Keylet,vector<RankGffEntryExonPair>* >*>::iterator i=alreadyin.begin();
		i!=alreadyin.end(); i++)
	{
		int numJnx=(*i).second->size();


		cerr<<numJnx<<" jnxs in "<<(*i).first<<endl;
		tjnx+=numJnx;


		map<Keylet, vector<RankGffEntryExonPair>* > *pJnxEntry=(*i).second;

		string curChr=(*i).first;




		for(map<Keylet, vector<RankGffEntryExonPair>* >::iterator j=pJnxEntry->begin();
			j!=pJnxEntry->end();
			j++)
		{

			if(StringUtil::split((*j).first.extra,"/",tmp).size()>3)
				boundOutCount++;

			jnxinfo<<setName<<":"<<(*j).first.JnxID<<"\t"<<curChr<<"\t"<<(*j).first.extra<<"\t";


			vector<RankGffEntryExonPair>* pGffEntries=(*j).second;
			jnxinfo<<pGffEntries->size()<<"\t";


			int strand=-1;
			bool ambig=false;

			for(vector<RankGffEntryExonPair>::iterator k=pGffEntries->begin();
				k!=pGffEntries->end();
				k++)
			{
				RankGffEntryExonPair& transcriptinfo=(*k);
				GffEntry& transcript=*transcriptinfo.gff;
				jnxinfo<<transcript.getExonRankFromIndex(transcriptinfo.rank)<<"\t"<<transcriptinfo.pExon->exonID<<"\t"<<transcriptinfo.pExon2->exonID<<"\t"<<transcriptinfo.gff->transcriptID<<"\t"<<transcript.name<<"\t"<<transcript.name2<<"\t"<<transcript.id<<"\t"<<transcript.strand<<"\t";

				if(strand==-1)
					strand=transcript.strand;

				if(strand!=transcript.strand)
					ambig=true;

			}

			jnxinfo<<((ambig)?"amb":"ok");
			jnxinfo<<endl;

			delete pGffEntries;
		}

		delete pJnxEntry;
	}


	delete[] leftCoord;
	delete[] rightCoord;*/

	cerr<<"Total jnx: "<<tjnx<<endl;
	//cerr<<"BoundOut Count= "<<boundOutCount<<endl;

	delete[] transferBuffer;
	GffEntry::resetEverything(); //free memory

	jnxinfo.close();
	snip.close();
}



#define DEBUG_ON

#ifdef DEBUG_ON
#define __DEBUG__ cerr<<"@"<<__LINE__<<endl
#else
#define __DEBUG__
#endif

//loci are strand sensitive, so no problem.

#define SNIP_PATH_MAX 10000


void getRecombineJnxome(int flank, string snipfile, string jnxinfofile,
		string chrDir, int dist, string gfffile, string locusfile,
		string setName, bool gffSrcFile, bool resetEverything) {

	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream jnxinfo(jnxinfofile.c_str()); //output file stream for jnx info file

	if (gfffile != "") {
		GffEntry::resetEverything();
		GffEntry::Loader loader(NULL);

		if (gffSrcFile) {
			cerr << "try to load " << endl;
			loader.loadGffFiles(gfffile);
		} else {
			loader.loadGffFile(gfffile, "", true);
		}
	} else {
		cerr << "using loaded Gffs" << endl;
	}

	if (locusfile != "") {
		GffEntry::Loader loader(NULL);
		loader.loadLoci(locusfile);
	} else {
		cerr << "using loaded loci" << endl;
	}

	int bufsize = flank*2 + 1;

	char *transferBuffer = new char[bufsize];



	map<string, int> names;
	int JnxID = 0;


	RandomAccessFile *raf = NULL; //the random access file to get the genomic sequence
	string curChr = ""; //the id of the current chr, used to decide to whether to open a new genomic seq file
	//vector<rj_state*> memPaths;
	
	map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::iterator j;
	map<GffEntry::ExonPtr, GffEntry::Locus*>::iterator ei;

	int locSize=GffEntry::loci.size();
	int procLoc=0;

	for (j=GffEntry::loci.begin();j!=GffEntry::loci.end();j++)
	{
		
		multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus=(*j).second;
		
	for (ei = chrMapLocus->begin(); ei != chrMapLocus->end(); ei++) {
		GffEntry::Locus& entry = *((*ei).second);

		procLoc++;
		if(procLoc%1000==1)
			cerr<<"processing locus "<<procLoc<<" of "<<locSize<<endl;

		map<string, int>::iterator i;

		map<Keylet, vector<KeyPair<rj_state*, rj_state*> >*> collisionCheck;
		map<Keylet, vector<KeyPair<rj_state*, rj_state*> >*>::iterator j;

		if (curChr != entry.chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{
			if (raf) //remember to close and delete the previous RandomAccessFile Object
			{
				raf->close();
				delete raf;
			}

			curChr = entry.chr; //and then open a new seq
			raf= new RandomAccessFile(chrDir + curChr + ".seq", bufsize, transferBuffer);
		}

		//construct all paths;

		vector<rj_state*> leftPaths;
		vector<rj_state*> rightPaths;

		queue<rj_state> Q;
		__DEBUG__;
		if(entry.exonCount>1000)
			continue;

		cerr<<procLoc<<",ec="<<entry.exonCount<<endl;
		//cal leftPaths;

		for (int x = 0; x < entry.exonCount; x++) {
			cerr<<"x="<<x<<",rj="<<rj_state::count<<endl;
			
			if(rj_state::count>SNIP_PATH_MAX)
			{
				cerr<<"abort constructing new rj_state; rj_state>SNIP_PATH_MAX:"<<rj_state::count<<">"<<SNIP_PATH_MAX<<endl;
					break;
			}	
			Q.push(rj_state(entry.strand, entry.exons[x]->end, x, flank));

			//go left;
			while (!Q.empty()) {
				rj_state& state = Q.front();
				int u = state.curExon;
				GffEntry::ExonPtr curExonPtr = entry.exons[u];

				int lenThisTime = MIN(state.remflank, curExonPtr->getLength());
				state.path.push_back(curExonPtr);

				KeyPair<int, int> bound(curExonPtr->end - lenThisTime,
						curExonPtr->end);

				state.blocks.push_back(bound);
				state.remflank -= lenThisTime;
				state.extra = StringUtil::str(bound.k1) + "/"
						+ StringUtil::str(bound.k2)
						+ ((state.extra != "") ? "/" : "") + state.extra;

				bool continueDown = false;

				if (state.remflank > 0) {
					for (int v = u - 1; v >= 0; v--) {

						if (entry.exons[v]->end + dist > entry.exons[u]->start)
							continue;

						continueDown = true;
						rj_state newState = state;//a new copy, a branch off

						newState.curExon = v;
						Q.push(newState); //queue the task;
					}
				}

				if (!continueDown) //either no more searchable or enough
				{
					rj_state* memState = new rj_state(state);
					leftPaths.push_back(memState);
				}
				Q.pop();
			}
		}
		__DEBUG__;
		//cal rightPaths;
		//Q.clear();
		
		
		
		for (int x = 0; x < entry.exonCount; x++) {

			Q.push(rj_state(entry.strand, entry.exons[x]->start, x, flank));
			//go right;
			if(rj_state::count>SNIP_PATH_MAX)
			{
				cerr<<"abort constructing new rj_state; rj_state>SNIP_PATH_MAX:"<<rj_state::count<<">"<<SNIP_PATH_MAX<<endl;
					break;
			}
			while (!Q.empty()) {
				rj_state& state = Q.front();
				int u = state.curExon;
				GffEntry::ExonPtr curExonPtr = entry.exons[u];

				int lenThisTime = MIN(state.remflank, curExonPtr->getLength());
				state.path.push_back(curExonPtr);

				KeyPair<int, int> bound(curExonPtr->start, curExonPtr->start
						+ lenThisTime);

				state.blocks.push_back(bound);
				state.remflank -= lenThisTime;
				state.extra += ((state.extra == "") ? string("") : string("/"))
						+ StringUtil::str(bound.k1) + "/" + StringUtil::str(
						bound.k2);

				bool continueDown = false;

				if (state.remflank > 0) {
					for (int v = u + 1; v < entry.exonCount; v++) {

						if (entry.exons[u]->end + dist > entry.exons[v]->start)
							continue;

						continueDown = true;
						rj_state newState = state;//a new copy, a branch off

						newState.curExon = v;

						Q.push(newState); //queue the task;
					}
				}

				if (!continueDown) //either no more searchable or enough
				{
					rj_state* memState = new rj_state(state);
					rightPaths.push_back(memState);
				}

				Q.pop();
			}
		}
		__DEBUG__;
		cerr<<leftPaths.size()<<"<><><>"<<rightPaths.size()<<endl;
		for (vector<rj_state*>::iterator x = leftPaths.begin(); x
				!= leftPaths.end(); x++) {
			for (vector<rj_state*>::iterator y = rightPaths.begin(); y
					!= rightPaths.end(); y++) {

				int gUp = (*x)->gBound;
				int gDown = (*y)->gBound;

				if (gUp + dist > gDown)
					continue;

				string extra = (*x)->extra + ":" + (*y)->extra;

				//strand-incensitive
				Keylet limits(gUp, gDown, extra);

				vector<KeyPair<rj_state*, rj_state*> >* group;

				if ((j = collisionCheck.find(limits)) == collisionCheck.end()) //if the region hasn't been outputed yet, do
				{

//					int id = collisionCheck.size();

					limits.JnxID = (JnxID++);

					if(limits.JnxID%1000==1)
						cerr<<"new jnx "<<limits.JnxID<<"ec="<<entry.exonCount<<"x="<<leftPaths.size()<<"y="<<rightPaths.size()<<endl;

					group = new vector<KeyPair<rj_state*, rj_state*> > ;
					collisionCheck.insert(map<Keylet, vector<KeyPair<rj_state*,
							rj_state*> >*>::value_type(limits, group));

					snip << ">" << curChr << "|" << setName << ":"
							<< limits.JnxID << endl;

					for (vector<KeyPair<int, int> >::reverse_iterator k =
							(*x)->blocks.rbegin(); k != (*x)->blocks.rend(); k++) {
						raf->transfer(snip, (*k).k1, (*k).k2);
						//snip<<"\t"<<(*k).k1<<","<<(*k).k2;
					}
					snip<<endl;
					for (vector<KeyPair<int, int> >::iterator k =
							(*y)->blocks.begin(); k != (*y)->blocks.end(); k++) {
						//snip<<"\t"<<(*k).k1<<","<<(*k).k2;
						raf->transfer(snip, (*k).k1, (*k).k2);
					}

					snip << endl; //write the sequence to snip file

				} else {

					group = (*j).second;
				}

				group->push_back(KeyPair<rj_state*, rj_state*> ((*x), (*y)));

			}
		}
		__DEBUG__;
		//free up unused paths


		for (map<Keylet, vector<KeyPair<rj_state*, rj_state*> >*>::iterator j =
				collisionCheck.begin(); j != collisionCheck.end(); j++) {

			//if (StringUtil::split((*j).first.extra, "/", tmp).size() > 3)
			//	boundOutCount++;

			jnxinfo << setName << ":" << (*j).first.JnxID << "\t" << curChr
					<< "\t" << (*j).first.extra << "\t";
			vector<KeyPair<rj_state*, rj_state*> >* exonPaths = (*j).second;

			jnxinfo << ((*exonPaths)[0]).k1->strand << "\t";
			jnxinfo<< ((*exonPaths)[0]).k1->blocks.size()<<"\t";
			jnxinfo<< ((*exonPaths)[0]).k2->blocks.size()<<"\t";

			jnxinfo << exonPaths->size() << "\t";

			for (vector<KeyPair<rj_state*, rj_state*> >::iterator z =
					exonPaths->begin(); z != exonPaths->end(); z++) {
				KeyPair<rj_state*, rj_state*> &kp = (*z);
				rj_state* leftState = kp.k1;
				rj_state* rightState = kp.k2;
				jnxinfo << leftState->path[0]->exonID << "\t";
				jnxinfo << rightState->path[0]->exonID << "\t";
				int ksize = leftState->path.size();
				jnxinfo << (ksize - 1) << "\t";
				for (int h = 1; h < ksize; h++) {
					jnxinfo << leftState->path[h]->exonID << "\t";
				}
				ksize = rightState->path.size();
				jnxinfo << (ksize - 1) << "\t";

				for (int h = 1; h < ksize; h++) {
					jnxinfo << rightState->path[h]->exonID << "\t";
				}



	//			int strand = -1;
				//bool ambig = false;



			}
			jnxinfo<<"ok";
			jnxinfo << endl;
			delete exonPaths;

		}
		cerr<<"before:"<<rj_state::count<<endl;
		for (vector<rj_state*>::iterator x = leftPaths.begin(); x
				!= leftPaths.end(); x++)
			delete (*x);
		for (vector<rj_state*>::iterator y = rightPaths.begin(); y
				!= rightPaths.end(); y++)
			delete (*y);
		cerr<<"after:"<<rj_state::count<<endl;
		__DEBUG__;
	}
	}

	if (raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
	{
		raf->close();
		delete raf;
	}
	cerr<<"done"<<endl;

	delete[] transferBuffer;
	GffEntry::resetEverything(); //free memory

	jnxinfo.close();
	snip.close();

}

//typedef GffEntry::Exon Exon;

void getdeNovoJnxome(int flank, int dist, string gfffile, string chrDir,string snipfile, string jnxinfofile, string setName, bool requireFuncFusion,int minIntronLength) {

	ifstream gff(gfffile.c_str()); //input file stream for gff file
	ofstream snip(snipfile.c_str()); //output file stream for snip file
	ofstream jnxinfo(jnxinfofile.c_str()); //output file stream for jnx info file

	KeyPair<int,int> *leftCoord=new KeyPair<int,int>[flank];
	KeyPair<int,int> *rightCoord=new KeyPair<int,int>[flank];
	int leftUsed;
	int rightUsed;

	vector<novoInfo> JnxInfoS;
	char line[GFF_BUFFER_SIZE]; //the buffer for reading lines from gff file

	map<string, map<Keylet, vector<RankGffEntryExonPair>* >* > alreadyin; //to remember which regions are already included (because multiple transcripts (of same gene) can have same promoter

	map<string,int> names;
	int JnxID=0;

//	int boundOutCount=0;
	map<Keylet,vector<RankGffEntryExonPair>* > *collisionCheck=NULL;

	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file


	GffEntry::resetEverything(); //reset the exonome

	if (!gff.good()) {

		cerr<<"Cannot open Gff"<<endl;
		return;

	}

	int lc=0;

	while (!gff.eof()) //read all entries from gff file
	{


		gff.getline(line, GFF_BUFFER_SIZE);


		++lc;

		 //if((++lc)%100==0)
		//	cerr<<"processing line "<<(lc)/*<<" "<<line*/<<endl;

		if(!strcmp(line,""))
			break;

		vector<string> spliton; //vector to carry the split tokens from each line of gff file
		StringUtil::split(string(line), string("\t"), spliton);

		if (spliton.size()<1)
			continue;
		if (spliton[0]=="#bin") //ignore commented lines from gff files
			continue;


		GffEntry &entry = *(GffEntry::createGffEntry(spliton)); //construct a GffEntry from the tokens


		map<string,int>::iterator i;
		map< Keylet, vector<RankGffEntryExonPair>* >::iterator j;

//		int seeks, seeks2, seeke, seeke2; //to remember the start and end of the region to get sequence


		if (curChr!=entry.chrom) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
		{

			curChr=entry.chrom; //and then open a new seq

			map<string, map<Keylet, vector<RankGffEntryExonPair>* >* >::iterator alr;

			if ((alr=alreadyin.find(curChr))==alreadyin.end()) {
				//no record of curChr
				collisionCheck=new map<Keylet, vector<RankGffEntryExonPair>* >;
				alreadyin.insert(map<string, map< Keylet, vector<RankGffEntryExonPair>* >* >::value_type(
						curChr, collisionCheck));


			}
			else
			{
				//cerr<<"***** "<<curChr<<" again!!"<<endl;
				collisionCheck=(*alr).second;
			}

		}

		for (int i=0; i<entry.exonCount-1; i++) {
			int gUp=entry.exons[i]->end; //1 index incl
			int gDown=entry.exons[i+1]->start; //0 index incl
			int lFlank=flank;
			int rFlank=flank;

			string extra;

			leftUsed=0;
			rightUsed=0;
			string lExtra="";
			string rExtra="";

			//care about upper
			int curI=i;

			while(lFlank>0)
			{
				int lenThisTime;
				KeyPair<int,int>& coord=leftCoord[leftUsed];
				if(curI>0)
					lenThisTime=MIN(lFlank,entry.exons[curI]->getLength());
				else
					lenThisTime=lFlank; //protrude into outside of TxStart!
				coord.k1=entry.exons[curI]->end-lenThisTime; //1-based end in exon => 0 based coord; no need to +1
				coord.k2=entry.exons[curI]->end; //1-based
				lFlank-=lenThisTime;
				curI--;
				leftUsed++;
				lExtra=StringUtil::str(coord.k1)+"/"+StringUtil::str(coord.k2)+((leftUsed==1)?"":"/")+lExtra;
			}

			curI=i+1; //start off from the exon to the right of the jnx

			while(rFlank>0)
			{
				int lenThisTime;
				KeyPair<int,int>& coord=rightCoord[rightUsed];
				if(curI<entry.exonCount-1)
					lenThisTime=MIN(rFlank,entry.exons[curI]->getLength());
				else
					lenThisTime=rFlank; //protrude into outside of TxEnd!!

				coord.k1=entry.exons[curI]->start; //0-based end in exon; no need to +1
				coord.k2=entry.exons[curI]->start+lenThisTime; //0-based to 1-based
				rFlank-=lenThisTime;
				curI++;
				rightUsed++;
				rExtra+=((rightUsed==1)?"":"/")+StringUtil::str(coord.k1)+"/"+StringUtil::str(coord.k2);

			}


			extra=lExtra+":"+rExtra;


			//strand-incensitive
			Keylet limits(gUp, gDown,extra);
			vector<RankGffEntryExonPair> *vgep;

			map< Keylet, vector<RankGffEntryExonPair>* >::iterator j;
			if ((j=collisionCheck->find(limits))==collisionCheck->end()) //if the region hasn't been outputed yet, do
			{
				//int id=collisionCheck->size();

				vgep=new vector<RankGffEntryExonPair>;
				limits.JnxID=(JnxID++);
				collisionCheck->insert(map<Keylet, vector<RankGffEntryExonPair>* >::value_type(limits,vgep));


			}else
			{

				vgep=(*j).second;
			}

			vgep->push_back(RankGffEntryExonPair(i,entry.exons[i],entry.exons[i+1],&entry));

		}

	}

	int novoIDStart=JnxID;
	//JnxID=0;

	cerr<<"total "<<lc<<" lines processed for known exons and transcripts"<<endl;
	cerr<<"Now Innovate! de Novo! ~~ Hakuna Matata!"<<endl;

	//map<Keylet,KeyPair<Exon*,Exon*> > novoInfo;

	//walk through Exonome! Exonome is actually ordered first by Chr, then by start, and then by end
	// and start gluing exons less than dist bp apart together! Hukuna Matata!
	curChr="";

	int novoExonCount=0;

	for (map<GffEntry::ExonPtr,vector<RankedGffEntryPair>*>::iterator ex=
			GffEntry::exonome.begin(); ex!=GffEntry::exonome.end(); ex++) {
		map<GffEntry::ExonPtr,vector<RankedGffEntryPair>*>::iterator cur=ex; //this points to the current exon;
		cur++; //cur points to the next downstream exon

		 GffEntry::Exon& eBait=*((*ex).first);
		const vector<RankedGffEntryPair>* assTBait=(*ex).second;

		RankedGffEntryPair rgepBait=(*assTBait)[0];

		//go right!

		while (cur!=GffEntry::exonome.end()) //hitting the outermost gate
		{

			 GffEntry::Exon& ePrey=*((*cur).first);
			const vector<RankedGffEntryPair>* assTranscripts=(*cur).second;
			if (curChr!=eBait.chr) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
			{

				if (raf) //remember to close and delete the previous RandomAccessFile Object
				{
					raf->close();
					delete raf;
				}

				curChr=eBait.chr; //and then open a new seq
				raf=new RandomAccessFile(chrDir+curChr+".seq");

				map<string, map<Keylet, vector<RankGffEntryExonPair>*>*>::iterator
						alr;

				if ((alr=alreadyin.find(curChr))==alreadyin.end()) {
					//no record of curChr //shouldn't happen I guess?
					collisionCheck=new map<Keylet, vector<RankGffEntryExonPair>*>;
					alreadyin.insert(map<string,
							map< Keylet, vector<RankGffEntryExonPair>*>* >::value_type(
							curChr, collisionCheck));

				} else {

					collisionCheck=(*alr).second;
				}

			}

			RankedGffEntryPair rgep=(*assTranscripts)[0]; //get the first associated transcript, for strand. assume same strand for all ass transcripts

			if (ePrey.chr!=eBait.chr || ePrey.start>eBait.end+dist) //hitting the wall!
			{
				break;
			}

			if (rgep.gff->strand!=rgepBait.gff->strand || eBait.end>ePrey.start-minIntronLength) // not the same strand, forget it, continue to the next Prey
			{
				cur++;
				continue;
				//go right!
			}

			int gUp=eBait.end; //1 index incl
			int gDown=ePrey.start; //0 index incl
		//	int lFlank=flank;
		//	int rFlank=flank;

			string extra;

			leftUsed=0;
			rightUsed=0;
			string lExtra="";
			string rExtra="";

			bool boundOut=false;
			if (eBait.end-eBait.start+1<flank) {
				boundOut=true;

			}

			if (ePrey.end-ePrey.start+1<flank) {
				boundOut=true;
			}

			//check for in frame
			novoInfo ninfo;
			ninfo.exons.k1=&eBait;
			ninfo.exons.k2=&ePrey;



			bool inFrame=false;
			char strand=rgepBait.gff->strand;
			if (strand==GffEntry::FORWARD) {
				if (eBait.frame==-1) {
					inFrame=(ePrey.frame==-1 || (rgepBait.gff->cdsStart>=ePrey.start && rgepBait.gff->cdsStart<=ePrey.end));
				} else if (eBait.frame==0) {
					int calFrame=(eBait.end-rgepBait.gff->cdsStart)%3;
					inFrame=(ePrey.frame==calFrame);
				} else {
					int calFrame=(eBait.end-eBait.start+eBait.frame)%3;
					inFrame=(ePrey.frame==calFrame);
				}

			} else {
				if (ePrey.frame==-1) {
					inFrame=(eBait.frame==-1 || (rgep.gff->cdsEnd>=eBait.start && rgep.gff->cdsEnd<=eBait.end));
				} else if (ePrey.frame==0) {
					int calFrame=(rgep.gff->cdsEnd-ePrey.start)%3;
					inFrame=(eBait.frame==calFrame);
				} else {
					int calFrame=(ePrey.end-ePrey.start+ePrey.frame)%3;
					inFrame=(eBait.frame==calFrame);
				}

			}

			ninfo.inFrame=inFrame;

			if (requireFuncFusion && !inFrame) {
				cur++;
				continue;
			}


			if (!boundOut) {
				//easy case:
				extra=StringUtil::str(gUp-flank)+"/"+StringUtil::str(gUp)+":"
						+StringUtil::str(gDown)+"/"
						+StringUtil::str(gDown+flank);

				Keylet limits(gUp, gDown, extra);

				vector<RankGffEntryExonPair> *vgep;

				map< Keylet, vector<RankGffEntryExonPair>*>::iterator j;

				if ((j=collisionCheck->find(limits))==collisionCheck->end()) {
					limits.JnxID=(JnxID++);

					snip<<">"<<setName<<":"<<limits.JnxID<<endl;
					JnxInfoS.push_back(ninfo);

					vgep=new vector<RankGffEntryExonPair>;

					raf->transfer(snip, gUp-flank, gUp);
					raf->transfer(snip, gDown, gDown+flank);

					snip<<endl;
					collisionCheck->insert(map<Keylet,
							vector<RankGffEntryExonPair>* >::value_type(limits,
							vgep));

					novoExonCount++;

					if (novoExonCount%100==1) {
						cerr<<novoExonCount<<" jnx(s) de novo invented"<<endl;
					}

				} else {
					//ignored for now
					//vgep=(*j).second;
				}

				//vgep->push_back(BiheadExoPair(&eBait,&ePrey));

			} else {
				//don't care about short exons yet

				//cerr<<"short exon begin proc"<<endl;
				vector<RankedGffEntryPair> &assT=*(eBait.assTranscripts);

				vector< vector<KeyPair<int,int> > > leftons;
				vector< vector<KeyPair<int,int> > > rightons;
				vector<string> leftExtras;
				vector<string> rightExtras;

				int nAssT=assT.size();
				for (int x=0; x<nAssT; x++) {
					RankedGffEntryPair &rgep=assT[x];
					int from=rgep.gff->getExonArrayIndex(rgep.rank);
					int lFlank=flank;
					string extra;
					vector<KeyPair<int,int> > path;
					//cerr<<"a1"<<endl;
					int g=0;

					while (lFlank>0 && from>=0) {
						const GffEntry::Exon *thisExon=rgep.gff->exons[from];

						int lenThisTime;
						if (from==0) {
							lenThisTime=lFlank;
						} else {
							lenThisTime=MIN(lFlank,thisExon->getLength());
						}

						KeyPair<int,int> B(thisExon->end-lenThisTime,
								thisExon->end);
						path.push_back(B);
						extra=StringUtil::str(B.k1)+"/"+StringUtil::str(B.k2)
								+(((g++)==0) ? "" : ("/"+extra));

						lFlank-=lenThisTime;

						from--;
					}
					leftons.push_back(path);
					leftExtras.push_back(extra);
				}

				//	cerr<<"a2"<<endl;
				vector<RankedGffEntryPair> &assTPrey=*(ePrey.assTranscripts);

				nAssT=assTPrey.size();
				//int g=0;

				for (int x=0; x<nAssT; x++) {
					RankedGffEntryPair &rgep=assTPrey[x];

					int from=rgep.gff->getExonArrayIndex(rgep.rank);

					int lFlank=flank;
					string extra;
					vector<KeyPair<int,int> > path;
					int g=0;
					while (lFlank>0 && from<rgep.gff->exonCount) {
						const GffEntry::Exon *thisExon=rgep.gff->exons[from];

						int lenThisTime;

						if (from==rgep.gff->exonCount-1) {
							lenThisTime=lFlank;
						} else {
							lenThisTime=MIN(lFlank,thisExon->getLength());
						}

						KeyPair<int,int> B(thisExon->start, thisExon->start
								+lenThisTime);

						path.push_back(B);
						extra+=(((g++)==0) ? "" : "/")+StringUtil::str(B.k1)
								+"/"+StringUtil::str(B.k2);

						lFlank-=lenThisTime;

						from++;
					}

					rightons.push_back(path);
					rightExtras.push_back(extra);
				}

				//cerr<<"a3"<<endl;

				int lefts=leftons.size();
				int rights=rightons.size();
				//all combination of leftons and rightons
				for (int x=0; x<lefts; x++)
					for (int y=0; y<rights; y++) {
						string &lExtra=leftExtras[x];
						string &rExtra=rightExtras[y];
						string extra=lExtra+":"+rExtra;

						vector<KeyPair<int,int> > &leftBounds=leftons[x];
						vector<KeyPair<int,int> > &rightBounds=rightons[y];

						vector<RankGffEntryExonPair> * vgep;

						Keylet limits(gUp, gDown, extra);

						map< Keylet, vector<RankGffEntryExonPair>*>::iterator j;

						if ((j=collisionCheck->find(limits))
								==collisionCheck->end()) //if the region hasn't been outputed yet, do
						{
							//int id=collisionCheck->size();

							vgep=new vector<RankGffEntryExonPair>;
							limits.JnxID=(JnxID++);
							collisionCheck->insert(map<Keylet,
									vector<RankGffEntryExonPair>* >::value_type(
									limits, vgep));

							snip<<">"<<setName<<":"<<limits.JnxID<<endl;
							JnxInfoS.push_back(ninfo);

							novoExonCount++;

							if (novoExonCount%100==1)
								cerr<<novoExonCount<<" jnx(s) de novo invented"
										<<endl;

							//cerr<<"smalljnx:"<<extra;

							int lbs=leftBounds.size();
							int rbs=rightBounds.size();

							for (int z=lbs-1; z>=0; z--) {
								raf->transfer(snip, leftBounds[z].k1,
										leftBounds[z].k2);
								//raf->transfer(cerr,leftBounds[z].k1,leftBounds[z].k2);
							}
							for (int z=0; z<rbs; z++) {
								raf->transfer(snip, rightBounds[z].k1,
										rightBounds[z].k2);
								//raf->transfer(cerr,rightBounds[z].k1,rightBounds[z].k2);

							}

							snip<<endl;

						} else {

						}

					}

			}

			cur++;
		}

	}

	if (raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
	{
		raf->close();
		delete raf;
	}

//	int tjnx=0;
	vector<string> tmp;

	cerr<<"finish writing fasta"<<endl;

	for (map<string, map<Keylet,vector<RankGffEntryExonPair>* >*>::iterator i=alreadyin.begin();
		i!=alreadyin.end(); i++)
	{
		//int numJnx=(*i).second->size();


		//cerr<<numJnx<<" jnxs in "<<(*i).first<<endl;
		//tjnx+=numJnx;


		map<Keylet, vector<RankGffEntryExonPair>* > *pJnxEntry=(*i).second;

		string curChr=(*i).first;




		for(map<Keylet, vector<RankGffEntryExonPair>* >::iterator j=pJnxEntry->begin();
			j!=pJnxEntry->end();
			j++)
		{

			vector<RankGffEntryExonPair>* pGffEntries=(*j).second;//pGffEntries=(*j).second;

			const Keylet& key=(*j).first;

			if(key.JnxID>=novoIDStart)
			{
				novoInfo &einfo=JnxInfoS[key.JnxID-novoIDStart];

				jnxinfo<<setName<<":"<<key.JnxID<<"\t"<<curChr<<":"<<key.extra<<"\t";
				jnxinfo<<(*(einfo.exons.k1->assTranscripts))[0].gff->strand<<"\t"<<((einfo.inFrame)?"inframe":"notinframe")<<"\t";
				jnxinfo<<einfo.exons.k1->exonID<<"\t"<<einfo.exons.k1->frame<<"\t"<<einfo.exons.k2->exonID<<"\t"<<einfo.exons.k2->frame;
				jnxinfo<<endl;
			}


			delete pGffEntries;
		}

		delete pJnxEntry;
	}

	delete[] leftCoord;
	delete[] rightCoord;

	//cerr<<"Total: "<<tjnx<<endl;
	//cerr<<"BoundOut Count= "<<boundOutCount<<endl;
	cerr<<"Total jnx de novo: "<<novoExonCount<<endl;
	GffEntry::resetEverything(); //free memory
	gff.close();
	jnxinfo.close();
	snip.close();
}

void getTS(int upstream,int downstream,string gfffile,string chrDir,string snipfile,string conservedblockfile,string chrLabel)
{

	ifstream gff(gfffile.c_str()); //input file stream for gff file
	ofstream snip(snipfile.c_str()); //output file stream for snip file
	char line[GFF_BUFFER_SIZE]; //the buffer for reading lines from gff file

	ifstream cbfin(S(conservedblockfile)); //the conserved block file stream (input)

	ConservedBlockQueue cbq(cbfin,snip); //the conserved block file reader (input from the conserved file stream, output to the snip file)

	set<KeyPair<int,int> > alreadyin; //to remember which regions are already included (because multiple transcripts (of same gene) can have same promoter
	map<string,int> names;

	RandomAccessFile *raf=NULL; //the random access file to get the genomic sequence
	string curChr=""; //the id of the current chr, used to decide to whether to open a new genomic seq file


	GffEntry::resetEverything(); //reset the exonome

	if(!gff.good())
		{
		curChr=chrLabel;
			cerr<<"no_info for "<<chrLabel<<endl;
			raf=new RandomAccessFile(chrDir+curChr+".seq");
			snip<<">\t"<<CHROMOSOME_NOINFO<<"\t"<<chrLabel;
			snip<<endl;
			raf->transfer(snip,0,-1);
			snip<<endl;
			cbq.output(0,INT_MAX);
			snip<<endl;
			snip.close();
			raf->close();
			delete raf;
			return;
		}

	while(!gff.eof()) //read all entries from gff file
	{

		gff.getline(line,GFF_BUFFER_SIZE);



		vector<string> spliton;  //vector to carry the split tokens from each line of gff file
		StringUtil::split(string(line),string("\t"),spliton);

		if(spliton.size()<1)
			continue;
		if(spliton[0]=="#bin") //ignore commented lines from gff files
			continue;

		GffEntry &entry = *(GffEntry::createGffEntry(spliton)); //construct a GffEntry from the tokens

		map<string,int>::iterator i;
		set< KeyPair<int,int> >::iterator j;

		int seeks,seeke; //to remember the start and end of the region to get sequence

		int gupstream;
		int gdownstream;

		//upstream and downstream (of promoters) are with respect to coding strand, so get read start and end of the genomic strand
		if(entry.strand==entry.REVERSE)
		{
			seeks=entry.txEnd-downstream;
			seeke=entry.txEnd+upstream;
			gupstream=downstream;
			gdownstream=upstream;
		}
		else
		{
			seeks=entry.txStart-upstream;
			seeke=entry.txStart+downstream;
			gupstream=upstream;
			gdownstream=downstream;
		}

		KeyPair<int,int> limits(seeks,seeke);



		if((j=alreadyin.find(limits))==alreadyin.end()) //if the region hasn't been outputed yet, do
		{
			//remember the frequency of the gene having alternative starts
			int c=1;
			if((i=names.find(entry.name2))!=names.end())
			{
				c=(*i).second;

				names.erase(i);
				names.insert(map<string,int>::value_type(entry.name2,++c));
			}
			else
			{
				names.insert(map<string,int>::value_type(entry.name2,1));
			}

			//snip<<">"<<"\t"<<entry.chrom<<"\t"<<limits<<"\t"<<entry
			//snip<<">"<<entry.name2<<"\t"<<c<<"\t"<<entry.strand<<"\t"<<limits<<endl;
			alreadyin.insert(set<KeyPair<int,int> >::value_type(limits));

			if(curChr!=entry.chrom) //if the current chrom to get seq from is not yet open, open it. else continue to use the current chrom
			{
				if(raf) //remember to close and delete the previous RandomAccessFile Object
				{
					raf->close();
					delete raf;
				}

				curChr=entry.chrom; //and then open a new seq
				raf=new RandomAccessFile(chrDir+curChr+".seq");

			}


			snip<<">\t"<<curChr<<"\t"<<seeks<<"\t"<<seeke<<"\t"<<(COMMONEINFO2)<<TYPE_TXSTART<<";"<<seeks<<";"<<seeke<<"\t";
			snip<<".\t."<<endl;
			raf->transfer(snip,seeks-1,seeke); //annotation starts from 1.
			//raf->transfer(snip,limits.k1-1,limits.k2); //get the genomic sequence limited by the region
			snip<<endl; //write the sequence to snip file
			cbq.output(seeks,seeke); //output the conserved block between mouse/human into snip file

		}
	}

	if(raf) //final clean up, close RandomAccessFile object, Gff object and the snip file
	{
		raf->close();
		delete raf;
	}
	GffEntry::resetEverything(); //free memory
	gff.close();
	snip.close();
}

void findUniqueFor(string annoSource,string snipfile)
{
	ofstream fil(snipfile.c_str());
	//vector<int> ids;	
	
	for (map<string,multimap<GffEntry::ExonPtr,GffEntry::Locus*>* >::iterator j=
			GffEntry::loci.begin(); j!=GffEntry::loci.end(); j++) {
		
		multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus=(*j).second;
		
		for (map<GffEntry::ExonPtr, GffEntry::Locus*>::iterator i=
				chrMapLocus->begin(); i!=chrMapLocus->end(); i++) {
			bool unique=true;

			GffEntry::Locus* pLoc=(*i).second;
			for (vector<GffEntry*>::iterator j=pLoc->transcripts.begin(); j
					!=pLoc->transcripts.end(); j++) {
				GffEntry*entry=(*j);

				if (entry->annoSource!=annoSource) {
					unique=false;
				}
			}

			if (unique) {
				for (vector<GffEntry*>::iterator j=pLoc->transcripts.begin(); j
						!=pLoc->transcripts.end(); j++) {
					fil<<(*j)->transcriptID<<endl;
				}
			}

		}
	}
	fil.close();
}

void splitGeneInfo(string infilename,string outfileprefix)
{
	ifstream gff(infilename.c_str());
	ofstream *output=NULL;

	char line[BUFFER_SIZE];

	string chrom;
	string curchrom("");


	set<string> opened_chrom;


	while(!gff.eof())

	{
		gff.getline(line,BUFFER_SIZE);

		vector<string> spliton;
		StringUtil::split(string(line),string("\t"),spliton);
		if(spliton.size()<1)
			continue;
		if(spliton[0]=="#bin")
			continue;


		chrom=spliton[2];
		if(curchrom!=chrom)
		{

			if(output)
			{
				cerr<<"closing "<<curchrom<<endl;
				output->close();
				delete output;
			}
			curchrom=chrom;
			cerr<<"writing to "<<curchrom<<endl;
			if(opened_chrom.find(curchrom)==opened_chrom.end())
			{
				opened_chrom.insert(set<string>::value_type(curchrom));
				output=new ofstream((outfileprefix+curchrom+".txt").c_str(),ios::out|ios::trunc);
			}
			else
			{
				output=new ofstream((outfileprefix+curchrom+".txt").c_str(),ios::out|ios::app);
			}
		}

		(*output)<<line<<endl;

	}

	if(output)
	{
		cerr<<"closing "<<curchrom<<endl;
		output->close();
		delete output;
	}

	gff.close();

}

int calculateConstituteExonFreq_inner_numOfReads(GffEntry::GBlock* block,GffEntry::GBlock::DIterator range)
{
	int count=0;

	int countMax=block->selexaMatches.size();

	for(GffEntry::GBlock::SIterator i=range.first;i!=range.second;i++)
	{
		count++;
		if(count>countMax)
		{
			cerr<<"count is bigger than max count"<<endl;
			die_exit("");

		}
	}
	return count;
}

inline pair<int,int> calculateConstituteExonFreq_inner_getChainResult(int readLength,map<GffEntry::GBlockPtr,BlockCountStruct>::iterator begin0,map<GffEntry::GBlockPtr,BlockCountStruct>::iterator end0,bool useUniquelyMappablePosition)
{
	map<GffEntry::GBlockPtr,BlockCountStruct>::iterator i;
	
	int pos=0;
	int count=0;
	for(i=begin0;i!=end0;i++)
	{

		GffEntry::GBlockPtr block=(*i).first;
		BlockCountStruct blockcount=i->second;
		
		if(blockcount.Bound!=block->getBound())
		{

			//Here
		//	cerr<<"aaa: block get"<<blockcount.getStart0()<<"-"<<blockcount.getEnd0()<<endl;
			count+=calculateConstituteExonFreq_inner_numOfReads(block,block->getSelexaMatchesWithinGenomicBound00(blockcount.getStart0(),blockcount.getEnd0()));
			pos+=block->getUniquelyMappablePositions(blockcount.Bound);
		//	cerr<<"aaa2"<<endl;

		}
		else
		{
			//cerr<<"ccc"<<endl;
			count+=block->selexaMatches.size();
		
			if(useUniquelyMappablePosition)
				pos+=block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR();
			else
				pos+=block->getNaiveMappablePos(); //the whole length is good
		}
	}
	

	GffEntry::GBlockPtr endblock=(*end0).first;
	BlockCountStruct blockcount=end0->second;

	
	if(blockcount.Bound!=endblock->getBound())
	{
		//Here

		KeyPair<int,int> useBound=::overlapBound(blockcount.Bound,::b00tob01(endblock->getBound00WithinGenomicPosForReadLength(readLength)));
		//cerr<<"bbb: block get"<<useBound.k1<<"-"<<useBound.k2<<endl;
		count+=calculateConstituteExonFreq_inner_numOfReads(endblock,endblock->getSelexaMatchesWithinGenomicBound00(useBound.k1,useBound.k2-1));
		pos+=endblock->getUniquelyMappablePositions(useBound);
		//cerr<<"bbb2"<<endl;
	}
	else
	{
		//cerr<<"ddd"<<endl;

		count+=calculateConstituteExonFreq_inner_numOfReads(endblock,endblock->getSelexaMatchesWithStartWithinRelPos(endblock->getBound00WithinRelPosForReadLength(readLength)));
		//cerr<<"ddd2"<<endl;
		if(useUniquelyMappablePosition)
			pos+=endblock->getUniquelyMappablePosL();
		else
			pos+=endblock->getNaiveMappablePos(readLength);
	}
	GffEntry::gTotalExonReadsUsed+=count;
	
	return pair<int,int>(count,pos);
}

/*void statEI(string chr,ostream& os,bool useUniquelyMappablePosition)
{
	
	
	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}
	
	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	
	if(useUniquelyMappablePosition && !GffEntry::isUniquelyMappablePositionLoaded())
	{
		die("Request to use Uniquely Mappable position but UMPs are not loaded");
	}

	//map<string, set<GBlockPtr> *> gblocks
	KeyPair<int,int> statChrExonic(0,0);
	KeyPair<int,int> statChrNonExonic(0,0);

	map<string,set<GffEntry::GBlockPtr>* >::iterator i=GffEntry::gblocks.find(chr);
	

 	//string chr=(*i).first;
	set<GffEntry::GBlockPtr>* sset=(*i).second;
	for(set<GffEntry::GBlockPtr>::iterator j=sset->begin();j!=sset->end();j++)
	{
			GffEntry::GBlockPtr block=(*j);
			if(block->isExonic())
			{
	
				statChrExonic.k1+=block->getSelexaMatchSize();
				if(useUniquelyMappablePosition)
					statChrExonic.k2+=block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR();
				else 
					statChrExonic.k2+=block->getLength();
				
			}else
			{
				statChrNonExonic.k1+=block->getSelexaMatchSize();
				if(useUniquelyMappablePosition)
					statChrNonExonic.k2+=block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR();
				else
					statChrNonExonic.k2+=block->getLength();
			}
			
			
			
	}
		
	//double readCountExonic=statChrExonic.k1;
	//double lengthExonic=statChrExonic.k2;
	//double readCountNonExonic=statChrNonExonic.k1;
	//double lengthNonExonic=statChrNonExonic.k2;	
	
	os<<chr<<"\t"<<statChrExonic.k1<<"\t"<<statChrExonic.k2<<"\t\t"<<statChrNonExonic.k1<<"\t"<<statChrNonExonic.k2<<endl;
	
}*/


void statEIReturnValues(string chr,ostream& os,bool useUniquelyMappablePosition,KeyPair<uint64_t,uint64_t>& statChrExonic,KeyPair<uint64_t,uint64_t>& statChrNonExonic,bool resetCount)
{
	
	if(resetCount)
	{
		statChrNonExonic.k1=0;
		statChrNonExonic.k2=0;
		statChrExonic.k1=0;
		statChrExonic.k2=0;
	}
	
	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}
	
	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	
	if(useUniquelyMappablePosition && !GffEntry::isUniquelyMappablePositionLoaded())
	{
		die("Request to use Uniquely Mappable position but UMPs are not loaded");
	}

	//map<string, set<GBlockPtr> *> gblocks
	//KeyPair<int,int> statChrExonic(0,0);
	//KeyPair<int,int> statChrNonExonic(0,0);

	map<string,set<GffEntry::GBlockPtr>* >::iterator i=GffEntry::gblocks.find(chr);
	
	if(i!=GffEntry::gblocks.end()){

 	//string chr=(*i).first;
	set<GffEntry::GBlockPtr>* sset=(*i).second;
	for(set<GffEntry::GBlockPtr>::iterator j=sset->begin();j!=sset->end();j++)
	{
			GffEntry::GBlockPtr block=(*j);
			if(block->isExonic())
			{
	
				statChrExonic.k1+=block->getSelexaMatchSize();
				if(useUniquelyMappablePosition)
					statChrExonic.k2+=block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR();
				else 
					statChrExonic.k2+=block->getLength();
				
			}else
			{
				statChrNonExonic.k1+=block->getSelexaMatchSize();
				if(useUniquelyMappablePosition)
					statChrNonExonic.k2+=block->getUniquelyMappablePosL()+block->getUniquelyMappablePosR();
				else
					statChrNonExonic.k2+=block->getLength();
			}
			
			
			
	}
	}
		
	//double readCountExonic=statChrExonic.k1;
	//double lengthExonic=statChrExonic.k2;
	//double readCountNonExonic=statChrNonExonic.k1;
	//double lengthNonExonic=statChrNonExonic.k2;	
	
	os<<chr<<"\t"<<statChrExonic.k1<<"\t"<<statChrExonic.k2<<"\t\t"<<statChrNonExonic.k1<<"\t"<<statChrNonExonic.k2<<endl;
	
}



void calculateConstituteExonFreq(string snipFile,string chr,double ratioRep,int readLength,bool ambigCheck,bool codingCheck,bool useUniquelyMappablePosition,double ignoreExonOfLowerPercentile,int ignoreExonOfFromMax,int ignoreExonOfFromMin)
{
	ofstream snip(snipFile.c_str());
	
	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}
	
	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}
	
	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}
	
	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	
	if(useUniquelyMappablePosition && !GffEntry::isUniquelyMappablePositionLoaded())
	{
		die("Request to use Uniquely Mappable position but UMPs are not loaded");
	}
	//now here comes something!!
	
	//for each Locus
	
	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;
	
	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);
	
	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}
	
	chrMapLocus=(*mi).second;
	
	
	//multimap<GffEntry::GBlockPtr,int> blockCountMap;
	
	
	/* algorithm get gblock of freq/totalTranscripts > ratioRep:	
	foreach locus l
		foreach transcript t in l
			foreach exon e in t
				foreach block b in e
					add b to blockCountMap <block,count>
					and increment count by 1
				rof
			rof
		rof
	rof
	*/
	
	int lc=0;
	int tc=chrMapLocus->size();
	


	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;
		


		map<GffEntry::GBlockPtr,BlockCountStruct> blockCountMap;


		GffEntry::Locus* locus=(*l).second;
		vector<GffEntry*>& transcripts=locus->transcripts;
		blockCountMap.clear();
		
		int rtranscriptCount=transcripts.size();
		
		////new
		int maxExonCount=0;
		int minExonCount=INT_MAX;
		
		
		
		//multimap<int,GffEntry*> orderedTranscripts;
		
		int *ecounts=new int[rtranscriptCount];
		//set<int> lengthsECount;
		
		int ei=0;
		
		//cerr<<"a"<<endl;

		for(vector<GffEntry*>::iterator t=transcripts.begin();t!=transcripts.end();t++)
		{
			int thisExonCount=(*t)->exonCount;
			if(thisExonCount>maxExonCount)
				maxExonCount=thisExonCount;
			
			if(thisExonCount<minExonCount)
				minExonCount=thisExonCount;


			//GffEntry* transcript=(*t);
			ecounts[ei++]=thisExonCount;//transcript->exonCount;
			
		}
		
		
		//cerr<<"b"<<endl;
		
		int LowerThresholdExonCount=INT_MIN;

		if(ignoreExonOfLowerPercentile>0){
			int precentileIndex=int(/*ceil*/((double)rtranscriptCount*ignoreExonOfLowerPercentile)-1);
			precentileIndex=MIN(precentileIndex,rtranscriptCount-1);
			precentileIndex=MAX(precentileIndex,0);
			nth_element(ecounts,ecounts+precentileIndex,ecounts+rtranscriptCount);
			LowerThresholdExonCount=ecounts[precentileIndex];
		}
		else if(ignoreExonOfFromMax>=0)
		{
			LowerThresholdExonCount=MAX(maxExonCount-ignoreExonOfFromMax,1);
		}
		else if(ignoreExonOfFromMin>=0)
		{
			LowerThresholdExonCount=minExonCount+ignoreExonOfFromMin;
		}

		//cerr<<"c"<<endl;
		
		
		delete[] ecounts;
		
		int transcriptCount=0;
		
		string _name=locus->getFirstName();
		cerr<<"processing locus "<<_name<<endl;
		cerr<<"lower threshold exoncount="<<LowerThresholdExonCount<<endl;




		for(vector<GffEntry*>::iterator t=transcripts.begin();t!=transcripts.end();t++)
		{
		
		
			
			GffEntry* transcript=(*t);
			
			
			////new
			if(transcript->exonCount<LowerThresholdExonCount)
				continue;
			
			transcriptCount++;
			////new
			cerr<<"transcript "<<transcript->name<<" of "<<transcriptCount<<" with "<<transcript->exonCount<<" exons"<<endl;
			
			for(int e=0;e<transcript->exonCount;e++)
			{
				
				GffEntry::ExonPtr exon=transcript->exons[e];
				if(!exon->blocks)
				{
					cerr<<">>>>>>>>>>>>>>>>>>>no block for exon"<<endl;
					continue;
				}

				//cerr<<"d"<<endl;
				for(vector<GffEntry::GBlockPtr>::iterator b=exon->blocks->begin();b!=exon->blocks->end();b++)
				{
					GffEntry::GBlockPtr block=(*b);
					cerr<<"block:"<<block->getStart1()<<"-"<<block->getEnd1()<<"\t";
					if(ambigCheck && block->getNaiveAmbiguity())
					{
						cerr<<"ambiguous"<<endl;
						continue;
					}
					
					KeyPair<int,int> bound=block->getBound();

					if(codingCheck)
					{
						if(block->end1<exon->maxCodingStart)
						{
							cerr<<"block end < coding start"<<endl;
							map<GffEntry::GBlockPtr,BlockCountStruct>::iterator i=blockCountMap.find(block);
							if(i!=blockCountMap.end())
								blockCountMap.erase(i);
							continue;
						}

						if(block->start0>=exon->minCodingEnd)
						{
							cerr<<"block start > coding end"<<endl;
							map<GffEntry::GBlockPtr,BlockCountStruct>::iterator i=blockCountMap.find(block);
							if(i!=blockCountMap.end())
								blockCountMap.erase(i);
							continue;
						}

						bound=KeyPair<int,int>(MAX(exon->maxCodingStart,block->start0),MIN(exon->minCodingEnd,block->end1));
						if(!::isBound01Valid(bound))
							continue;
					}


					//cerr<<"e"<<endl;
					map<GffEntry::GBlockPtr,BlockCountStruct>::iterator i=blockCountMap.find(block);
					if(i==blockCountMap.end())
					{
						cerr<<"inserted"<<endl;
						blockCountMap.insert(map<GffEntry::GBlockPtr,BlockCountStruct>::value_type(block,BlockCountStruct(1,bound)));
					}
					else
					{
						cerr<<"incremented"<<endl;
						i->second++;
						i->second.intersect(bound);
					}

				}//rof b
			}//rof e
		}//rof t
		

		//now destroy those with invalid bound

		vector<map<GffEntry::GBlockPtr,BlockCountStruct>::iterator> toErase;

		for(map<GffEntry::GBlockPtr,BlockCountStruct>::iterator i=blockCountMap.begin();i!=blockCountMap.end();i++)
		{
			if(!isBound01Valid(i->second.Bound))
				toErase.push_back(i);
		}

		for(vector<map<GffEntry::GBlockPtr,BlockCountStruct>::iterator>::iterator i=toErase.begin();i!=toErase.end();i++)
		{
			blockCountMap.erase(*i);
		}

		cerr<<"count map created"<<endl;
		//now go through all blocks in the BlockCountMap
		//find contaguous pieces of blocks with count/transcriptTotal>ratioRep
		//get number of possible positions for each piece and the associated read counts
		//add to the records the pos and the reads
		
		map<GffEntry::GBlockPtr,BlockCountStruct>::iterator lMarker;
		
		string _chrom=locus->chr;
		char _strand=locus->strand;

		int _start=INT_MAX;
		int _end=INT_MIN;
		vector<int> _starts;
		vector<int> _sizes;




		snip<<locus->chr<<"\t";
		snip<<locus->strand<<"\t";
		
		/*for(set<string>::iterator i=locus->names.begin();i!=locus->names.end();i++)
		{
			snip<<(*i)<<",";
		}*/ //changed to below
		
		//set<string>::iterator nameI=locus->names.begin();
		//if(nameI!=locus->names.end())
		//{
		snip<<_name;
		//}
		
		snip<<"\t";
		

		
		//output basic info for block
		
		//find the first eligible block
		
		for(lMarker=blockCountMap.begin();lMarker!=blockCountMap.end();lMarker++)
		{
			//int count=(*lMarker).second;
			GffEntry::GBlockPtr pbl=(*lMarker).first;
			BlockCountStruct blockcount=lMarker->second;
			int count=blockcount.getCount();

			//
			if((double)count/transcriptCount>=ratioRep)
			{
				//cerr<<" accepted"<<endl;
				//_start=pbl->start0;
				_start=blockcount.getStart0();
				//_end=_start;
				//_starts.push_back(0); //pbl->getStart0()-_start
				//_sizes.push_back(blockcount.getLength()); //pbl->getLength()
				//_scores.push_back(((double)count/transcriptCount)*1000);
				//_end=blockcount.getEnd1();//pbl->getEnd1();

				break;
			}else
			{
				cerr<<"block "<<locus->chr<<":"<<pbl->getStart1()<<"-"<<pbl->getEnd1()<<"\t"<<"representing "<<count<<"/"<<transcriptCount<<"("<<setprecision(3)<<((double)count/transcriptCount*100)<<"%) of transcript ,clipped as "<<locus->chr<<":"<<blockcount.getStart1()<<"-"<<blockcount.getEnd1();

				cerr<<" ignored: low representative"<<endl;
			}

		}
		
		//lMarker reached end
		if(lMarker==blockCountMap.end())
		{	
			cerr<<"No Blocks are good for locus "<<_name<<endl;
			snip<<"NBG\t;";
			snip<<endl;
			continue;
		}
		





		//an imaginative end
		map<GffEntry::GBlockPtr,BlockCountStruct>::iterator gi;
		gi=lMarker;
		
		int prevEnd=gi->second.getStart0();  /////gi->first->start0;
		
		int readcount=0;
		int pos=0;
		
		string blockIResult;
		int nBlockCounted=0;
		
		
		map<GffEntry::GBlockPtr,BlockCountStruct >::iterator prevGI;
		
	//	int bn=0;
		bool prevIsMinor=false;
		
		for(;gi!=blockCountMap.end();gi++)
		{
			//bn++;
//			cerr<<"a"<<endl;
			

//			cerr<<"b"<<endl;
			GffEntry::GBlockPtr pbl=gi->first;
			BlockCountStruct blockcount=gi->second;

			int count=blockcount.getCount();
//			cerr<<"c"<<endl;			
			//moved  bool prevIsMinor=false;
			
			cerr<<"block "<<locus->chr<<":"<<pbl->getStart1()<<"-"<<pbl->getEnd1()<<"\t"<<"representing "<<count<<"/"<<transcriptCount<<"("<<setprecision(3)<<((double)count/transcriptCount*100)<<"%) of transcript ,clipped as "<<locus->chr<<":"<<blockcount.getStart1()<<"-"<<blockcount.getEnd1();

			
			bool isMinor=((double)count/transcriptCount<ratioRep);
			
			
					
			
			if(!isMinor)
			{
				cerr<<" accepted"<<endl;

				_starts.push_back(blockcount.getStart0()-_start); //pbl->getStart0()-_start
				_sizes.push_back(blockcount.getLength()); //pbl->getLength()
				//_scores.push_back(((double)count/transcriptCount)*1000);
				_end=blockcount.getEnd1();  //pbl->getEnd1();
			}
			else
			{
				cerr<<" ignored: low representative"<<endl;
			}
			
			if(prevEnd!=blockcount.getStart0() || isMinor) //if jumped (like: across exon) or not enough count, break the chain . was: pbl->start0
			{
	//			cerr<<"c1"<<endl;				
				if (!prevIsMinor) {
					//cerr<<"f1.1"<<endl;
					pair<int,int> result=
							calculateConstituteExonFreq_inner_getChainResult(
									readLength, lMarker, prevGI, useUniquelyMappablePosition);
					//cerr<<"f1.2"<<endl;
					readcount+=result.first;
					pos+=result.second;

					{//outresultlocal
						nBlockCounted++;
//						cerr<<"c1.3"<<endl;							
						blockIResult+=StringUtil::str(lMarker->second.getStart1());  //lMarker->first->getStart1()
						blockIResult+="-";
//						cerr<<"c1.4"<<endl;							
						blockIResult+=StringUtil::str(prevGI->second.getEnd1()); ///(*prevGI).first->getEnd1()
						blockIResult+="\t";
						blockIResult+=StringUtil::str(result.first);
						blockIResult+="\t";
						blockIResult+=StringUtil::str(result.second);
						blockIResult+="\t";
//						cerr<<"c1.5"<<endl;							
					}
					
					//cerr<<"f1.e"<<endl;
				}
//				cerr<<"c2"<<endl;					
				lMarker=gi;
//				cerr<<"c3"<<endl;					
				if(isMinor)
					lMarker++;
//				cerr<<"c4"<<endl;				
				if(lMarker==blockCountMap.end())
					break;
			}
//			cerr<<"d"<<endl;				
			prevIsMinor=isMinor;
			
			
			prevGI=gi;
			
			prevEnd=blockcount.getEnd1(); //pbl->end1;
//			cerr<<"e"<<endl;				
			if(isMinor)
				prevGI++;
//			cerr<<"f"<<endl;				
		}
		
		if (lMarker!=blockCountMap.end()) {

			//cerr<<"g"<<endl;
			pair<int,int> result=
					calculateConstituteExonFreq_inner_getChainResult(
							readLength, lMarker, prevGI, useUniquelyMappablePosition);

			//cerr<<"h"<<endl;
			readcount+=result.first;
			pos+=result.second;


			//cerr<<(*lMarker).first->chr<<"\t"<<(*lMarker).first->start0<<"\t"<<(*lMarker).first->end1<<"\t"<<_name<<"\t"<<0<<"\t"<<_strand<<"\t"<<_start<<"\t"<<_end<<"\t"<<"0,0,0"<<"\t"<<_starts.size()<<"\t";


			{//outresultlocal
				nBlockCounted++;
				blockIResult+=StringUtil::str(lMarker->second.getStart1());
				blockIResult+="-";
				blockIResult+=StringUtil::str(prevGI->second.getEnd1());
				blockIResult+="\t";
				blockIResult+=StringUtil::str(result.first);
				blockIResult+="\t";
				blockIResult+=StringUtil::str(result.second);
				blockIResult+="\t";
			}
		}
		
		{//out global Result
		
		snip<<readcount;
		snip<<"\t";
		snip<<pos;
		snip<<"\t";
		
		cerr<<"result: "<<"\t";
		cerr<<readcount;
		cerr<<"\t";
		cerr<<pos;
		cerr<<endl;
						
		}
		
		
		snip<<nBlockCounted<<"\t";
		
		snip<<blockIResult;
		
		snip<<";"<<endl;
		
		
		if(_end==_start)
			continue;

		cout<<_chrom<<"\t"<<_start<<"\t"<<_end<<"\t"<<_name<<"\t"<<0<<"\t"<<_strand<<"\t"<<_start<<"\t"<<_end<<"\t"<<"0,0,0"<<"\t"<<_starts.size()<<"\t";
		for(vector<int>::iterator kk=_sizes.begin();kk!=_sizes.end();kk++)
		{
			cout<<StringUtil::str(*kk)<<",";
		}
		cout<<"\t";

		for(vector<int>::iterator kk=_starts.begin();kk!=_starts.end();kk++)
		{
			cout<<StringUtil::str(*kk)<<",";
		}
		cout<<endl;

	}//rof l
	
	
	
	
	///////here
	snip.close();
	

	//clean up
	cerr<<"cleaning up after chr"<<chr<<endl;
	GffEntry::resetGBlocks(chr);
	//GffEntry::resetJnxTags(chr);
	cerr<<"done clean up"<<endl;

}



class SExonTraversalThread
{
public:
	GffEntry::ExonPtr u;
	int distance;
	vector<GffEntry::ExonPtr> exons;
	vector<GffEntry::Jnx*> jnxs;
	SExonTraversalThread(GffEntry::ExonPtr _u,int _distance): u(_u),distance(_distance){}
};


typedef multimap<int,SExonTraversalThread>::iterator ETSIterator;
typedef pair<ETSIterator, ETSIterator > ETDIterator;




void calculateAEP(Splidar_OpFlag op,int minDistanceA,int maxDistanceA,int minDistanceB,int maxDistanceB,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{

	cerr<<"calling calculateAEP of AEP("<<minDistanceA<<","<<maxDistanceA<<";"<<minDistanceB<<","<<maxDistanceB<<")"<<endl;



	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=l->second;
		string locusName=locus->getFirstName();//(*locus->names.begin());

		//replaced by new exon group assignment
		locus->assignExonBoundGroups(jobID);


		locus->printExonBounds(cerr);


		set<string> isoBoundRegistry; //moved out of loop 1/31/2010 night

		for(map<string,NESuperGroup*>::iterator i=locus->rightSuperGroups.begin();i!=locus->rightSuperGroups.end();i++)
		{
				//do something with range;


				NESuperGroup* sgroup=i->second;
				cerr<<"passing right super group "<<sgroup->getID()<<endl;


				//AEP_SplidarGraph_REGJ(Splidar_OpFlag _op,int _minDistanceA,int _maxDistanceA,int _minDistanceB,int _maxDistanceB,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NESuperGroupPtr  _sgroup, GffEntry::Locus* _locus,
				//string _locusName, bool _checkLocusName, int _readLength, string _eventTypeName)
				AEP_SplidarGraph_REGJ sgraph(op,minDistanceA,maxDistanceA,minDistanceB,maxDistanceB,fout,fSeqOut,rafIn,sgroup,locus,locusName,true,readLength,isoBoundRegistry,op.eventTypeName);
				sgraph.enterGraphLoop();
		}





	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}


void calculateA53SS(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{

	cerr<<"calling calculateA53SS"<<endl;



	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=l->second;
		string locusName=locus->getFirstName();

		//replaced by new exon group assignment

		locus->assignExonBoundGroups(jobID);


		locus->printExonBounds(cerr);

		set<string> isoBoundRegistry;

		for(map<string,NESuperGroup*>::iterator i=locus->rightSuperGroups.begin();i!=locus->rightSuperGroups.end();i++)
		{
				//do something with range;


				NESuperGroup* sgroup=i->second;
				cerr<<"passing right super group "<<sgroup->getID()<<endl;


				//AEP_SplidarGraph_REGJ(Splidar_OpFlag _op,int _minDistanceA,int _maxDistanceA,int _minDistanceB,int _maxDistanceB,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NESuperGroupPtr  _sgroup, GffEntry::Locus* _locus,
				//string _locusName, bool _checkLocusName, int _readLength, string _eventTypeName)
				//A53SS_Splidar_REGJ(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,NESuperGroupPtr  _sgroup,string _locusName,int _readLength):op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)

				A53SS_Splidar_REGJ sgraph(op,fout,fSeqOut,rafIn,sgroup,locus,locusName,readLength,isoBoundRegistry); //check locus name?!~



		}



	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}


void calculateATE(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{

	cerr<<"calling calculateATE"<<endl;



	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=l->second;
		string locusName=locus->getFirstName();

		//replaced by new exon group assignment

		locus->assignExonBoundGroups(jobID);


		locus->printExonBounds(cerr);

		set<string> isoBoundRegistry;

		//ATE_Splidar_REGJ(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry)
		ATE_Splidar_REGJ sgraph(op,fout,fSeqOut,rafIn,locus,locusName,readLength,isoBoundRegistry); //check locus name?!~







	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}

void calculateRI(Splidar_OpFlag op,bool parentsNeedToBeSplicedEitherSide,bool parentsNeedToBeSplicedOnBothSides, bool parentsNeedToBeSplicedBothSidesInSameTranscript,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn,bool retainedIntronsHaveToBeFreeOfExons,bool useUnambiguousRegionsOfRetainedIntron)
{

	cerr<<"calling calculateRI"<<endl;



	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=l->second;
		string locusName=locus->getFirstName();

		//replaced by new exon group assignment

		locus->assignExonBoundGroups(jobID);


		locus->printExonBounds(cerr);

		set<string> isoBoundRegistry;

		//ATE_Splidar_REGJ(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry)

		//RI_Splidar_REG(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry,bool parentsNeedToBeSpliced=true): op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)

		RI_Splidar_REG sgraph(op,fout,fSeqOut,rafIn,locus,locusName,readLength,isoBoundRegistry, parentsNeedToBeSplicedEitherSide, parentsNeedToBeSplicedOnBothSides,  parentsNeedToBeSplicedBothSidesInSameTranscript, retainedIntronsHaveToBeFreeOfExons, useUnambiguousRegionsOfRetainedIntron); //check locus name?!~







	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}
/*


void calculateSpliceMMGraph(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{


	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;




	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());


		cerr<<"locus:"<<locusName<<endl;

		cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		DenGraph dengraph(locus,GffEntry::JnxTagSet::Union);

		dengraph.printDenObjsinBEDFormat(cout);



		if(locus->strand==GffEntry::REVERSE)
		{
			cerr<<"graph before inversion"<<endl;
			cerr<<dengraph;
			dengraph.inverseGraphInPlace();

		}

		cerr<<"graph before collapse"<<endl;
		cerr<<dengraph;

		dengraph.printDenGraphDFSInSIFFormat(cerr,false,string("#")+locusName+".SIF ","goto");

		dengraph.collapseGraphInPlace();
		dengraph.printDenGraphDFSInSIFFormat(cerr,true,string("#")+locusName+".CSIF ","goto");

		cerr<<"graph after collapse"<<endl;
		cerr<<dengraph;



	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}
*/

void writeFastQEntry(ostream & os,string seqName,string seq,char quality='Z')
{
	int seqlength=seq.length();
	string qualityString(seqlength,quality);
	os<<"@"<<seqName<<endl;
	os<<seq<<endl;
	os<<"+"<<seqName<<endl;
	os<<qualityString<<endl;
	
}


class SimSeqInfo
{
public:
	string seq;
	string seqName;
	char strand;
	
	SimSeqInfo(){}
	SimSeqInfo(string _seqName,string _seq, char _strand):seq(_seq),seqName(_seqName),strand(_strand)
	{}
	string getPrefix(int length)
	{
		return seq.substr(0,length);
	}
	
};

ostream& operator <<(ostream& os, const SimSeqInfo& ssi)
{
	//cerr<<"out "<<ssi.seq<<"\t"<<ssi.seqName<<"\t"<<ssi.strand<<endl;
	os<<ssi.strand<<"\t"<<ssi.seqName<<"\t"<<ssi.seq;
	return os;
}

void simulateReadExonicReadsHashByPrefix(string snipfileprefix,string snipfilesuffix,int readLength,string chrDir,int prefixLength)
{
	//for each block
	//from position 0 of that block to l-39
	//write the fastq file
	if(!GffEntry::isGBlocksLoaded())
	{
		die("Gblock not loaded");
	}
	
	RandomAccessFile *raf=NULL;
	
	
	typedef buffered_ofstream<SimSeqInfo,vector<SimSeqInfo>,vector<SimSeqInfo>::iterator> buffed_ofstream_ssi;
	map<string,buffed_ofstream_ssi* > bufferedstreams;
	
	cerr<<"simulate exonic reads for chrs "<<GffEntry::gblocks.size()<<endl;
	
	int n=0;
	
	for(map<string, set<GffEntry::GBlockPtr> *> ::iterator i=GffEntry::gblocks.begin();i!=GffEntry::gblocks.end();i++)
	{	n++;
		if(n>=11 && n<=49)
		{
			
		}
		else
		{
			continue;
			
		}
		
		cerr<<n;
		string curChr=(*i).first;
		set<GffEntry::GBlockPtr>* sset=(*i).second;
		cerr<<">processing gblocks on chr "<<curChr<<endl;
		raf=new RandomAccessFile(chrDir+curChr+".seq");
		
		
		for(set<GffEntry::GBlockPtr>::iterator j=sset->begin();j!=sset->end();j++)
		{
			GffEntry::GBlockPtr block=(*j);
			
			if(block->strand==GffEntry::UNKNOWN)
				continue;
			
			int gblockAbsStart0=block->getStart0();
			int gblockLength=block->getLength();
			for(int k=0;k<=gblockLength-readLength;k++) //k is relative to start of gblock
			{
				
				int genomPos=gblockAbsStart0+k;
				string seq=StringUtil::toUpper(raf->get(genomPos,genomPos+readLength));
				
				string seqName=StringUtil::str(block->ID)+":"+curChr+":"+StringUtil::str(genomPos);
				
				for(int l=0;l<2;l++)
				{
					
					
					char strand=((l==0)?'F':'R');
					
					if(l==1)
						seq=reverse_complement(seq);
				
					SimSeqInfo ssi(seqName,seq,strand);
					string prefix=ssi.getPrefix(prefixLength);
					buffed_ofstream_ssi* bs=NULL;
					map<string,buffed_ofstream_ssi* > ::iterator ssii=bufferedstreams.find(prefix);
					if(ssii==bufferedstreams.end())
					{
						bs=new buffed_ofstream_ssi(snipfileprefix+prefix+snipfilesuffix);
						bufferedstreams.insert(map<string,buffed_ofstream_ssi *>::value_type(prefix,bs));
					}else
					{
						bs=(*ssii).second;
					}
				
					bs->push(ssi);
				}
			}
			
			
		}
		
		raf->close();

		delete raf;	
	}
	
	
	for(map<string,buffed_ofstream_ssi* > ::iterator ssii=bufferedstreams.begin();
			ssii!=bufferedstreams.end();
			ssii++)
	{
		buffed_ofstream_ssi *bs=(*ssii).second;
		bs->flush();
		delete bs;
	}
	
	
}
void simulateReadExonicReadsHashByChr(string snipfileprefix,string snipfilesuffix,int readLength,string chrDir)
{
	//for each block
	//from position 0 of that block to l-39
	//write the fastq file
	if(!GffEntry::isGBlocksLoaded())
	{
		die("Gblock not loaded");
	}
	
	RandomAccessFile *raf=NULL;
	ofstream *snip=NULL;
	
	cerr<<"simulate exonic reads for chrs "<<GffEntry::gblocks.size()<<endl;
	
	int n=0;
	
	for(map<string, set<GffEntry::GBlockPtr> *> ::iterator i=GffEntry::gblocks.begin();i!=GffEntry::gblocks.end();i++)
	{	n++;
		cerr<<n;
		string curChr=(*i).first;
		set<GffEntry::GBlockPtr>* sset=(*i).second;
		cerr<<"processing gblocks on chr "<<curChr<<endl;
		raf=new RandomAccessFile(chrDir+curChr+".seq");
		snip=new ofstream((snipfileprefix+curChr+snipfilesuffix).c_str());
		
		for(set<GffEntry::GBlockPtr>::iterator j=sset->begin();j!=sset->end();j++)
		{
			GffEntry::GBlockPtr block=(*j);
			
			if(block->strand==GffEntry::UNKNOWN)
				continue;
			
			int gblockAbsStart0=block->getStart0();
			int gblockLength=block->getLength();
			for(int k=0;k<=gblockLength-readLength;k++) //k is relative to start of gblock
			{
				
				int genomPos=gblockAbsStart0+k;
				string seqName=StringUtil::str(block->ID)+":"+curChr+":"+StringUtil::str(genomPos);
				writeFastQEntry(*snip,seqName,raf->get(genomPos,genomPos+readLength),'Z');
			}
			
			
		}
		
		
		raf->close();
		snip->close();
		delete raf;
		delete snip;
		
	}
}
