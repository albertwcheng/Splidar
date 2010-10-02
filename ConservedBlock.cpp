#include "snipIncludes.h"


int getActualPosGappedSeq(const char*seq, const int* sorted_pos, int *out_pos, int num_pos, char gap)
{
	int apos=0;
	int rpos=0;
	int cur_task=0;
	
	for(const char*cur=seq;*cur!='\0';cur++,apos++)
	{
		if(*cur==gap)
		{}
		else
		{

			while(rpos==*sorted_pos) //if contaguos identical request
			{
					*out_pos=apos;
					
					cur_task++;
					sorted_pos++;
					out_pos++;
			}
	
			if(cur_task==num_pos)
			{
				return cur_task;
			}
			rpos++;
		}
		
	}
	
	return cur_task; //out of range
	
}

BlockHeader::BlockHeader() {

}
BlockHeader::BlockHeader(int _chrom, int _gstart, int _gend) :
	chrom(_chrom), gstart(_gstart), gend(_gend) {

}
void BlockHeader::addName(string name) {
	names.push_back(name);
}
void BlockHeader::addName2(string name2) {
	names2.push_back(name2);
}
void BlockHeader::addExtData(string key, string value) {
	extData.insert(map<string, string>::value_type(key, value));
}
void BlockHeader::addTypeDesc(string key, string value) {
	typeDesc.insert(map<string, string>::value_type(key, value));
}

ostream& operator << (ostream& os, BlockHeader& bh)
{
	os<<bh.chrom
	  <<" "
	  <<bh.gstart
	  <<" "
	  <<bh.gend
	  <<" "
	  <<StringUtil::formKeyValueString(bh.typeDesc)
	  <<" "
	  <<StringUtil::formArrayString(bh.ids);
	 return os;
}


bool ConservedBlock::overlaps(int start, int end) {
	return (start<=getPrimaryStart() && end>=getPrimaryStart()) || (start
			<=getPrimaryEnd() && end>=getPrimaryEnd()) || (start
			>=getPrimaryStart() && start<=getPrimaryEnd()) || (end
			>=getPrimaryStart() && end<=getPrimaryEnd());
}
void ConservedBlock::init(istream& fin) {
	char buf[100];
	 int pos=0;
	 
	 //skip header if there is one
	 while(!fin.eof()){
	 
	 fin.getline(buf,100);
	 if(buf[0]=='#')
		 pos=fin.tellg();
	 else
		 break;
	 //fin.getline(buf,100);	
	 }
	 
	 fin.seekg(pos,ios::beg);
	 int c;
	 fin>>c;
	 if(c!=0)
		 cerr<<"Error reading header of conserved block"<<endl;
	 
	 //string c;
	 //fin>>c;
	 //cerr<<c;
	 //fin>>c;
	 //cerr<<c;
	 fin.seekg(pos,ios::beg);
}
int ConservedBlock::getAlignmentNumber() const {
	return StringUtil::atoi(value[0]);
}
const string& ConservedBlock::getPrimaryChr() const {
	return value[1];
}
int ConservedBlock::getPrimaryStart() const {
	return StringUtil::atoi(value[2]);
}
int ConservedBlock::getPrimaryEnd() const {
	return StringUtil::atoi(value[3]);
}
const string& ConservedBlock::getSecondaryChr() const {
	return value[4];
}
int ConservedBlock::getSecondaryStart() const {
	return StringUtil::atoi(value[5]);
}

int ConservedBlock::getSecondaryEnd() const {
	return StringUtil::atoi(value[6]);
}
char ConservedBlock::getStrand() const {
	return value[7][0];
}
int ConservedBlock::getBLASTZScore() const {
	return StringUtil::atoi(value[8]);
}
const string& ConservedBlock::getPrimarySeq() const {
	return value[9];
}
const string& ConservedBlock::getSecondarySeq() const {
	return value[10];
}


void ConservedBlockQueue::resetIterator() {
	it=vecs.begin();
}
ConservedBlockQueue::ConservedBlockQueue(istream& _fin, ostream& _fout) :
	fin(_fin), fout(_fout) {
	if(!fin.good())
	{
		it=vecs.end();
		cerr<<"cannot open conserved block file"<<endl;
		return;
	}
	cb.init(_fin);
	
	while (!fin.eof()) {
		if (fin.eof())
			break;
		fin>>cb;
		//cerr<<"finish1"<<endl;
		vecs.push_back(cb);
	}
	cerr<<"done reading cb"<<endl;
	resetIterator();
}

void ConservedBlockQueue::output(int gstart, int gend) //chrom assumed same, gend-inclusive
{
	//fout<<">>>>>>"<<gstart<<","<<gend<<endl;
	
	//cerr<<">>>>>>"<<gstart<<","<<gend<<endl;
	//leftsearch
	if(vecs.size()<1)
		return;
	
	if(it==vecs.end())
		return;
	
	vector<ConservedBlock>::iterator i;
	vector<ConservedBlock>::iterator newpos=it;
	vector<ConservedBlock>::iterator leftmost=vecs.end();
	vector<ConservedBlock>::iterator rightmost=vecs.end();

	i=it;
	//bool productive=false;
	bool prodl=false;
	bool prodr=false;

	if ((*i).getPrimaryEnd()>gstart) {
		//cerr<<"::"<<(*i).getPrimaryEnd()<<","<<gstart<<endl;
		while (!(*i).overlaps(gstart, gend) && (*i).getPrimaryEnd()>gstart && i!=vecs.begin()) {
			i--;
			//cerr<<"L"<<(*i).getPrimaryStart()<<","<<(*i).getPrimaryEnd()<<endl;
		}
		while ((*i).overlaps(gstart, gend)) {
			prodl=true;
			//fout<<(*i);
			//cerr<<">>"<<(*i).getPrimaryStart()<<","<<(*i).getPrimaryEnd()<<endl;
			if (rightmost==vecs.end()) {
				rightmost=i;
			}

			leftmost=i;

			if (i==vecs.begin()) {
				break;

			}

			i--;
		}

	}

	if (prodl)
		newpos=(i==vecs.end()) ? i-1 : i;

	i=it;
	if ((*i).getPrimaryStart()<gend) {
		while ( i!=vecs.end() && !(*i).overlaps(gstart, gend) && (*i).getPrimaryStart()<gend ) {
		//	cerr<<"R";
			i++;
		}

		if (i!=vecs.end()) {

			while ((*i).overlaps(gstart, gend)) {
				prodr=true;
				//fout<<(*i);
			//	cerr<<">>"<<(*i).getPrimaryStart()<<","<<(*i).getPrimaryEnd()<<endl;
				if (leftmost==vecs.end()) {
					leftmost=i;
				}

				rightmost=i;

				i++;
				if (i==vecs.end()) {
					break;
				}

			}
		}
	}

	if (prodr) {
		newpos=(i==vecs.end()) ? i-1 : i;
	}
	if (prodr || prodl) {
		//cerr<<"cbq Conserved Blocks found"<<endl;
		for (i=leftmost; i!=rightmost+1; i++) {
			int rstart, rend;
			int estart, eend;
			estart=(*i).getPrimaryStart();
			eend=(*i).getPrimaryEnd();

			int wrstart=estart;
			int wrend=eend;

			if (gstart>estart) {
				rstart=gstart-estart; //read start rel to seq from 0 ungapped
				wrstart=gstart;
			} else {
				rstart=0;
			}

			if (gend<eend) {
				rend=gend-estart;//read end rel to seq from 0 ungapped
				wrend=gend;

			} else {

				rend=eend-estart;
			}

			const string& pseq=(*i).getPrimarySeq();
			const string& sseq=(*i).getSecondarySeq();
			int rpos[2];
			int apos[2];

			rpos[0]=rstart;
			rpos[1]=rend;
::			getActualPosGappedSeq(S(pseq),rpos,apos,2);

			int tstart=apos[0];
			int tend=apos[1];
			int tlen=tend-tstart+1; //inclusive

			string subpseq=pseq.substr(tstart,tlen);
			string subsseq=sseq.substr(tstart,tlen);

			//fout<<(*i);
			fout<<=(*i);
			fout<<"\t";
			fout<<wrstart<<"\t"<<wrend<<endl;

			fout<<subpseq<<endl;
			fout<<subsseq<<endl;

			//cerr<<"cbq>>"<<(*i).getPrimaryStart()<<","<<(*i).getPrimaryEnd()<<endl;

		}
	}

	it=newpos;

}





ostream& operator <<= (ostream& os,const ConservedBlock& cb)
{
	os<<">> ";
		for(int i=1;i<=8;i++)
			os<<cb.value[i]<<" ";
	
		return os;
		//os<<endl;
}

istream & operator >> (istream& is,ConservedBlock& cb)
{
	/*string dummy;
	cb.value[0]="#";
	while(cb.value[0][0]=='#')
	{
		is>>cb.value[0];
	}*/
	for(int i=0;i<11;i++)
	{
		is>>cb.value[i];
		if(is.eof())
			return is;
	}
	//is>>dummy;
	return is;
}


ostream & operator << (ostream& os,const ConservedBlock& cb)
{
	os<<">> ";
	for(int i=1;i<=8;i++)
		os<<cb.value[i]<<" ";
	
	os<<endl;
	os<<cb.value[9]<<endl;
	os<<cb.value[10]<<endl;
	return os;
}




