#ifndef CONSERVEDBLOCK_H_
#define CONSERVEDBLOCK_H_

class BlockHeader
{
public:
	
	int chrom;
	int gstart;
	int gend;
	map<string,string> typeDesc;
	vector<string> ids;
	vector<string> names;
	vector<string> names2;
	map<string,string> extData;
	BlockHeader();
	BlockHeader(int _chrom,int _gstart,int _gend);
	void addName(string name);
	void addName2(string name2);
	void addExtData(string key,string value);
	void addTypeDesc(string key,string value);
	
	
};


class ConservedBlock
{
private:
	string value[11];
	

public:
	bool overlaps(int start,int end);
	void init(istream& fin);
	int getAlignmentNumber() const;
	const string& getPrimaryChr() const;
	int getPrimaryStart() const;
	int getPrimaryEnd() const;
	const string& getSecondaryChr() const;
	int getSecondaryStart() const;
	int getSecondaryEnd() const;
	char getStrand() const;
	int getBLASTZScore() const;
	const string& getPrimarySeq() const;
	const string& getSecondarySeq() const;
	friend istream & operator >> (istream& is,ConservedBlock& cb);
	friend ostream & operator << (ostream& os,const ConservedBlock& cb);
	friend ostream & operator <<= (ostream& os,const ConservedBlock& cb);
};


class ConservedBlockQueue
{

public:
	vector<ConservedBlock> vecs;
	ConservedBlock cb;
	ostream &fout;
	istream &fin;
	vector<ConservedBlock>::iterator it;
	void resetIterator();
	ConservedBlockQueue(istream& _fin,ostream& _fout);
	void output(int gstart,int gend);
	
};

istream & operator >> (istream& is,ConservedBlock& cb);
ostream & operator << (ostream& os,const ConservedBlock& cb);

#endif /*CONSERVEDBLOCK_H_*/
