#ifndef RANDOMACCESSFILE_H_
#define RANDOMACCESSFILE_H_
#include <ostream>
#include <string>
#include <fstream>

#include <string>
using namespace std;

#define READLENGTH 1024
#define READLENGTHP	1025
#ifndef MIN
#define MIN(x,y) (((x)<(y))?(x):(y))
#endif
class RandomAccessFile: public ifstream
{
public:
	char *rafbuf;
	int rafbuflen;
	bool externBuffer;
	//string get(int start,int endEx);
	void transfer(ostream&os,int start,int endEx);
	string get(int start,int endEx);
	RandomAccessFile(string filename,int bufsize=READLENGTHP, char *_rafbuf=NULL);
	~RandomAccessFile();
};


#endif /*RANDOMACCESSFILE_H_*/
