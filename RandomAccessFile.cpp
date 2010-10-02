#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;
#include "StringUtil.h"
#include <sstream>
#include "RandomAccessFile.h"

#include <limits.h>



/*string RandomAccessFile::get(int start, int endEx) {

	if(start<0)
		start=0;

	if(endEx<0)
		endEx=INT_MAX;

	int lenRead=0;
	int lenToRead;

	seekg(start, ios::beg);



	string returna;

	while(start+lenRead<endEx && !eof())
	{

		lenToRead=MIN(buflen,endEx-start-lenRead);

		read(buf, lenToRead);

		buf[lenToRead+1]='\0';
		returna+=buf;

		lenRead+=lenToRead;
	}

	return returna;
}*/

RandomAccessFile::RandomAccessFile(string filename, int bufsize, char *_rafbuf) :
	ifstream(filename.c_str()), rafbuflen(bufsize) {
	if(_rafbuf)
	{
		rafbuf=_rafbuf;
		externBuffer=true;
	}
	else
	{
		externBuffer=false;
		rafbuf=new char[rafbuflen];
	}

}

RandomAccessFile::~RandomAccessFile()
{
	close();
	if(!externBuffer)
		delete[] rafbuf;
}

string RandomAccessFile::get(int start,int endEx)
{
	ostringstream ss;
	this->transfer(ss,start,endEx);
	return ss.str();
}

void RandomAccessFile::transfer(ostream&os,int start,int endEx)
{
	if(start<0)
		start=0;

	if(endEx<0)
		endEx=INT_MAX;

	int lenRead=0;
	int lenToRead;

	seekg(start, ios::beg);


	while(start+lenRead<endEx && !eof())
	{

		lenToRead=MIN(rafbuflen-1,endEx-start-lenRead);

		read(rafbuf, lenToRead);


		int thisRead=gcount();

		rafbuf[thisRead]='\0';

		os<<rafbuf;
		lenRead+=thisRead;
	}

}


