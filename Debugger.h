#ifndef DEBUGGER_H_
#define DEBUGGER_H_

#define GLOBAL_DEBUG_LEVEL 2 ////////set global debug level!


#include <iostream>
#include <time.h>

using namespace std;

class Timer
{
private:
	clock_t _start;
	clock_t _end;
public:
	inline void start()
	{
		_start=clock();
	}
	inline void end()
	{
		_end=clock();
	}
	inline operator double()
	{
		return ((double)(_end-_start))/CLOCKS_PER_SEC;
	}
};


/*class ControlledStream: public ostream
{
public:
	int level;
	ostream& dst;
	int globalDebugLevel;
	inline ControlledStream(ostream& _dst,int _level, int _globalDebugLevel):dst(_dst),level(_level),globalDebugLevel(_globalDebugLevel){}
};

template<typename T>
ControlledStream& operator << (ControlledStream& os,const T& obj)
{
	if (os.level<=os.globalDebugLevel)
	{
		cerr<<"aaaa";
		os.dst<<obj;
	}
	return os;
}



extern ControlledStream cerr1; //will appear in final build, final private build and ultimate debug build, basically the same as cerr
extern ControlledStream cerr2; //will appear in final private build and ultimate debug build
extern ControlledStream cerr3; //will appear in ultimate debug build, for line markers, e.g.
*/




#define CONTROLLED_STREAM(OSTREAM,OSCHAIN,LEVEL,GLOBALLEVEL) if(LEVEL<=GLOBALLEVEL) OSTREAM<<OSCHAIN;

#define CONTROLLED_STDERR_STREAM(OSCHAIN,LEVEL,GLOBALLEVEL) CONTROLLED_STREAM(cerr,OSCHAIN,LEVEL,GLOBALLEVEL)

#define STDERR1(OSCHAIN,GLOBALLEVEL) CONTROLLED_STDERR_STREAM(OSCHAIN,1,GLOBALLEVEL)
#define STDERR2(OSCHAIN,GLOBALLEVEL) CONTROLLED_STDERR_STREAM(OSCHAIN,2,GLOBALLEVEL)
#define STDERR3(OSCHAIN,GLOBALLEVEL) CONTROLLED_STDERR_STREAM(OSCHAIN,3,GLOBALLEVEL)

#define DEBUG1(OSCHAIN) STDERR1(OSCHAIN,GLOBAL_DEBUG_LEVEL)
#define DEBUG2(OSCHAIN) STDERR2(OSCHAIN,GLOBAL_DEBUG_LEVEL)
#define DEBUG3(OSCHAIN) STDERR3(OSCHAIN,GLOBAL_DEBUG_LEVEL)





#endif /*DEBUGGER_H_*/
