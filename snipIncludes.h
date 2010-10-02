/**
 * \file snipIncludes.h
 * The general includes for the Snip program
 * \author W. Albert Cheng (awcheng@mit.edu)
 */

#ifndef SNIPINCLUDES_H_
#define SNIPINCLUDES_H_

#include <string>
using namespace std;

class Splidar_OpFlag
{
	public:
	
	bool getSequence;
	bool getCount;
	bool egStringOutputCommonFlankingAsCoord;
	bool commonFlankingCobound;
	string excelHyperLinkPrefix;
	string excelHyperLinkSuffix;
	int seqGetI5Len;
	int seqGetI3Len;
	int seqGetE5Len;
	int seqGetE3Len;
	string eventTypeName;
	bool avoidIsoBound;
	inline Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len,bool _avoidIsoBound=true):
		eventTypeName(_eventTypeName),getSequence(_getSequence),getCount(_getCount),excelHyperLinkPrefix(_excelHyperLinkPrefix),excelHyperLinkSuffix(_excelHyperLinkSuffix),egStringOutputCommonFlankingAsCoord(_egStringOutputCommonFlankingAsCoord),commonFlankingCobound(_commonFlankingCobound),
		 seqGetI5Len(_seqGetI5Len), seqGetI3Len(_seqGetI3Len), seqGetE5Len(_seqGetE5Len), seqGetE3Len(_seqGetE3Len),avoidIsoBound(_avoidIsoBound){}
	const static Splidar_OpFlag GET_ONLY_SEQUENCE;
	const static Splidar_OpFlag GET_ONLY_COUNT;
	const static Splidar_OpFlag GET_BOTH_COUNT_AND_SEQUENCE;
};

#include "Debugger.h"

#include <iostream>
#include <fstream>

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <list>
#include <queue>
#include "Commonf.h"
#include "StringUtil.h"
#include "GffEntry.h"
#include "NucleicBit.h"
#include "FastaFile.h"
#include "KeyedPosition.h"
#include "RandomAccessFile.h"
#include "ConservedBlock.h"
#include "SpliceGraph.h"
#include "SplidarGraph.h"
#include "snip.h"

#include "Vector2Array.h"
#include "SuperConfigFile.h"
#include "PyConfigFile.h"

//#include "SpliceMMGraph.h"

/**
 * \def S(x)
 * conversion of std::string to const char*
 */
#define S(x) (x).c_str()

//////#define BUFFER_SIZE 1000000  this may be the problem!

#define THRESHOLD_TRIM_BUFFER_TIME 10
#define THRESHOLD_TRIM_BUFFER_COMP 10240

#endif /*SNIPINCLUDES_H_*/
