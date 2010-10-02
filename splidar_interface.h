/*
 * splidar_interface.h
 *
 *  Created on: Jan 22, 2010
 *      Author: albertwcheng
 */

#ifndef SPLIDAR_INTERFACE_H_
#define SPLIDAR_INTERFACE_H_

#include "AEP.h"
#include "A53SS.h"
#include "ATE.h"
#include "RI.h"
#include "A3UTR.h"
#include "A3UTRMISO.h"
#include "Expression.h"
#include "statEI.h"

class SplidarSubProgram{
public:
	virtual void operator()( PyConfigFile& params)=0;
	virtual ~SplidarSubProgram(){}
	struct splidar_exception{
		string msg;
		inline splidar_exception(string _msg):msg(_msg){}
	};
};

class SSP_Echo:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{
		try{
		cerr<<params.getString("msg")<<endl;
		}catch(PyConfigFile::key_not_found & knf)
		{
			cerr<<"key not found. please supply msg to echo"<<endl;
		}
	}
};

class SSP_Add:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{
		cerr<<(params.getLongDirect("x")+params.getLongDirect("y"))<<endl;
	}
};

class SplidarSubProgramManager{
public:
	map<string,SplidarSubProgram*> programs;
	void registerProgram(string name, SplidarSubProgram* _program)
	{
		programs.insert(map<string,SplidarSubProgram*>::value_type(name,_program));
	}

	bool call(string name, PyConfigFile& params)
	{
		//params.listAllItemsAsSuperString(cerr);

		map<string,SplidarSubProgram*>::iterator i=programs.find(name);
		if(i==programs.end())
		{
			cerr<<"program "<<name<<" not found"<<endl;
			cerr<<"Available Programs:"<<endl;
			for(i=programs.begin();i!=programs.end();i++)
			{
				cerr<<i->first<<endl;
			}
			return false;
		}

		SplidarSubProgram& programToCall=*(i->second);
		//cerr<<"a"<<endl;
		cerr<<"calling program"<<endl;
		programToCall(params);
		cerr<<"program exited"<<endl;
		//cerr<<"b"<<endl;
		return true;
	}
	bool call( PyConfigFile& params)
	{
		string name;
		try{
			 name=params.getString("program");
			 bool returnvalue=call(name,params);
			 return returnvalue;
		}catch(PyConfigFile::key_not_found& knf)
		{
			cerr<<"program name not supplied. Abort"<<endl;
			return false;
		}/*catch(PyConfigFile::python_error& python_error)
		{
			cerr<<"python error occurred. Abort:"<<python_error.errlog<<endl;
			return false;
		}*/




		return false;
	}
	bool call(string cmdfile)
	{

		PyConfigFile conf(cmdfile);

		return call(conf);
	}
	bool call(int argc,const char**argv)
	{
		PyConfigFile conf(argc,argv);
		return call(conf);
	}

};



//this is the new template!
class SSP_SE_AEP:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	//string mapviewPrefix=params.readSuperString("mapviewPrefix");

	//string outputPrefix=params.readSuperString("outputPrefix");


	string jnxInfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//AEP Specific



	int minDistanceA=params.readInt("minDistanceA");
	int maxDistanceA=params.readInt("maxDistanceA");
	int minDistanceB=params.readInt("minDistanceB");
	int maxDistanceB=params.readInt("maxDistanceB");







	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_calAEPGeneric
		( op,
		 chromLabels,
		  minDistanceA,
		  maxDistanceA,
		  minDistanceB,
		  maxDistanceB,
		 &genomeReadFiles,
		 &jnxReadFiles,
		 readMapFormatCode,
		  chrDir,
		 &outputCountFiles,
		 &outputSeqFiles,
		  transcriptome,
		  jnxInfoFile,
		  minSpan,
		  readLength,
		  useUniquelyMappablePosition,
		  uniquelyMappableChrRef,
		  uniquelyMappableCompactPositionFilePrefix,
		  GffFileName,
		  GBlockFile,
		  lociFile);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};


//this is the new template!
class SSP_A53SS:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	//string mapviewPrefix=params.readSuperString("mapviewPrefix");

	//string outputPrefix=params.readSuperString("outputPrefix");


	string jnxInfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//AEP Specific










	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_calA53SSGeneric
		( op,
		 chromLabels,
		 &genomeReadFiles,
		 &jnxReadFiles,
		 readMapFormatCode,
		  chrDir,
		 &outputCountFiles,
		 &outputSeqFiles,
		  transcriptome,
		  jnxInfoFile,
		  minSpan,
		  readLength,
		  useUniquelyMappablePosition,
		  uniquelyMappableChrRef,
		  uniquelyMappableCompactPositionFilePrefix,
		  GffFileName,
		  GBlockFile,
		  lociFile);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};


class SSP_ATE:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	//string mapviewPrefix=params.readSuperString("mapviewPrefix");

	//string outputPrefix=params.readSuperString("outputPrefix");


	string jnxInfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//AEP Specific










	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_calATEGeneric
		( op,
		 chromLabels,
		 &genomeReadFiles,
		 &jnxReadFiles,
		 readMapFormatCode,
		  chrDir,
		 &outputCountFiles,
		 &outputSeqFiles,
		  transcriptome,
		  jnxInfoFile,
		  minSpan,
		  readLength,
		  useUniquelyMappablePosition,
		  uniquelyMappableChrRef,
		  uniquelyMappableCompactPositionFilePrefix,
		  GffFileName,
		  GBlockFile,
		  lociFile);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};


class SSP_RI:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	//string mapviewPrefix=params.readSuperString("mapviewPrefix");

	//string outputPrefix=params.readSuperString("outputPrefix");


	string jnxInfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//RI Specific
	bool parentsNeedToBeSplicedEitherSide=params.readBool("parentsNeedToBeSplicedEitherSide");
	bool parentsNeedToBeSplicedOnBothSides=params.readBool("parentsNeedToBeSplicedOnBothSides");
	bool parentsNeedToBeSplicedBothSidesInSameTranscript=params.readBool("parentsNeedToBeSplicedBothSidesInSameTranscript");

	bool retainedIntronsHaveToBeFreeOfExons=params.readBool("retainedIntronsHaveToBeFreeOfExons");
	bool useUnambiguousRegionsOfRetainedIntron=params.readBool("useUnambiguousRegionsOfRetainedIntron");







	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_calRIGeneric
		( op,
				parentsNeedToBeSplicedEitherSide,
				parentsNeedToBeSplicedOnBothSides,
				parentsNeedToBeSplicedBothSidesInSameTranscript,
		 chromLabels,
		 &genomeReadFiles,
		 &jnxReadFiles,
		 readMapFormatCode,
		  chrDir,
		 &outputCountFiles,
		 &outputSeqFiles,
		  transcriptome,
		  jnxInfoFile,
		  minSpan,
		  readLength,
		  useUniquelyMappablePosition,
		  uniquelyMappableChrRef,
		  uniquelyMappableCompactPositionFilePrefix,
		  GffFileName,
		  GBlockFile,
		  lociFile,
		  retainedIntronsHaveToBeFreeOfExons,
		  useUnambiguousRegionsOfRetainedIntron
		  );


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_eA3UTR:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	//int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");



	//string jnxInfoFile=params.readSuperString("jnxinfoFile");
	//string GffFileName=params.readSuperString("GffFileName");
	string A3UTRBlockFile=params.readSuperString("A3UTRBlockFile");
	//string lociFile=params.readSuperString("lociFile");
	//string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//e3UTR Specific










	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_expressA3UTRGeneric
			( op,
			chromLabels,
			A3UTRBlockFile,
			&genomeReadFiles,
			readMapFormatCode,
			 chrDir,
			&outputCountFiles,
			&outputSeqFiles,
			 transcriptome,
			 readLength,
			 useUniquelyMappablePosition,
			 uniquelyMappableChrRef,
			 uniquelyMappableCompactPositionFilePrefix
	);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_A3UTRMISO:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");

	//int numChrom=params.readIntWithDefaultValue("numChrom",chromlabels.size());
	//int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");



	//string jnxInfoFile=params.readSuperString("jnxinfoFile");
	//string GffFileName=params.readSuperString("GffFileName");
	string A3UTRBlockFile=params.readSuperString("A3UTRBlockFile");
	//string lociFile=params.readSuperString("lociFile");
	//string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");


	string eventName=params.readSuperString("eventName");

	bool wantSequence=params.readBool("wantSequence");
	bool wantCounts=params.readBool("wantCounts");
	string excelHyperLinkPrefix=params.readSuperString("excelHyperlinkPrefix");
	string excelHyperLinkSuffix=params.readSuperString("excelHyperlinkSuffix");
	int seqGetI3Len=params.readInt("seqGetI3Len");
	int seqGetI5Len=params.readInt("seqGetI5Len");
	int seqGetE5Len=params.readInt("seqGetE5Len");
	int seqGetE3Len=params.readInt("seqGetE3Len");

	//AEP Spec: change to default value in other application
	bool egstringOutputFlankingAsCoord=params.readBool("EGStringOutputFlankingAsCoord");
	bool requiresFlankingCobound=params.readBool("requiresFlankingCobound");

	//Splidar_OpFlag(const string& _eventTypeName,bool _getSequence,bool _getCount,string _excelHyperLinkPrefix,string _excelHyperLinkSuffix,bool _egStringOutputCommonFlankingAsCoord,bool _commonFlankingCobound,int _seqGetI5Len,int _seqGetI3Len,int _seqGetE5Len,int _seqGetE3Len)
	Splidar_OpFlag op(eventName,wantSequence,wantCounts,excelHyperLinkPrefix,excelHyperLinkSuffix,egstringOutputFlankingAsCoord,requiresFlankingCobound,seqGetI5Len,seqGetI3Len,seqGetE5Len,seqGetE3Len,true);

	//map<string,PyObject*> genomeMaps1=params.getStringDictOfPyObjects("genomeMaps");
	//map<string,PyObject*> jnxMaps1=params.getStringDictOfPyObjects("jnxMaps");


	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	params.getStringMapOfStringList("jnxReadFiles",jnxReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");
	outputSeqFiles=params.getStringDictOfStrings("outputSeqFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);


	//e3UTR Specific










	//splidar_calSkippedExon2Generic(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);
	/*splidar_calAEPGeneric
		(Splidar_OpFlag op,
		 const vector<string>& chromLabels,
		 int minDistanceA,
		 int maxDistanceA,
		 int minDistanceB,
		 int maxDistanceB,
		 const map<string, vector<string> >* genomeReadFiles,
		 const map<string, vector<string> >* jnxReadFiles,
		 int readMapFormat,
		 string chrDir,
		 const map<string, string> *outputCountFiles,
		 const map<string, string> *outputSeqFiles,
		 string transcriptome,
		 string JnxInfoFile,
		 int minSpan,
		 int readLength,
		 bool useUniquelyMappablePosition,
		 string uniquelyMappableChrRef,
		 string uniquelyMappableCompactPositionFilePrefix,
		 string GffFileName,
		 string GBlockFile,
		 string lociFile)*/

	splidar_A3UTRMISOGeneric
			( op,
			chromLabels,
			A3UTRBlockFile,
			&genomeReadFiles,
			readMapFormatCode,
			 chrDir,
			&outputCountFiles,
			&outputSeqFiles,
			 transcriptome,
			 readLength,
			 useUniquelyMappablePosition,
			 uniquelyMappableChrRef,
			 uniquelyMappableCompactPositionFilePrefix
	);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};






class SSP_PrepAnnotation:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{
		GffEntry::resetEverything();

			params.listAllItemsAsSuperString(cerr);
			//cerr<<"after"<<endl;


			try
			{
			string transcriptome=params.readSuperString("transcriptome");


			string locusannopath=params.readSuperString("locusannopath");

			string GffFileName=params.readSuperString("GffFileName");
			string GBlockFile=params.readSuperString("GBlockFile");
			string lociFile=params.readSuperString("lociFile");

			string gffLabel=params.readSuperStringWithDefaultValue("gffLabel","");
			bool linkViaName2=params.readBoolWithDefaultValue("linkViaName2",true); //true?
			bool linkViaExon=params.readBoolWithDefaultValue("linkViaExon",false); //false?
			bool linkViaExonGroup=params.readBoolWithDefaultValue("linkViaExonGroup",false); //false?
			bool gffSrcFile=params.readBoolWithDefaultValue("gffSrcFile",true);
			bool getKnownJnxom=params.readBoolWithDefaultValue("getKnownJnxome",false);

			string locusOutput=params.readSuperStringWithDefaultValue("locusOutputFileName",locusannopath+"/"+transcriptome+".n.loci");
			string exonGroupOutput=params.readSuperStringWithDefaultValue("exonGroupOutputFileName",locusannopath+"/"+transcriptome+".n.exongrp");

			string gblockOutput=params.readSuperStringWithDefaultValue("GBlockOutputFilename",locusannopath+"/"+transcriptome+".b");

			cerr<<"computing loci from trxome "<<transcriptome<<endl;
			//getLoci(string gfffile,string snipfile, string siblingVectorFile ,bool gffSrcFile,string gffLabel,bool linkViaName2, bool linkViaExon, bool linkViaExonGroup ,bool writeExonXref,bool writeTranscriptXref,bool resetEverything)

			getLoci(GffFileName,locusOutput,exonGroupOutput,gffSrcFile,gffLabel,linkViaName2,linkViaExon,linkViaExonGroup);

			cerr<<"writing gblocks"<<endl;

			getGBlock(gblockOutput, GffFileName, transcriptome);


			if(getKnownJnxom)
			{
				int readLength=params.readInt("readLength");
				int minSpanForJnx=params.readIntWithDefaultValue("minSpanForJnx",1);
				string outputFastaFile=params.readSuperString("outputJnxomeFastaFile");
				string outputInfoFile=params.readSuperString("outputJnxomeInfoFile");
				string chrDir=params.readSuperString("chrDir");


				cerr<<"getting known jnxome for "<<transcriptome<<endl;

				//getJnxome(35,HUMAN_OTHER_ANNOS+transcriptome+".rg",HUMAN_CHROM,string("/mit/awcheng/")+"hg18."+transcriptome+"2.jnx70.fasta",string("/mit/awcheng/")+"hg18."+transcriptome+"2.jnx70.info",transcriptome+"Jnx");

				/*human EMT readlength=39*/ //getKnownJnxome(39, 1, string("/scratch/jnx6")+"hg18."+transcriptome+".n.kj.jnx70.v6_r39m1.fasta", string("/scratch/jnx6")+"hg18."+transcriptome+".n.kj.jnx70.v6_r39m1.info" , HUMAN_CHROM,  string("source")+transcriptome+".txt", transcriptome);
				/*mouse Bill Wong Red Blood Cell readlength=36*/

				/*void getKnownJnxome(int readLength, int minExonSpan, string snipfile,
						string jnxinfofile, string chrDir, string gfffile, string setName,
						bool gffSrcFile, bool resetEverything)*/

				getKnownJnxome(readLength, minSpanForJnx, outputFastaFile, outputInfoFile , chrDir ,  GffFileName, transcriptome);

			}



			}catch(PyConfigFile::key_not_found& knfException)
			{
				cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
				//die("");
				//throw knfException;
			}catch(splidar_exception& splidarException)
			{
				cerr<<"splidar exception. Abort"<<endl;
				//throw splidarException;
			}



			GffEntry::resetEverything();
			cerr<<"<Done>"<<endl;

	}
};


class SSP_Expression:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");


	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");



	map<string, vector<string> > genomeReadFiles;
	map<string, vector<string> > jnxReadFiles;
	map<string, string> outputCountFiles;
	map<string, string> outputSeqFiles;

	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);
	outputCountFiles=params.getStringDictOfStrings("outputCountFiles");


	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);



	//expression specific
	double constitutiveThreshold=params.getFloatDirect("constitutiveThreshold");
	bool ambigCheck=params.readBool("ambigCheck");
	bool codingCheck=params.readBool("codingCheck");

	double ignoreExonOfLowerPercentile=params.getFloatDirectWithDefaultValue("ignoreExonOfLowerPercentile",-1);
	int ignoreExonOfFromMax=params.readIntWithDefaultValue("ignoreExonOfFromMax",-1);
	int ignoreExonOfFromMin=params.readIntWithDefaultValue("ignoreExonOfFromMin",-1);
	/*void splidar_expressionAnalysisGeneric(
			const vector<string>&chromLabels,
			const map<string, vector<string> >* genomeReadFiles,
			int readMapFormat,
			const map<string,string> *outputCountFiles,
			string transcriptome,
			int minSpan,
			int readLength,
			double constitutiveThreshold,
			bool ambigCheck,
			bool codingCheck,
			bool useUniquelyMappablePosition,
			string uniquelyMappableChrRef,
			string uniquelyMappableCompactPositionFile,
			string GffFileName,
			string GBlockFile,
			string lociFile)*/


	splidar_expressionAnalysisGeneric
		(
		 chromLabels,
		 &genomeReadFiles,
		 readMapFormatCode,
		 &outputCountFiles,
		  transcriptome,
		  minSpan,
		  readLength,
		  constitutiveThreshold,
		  ambigCheck,
		  codingCheck,
		  useUniquelyMappablePosition,
		  uniquelyMappableChrRef,
		  uniquelyMappableCompactPositionFilePrefix,
		  ignoreExonOfLowerPercentile,
		  ignoreExonOfFromMax,
		  ignoreExonOfFromMin,
		  GffFileName,
		  GBlockFile,
		  lociFile);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_StatEI:public SplidarSubProgram{
public:
	void operator()(PyConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);
	//cerr<<"after"<<endl;


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	vector<string> chromLabels=params.readArray("chromlabels");


	int minSpan=params.readInt("minSpan");
	int readLength=params.readInt("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFilePrefix=params.readSuperString("uniquelyMappableCompactPositionFile");

	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");



	map<string, vector<string> > genomeReadFiles;



	params.getStringMapOfStringList("genomeReadFiles",genomeReadFiles);



	string readMapFormat=StringUtil::toUpper(params.readSuperString("readMapFormat")); //MAPVIEW,SAM
	int readMapFormatCode=::getReadMapFormatCode(readMapFormat);





	splidar_calStatEIGeneric(
			sampleName,
			 chromLabels,
			& genomeReadFiles,
			readMapFormatCode,
			 transcriptome,
			 minSpan,
			 readLength,
			 useUniquelyMappablePosition,
			 uniquelyMappableChrRef,
			 uniquelyMappableCompactPositionFilePrefix,
			 GffFileName,
			 GBlockFile,
			 lociFile);


	}catch(PyConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<". Abort"<<endl;
		//die("");
		//throw knfException;
	}catch(splidar_exception& splidarException)
	{
		cerr<<"splidar exception. Abort"<<endl;
		//throw splidarException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};


#endif /* SPLIDAR_INTERFACE_H_ */
