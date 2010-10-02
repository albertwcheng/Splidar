/*
 * Expression.h
 *
 *  Created on: Feb 25, 2010
 *      Author: awcheng
 */

#ifndef EXPRESSION_H_
#define EXPRESSION_H_


/*
 *
 * splidar_calA53SSPerChromGeneric(
 * Splidar_OpFlag op,
 * string chr,
 * RandomAccessFile*fChrSeq,
 * const vector<string>& fSelexa,
 * const vector<string>& fJnxSelexa,
 * int readMapFormat,
 * string snipFile,
 * string seqoutfile,
 * string transcriptome,
 * string jnxinfoFile,
 * int minSpan,
 * int readLength,
 * bool useUniquelyMappablePosition,
 * string uniquelyMappableChrRef,
 * string uniquelyMappableCompactPositionFilePrefix,
 * string GffFileName,
 * string GBlockFile,
 * string lociFile)
 *
 *
 *
 *
 */
void splidar_expressionAnalysisPerChromGeneric(
		string chr,
		const vector<string>& fSelexa,
		int readMapFormat,
		string resultOutFile,
		string transcriptome,
		int minSpan,
		int readLength,
		double constitutiveThreshold,
		bool ambigCheck,
		bool codingCheck,
		bool useUniquelyMappablePosition,
		string uniquelyMappableChrRef,
		string uniquelyMappableCompactPositionFile,
		double ignoreExonOfLowerPercentile,
		int ignoreExonOfFromMax,
		int ignoreExonOfFromMin,
		string GffFileName,
		string GBlockFile,
		string lociFile)
{
	//transcriptome,jnxinfoFile,readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile
	//GffFileName
	//GBlockFile
	//lociFile

	splidar_loadGffEntryGeneric(transcriptome,GffFileName); //splidar_loadGffEntryGeneric(string transcriptome,string GffFileName)
	splidar_loadGBlocksGeneric(transcriptome,GBlockFile); //splidar_loadGBlocksGeneric(string transcriptome,string GBlockFile)
	splidar_loadLociGeneric(transcriptome,lociFile); //splidar_loadLociGeneric(string transcriptome,string lociFile)

	if(useUniquelyMappablePosition)
	{
		//splidar_loadUniquelyMappablePositionGeneric(string transcriptome,string jnxinfoFile,int readLength,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile)
		splidar_loadUniquelyMappablePositionGeneric(transcriptome,""/*no need to load jnx*/,minSpan,readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile,true,false);
	}

	for(vector<string>::const_iterator fSelexaI=fSelexa.begin();fSelexaI!=fSelexa.end();fSelexaI++)
	{
		cerr<<"loading selexa file "<<*fSelexaI<<endl;
		splidar_loadSelexaMatchesGeneric(*fSelexaI,transcriptome,readMapFormat);
	}

	cerr<<"finish loading selexa matches"<<endl;

	//splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);
	//now everything is ready.
	//calculateConstituteExonFreq(string snipFile,string chr,double ratioRep,int readLength,bool ambigCheck,bool codingCheck,bool useUniquelyMappablePosition)
	calculateConstituteExonFreq(resultOutFile,chr,constitutiveThreshold,readLength,ambigCheck,codingCheck,useUniquelyMappablePosition, ignoreExonOfLowerPercentile, ignoreExonOfFromMax, ignoreExonOfFromMin);

	cerr<<"cleaning up after chr"<<chr<<endl;
	GffEntry::resetGBlocks(chr);
	cerr<<"done clean up"<<endl;

}

/*
 * Splidar_OpFlag op,
	 const vector<string>& chromLabels,
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
	 string lociFile)
 *
 *
 */


void splidar_expressionAnalysisGeneric(
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
		double ignoreExonOfLowerPercentile,
		int ignoreExonOfFromMax,
		int ignoreExonOfFromMin,
		string GffFileName,
		string GBlockFile,
		string lociFile)
{
	GffEntry::gTotalExonReadsUsed=0;

	int numChrom=chromLabels.size();

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];

		string outputCountFilePerChrom;

		cerr<<"cal expression "<<i<<":"<<chr<<endl;

	    map<string,vector<string> >::const_iterator genomeReadFilesPerChromI=genomeReadFiles->find(chr);

		if(genomeReadFilesPerChromI==genomeReadFiles->end() || genomeReadFilesPerChromI->second.size()==0)
		{
			cerr<<"genome read files for chromosome "<<chr<<" is not specified";
			continue;
		}

		map<string,string>::const_iterator outputCountPerChromI=outputCountFiles->find(chr);

		if(outputCountPerChromI==outputCountFiles->end())
		{
			cerr<<"requiring count output but filename not specified"<<endl;
			die("");
		}

		outputCountFilePerChrom=outputCountPerChromI->second;
		cerr<<"output count file to "<<outputCountFilePerChrom<<endl;


		//splidar_expressionAnalysisPerChromGeneric(string chr,string fSelexa,string resultOutFile,string transcriptome,string jnxinfoFile,int readLength,double constitutiveThreshold,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)
		//outputPrefix+"."+chr+".count"
		splidar_expressionAnalysisPerChromGeneric(chr,genomeReadFilesPerChromI->second, readMapFormat,outputCountFilePerChrom,transcriptome,minSpan,readLength,constitutiveThreshold,ambigCheck,codingCheck,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile, ignoreExonOfLowerPercentile,ignoreExonOfFromMax,
				 ignoreExonOfFromMin,GffFileName, GBlockFile,lociFile);
	}

	cerr<<"Total Exon Reads Used="<<GffEntry::gTotalExonReadsUsed<<endl;
}

#endif /* EXPRESSION_H_ */
