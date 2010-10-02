/*
 * statEI.h
 *
 *  Created on: Mar 3, 2010
 *      Author: awcheng
 */

#ifndef STATEI_H_
#define STATEI_H_

/*
 *
 *
 *
 * void splidar_expressionAnalysisPerChromGeneric(
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

 */
void splidar_statEIGenericPerChrom(
		string chr,
		const vector<string>& fSelexa,
		int readMapFormat,
		string transcriptome,
		int minSpan,
		int readLength,
		bool useUniquelyMappablePosition,
		string uniquelyMappableChrRef,
		string uniquelyMappableCompactPositionFile,
		string GffFileName,string GBlockFile,
		string lociFile,
		KeyPair<uint64_t,uint64_t>& totalExonic,
		KeyPair<uint64_t,uint64_t>& totalNonExonic)
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
		splidar_loadUniquelyMappablePositionGeneric(transcriptome,"",minSpan,readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile,true,false);
	}


	for(vector<string>::const_iterator fSelexaI=fSelexa.begin();fSelexaI!=fSelexa.end();fSelexaI++)
	{
		cerr<<"loading selexa file "<<*fSelexaI<<endl;
		splidar_loadSelexaMatchesGeneric(*fSelexaI,transcriptome,readMapFormat);
	}

	cerr<<"finish loading selexa matches"<<endl;


	//now everything is ready.
	KeyPair<uint64_t,uint64_t> statChrExonic(0,0);
	KeyPair<uint64_t,uint64_t> statChrNonExonic(0,0);

	//statEIReturnValues(string chr,ostream& os,bool useUniquelyMappablePosition,KeyPair<int,int>& statChrExonic,KeyPair<int,int>& statChrNonExonic,bool resetCount)
	statEIReturnValues(chr,cout, useUniquelyMappablePosition,statChrExonic,statChrNonExonic,true);

	totalExonic.k1+=statChrExonic.k1;
	totalExonic.k2+=statChrExonic.k2;
	totalNonExonic.k1+=statChrNonExonic.k1;
	totalNonExonic.k2+=statChrNonExonic.k2;

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

void splidar_calStatEIGeneric(
		string sample,
		const vector<string>& chromLabels,
		const map<string, vector<string> >* genomeReadFiles,
		int readMapFormat,
		string transcriptome,
		int minSpan,
		int readLength,
		bool useUniquelyMappablePosition,
		string uniquelyMappableChrRef,
		string uniquelyMappableCompactPositionFile,
		string GffFileName,
		string GBlockFile,
		string lociFile)
{
	//string sample="Sam2";
	cout<<"#"<<sample<<endl;
	cout<<"#\tchr\tExonicReads\tExonicLen\t\tNonExonicReads\tNonExonicLen"<<endl;

	KeyPair<uint64_t,uint64_t> totalExonic(0,0);
	KeyPair<uint64_t,uint64_t> totalNonExonic(0,0);

	int numChrom=chromLabels.size();

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];

	    map<string,vector<string> >::const_iterator genomeReadFilesPerChromI=genomeReadFiles->find(chr);

		if(genomeReadFilesPerChromI==genomeReadFiles->end() || genomeReadFilesPerChromI->second.size()==0)
		{
			cerr<<"genome read files for chromosome "<<chr<<" is not specified";
			continue;
		}


		cout<<"#\t";
		splidar_statEIGenericPerChrom(chr,genomeReadFilesPerChromI->second,readMapFormat,transcriptome, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef,  uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile, totalExonic,totalNonExonic);
	}

	cout<<endl;

	cout<<"totalExonicReads="<<totalExonic.k1<<endl;
	cout<<"totalExonicPos="<<totalExonic.k2<<endl;
	cout<<"totalNonExonicReads="<<totalNonExonic.k1<<endl;
	cout<<"totalNonExonicPos="<<totalNonExonic.k2<<endl;



}


#endif /* STATEI_H_ */
