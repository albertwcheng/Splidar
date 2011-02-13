/*
 * RI.h
 *
 *  Created on: Feb 4, 2010
 *      Author: awcheng
 */

#ifndef RI_H_
#define RI_H_

void splidar_calRIPerChromGeneric(Splidar_OpFlag op,bool parentsNeedToBeSplicedEitherSide,bool parentsNeedToBeSplicedOnBothSides, bool parentsNeedToBeSplicedBothSidesInSameTranscript,string chr,RandomAccessFile*fChrSeq,const vector<string>& fSelexa,const vector<string>& fJnxSelexa,int readMapFormat,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile,bool retainedIntronsHaveToBeFreeOfExons,bool useUnambiguousRegionsOfRetainedIntron)
{


	splidar_loadGffEntryGeneric(transcriptome,GffFileName); //splidar_loadGffEntryGeneric(string transcriptome,string GffFileName)
	splidar_loadGBlocksGeneric(transcriptome,GBlockFile); //splidar_loadGBlocksGeneric(string transcriptome,string GBlockFile)
	splidar_loadLociGeneric(transcriptome,lociFile); //splidar_loadLociGeneric(string transcriptome,string lociFile)

	if(useUniquelyMappablePosition)
	{
		//splidar_loadUniquelyMappablePositionGeneric(string transcriptome,string jnxinfoFile,int readLength,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile)
		splidar_loadUniquelyMappablePositionGeneric(transcriptome,jnxinfoFile,minSpan,readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,true,true); //toJnxBlock also!!
	}

	//splidar_loadJnxInfo(jnxinfoFile); //this is dor yu and old, even if use, should use the generic version

	cerr<<"finished loading jnxinfo and uniquely mappable position (if applicable)"<<endl;

	if(op.getCount && snipFile.length()>0)
	{
		cerr<<"load selexa matches"<<endl;

		for(vector<string>::const_iterator fJnxSelexaI=fJnxSelexa.begin();fJnxSelexaI!=fJnxSelexa.end();fJnxSelexaI++)
		{
			cerr<<"loading jnx selexa file "<<*fJnxSelexaI<<endl;
			splidar_loadJnxSelexaMatchesGeneric(*fJnxSelexaI,readMapFormat);
		}

		cerr<<"finished loading jnx selexa"<<endl;

		for(vector<string>::const_iterator fSelexaI=fSelexa.begin();fSelexaI!=fSelexa.end();fSelexaI++)
		{
			cerr<<"loading selexa file "<<*fSelexaI<<endl;
			splidar_loadSelexaMatchesGeneric(*fSelexaI,transcriptome,readMapFormat);
		}

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateRI(op, parentsNeedToBeSplicedEitherSide, parentsNeedToBeSplicedOnBothSides,  parentsNeedToBeSplicedBothSidesInSameTranscript,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq, retainedIntronsHaveToBeFreeOfExons, useUnambiguousRegionsOfRetainedIntron);

	//clean up
	cerr<<"cleaning up after chr"<<chr<<endl;
	GffEntry::resetGBlocks(chr);
	GffEntry::resetJnxTags(chr);
	cerr<<"done clean up"<<endl;
}


void splidar_calRIGeneric
	(Splidar_OpFlag op,
	bool parentsNeedToBeSplicedEitherSide,
	bool parentsNeedToBeSplicedOnBothSides,
	bool parentsNeedToBeSplicedBothSidesInSameTranscript,
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
	 string lociFile,
	 bool retainedIntronsHaveToBeFreeOfExons,
	 bool useUnambiguousRegionsOfRetainedIntron)

{
	int numChrom=chromLabels.size();

	//check info here:
	if(!genomeReadFiles || !jnxReadFiles)
	{
		cerr<<"genomeReadFiles and/or jnxReadFiles not specified"<<endl;
		die("");
	}

	if(op.getSequence && !outputSeqFiles)
	{
		cerr<<"requiring sequence output but outputSeqFiles not specified"<<endl;
		die("");
	}

	if(op.getCount && !outputCountFiles)
	{
		cerr<<"requiring count output but outputCountFiles not specified"<<endl;
		die("");
	}

	for(int i=0;i<numChrom;i++)
	{
		string chr=chromLabels[i];



		cerr<<"cal Alt 5'/3' Splice Sites for"<<i<<":"<<chr<<endl;


	    map<string,vector<string> >::const_iterator genomeReadFilesPerChromI=genomeReadFiles->find(chr);

		if(genomeReadFilesPerChromI==genomeReadFiles->end())// || genomeReadFilesPerChromI->second.size()==0)
		{
			cerr<<"genome read files for chromosome "<<chr<<" is not specified";
			continue;
		}

	     map<string,vector<string> >::const_iterator jnxReadFilesPerChromI=jnxReadFiles->find(chr);

		if(jnxReadFilesPerChromI==jnxReadFiles->end())// || jnxReadFilesPerChromI->second.size()==0)
		{
			cerr<<"jnx read files for chromosome "<<chr<<" is not specified";
			continue;

		}

		string outputCountFilePerChrom;
		string outputSeqFilePerChrom;







		RandomAccessFile *fRAF=NULL;

		if(op.getCount)
		{
			 map<string,string>::const_iterator outputCountPerChromI=outputCountFiles->find(chr);

			if(outputCountPerChromI==outputCountFiles->end())
			{
				cerr<<"requiring count output but filename not specified"<<endl;
				die("");
			}

			outputCountFilePerChrom=outputCountPerChromI->second;
			cerr<<"output count file to "<<outputCountFilePerChrom<<endl;

		}

		if(op.getSequence)
		{
			 map<string,string>::const_iterator outputSeqPerChromI=outputSeqFiles->find(chr);

			if(outputSeqPerChromI==outputSeqFiles->end())
			{
				cerr<<"requiring seq output but filename not specified"<<endl;
				die("");
			}

			outputSeqFilePerChrom=outputSeqPerChromI->second;
			cerr<<"output seq file to "<<outputSeqFilePerChrom<<endl;

			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		splidar_calRIPerChromGeneric(op, parentsNeedToBeSplicedEitherSide, parentsNeedToBeSplicedOnBothSides,  parentsNeedToBeSplicedBothSidesInSameTranscript,chr,fRAF,genomeReadFilesPerChromI->second,jnxReadFilesPerChromI->second,readMapFormat,outputCountFilePerChrom,outputSeqFilePerChrom,transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile, retainedIntronsHaveToBeFreeOfExons, useUnambiguousRegionsOfRetainedIntron);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}


#endif /* RI_H_ */
