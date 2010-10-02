/*
 * snipMain_backup.h
 *
 *  Created on: Jan 22, 2010
 *      Author: albertwcheng
 */

#ifndef SNIPMAIN_BACKUP_H_
#define SNIPMAIN_BACKUP_H_

//PENDING: make this generic or remove
bool validateJnxomeNaive(string jnxinfoFile,string transcriptome)
{


	//splidar_loadGffEntry(transcriptome);
	//splidar_loadGBlocks(transcriptome);
	//splidar_loadLoci(transcriptome);

	splidar_loadJnxInfo(jnxinfoFile);

	//check everything exon-exon is represented;

	for( map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator i=GffEntry::loci.begin();i!=GffEntry::loci.end();i++)
	{
		multimap<GffEntry::ExonPtr,GffEntry::Locus*>* m2=i->second;
		for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator j=m2->begin();j!=m2->end();j++)
		{
			GffEntry::Locus* locus=j->second;

			for(vector<GffEntry*>::iterator k=locus->transcripts.begin();k!=locus->transcripts.end();k++)
			{
				GffEntry* transcript=*k;
				for(int i=0;i<transcript->exonCount-1;i++)
				{
					GffEntry::ExonPtr leftExon=transcript->exons[i];
					GffEntry::ExonPtr rightExon=transcript->exons[i+1];
					if(!leftExon->outJnxs)
					{
						cerr<<"JNXError: leftExon don't branch out to right Exon: outJnxs==NULL:"<<locus->getFirstName()<<":"<<leftExon->exonID<<":"<<leftExon->getBound()<<" to "<<rightExon->exonID<<":"<<rightExon->getBound()<<endl;
						return false;
					}

					if(leftExon->outJnxs->find(rightExon)==leftExon->outJnxs->end())
					{
						cerr<<"JNXError: leftExon don't connect to right Exon: outJnxs->find(rightExon)==end():"<<locus->getFirstName()<<":"<<leftExon->exonID<<":"<<leftExon->getBound()<<" to "<<rightExon->exonID<<":"<<rightExon->getBound()<<endl;
						return false;
					}
				}
			}
		}
	}

	return true;

}


void splidar_calMutuallyExclusiveExonPerChromGeneric(Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)
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

		splidar_loadJnxSelexaMatchesGeneric(fJnxSelexa);

		cerr<<"finished loading jnx selexa"<<endl;

		splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateMXE(op,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}


void splidar_calMutuallyExclusiveExonGeneric(Splidar_OpFlag op,const string* chromLabels, int numChrom,string mapviewPrefix,string chrDir,string outputPrefix,string transcriptome,string JnxInfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

{


	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		splidar_calMutuallyExclusiveExonPerChromGeneric(op,chr,fRAF,mapviewPrefix+"."+chr+".g.mapview.s.q10",mapviewPrefix+"."+chr+".j.mapview.s.q10",outputPrefix+"."+chr+".exinc.xls",outputPrefix+"."+chr+".exinc.seq.xls",transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}

void splidar_calRetainedIntronPerChromGeneric(Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)
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

		splidar_loadJnxSelexaMatchesGeneric(fJnxSelexa);

		cerr<<"finished loading jnx selexa"<<endl;

		splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateRI(op,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}


void splidar_calRetainedIntronGeneric(Splidar_OpFlag op,const string* chromLabels, int numChrom,string mapviewPrefix,string chrDir,string outputPrefix,string transcriptome,string JnxInfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

{


	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		splidar_calRetainedIntronPerChromGeneric(op,chr,fRAF,mapviewPrefix+"."+chr+".g.mapview.s.q10",mapviewPrefix+"."+chr+".j.mapview.s.q10",outputPrefix+"."+chr+".exinc.xls",outputPrefix+"."+chr+".exinc.seq.xls",transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}

void splidar_calAlternativeTerminalExonPerChromGeneric(Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)
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

		splidar_loadJnxSelexaMatchesGeneric(fJnxSelexa);

		cerr<<"finished loading jnx selexa"<<endl;

		splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateATE(op,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}


void splidar_calAlternativeTerminalExonGeneric(Splidar_OpFlag op,const string* chromLabels, int numChrom,string mapviewPrefix,string chrDir,string outputPrefix,string transcriptome,string JnxInfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

{


	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		splidar_calAlternativeTerminalExonPerChromGeneric(op,chr,fRAF,mapviewPrefix+"."+chr+".g.mapview.s.q10",mapviewPrefix+"."+chr+".j.mapview.s.q10",outputPrefix+"."+chr+".exinc.xls",outputPrefix+"."+chr+".exinc.seq.xls",transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}


void splidar_calAlternative53SpliceSitePerChromGeneric(Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)
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

		splidar_loadJnxSelexaMatchesGeneric(fJnxSelexa);

		cerr<<"finished loading jnx selexa"<<endl;

		splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateA53SS(op,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}


void splidar_calAlternative53SpliceSiteGeneric(Splidar_OpFlag op,const string* chromLabels, int numChrom,string mapviewPrefix,string chrDir,string outputPrefix,string transcriptome,string JnxInfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

{


	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		splidar_calAlternative53SpliceSitePerChromGeneric(op,chr,fRAF,mapviewPrefix+"."+chr+".g.mapview.s.q10",mapviewPrefix+"."+chr+".j.mapview.s.q10",outputPrefix+"."+chr+".exinc.xls",outputPrefix+"."+chr+".exinc.seq.xls",transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}


void splidar_calSkippedExonNPerChrom(Splidar_OpFlag op,int minDistance,int maxDistance,string snipFile,string seqoutfile,string chr,string fSelexa,string fJnxSelexa, string jnxinfoFile,RandomAccessFile*fChrSeq=NULL,bool useUniquelyMappablePosition=true,string transcriptome=DEFAULT_TRANSCRIPTOME,int readLength=DEFAULT_READLENGTH)
{


	splidar_loadGffEntry(transcriptome);
	splidar_loadGBlocks(transcriptome);
	splidar_loadLoci(transcriptome);

//	GffEntry::indexLociByName();

	if(useUniquelyMappablePosition)
		::splidar_loadUniquelyMappablePosition();

	splidar_loadJnxInfo(jnxinfoFile);

	cerr<<"finished loading jnxinfo"<<endl;
	if(op.getCount && snipFile.length()>0)
	{
		cerr<<"load selexa matches"<<endl;
		splidar_loadJnxSelexaMatches(fJnxSelexa);
		cerr<<"finished loading jnx selexa"<<endl;
		splidar_loadSelexaMatches(fSelexa,transcriptome);
		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateSEn(op,minDistance,maxDistance,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}




void splidar_calSkippedExonN(Splidar_OpFlag op,int minDistance,int maxDistance,string sample,string chrDir=HUMAN_CHROM,string JnxInfoFile=DEFAULT_JUNCTIONINFO_FILE,bool useUniquelyMappablePosition=true,string transcriptome=DEFAULT_TRANSCRIPTOME,int readLength=DEFAULT_READLENGTH)

{


	for(int i=1;i<2/*CHROM_NUMBER*/;i++)
	{

		string chr=hg18_chromlabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

	//HUMAN_MAPRR+sample+"."+chr+".g.mapview.s.q10",HUMAN_MAPRR_ANALYSIS_OUT+sample+"."+chr+".count"
		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}
		splidar_calSkippedExonNPerChrom(op,minDistance,maxDistance,HUMAN_SPLICING_ANALYSIS_OUT+sample+"."+chr+".exinc.xls",HUMAN_SPLICING_ANALYSIS_OUT+sample+"."+chr+".exinc.seq.xls",chr,HUMAN_MAPRR+sample+"."+chr+".g.mapview.s.q10",HUMAN_MAPRR+sample+"."+chr+".j.mapview.s.q10",JnxInfoFile,fRAF,useUniquelyMappablePosition,transcriptome);
		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}

void splidar_expressionAnalysisPerChromGeneric(string chr,string fSelexa,string resultOutFile,string transcriptome,int minSpan,int readLength,double constitutiveThreshold,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)
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


	splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);
	//now everything is ready.
	calculateConstituteExonFreq(resultOutFile,chr,constitutiveThreshold,readLength,useUniquelyMappablePosition);

}

//this is generic already
void spliceMMGraph_calSpliceMMPerChrom(Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)
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

		splidar_loadJnxSelexaMatchesGeneric(fJnxSelexa);

		cerr<<"finished loading jnx selexa"<<endl;

		splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

		cerr<<"finish loading selexa matches"<<endl;

	}


	calculateSpliceMMGraph(op,snipFile,seqoutfile,chr,readLength,useUniquelyMappablePosition,fChrSeq);


	//now everything is ready.
	//calculateConstituteExonFreq(resultOutFile,chr,0.5,39);
}


void spliceMMGraph_calSpliceMM(Splidar_OpFlag op,const string* chromLabels, int numChrom,string mapviewPrefix,string chrDir,string outputPrefix,string transcriptome,string JnxInfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

{

	DenObj::printBEDHeader(cout);

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromLabels[i];
		cerr<<"cal Skipped Exon for"<<i<<":"<<chr<<endl;

		RandomAccessFile *fRAF=NULL;
		if(op.getSequence)
		{
			cerr<<"open chr sequence file"<<endl;
			fRAF=new RandomAccessFile(chrDir+chr+".seq");
		}

		//Splidar_OpFlag op,string chr,RandomAccessFile*fChrSeq,string fSelexa,string fJnxSelexa,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile

		spliceMMGraph_calSpliceMMPerChrom(op,chr,fRAF,mapviewPrefix+"."+chr+".g.mapview.s.q10",mapviewPrefix+"."+chr+".j.mapview.s.q10",outputPrefix+"."+chr+".exinc.xls",outputPrefix+"."+chr+".exinc.seq.xls",transcriptome,JnxInfoFile,minSpan,readLength,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,GffFileName,GBlockFile,lociFile);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
		}

	}
}






void splidar_statEIGenericPerChrom(string chr,string fSelexa,string transcriptome,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile, KeyPair<int,int>& totalExonic, KeyPair<int,int>& totalNonExonic)
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


	splidar_loadSelexaMatchesGeneric(fSelexa,transcriptome);

	//now everything is ready.
	KeyPair<int,int> statChrExonic(0,0);
	KeyPair<int,int> statChrNonExonic(0,0);

	//statEIReturnValues(string chr,ostream& os,bool useUniquelyMappablePosition,KeyPair<int,int>& statChrExonic,KeyPair<int,int>& statChrNonExonic,bool resetCount)
	statEIReturnValues(chr,cout, useUniquelyMappablePosition,statChrExonic,statChrNonExonic,true);

	totalExonic.k1+=statChrExonic.k1;
	totalExonic.k2+=statChrExonic.k2;
	totalNonExonic.k1+=statChrNonExonic.k1;
	totalNonExonic.k2+=statChrNonExonic.k2;
}

void splidar_calStatEIGeneric(string sample,const string* chromlabels,int numChrom,string mapviewPrefix,string transcriptome,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)
{
	//string sample="Sam2";
	cout<<"#"<<sample<<endl;
	cout<<"#\tchr\tExonicReads\tExonicLen\t\tNonExonicReads\tNonExonicLen"<<endl;

	KeyPair<int,int> totalExonic(0,0);
	KeyPair<int,int> totalNonExonic(0,0);

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromlabels[i];

		cout<<"#\t";
		splidar_statEIGenericPerChrom(chr,mapviewPrefix+"."+chr+".g.mapview.s.q10",transcriptome, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef,  uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile, totalExonic,totalNonExonic);
	}

	cout<<endl;

	cout<<"totalExonicReads="<<totalExonic.k1<<endl;
	cout<<"totalExonicPos="<<totalExonic.k2<<endl;
	cout<<"totalNonExonicReads="<<totalNonExonic.k1<<endl;
	cout<<"totalNonExonicPos="<<totalNonExonic.k2<<endl;



}

void splidar_expressionAnalysisGeneric(const string* chromlabels,int numChrom,string mapviewPrefix,string outputPrefix,string transcriptome,int minSpan,int readLength,double constitutiveThreshold,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)
{
	GffEntry::gTotalExonReadsUsed=0;

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromlabels[i];
		cerr<<"cal expression "<<i<<":"<<chr<<endl;
		//splidar_expressionAnalysisPerChromGeneric(string chr,string fSelexa,string resultOutFile,string transcriptome,string jnxinfoFile,int readLength,double constitutiveThreshold,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)

		splidar_expressionAnalysisPerChromGeneric(chr,mapviewPrefix+"."+chr+".g.mapview.s.q10",outputPrefix+"."+chr+".count",transcriptome,minSpan,readLength,constitutiveThreshold,useUniquelyMappablePosition,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile, GffFileName, GBlockFile,lociFile);
	}

	cerr<<"Total Exon Reads Used="<<GffEntry::gTotalExonReadsUsed<<endl;
}

void splidar_findExonGroupForLocus(string locusName,string transcriptome)
{


	if(!GffEntry::isGffLoaded())
		die("while finding exon group for locus: Gff not loaded: Abort");

	if(!GffEntry::isGBlocksLoaded())
		die("while finding exon group for locus: Gff not loaded: Abort");

	if(!GffEntry::isLocusLoaded())
		die("while finding exon group for locus: locus not loaded: Abort");
	//splidar_loadGffEntry(transcriptome);
	//splidar_loadGBlocks(transcriptome);
	//splidar_loadLoci(transcriptome);

	GffEntry::indexLociByName();

	cerr<<" indexing done "<<endl;

	GffEntry::Locus* locus=
	GffEntry::findLocusByName(locusName);


	if(!locus)
	{
		cerr<<"locus not found"<<endl;
		return;
	}

	cerr<<"locus "<<locusName<<"\t"<<locus->chr<<":"<<locus->getHeadExon()->getStart1()<<"-"<<locus->getTailExon()->getEnd1()<<endl;
	//string wait;
	//cin>>wait;
	cerr<<" assign Start " <<endl;
	::assignExonGroupPerLocus(locus,-100);

	cerr<<" printing "<<endl;
	::printExonTree(cerr,locus);


}




//pending: make generic or remove
void splidar_findExonGroup_AllLoci(string transcriptome)
{


	if(!GffEntry::isGffLoaded())
		die("while finding exon group for all loci: Gff not loaded: Abort");
	//splidar_loadGffEntry(transcriptome);
	//splidar_loadGBlocks(transcriptome);
	//splidar_loadLoci(transcriptome);
	if(!GffEntry::isGBlocksLoaded())
		die("while finding exon group for all loci: gblocks not loaded. Abort");

	if(!GffEntry::isLocusLoaded())
		die("while finding exon group for all loci: loci not loaded. Abort");


	GffEntry::indexLociByName();

	cerr<<" indexing done "<<endl;


	int jobID=-1;

	for(multimap<string,GffEntry::Locus*>::iterator i=GffEntry::nameIndexedLoci.begin();i!=GffEntry::nameIndexedLoci.end();i++)
	{
		jobID--;
		string locusName=(*i).first;
		GffEntry::Locus* locus= (*i).second;
		cerr<<"locus "<<locusName<<"\t"<<locus->chr<<":"<<locus->getHeadExon()->getStart1()<<"-"<<locus->getTailExon()->getEnd1()<<endl;
		::assignExonGroupPerLocus(locus,jobID);
		::printExonTree(cerr,locus);
	}
	//string wait;
	//cin>>wait;
	//cerr<<" assign Start " <<endl;


	//cerr<<" printing "<<endl;
	//::printExonTree(cerr,locus);


}

void splidar_expressionBedPerChromGeneric(string chr,string transcriptome,double constitutiveThreshold,string GffFileName,string GBlockFile,string lociFile)
{
	//transcriptome,jnxinfoFile,readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile
	//GffFileName
	//GBlockFile
	//lociFile

	splidar_loadGffEntryGeneric(transcriptome,GffFileName); //splidar_loadGffEntryGeneric(string transcriptome,string GffFileName)
	splidar_loadGBlocksGeneric(transcriptome,GBlockFile); //splidar_loadGBlocksGeneric(string transcriptome,string GBlockFile)
	splidar_loadLociGeneric(transcriptome,lociFile); //splidar_loadLociGeneric(string transcriptome,string lociFile)


	//now everything is ready.
	getConstitutiveExonBed( chr, constitutiveThreshold, true);

}

void splidar_expressionBedGeneric(const string* chromlabels,int numChrom,string transcriptome,double constitutiveThreshold,string GffFileName,string GBlockFile,string lociFile)
{
	GffEntry::gTotalExonReadsUsed=0;

	cout<<"track name=expressionBlocks_"<<transcriptome<<" description=\"blocks used for calculating expression using annotation "<<transcriptome <<" and requiring "<<constitutiveThreshold<<" of transcripts passing block\" useScore=0 visibility=full"<<endl;

	for(int i=0;i<numChrom;i++)
	{

		string chr=chromlabels[i];
		cerr<<"out BED: "<<i<<":"<<chr<<endl;
		//splidar_expressionAnalysisPerChromGeneric(string chr,string fSelexa,string resultOutFile,string transcriptome,string jnxinfoFile,int readLength,double constitutiveThreshold,bool useUniquelyMappablePosition,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,string GffFileName,string GBlockFile,string lociFile)

		splidar_expressionBedPerChromGeneric(chr,transcriptome,constitutiveThreshold, GffFileName, GBlockFile,lociFile);
	}

	cerr<<"Total Exon Reads Used="<<GffEntry::gTotalExonReadsUsed<<endl;
}

#endif /* SNIPMAIN_BACKUP_H_ */
