/**
 * \file snipMain.cpp
 * The entry point of the main program of Snip
 * \author W. Albert Cheng (awcheng@mit.edu)
 */

#include "snipIncludes.h"
#include "SplidarGraph.h"
//#include "BayesPsi3UTR.h"
#include "CollapseWig.h"
#include "scanStatistics.h"



#define NUM_T 10

class Lis:public GffEntry::TranscriptLoadingListener
{
public:
	void transcriptLoaded(GffEntry& entry,string annoSource)
	{
		cerr<<entry.name<<" with ID="<<entry.transcriptID<<" loaded from "<<annoSource<<endl;
	}
};



void splidar_loadGffEntryGeneric(string transcriptome,string GffFileName)
{
	if(GffEntry::isGffLoaded())
		return;
	//string transcriptome="RGuKEA";
	cerr<<"load GffEntry "<<transcriptome<<endl;
	GffEntry::Loader loader(NULL);
	loader.loadGffFiles(GffFileName);

}


void splidar_loadLociGeneric(string transcriptome,string lociFile)
{

	if(GffEntry::isLocusLoaded())
		return;
	
	//splidar_loadGffEntry();
	
	
//	string transcriptome="RGuKEA";
	cerr<<"load loci for "<<transcriptome<<endl;
	GffEntry::Loader loader(NULL);
	loader.loadLoci(lociFile);
	cerr<<"finish loading loci"<<endl;
	//not require inframe
	//getLoci(string("source")+transcriptome+".txt",HUMAN_LOCUS_ANNOS+transcriptome+".loci");

}



void splidar_loadGBlocksGeneric(string transcriptome,string GBlockFile)
{
	if(GffEntry::isGBlocksLoaded())
		return;
	
	//splidar_loadGffEntry();

	
	cerr<<"load gblocks for "<<transcriptome<<endl;
	GffEntry::Loader loader(NULL);
	loader.loadGBlocks(/*HUMAN_LOCUS_ANNOS+transcriptome+".b"*/GBlockFile);
	cerr<<"finish gblocks loading"<<endl;
	
}




void splidar_loadSelexaMatchesGeneric(string filename,string transcriptome,int format)
{
	
	//splidar_loadGffEntry();
	//splidar_loadGBlocks(transcriptome);
	//cerr<<"lsm a"<<endl;

	if(!GffEntry::isGffLoaded())
		cerr<<("warning: while loading selexa matches: gff not loaded");

	//cerr<<"lsm b"<<endl;

	if(!GffEntry::isGBlocksLoaded())
		die("while loading selexa matches: gblocks not loaded: Abort");

	cerr<<"loading selexa reads from "<<filename<<" for transcriptome "<<transcriptome<<endl;
	GffEntry::Loader loader(NULL);
	loader.loadSelexaMatches(filename,format);
	cerr<<"finish loading selexa"<<endl;
}


typedef vector<GffEntry::GBlockPtr>::iterator I1;
typedef vector<GffEntry::GBlockPtr>::reverse_iterator RI1;
typedef pair<I1,I1> I2;



void printVR(I1 firstBlock,I1 i,I1 j,bool printEndl=true)
{
	
	int dead=6;
	int counter=0;
	
	for(I1 x=i; x!=j;x++)
	{
		if(++counter==dead)
				break;
		cerr<<(x-firstBlock)<<" ";
	}
		
	if(printEndl)
		cerr<<endl;
}

void printVR2(I1 i,I1 j,bool printEndl=true)
{
	
	int dead=60000	;
	int counter=0;
	
	for(I1 x=i; x!=j;x++)
	{
		if(++counter==dead)
				break;
		cerr<<(*x)->getStart1()<<"-"<<(*x)->getEnd1()<<" ";
	}
		
	if(printEndl)
		cerr<<endl;
}


ostream& printVR3(I1 i)
{		
	cerr<<(*i)->getStart1()<<"-"<<(*i)->getEnd1();	
	return cerr;
}

ostream& printVR3(RI1 i)
{		
	cerr<<(*i)->getStart1()<<"-"<<(*i)->getEnd1();	
	return cerr;
}


//#include "A3UTR.h"







void splidar_loadJnxSelexaMatchesGeneric(string fJnxSelexa,int format)
{
	//splidar_loadJnxSelexaMatches(fJnxSelexa,format); //the original version is generic enough. so directly call it
	GffEntry::Loader loader(NULL);
	loader.loadJnxSelexaMatches(fJnxSelexa,format);
}



void splidar_loadJnxInfoGeneric(string filename,int minSpan,int readLength)
{
	GffEntry::Loader loader(NULL);
	loader.loadJnxInfo(filename,minSpan,readLength);
}




void splidar_loadUniquelyMappablePositionGeneric(string transcriptome,string jnxinfoFile,int minSpan,int readLength,string uniquelyMappableChrRef, string uniquelyMappableCompactPositionFile,bool toGBlock,bool toJnxBlock)
{
	//splidar_loadGffEntry(transcriptome);
	//splidar_loadGBlocks(transcriptome);
	if(toJnxBlock)
		splidar_loadJnxInfoGeneric(jnxinfoFile,minSpan,readLength);
	
	GffEntry::Loader loader(NULL);
	loader.loadUniquelyMappablePosition(readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFile,COMPACTPOSITIONSUFFIX,toGBlock,toJnxBlock);
}















void printFGHelp(const char* programname)
{
	cerr<<programname<<endl;
	cerr<<"-b fasta output_prefix output_suffix prefixLength outchrRef readLength"<<endl;
	cerr<<"\tGenerete binary files of simulated reads binned by prefix"<<endl;

	cerr<<"-b2 fasta output_prefix output_suffix prefixLength outchrRef readLength prefix2"<<endl;
	cerr<<"\tGenerete binary files of simulated reads binned by prefix (only those with prefix2)"<<endl;
	
	
	cerr<<"-f binaryUnfold binaryFold threshold"<<endl;
	cerr<<"\tFold a file by removing reads redundant for threshold time(s)"<<endl;
	
	cerr<<"-f2 binaryUnfold binaryFold threshold outputformat=[k,p]"<<endl;
	cerr<<"\tFold a file by removing reads redundant for threshold time(s) by std::sort"<<endl;

	cerr<<"-pk binary"<<endl;
	cerr<<"\tPrint the content of a binary KEYEDPOSITION file"<<endl;

	cerr<<"-pp binary"<<endl;
	cerr<<"\tPrint the content of a binary PLAIN POSITION file"<<endl;

	cerr<<"-pc binary"<<endl;
	cerr<<"\tPrint the content of a binary COMPACT POSITION file"<<endl;
	
	cerr<<"-c chrRef outputPrefix outputSuffix filebin formatBin=[k,p]"<<endl;
	cerr<<"\tPartition the simulated reads by chr. k=keyed, p=non-keyed bin file input"<<endl;
	
	cerr<<"-s binPosition sortedOutput"<<endl;
	cerr<<"\tsort the simulated reads by position [encoded as PLAIN POSITION file]"<<endl;
	
	cerr<<"-sc binPosition sortedOutput"<<endl;
	cerr<<"\tsort the simulated reads by position [input as POSITION file, output encoded as COMPACT POSITION file]"<<endl;
	
}

void foldGenomics_generateBinary2(int argc,const char** argv)
{
	if(argc<9)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],StringUtil::atoi(argv[7]));
	encoder.transferFromFastaFile(argv[2],argv[8]);
	
}
void foldGenomics_generateBinary(int argc,const char** argv)
{
	if(argc<8)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersEncoder encoder(argv[3],argv[4],StringUtil::atoi(argv[5]),argv[6],StringUtil::atoi(argv[7]));
	encoder.transferFromFastaFile(argv[2]);

}


void foldGenomics_fold_stdsort(int argc,const char** argv)
{
	if(argc<6)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);
	int numEntries=decoder.numEntriesPending(KEYEDPOSITION);

	KeyedPosition* kps=new KeyedPosition[numEntries];

	ofstream ffileOut(argv[3],ios::binary|ios::out);

	int nEntries=0;
	int fEntries=0;

	int thr=StringUtil::atoi(argv[4]);
	
	
	int i=0;
	
	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		kps[i++]=decoder.getKeyedPosition();

		nEntries++;

	}
	
	cerr<<nEntries<<" of sim reads read"<<endl;
	
	//cerr<<"*** Now fold ***"<<endl;
	//sort entries now

	std::sort(kps,kps+numEntries);

	bool outputPlainPosition=!strcmp(argv[5],"p");

	//now entries happening for a particular number of time and output

	KeyedPosition prePos=kps[0];
	int freq=1;
	
	for(i=1;i<numEntries;i++)
	{
		const KeyedPosition& thisPos=kps[i];

		if(thisPos==prePos)
		{
			freq++;
		}
		else
		{
			if(freq<=thr)
			{
				if(outputPlainPosition)
					ffileOut<<(Position)prePos;
				else
					ffileOut<<prePos;
				fEntries++;
			}

			prePos=thisPos;
			freq=1;
		}
	}

	if(freq<=thr)
	{
		if(outputPlainPosition)
			ffileOut<<(Position)prePos;
		else
			ffileOut<<prePos;
		fEntries++;
	}

	delete[] kps;


	cout<<nEntries<<" folded to "<<fEntries<<" with at most "<<thr<<" ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;

	ffileOut.close();

}

void foldGenomics_fold(int argc,const char** argv)
{
	if(argc<5)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);
	
	map<KeyedPosition,uint> kpf;
	typedef map<KeyedPosition,uint>::iterator I;
	typedef pair<I,bool> IStat;
	typedef map<KeyedPosition,uint>::value_type V;
	
	ofstream ffileOut(argv[3],ios::binary|ios::out);

	int nEntries=0;
	int fEntries=0;

	unsigned int thr=StringUtil::atoi(argv[4]); //changed unsigned-signed

	while(decoder.readEntry())
	{
		//decoder.kpos.printAsText(cerr);
		//cerr<<endl;
		IStat stat=kpf.insert(V(decoder.kpos,1));
		if(!stat.second) //more than once!!
		{
			stat.first->second++;
		}

		nEntries++;

	}

	cerr<<nEntries<<" of sim reads read"<<endl;

	//cerr<<"*** Now fold ***"<<endl;
	for(I i=kpf.begin();i!=kpf.end();i++)
	{
		if(i->second<=thr)
		{
			fEntries++;

			ffileOut<<(i->first);

			//i->first.printAsText(cerr);
			//cerr<<endl;
		}
	}
	cout<<nEntries<<" folded to "<<fEntries<<" with at most "<<thr<<" ocurrence(s)"<<endl;
	cerr<<fEntries<<" of sim reads after fold"<<endl;
	
	ffileOut.close();

}

void foldGenomics_print(int argc,const char** argv,int formatBin)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}

	GenomeNmersDecoder decoder(argv[2]);



	int nEntries=0;

	
	while(decoder.readEntry(formatBin))
	{
		switch(formatBin)
		{

		case KEYEDPOSITION:
			decoder.getKeyedPosition().printAsText(cout);
			break;
		case COMPACTPOSITION:
			decoder.getCompactPosition().printAsText(cout);
			break;
		default:
			decoder.getPosition().printAsText(cout);
		}
		cout<<endl;

		nEntries++;

	}
	
	cerr<<nEntries<<" of entries read and printed"<<endl;

}


void foldGenomics_partitionChr(int argc,const char**argv)
{
	if(argc<7)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionChrPartitioner pcp(argv[2],argv[3],argv[4]);
	pcp.partition(argv[5],(!strcmp(argv[6],"k")?KEYEDPOSITION:POSITION));
}

void foldGenomics_sort(int argc,const char**argv)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionSorter::sort(argv[2],argv[3]);
}

void foldGenomics_sortCompact(int argc,const char**argv)
{
	if(argc<3)
	{
		printFGHelp(argv[0]);
		return;
	}
	
	PositionSorter::sortCompact(argv[2],argv[3]);
}

int foldGenomics_main(int argc,const char** argv)
{
		cerr<<"[Fold Genomics Olympics 3.1b]"<<endl;
		cerr<<"[Built:"<<__DATE__<<" "<<__TIME__<<"]"<<endl;
		if(argc<2 || !strcmp(argv[1],"-h"))
		{
			printFGHelp(argv[0]);
		}else if(!strcmp(argv[1],"-b"))
		{
			foldGenomics_generateBinary(argc,argv);
		}
		else if(!strcmp(argv[1],"-b2"))
		{
				foldGenomics_generateBinary2(argc,argv);
		}
		else if(!strcmp(argv[1],"-f"))
		{
			foldGenomics_fold(argc,argv);
		}else if(!strcmp(argv[1],"-pk"))
		{
			foldGenomics_print(argc,argv,KEYEDPOSITION);
		}else if(!strcmp(argv[1],"-pp"))
		{
			foldGenomics_print(argc,argv,POSITION);
		}else if(!strcmp(argv[1],"-pc"))
		{
			foldGenomics_print(argc,argv,COMPACTPOSITION);
		}
		else if(!strcmp(argv[1],"-c"))
		{
			foldGenomics_partitionChr(argc,argv);
		}else if(!strcmp(argv[1],"-s"))
		{
			foldGenomics_sort(argc,argv);
		}
		else if(!strcmp(argv[1],"-sc"))
		{
			foldGenomics_sortCompact(argc,argv);
		}
		else if(!strcmp(argv[1],"-f2"))
		{
			foldGenomics_fold_stdsort(argc,argv);
		}
		cerr<<"<Done>"<<endl;
		return 0;
}

//PENDING: Make generic or remove
int bill_collapseWig(int argc,const char** argv)
{
	//	CollapseWig(double _scoreToAcceptBlock,double _scoreDensityToAcceptBlock, int _distanceToCollapseBlocks, int _minBlockSize)
	CollapseWig cw(3.0,0); //seed
//	cw.addScoreRange(10,1,0);
//	cw.addScoreRange(40,1,0);
	cw.addScoreRange(70,1,0);
//	cw.addScoreRange(100,1.5,0);
	cw.addScoreRange(130,1.8,0);
	cw.readAndCollapseWigs("/net/coldfact/data/awcheng/Bill/Blood/maps.rr/chromWig/mergedWig/chrWiseMerged.wig");
	cw.readAndCollapseWigs("/net/coldfact/data/awcheng/Bill/Blood/maps.rr/chromWig/mergedWig/chrWiseMerged.wig");
	
	cw.finalize(60);
	cw.writeBed(cout,"BLOODTRANSCRIPTv1.");
	return 1;
}

//PENDING: make generic or remove
int test_ScanStat(int argc,const char ** argv)
{
/*	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(4, 1, 10))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(4, 1, 10))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(4,0.1,10,100))<<endl;
	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(5, 12*0.25, int(1/0.25)))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(5, 12*0.25, int(1/0.25)))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(5,0.12,25,100))<<endl;
	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(5, (float(8)/6), 6))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(5, (float(8)/6), 6))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(5,0.08 ,100.0/6,100))<<endl;
	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(41, 4800*(float(1)/240), 240))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(41, 4800*(float(1)/240), 240))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(41,48 ,100.0/240,100))<<endl;
	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(11, 33*(float(1)/6), 6))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(11, 33*(float(1)/6), 6))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(11,0.33,100.0/6,100))<<endl;
	cerr<<float_approx(SS_Poisson_Prob_Approx3_3(6, 342*(float(1)/674.51), 674.51))<<endl;
	cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(6, 342*(float(1)/674.51), 674.51))<<endl;
	cerr<<(::SS_Poisson_Prob_Alm_Approx(6,3.42,100.0/674.51,100))<<endl;		*/
	//L= genomic uniq pos = 2177505976
	//uniq reads = 14646824
	//(w)indow size = 100
	//lambda (per base) = 14646824/2177505976
	//psi= w = lambda * w
	
/*	for(int k=4;k<=20;k++)
	{
		cerr<<float_approx(::SS_Poisson_Prob_Naus_Approx(k, 1, 10))<<endl;
	}*/
	
	cl_R window=StringUtil::atoi(argv[1]);
	cl_R uniqPos="2177505976";
	cl_R uniqRead="14646824";
	cl_R lambda=uniqRead/uniqPos;
	cl_R psi=lambda*window;
	cl_I L=round1(uniqPos/window);
	cl_R Prob=1.0;
	int k=StringUtil::atoi(argv[2]);
	double thresold=StringUtil::atof(argv[3]);
	if(k==-1)
	{
		k=int(float_approx(psi))+1;
	}
	while(Prob>=thresold)
	{
		k++;
		//cerr<<k<<"\t"<<float_approx(psi)<<"\t"<<L<<"\t";
		cerr<<"k="<<k<<"\t"<<"lambda="<<float_approx(lambda)<<"\t"<<"window="<<window<<"\t"<<"uniqPos="<<uniqPos<<"\t"<<"Prob="<<"\t";
	    Prob=(::___SS_Poisson_Prob_Alm_Approx(k,(lambda),(window),(uniqPos)));//(::SS_Poisson_Prob_Naus_Approx(k,psi,L));//float_approx(::SS_Poisson_Prob_Approx3_3(k,psi,L));//(::SS_Poisson_Prob_Alm_Approx(k,lambda,window,uniqPos));
		//Prob=(::SS_Poisson_Prob_Naus_Approx(k,psi,L));
		cerr<<Prob<<endl;
	}


	return 1;
}



#include "splidar_interface.h"


/*class SSP_MXE:public SplidarSubProgram{
public:
	void operator()(SuperConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	Vector2Array<string> chromlabels=params.readArray("chromlabels");

	int numChrom=params.read<int>("numChrom",chromlabels.size());
	int minSpan=params.read<int>("minSpan");
	int readLength=params.read<int>("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFile=params.readSuperString("uniquelyMappableCompactPositionFile");

	string mapviewPrefix=params.readSuperString("mapviewPrefix");

	string outputPrefix=params.readSuperString("outputPrefix");


	string jnxinfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");

	Splidar_OpFlag op(params.readBool("wantSequence"),params.readBool("wantCounts"));

	splidar_calMutuallyExclusiveExonGeneric(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);



	}catch(ConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<endl;
		throw knfException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_RI:public SplidarSubProgram{
public:
	void operator()(SuperConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	Vector2Array<string> chromlabels=params.readArray("chromlabels");

	int numChrom=params.read<int>("numChrom",chromlabels.size());
	int minSpan=params.read<int>("minSpan");
	int readLength=params.read<int>("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFile=params.readSuperString("uniquelyMappableCompactPositionFile");

	string mapviewPrefix=params.readSuperString("mapviewPrefix");

	string outputPrefix=params.readSuperString("outputPrefix");


	string jnxinfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");

	Splidar_OpFlag op(params.readBool("wantSequence"),params.readBool("wantCounts"));

	splidar_calRetainedIntronGeneric(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);



	}catch(ConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<endl;
		throw knfException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_A53SS:public SplidarSubProgram{
public:
	void operator()(SuperConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	Vector2Array<string> chromlabels=params.readArray("chromlabels");

	int numChrom=params.read<int>("numChrom",chromlabels.size());
	int minSpan=params.read<int>("minSpan");
	int readLength=params.read<int>("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFile=params.readSuperString("uniquelyMappableCompactPositionFile");

	string mapviewPrefix=params.readSuperString("mapviewPrefix");

	string outputPrefix=params.readSuperString("outputPrefix");


	string jnxinfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");

	Splidar_OpFlag op(params.readBool("wantSequence"),params.readBool("wantCounts"));

	splidar_calAlternative53SpliceSiteGeneric(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);



	}catch(ConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<endl;
		throw knfException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};

class SSP_ATE:public SplidarSubProgram{
public:
	void operator()(SuperConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	Vector2Array<string> chromlabels=params.readArray("chromlabels");

	int numChrom=params.read<int>("numChrom",chromlabels.size());
	int minSpan=params.read<int>("minSpan");
	int readLength=params.read<int>("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFile=params.readSuperString("uniquelyMappableCompactPositionFile");

	string mapviewPrefix=params.readSuperString("mapviewPrefix");

	string outputPrefix=params.readSuperString("outputPrefix");


	string jnxinfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");

	Splidar_OpFlag op(params.readBool("wantSequence"),params.readBool("wantCounts"));

	splidar_calAlternativeTerminalExonGeneric(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);



	}catch(ConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<endl;
		throw knfException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};




*/

/*class ListParam_program:public SplidarSubProgram{
public:
	void operator()(PyConfigFile & params)
	{


		params.listAllVariables(cerr);


	}

};*/

/*
class SpliceMMGraph_program:public SplidarSubProgram{
public:
	void operator()(SuperConfigFile& params)
	{

	GffEntry::resetEverything();

	params.listAllItemsAsSuperString(cerr);


	try
	{
	string transcriptome=params.readSuperString("transcriptome");


	string locusannopath=params.readSuperString("locusannopath");
	string sampleName=params.readSuperString("sampleName");
	Vector2Array<string> chromlabels=params.readArray("chromlabels");

	int numChrom=params.read<int>("numChrom",chromlabels.size());
	int minSpan=params.read<int>("minSpan");
	int readLength=params.read<int>("readLength");

	bool useUniquelyMappablePosition=params.readBool("useUniquelyMappablePosition");



	string uniquelyMappableChrRef=params.readSuperString("uniquelyMappableChrRef");
	string uniquelyMappableCompactPositionFile=params.readSuperString("uniquelyMappableCompactPositionFile");

	string mapviewPrefix=params.readSuperString("mapviewPrefix");

	string outputPrefix=params.readSuperString("outputPrefix");


	string jnxinfoFile=params.readSuperString("jnxinfoFile");
	string GffFileName=params.readSuperString("GffFileName");
	string GBlockFile=params.readSuperString("GBlockFile");
	string lociFile=params.readSuperString("lociFile");
	string exonGroupFile=params.readSuperString("exonGroupFile");

	string chrDir=params.readSuperString("chrDir");

	Splidar_OpFlag op(params.readBool("wantSequence"),params.readBool("wantCounts"));

	spliceMMGraph_calSpliceMM(op,chromlabels,  numChrom, mapviewPrefix, chrDir, outputPrefix, transcriptome, jnxinfoFile, minSpan, readLength, useUniquelyMappablePosition, uniquelyMappableChrRef, uniquelyMappableCompactPositionFile, GffFileName, GBlockFile, lociFile);



	}catch(ConfigFile::key_not_found& knfException)
	{
		cerr<<"unknown key "<<knfException.key<<endl;
		throw knfException;
	}



	GffEntry::resetEverything();
	cerr<<"<Done>"<<endl;

	}

};
*/

void printUsage(int argc,const char**argv)
{
		cerr<<"Usage:\t"<<argv[0]<<" commandFile"<<endl;
		cerr<<"\t"<<argv[0]<<"-"<<" -key values..."<<endl;
}



int snip_main(int argc,const char**argv)
{
		cerr<<"[Splidar Vodka 5.1 PyPlus]"<<endl;
		cerr<<"[Built:"<<__DATE__<<" "<<__TIME__<<"]"<<endl;
		cerr<<"called with args:"<<endl;
		for(int iarg=0;iarg<argc;iarg++)
		{
			cerr<<iarg<<"\t"<<argv[iarg]<<endl;
		}
		cerr<<"***********"<<endl;;

		cerr<<"init Python"<<endl;
		PyConfigFile::initPython();
		cerr<<"done init Python"<<endl;
		/*try
		{
		PyConfigFile pcf;
		pcf.loadFromArgs(argc,argv);
		cerr<<"finish loading"<<endl;
		//cerr<<pcf.getString("everythingPrefix")<<endl;
		pcf.listAllVariables(cerr);
		cerr<<"Done"<<endl;
		pcf.finalize();
		}catch(PyConfigFile::file_not_found& e)
		{
			cerr<<"cannot load "<<e.filename<<endl;
		}
		catch(PyConfigFile::key_not_found& e)
		{
			cerr<<"cannot find object "<<e.key<<endl;
		}
		return 1;*/



		//return test_splidar_subprograms();

		//return test_config();

		try{

		SplidarSubProgramManager manager;



		SSP_Echo Echo;
		manager.registerProgram("Echo",&Echo);


		//The Splicing Programs:
		SSP_SE_AEP Splidar_SE_AEP;
		SSP_A53SS Splidar_A53SS;
		//SSP_MXE Splidar_MXE;
		SSP_RI Splidar_RI;
		SSP_ATE Splidar_ATE;
		SSP_eA3UTR Splidar_eA3UTR;
		SSP_A3UTRMISO Splidar_A3UTRMISO;
		SSP_Expression Splidar_Expression;

		SSP_PrepAnnotation Splidar_prepAnnotation;
		SSP_StatEI Splidar_StatEI;


		////SSP_A53SS Splidar_A53SS;
		//SpliceMMGraph_program smm_program;

		manager.registerProgram("AEP",&Splidar_SE_AEP);
		manager.registerProgram("A53SS",&Splidar_A53SS);
		//manager.registerProgram("MXE",&Splidar_MXE);
		manager.registerProgram("RI",&Splidar_RI);
		manager.registerProgram("ATE",&Splidar_ATE);
		manager.registerProgram("A3UTR",&Splidar_eA3UTR);
		manager.registerProgram("A3UTRMISO",&Splidar_A3UTRMISO);
		manager.registerProgram("EXP",&Splidar_Expression);
		manager.registerProgram("PrepAnno",&Splidar_prepAnnotation);
		manager.registerProgram("StatEI",&Splidar_StatEI);

		//manager.registerProgram("A53SS",&Splidar_A53SS);
		//manager.registerProgram("SMM",&smm_program);

		if(argc<2)
		{
			printUsage(argc,argv);
			return 1;
		}

		if(argc==2)
			return manager.call(argv[1]);

		if(argc>2)
			return manager.call(argc,argv);

		}catch(PyConfigFile::file_not_found& fnfException)
		{
			cerr<<"file not found "<<fnfException.filename<<endl;
			//throw fnfException;

		}
		catch(PyConfigFile::key_not_found &knfException)
		{
			cerr<<"key not found thrown: Abort"<<endl;
		}
		catch(PyConfigFile::python_error& python_error)
		{
					cerr<<"python error occurre. Abort:"<<python_error.errlog<<endl;

		}
		catch(SplidarSubProgram::splidar_exception& splidar_e)
		{
			cerr<<"splidar exception thrown: "<<splidar_e.msg<<endl;
			//dier("",1);
		}

		cerr<<"start deinit python"<<endl;
		PyConfigFile::deinit_Python();

		cerr<<"<done done>"<<endl;

		return 1;
}



int main(int argc, const char **argv)
{

	/*int tryLongMax=INT_MAX;
	cerr<<tryLongMax<<endl;
	cerr<<INT_MAX<<endl;

	string str;
	ifstream fin("tryConsumeLine.txt");
	fin>>str;
	cerr<<"["<<str<<"]"<<endl;
	while(!fin.eof())
	{
	str=::consumeStreamUntilNewLine(fin);
	cerr<<"["<<str<<"]"<<endl;
	}
	fin.close();

	return 1;*/

	//return foldGenomics_main(argc,argv);
	return snip_main(argc,argv);
}


