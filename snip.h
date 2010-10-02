#ifndef SNIP_H_
#define SNIP_H_




#include "snipIncludes.h"

#include "GffEntry.h"
//#include "SplidarGraph.h"

#define GFF_BUFFER_SIZE 10000000
#define NUL_FILE "/dev/null"

const string CHROMOSOME_NOINFO("CHROMOSOME_NOINFO");
void getUltimateNonExonic(string gfffile,string chrDir,string snipfile,string conservedblockfile,string chrLabel="");
void splitGeneInfo(string infilename,string outfileprefix);
void getTS(int upstream,int downstream,string gfffile,string chrDir,string snipfile,string conservedblockfile,string chrLabel="");

void getKnownJnxome(int readLength, int minExonSpan, string snipfile, string jnxinfofile , string chrDir,  string gfffile="",  string setName="jnxome",bool gffSrcFile=true,bool resetEverything=false);
//(int flank, string gfffile, string chrDir,string snipfile, string jnxinfofile, string setName="jnxome");


void getdeNovoJnxome(int flank, int dist, string gfffile, string chrDir,string snipfile, string jnxinfofile, string setName="dnjnxome", bool requireFuncFusion=false, int minIntronLength=50);
void getUltimateGBlockAnnotation(string gfffile,string snipfile);
void getLoci(string gfffile="",string snipfile=NUL_FILE,string siblingVectorFile=NUL_FILE,bool gffSrcFile=true,string gffLabel="",bool linkViaName2=true, bool linkViaExon=true, bool linkViaExonGroup=true ,bool writeExonXref=true,bool writeTranscriptXref=true,bool resetEverything=false);

void getRecombineJnxome(int flank, string snipfile, string jnxinfofile , string chrDir,  int dist =5,string gfffile="", string locusfile="",  string setName="jnxome",bool gffSrcFile=true,bool resetEverything=false);
void writeExonXRef(string dstFile,string gfffile="", bool gffSrcFile=true, bool resetEverything=false);
void writeTranscriptXRef(string dstFile,string gfffile="", bool gffSrcFile=true, bool resetEverything=false);


void findUniqueFor(string annoSource,string snipfile);

void getGBlock(string snipfile, string gfffile, string setName, bool gffSrcFile=true, bool resetEverything=false);
//calculateConstituteExonFreq(string snipFile,string chr,double ratioRep,int readLength,bool ambigCheck,bool codingCheck,bool useUniquelyMappablePosition)
void calculateConstituteExonFreq(string snip,string chr,double ratioRep,int readLength,bool ambigCheck,bool codingCheck,bool useUniqueMappablePosition,double ignoreExonOfLowerPercentile,int ignoreExonOfFromMax,int ignoreExonOfFromMin);
void  statEI(string chr,ostream& os,bool useUniquelyMappablePositions=true);
void statEIReturnValues(string chr,ostream& os,bool useUniquelyMappablePosition,KeyPair<uint64_t,uint64_t>& statChrExonic,KeyPair<uint64_t,uint64_t>& statChrNonExonic,bool resetCount=true);

void simulateReadExonicReadsHashByChr(string snipfileprefix,string snipfilesuffix,int readLength,string chrDir);
void simulateReadExonicReadsHashByPrefix(string snipfileprefix,string snipfilesuffix,int readLength,string chrDir,int prefixLength=6);


/*void assignExonGroupPerLocus(GffEntry::Locus* locus);

void printExonTree(ostream& os,GffEntry::Locus* locus,string ident="\t");
void printExonTree(ostream& os,NExonGroup* root,string ident="\t",string prefix="");
*/
//void printExonTree(ostream& os,NExonGroup* root,string ident="\t",string prefix="",string chr="",string geneName="");
//void printExonTree(ostream& os,GffEntry::Locus* locus,string ident="\t");
//void assignExonGroupPerLocus(GffEntry::Locus* locus,int jobID,bool validateAssignment=true);
//void calculateSE(string snipFile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);


void calculateAEP(Splidar_OpFlag op,int minDistanceA,int maxDistanceA,int minDistanceB,int maxDistanceB,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn);
void calculateA53SS(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn);
void calculateATE(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn);
void calculateRI(Splidar_OpFlag op,bool parentsNeedToBeSplicedEitherSide,bool parentsNeedToBeSplicedOnBothSides, bool parentsNeedToBeSplicedBothSidesInSameTranscript,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn,bool retainedIntronsHaveToBeFreeOfExons,bool useUnambiguousRegionsOfRetainedIntron);

//void calculateSEn(Splidar_OpFlag op,int minDistance,int maxDistance,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);
//void calculateSE2(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);
//void calculateMXE(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);
//void calculateA53SS(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);
//void calculateRI(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);
//void calculateATE(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);

void calculateSpliceMMGraph(Splidar_OpFlag op,string snipFile,string seqfileOut,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile *rafIn=NULL);


int calculateConstituteExonFreq_inner_numOfReads(GffEntry::GBlock* block,GffEntry::GBlock::DIterator range);

void getConstitutiveExonBed(string chr,double ratioRep,bool ambigCheck);



#endif /*SNIP_H_*/
