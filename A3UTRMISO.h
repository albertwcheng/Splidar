/*
 * A3UTRMISO.h
 *
 *  Created on: Jul 2, 2010
 *      Author: awcheng
 */

#ifndef A3UTRMISO_H_
#define A3UTRMISO_H_

#include <deque>
#include "BayesPsi3UTR.h"

inline string makemisobinarystring3UTR(int howmany0s,int howmany1s)
{
	string binstr;

	for(int i=0;i<howmany0s;i++)
		binstr+="0";

	for(int i=0;i<howmany1s;i++)
		binstr+="1";

	return binstr;

}

//put it in the splidar_interface.h
inline void splidar_expressA3UTRMISOPerChrom_output(Splidar_OpFlag op,ostream *os, ostream* osSeq, RandomAccessFile* raf, const string& locusName,string& chr,char strand,I2 & w,int utrGroupIndex,int readLength, bool useUniquelyMappablePosition)
{


	deque<string> jnxstringsnips; //field 4
	deque<string> blockboundsnips; //field 7


	vector<KeyPair<int,int> > blockden;

	//now construct those:

	GffEntry::GBlockPtr firstBlock;
	GffEntry::GBlockPtr lastBlock;

//	cerr<<"a"<<endl;

	//int numFragments=w.second-w.first;
	KeyPair<int,int> contriNext(0,0);

	for(I1 block=w.first;block!=w.second;block++)
	{



		GffEntry::GBlockPtr blockPtr=*block;

		lastBlock=blockPtr; //extend this until the last one;

		if(block==w.first)
		{
			firstBlock=blockPtr;
			jnxstringsnips.push_back(chr+":"+StringUtil::str(blockPtr->getStart1())+strand);
		}
		jnxstringsnips.push_back(chr+":"+StringUtil::str(blockPtr->getEnd1())+strand);

		if(strand==GffEntry::FORWARD)
		{
			blockboundsnips.push_back(StringUtil::str(blockPtr->getStart1())+"-"+StringUtil::str(blockPtr->getEnd1()));
		}
		else
		{
			blockboundsnips.push_front(StringUtil::str(blockPtr->getStart1())+"-"+StringUtil::str(blockPtr->getEnd1()));
		}

		if(os)
		{

			KeyPair<int,int> thisDen(0,0);
			KeyPair<int,int> thisDenProtrude(0,0);

			thisDen=blockPtr->getDensity(readLength,useUniquelyMappablePosition,false,true);

			if(w.second-block>1)
				thisDenProtrude=blockPtr->getDensity(readLength,useUniquelyMappablePosition,true,/*ambigCheck*/true);
			else
				thisDenProtrude=thisDen;

			if(strand==GffEntry::FORWARD)
			{
				KeyPair<int,int> thisRealDen(thisDen.k1+contriNext.k1,thisDen.k2+contriNext.k2); //add up len and pos from previous round of protruded region

				blockden.push_back(thisRealDen);

				contriNext.k1=thisDenProtrude.k1-thisDen.k1;
				contriNext.k2=thisDenProtrude.k2-thisDen.k2;
			}
			else
			{

				blockden.push_back(thisDenProtrude);

			}


		}

	}
	//cerr<<"b"<<endl;

	int numFragments=blockboundsnips.size();

	string jnxstring=StringUtil::join<deque<string>, deque<string>::const_iterator>(jnxstringsnips,",");
	string blockbound=StringUtil::join<deque<string>, deque<string>::const_iterator>(blockboundsnips,"/");

	if(os)
	{

		deque<string> misostringsnips; //field 10
		deque<string> SMMPsi;
		deque<string> blockdensnips;

		int numBlocksInGroup=numFragments;//blockden.size();

		KeyPair<int,int> firstCoreDen(0,0);

		KeyPair<int,int> cumDen(0,0);
		KeyPair<int,int> prevDen(0,0);

		bool NotANumber=false;

		double current1Psi=1.0;



		for(int bi=0;bi<numBlocksInGroup;bi++)
		{


			string misobinarystring=makemisobinarystring3UTR(bi,numBlocksInGroup-bi);  //1111  0111 0011 0001

			KeyPair<int,int> thisDen;

			if(strand==GffEntry::FORWARD)
			{
				thisDen=blockden[bi];

			}
			else
			{
				thisDen=blockden[numBlocksInGroup-bi-1]; //from end
			}

			blockdensnips.push_back(StringUtil::str(thisDen.k2)+":"+StringUtil::str(thisDen.k1));

			prevDen=thisDen;
			cumDen.k1+=thisDen.k1;
			cumDen.k2+=thisDen.k2;

			misostringsnips.push_back(misobinarystring+":"+StringUtil::str(cumDen.k2)+":"+StringUtil::str(cumDen.k1));

			//::getBayesPsi_3UTR(extendedDen.k2,extendedDen.k1,coreDen.k2,coreDen.k1);
			if(bi==0)
			{
				firstCoreDen=thisDen;
			}
			else //(bi>0)
			{
				if(NotANumber)
				{
					SMMPsi.push_back("nan");
				}
				else
				{
					double bayesPsi=::getBayesPsi_3UTR(thisDen.k2,thisDen.k1,prevDen.k2,prevDen.k1);
					if(bi==1) //the first pair of blocks
					{
						if(bayesPsi<0)
						{
							SMMPsi.push_back("nan");
							//SMMPsi.push_back("nan");
							NotANumber=true;
						}
						else
						{
							SMMPsi.push_back(StringUtil::str(1.0-bayesPsi));
							//SMMPsi.push_back(StringUtil::str(bayesPsi));
							current1Psi*=bayesPsi;
						}
					}
					else
					{
						if(bayesPsi<0)
						{
							//SMMPsi.push_back("nan");
							//the last one gives the current1Psi (assuming 1.0)
							SMMPsi.push_back(StringUtil::str(current1Psi));

							NotANumber=true;
						}
						else
						{
							double curExitPsi=current1Psi*(1.0-bayesPsi);
							SMMPsi.push_back(StringUtil::str(curExitPsi));
							current1Psi*=bayesPsi;
						}
					}

				}
			}






		}
		//cerr<<"c"<<endl;

		//the end
		if(NotANumber)
		{
			SMMPsi.push_back("nan");
		}
		else
		{
			SMMPsi.push_back(StringUtil::str(current1Psi));
		}
		//cerr<<"d"<<endl;

			//cerr<<"o4"<<endl;
		(*os)<<"Alt3UTR"<<"\t";
		(*os)<<locusName<<"\t";
		(*os)<<locusName<<":"<<utrGroupIndex<<"\t"; //egstring
		(*os)<<jnxstring<<"\t";
		(*os)<<chr<<"\t";
		(*os)<<strand<<"\t";
		(*os)<<numFragments<<"\t";
		(*os)<<blockbound<<"\t";
		(*os)<<firstBlock->getStart1()<<"\t";
		(*os)<<lastBlock->getEnd1()<<"\t";
		(*os)<<op.excelHyperLinkPrefix<<chr<<":"<<firstBlock->getStart1()<<"-"<<lastBlock->getEnd1()<<op.excelHyperLinkSuffix<<"\t";
		(*os)<<firstCoreDen.k2<<"\t"<<firstCoreDen.k1<<"\t";
		(*os)<<cumDen.k2<<"\t"<<cumDen.k1<<"\t";

		(*os)<<StringUtil::join<deque<string>, deque<string>::const_iterator>(blockdensnips,";")<<"\t";
		(*os)<<StringUtil::join<deque<string>, deque<string>::const_iterator>(misostringsnips,";")<<"\t";
		(*os)<<StringUtil::join<deque<string>, deque<string>::const_iterator>(SMMPsi,";");
		(*os)<<endl;

	}

	//cerr<<"e"<<endl;


	if(osSeq && raf)
	{
		(*osSeq)<<"Alt3UTR"<<"\t";
		(*osSeq)<<locusName<<"\t";
		(*osSeq)<<locusName<<":"<<utrGroupIndex<<"\t"; //egstring
		(*osSeq)<<jnxstring<<"\t";
		(*osSeq)<<chr<<"\t";
		(*osSeq)<<strand<<"\t";
		(*osSeq)<<numFragments<<"\t";
		(*osSeq)<<blockbound<<"\t";
		(*osSeq)<<firstBlock->getStart1()<<"\t";
		(*osSeq)<<lastBlock->getEnd1()<<"\t";
		(*osSeq)<<op.excelHyperLinkPrefix<<chr<<":"<<firstBlock->getStart1()<<"-"<<lastBlock->getEnd1()<<op.excelHyperLinkSuffix<<"\t";



		deque<string> sequences;

	//	cerr<<"f"<<endl;



		for(I1 block=w.first;block!=w.second;block++){

			//cerr<<"fget "<<(*block)->getStart0()<<","<<(*block)->getEnd1()<<endl;
			string gseq=raf->get((*block)->getStart0(),(*block)->getEnd1());

			if(strand==GffEntry::FORWARD)
			{
				sequences.push_back(gseq);
			}
			else
			{
			//	cerr<<"reva"<<endl;
				sequences.push_front(::reverse_complement(gseq));
			//	cerr<<"revb"<<endl;
			}


		}
		//cerr<<"g"<<endl;

		(*osSeq)<<StringUtil::join<deque<string>, deque<string>::const_iterator>(sequences,"|");

		(*osSeq)<<endl;

	}
	//cerr<<"h"<<endl;

	//cerr<<"o8"<<endl;

}

//splidar_calAEPPerChromGeneric(Splidar_OpFlag op,int minDistanceA,int maxDistanceA,int minDistanceB,int maxDistanceB,string chr,RandomAccessFile*fChrSeq,const vector<string>& fSelexa,const vector<string>& fJnxSelexa,int readMapFormat,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

void splidar_A3UTRMISOPerChromGeneric(Splidar_OpFlag op,string chr,const vector<string>& fSelexa,int readMapFormat,ostream* os,ostream* osSeq,RandomAccessFile *fRAF,string transcriptome,int readLength,bool useUniquelyMappablePosition)
{


	GffEntry::Loader loader(NULL);
	if(op.getCount && os) //load the reads for this specific chromosome.
	{
		for(vector<string>::const_iterator fSelexaI=fSelexa.begin();fSelexaI!=fSelexa.end();fSelexaI++)
		{
			cerr<<"loading selexa file "<<*fSelexaI<<endl;
			splidar_loadSelexaMatchesGeneric(*fSelexaI,transcriptome,readMapFormat);
		}
	}

	map<string,vector<GffEntry::GBlockPtr>* > genesUTR;

	map<string,set<GffEntry::GBlockPtr>* >::iterator fgbi=GffEntry::gblocks.find(chr);
	if(fgbi==GffEntry::gblocks.end())
	{
		cerr<<"chr not found!"<<endl;
		return;
	}


	set<GffEntry::GBlockPtr>* gblockChr=fgbi->second;
	cerr<<"finding alternative 3'UTR in chr "<<chr<<endl;

	//make a dictionary of genesUTR  with name-> [ GBlocks... ]
	//the gblocks can contain multiple groups of contaguous blocks

	for(set<GffEntry::GBlockPtr>::iterator j=gblockChr->begin();j!=gblockChr->end();j++)
	{
		GffEntry::GBlockPtr utrblock=(*j);
		map<string,vector<GffEntry::GBlockPtr>* >::iterator k=genesUTR.find(utrblock->name);
		vector<GffEntry::GBlockPtr>* utrs;
		if(k==genesUTR.end())
		{
			utrs=new vector<GffEntry::GBlockPtr>;
			genesUTR.insert(map<string,vector<GffEntry::GBlockPtr>* >::value_type(utrblock->name,utrs));

		}else
			utrs=k->second;

		utrs->push_back(utrblock);


	}



	cerr<<"now start!"<<endl;



	for(map<string,vector<GffEntry::GBlockPtr>* >::iterator i=genesUTR.begin();i!=genesUTR.end();i++)
	{

		//for each gene, the UTR regions of that gene

		string geneName=i->first;
		vector<GffEntry::GBlockPtr> * utrRegions=i->second;
		vector<I2> UTRgroups;

		char strand=GffEntry::BOTH_STRANDS;

		//cerr<<"a"<<endl;

		I2 bound;

		if(utrRegions->size()<2){
			delete utrRegions;
			continue;
		}

		//cerr<<"b"<<endl;

		I1 firstBlock=utrRegions->begin();

		bound.first=firstBlock;
		bound.second=bound.first;
		GffEntry::GBlockPtr prevBlock=(*bound.first);

		strand=prevBlock->strand;
		string chr=prevBlock->chr;



		for(I1 i1=utrRegions->begin();i1!=utrRegions->end();i1++)
		{
			GffEntry::GBlockPtr curBlock=*i1;

			if (curBlock->strand!=strand)
			{
				cerr<<"serious error: strand incompat"<<endl;
			}

			if(curBlock->getStart0()<=prevBlock->getEnd1())
			{
				//extend
				bound.second=i1;

			}
			else
			{
				//contaguous blocks broken

				//if(bound.first!=bound.second) //check not a singleton block //NEW do not check: singleton block accepted
				//{
					//two contaguous blocks qualify a UTR group

					bound.second++;
					UTRgroups.push_back(bound);
				//}

				bound.first=i1;
				bound.second=bound.first;

			}

			prevBlock=curBlock;
		}

		//cerr<<"c"<<endl;

		//for the left over groups
		//if(bound.first!=bound.second) //NEW dot not check singleton block: singleton block accepted
		//{
			//two contaguous blocks qualify a UTR group

			bound.second++;
			UTRgroups.push_back(bound);
		//}

		//cerr<<"d"<<endl;
		//cerr<<"num of UTR="<<utrRegions->size()<<endl;
		//cerr<<"UTRs=";
		//printVR2(utrRegions->begin(),utrRegions->end());

		//now process this gene!

		//for each UTR group of contaguous blocks
		//w is I2, which specifies the range of blocks for the UTR group
		int utrGroupIndex=0;
		for(vector<I2>::iterator w=UTRgroups.begin();w!=UTRgroups.end();w++)
		{
			utrGroupIndex++;

			//splidar_expressA3UTRMISOPerChrom_output(op,os, osSeq, fRAF, geneName, chr, strand,core, extended,utrGroupIndex,yindex,readLength);

			if(strand!=GffEntry::FORWARD && strand!=GffEntry::REVERSE)
			{
				cerr<<"serious error!: strand is neither forward nor reverse!!";
				continue;
			}

			//put in op, os, osSeq, fRAF, geneName, chr, strand, w->first, w->second, utrGroupIndex, readLength
			splidar_expressA3UTRMISOPerChrom_output(op,os, osSeq, fRAF, geneName, chr, strand,*w ,utrGroupIndex,readLength, useUniquelyMappablePosition);


		}

		//cerr<<"e"<<endl;
		//now release memory
		delete utrRegions;
	}

	//clean up
	cerr<<"cleaning up after chr"<<chr<<endl;
	GffEntry::resetGBlocks(chr);
	GffEntry::resetJnxTags(chr);
	cerr<<"done clean up"<<endl;

}


void splidar_A3UTRMISOGeneric
		(Splidar_OpFlag op,
		const vector<string>&chromLabels,
		string UTRblocks,
		const map<string, vector<string> >* genomeReadFiles,
		int readMapFormat,
		string chrDir,
		const map<string, string> *outputCountFiles,
		const map<string, string> *outputSeqFiles,
		string transcriptome,
		int readLength,
		bool useUniquelyMappablePosition,
		string uniquelyMappableChrRef,
		string uniquelyMappableCompactPositionFilePrefix
)
{
	GffEntry::Loader loader(NULL);
	loader.loadGBlocks(UTRblocks);
	cerr<<"finish 3UTR block loading"<<endl;



	//go each UTR block, add to the genesUTR record

	//static map<string, set<GBlockPtr> *> gblocks; //chr -> start-ordered Gblocks Ptr


	//cerr<<"end associating genes with 3'UTR"<<endl;

	//cerr<<"finish loading selexa Matches"<<endl;

	//splidar_loadJnxInfo(jnxinfoFile);
	if(useUniquelyMappablePosition)
	{
		GffEntry::Loader loaderUMP(NULL);
		loaderUMP.loadUniquelyMappablePosition(readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,COMPACTPOSITIONSUFFIX,true,false);
	}


	for(unsigned int i=0;i<chromLabels.size();i++)
	{



		string chr=chromLabels[i];

		cerr<<"cal 3UTR "<<i<<":"<<chr<<endl;
		RandomAccessFile *fRAF=NULL;
		ofstream *foutSeq=NULL;
		ofstream *fout=NULL;

	    map<string,vector<string> >::const_iterator genomeReadFilesPerChromI=genomeReadFiles->find(chr);

		if(genomeReadFilesPerChromI==genomeReadFiles->end())// || genomeReadFilesPerChromI->second.size()==0)
		{
			cerr<<"genome read files for chromosome "<<chr<<" is not specified";
			continue;
		}


		if(op.getSequence)
		{
			fRAF=new RandomAccessFile(chrDir+chr+".seq");

			 map<string,string>::const_iterator outputSeqPerChromI=outputSeqFiles->find(chr);

			if(outputSeqPerChromI==outputSeqFiles->end())
			{
				cerr<<"requiring seq output but filename not specified"<<endl;
				die("");
			}

			string outputSeqFilePerChrom=outputSeqPerChromI->second;
			cerr<<"output seq file to "<<outputSeqFilePerChrom<<endl;

			foutSeq=new ofstream(outputSeqFilePerChrom.c_str());
		}


		if(op.getCount)
		{

			 map<string,string>::const_iterator outputCountPerChromI=outputCountFiles->find(chr);

			if(outputCountPerChromI==outputCountFiles->end())
			{
				cerr<<"requiring count output but filename not specified"<<endl;
				die("");
			}

			string outputCountFilePerChrom=outputCountPerChromI->second;
			cerr<<"output count file to "<<outputCountFilePerChrom<<endl;


			fout=new ofstream(outputCountFilePerChrom.c_str());
			//cerr<<"open for writing:"<<(HUMAN_ALT3UTR_ANALYSIS_OUT+sample+"."+chr+".exinc.xls")<<endl;


		}

		//splidar_expressA3UTRPerChromGeneric(Splidar_OpFlag op,string chr,const vector<string>& selexaFileName,int readMapFormat,ostream* os,ostream* osSeq,RandomAccessFile *fRAF,string transcriptome,int readLength)
		splidar_A3UTRMISOPerChromGeneric(op,chr,genomeReadFilesPerChromI->second,readMapFormat,fout,foutSeq,fRAF,transcriptome,readLength, useUniquelyMappablePosition);

		if(fRAF)
		{
			fRAF->close();
			delete fRAF;
			foutSeq->close();
			delete foutSeq;

		}

		if(fout)
		{

			fout->close();
			delete fout;
		}


		cerr<<"cleaning up after chr"<<chr<<endl;
		GffEntry::resetGBlocks(chr);
		cerr<<"done cleaning up"<<endl;
	}


}

#endif /* A3UTRMISO_H_ */
