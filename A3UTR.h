/*
 * A3UTR.h
 *
 *  Created on: Jan 21, 2010
 *      Author: albertwcheng
 */

#ifndef A3UTR_H_
#define A3UTR_H_

#include "BayesPsi3UTR.h"

inline void splidar_expressA3UTRPerChrom_output(Splidar_OpFlag op,ostream *os, ostream* osSeq, RandomAccessFile* raf, const string& locusName,string& chr,char strand,I2& core,I2& extended,int utrGroupIndex,int yindex,int readLength)
{




	KeyPair<int,int> boundExtendedRegion((*extended.first)->getStart1(),(*(extended.second-1))->getEnd1());
	//cerr<<"o5"<<endl;
	KeyPair<int,int> boundCoreRegion((*core.first)->getStart1(),(*(core.second-1))->getEnd1());
	//cerr<<"o6"<<endl;

	KeyPair<int,int> extendedUTRBound;
	KeyPair<int,int> coreUTRBound;

	if(strand==GffEntry::FORWARD)
	{
		extendedUTRBound=KeyPair<int,int> (boundCoreRegion.k1,boundExtendedRegion.k2);
		coreUTRBound=boundCoreRegion;
	}
	else
	{
		extendedUTRBound=KeyPair<int,int> (boundExtendedRegion.k1,boundCoreRegion.k2);
		coreUTRBound=boundCoreRegion;

	}
	//cerr<<"o7"<<endl;

	string jnxstring;

	if(strand==GffEntry::FORWARD)
		jnxstring=chr+":"+StringUtil::str(coreUTRBound.k1)+strand+","+chr+":"+StringUtil::str(coreUTRBound.k2)+strand+","+chr+":"+StringUtil::str(extendedUTRBound.k2)+strand;
	else
		jnxstring=chr+":"+StringUtil::str(extendedUTRBound.k1)+strand+","+chr+":"+StringUtil::str(coreUTRBound.k1)+strand+","+chr+":"+StringUtil::str(coreUTRBound.k2)+strand;




	//const char*linkstringpre="=HYPERLINK(http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
	//const char *linkstringsuf=",UCSC Browse)";

	if(os)
	{
		//cerr<<"o1"<<endl;
		KeyPair<int,int> extendedDen=GffEntry::GBlock::getDensityOfContigBlocks<I1,RI1>(extended,readLength);
		//cerr<<"o2"<<endl;
		KeyPair<int,int> coreDen=GffEntry::GBlock::getDensityOfContigBlocks<I1,RI1>(core,readLength);
		//cerr<<"o3"<<endl;




		double bayesPsi;//=::getBayesPsi_3UTR(extendedDen.k2,extendedDen.k1,coreDen.k2,coreDen.k1);
		if(extendedDen.k2==0 || coreDen.k2==0)
		{
			bayesPsi=-1;
		}else
		{
			bayesPsi=::getBayesPsi_3UTR(extendedDen.k2,extendedDen.k1,coreDen.k2,coreDen.k1);
		}

			//cerr<<"o4"<<endl;
		(*os)<<"Alt3UTR"<<"\t";
		(*os)<<locusName<<"\t";
		(*os)<<locusName<<":"<<utrGroupIndex<<":"<<yindex<<"\t"; //egstring
		(*os)<<jnxstring<<"\t";
		(*os)<<chr<<"\t";
		(*os)<<strand<<"\t";
		(*os)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"/"<<(coreUTRBound.k1+1)<<"-"<<coreUTRBound.k2<<"\t";
		(*os)<<"1"<<"\t"; //the cobound flag
		(*os)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"/"<<(coreUTRBound.k1+1)<<"-"<<coreUTRBound.k2<<"\t"; //for the bound exonic slot
		(*os)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"\t";//incSpecString
		(*os)<<""<<"\t"; //excSpecString

		(*os)<<op.excelHyperLinkPrefix<<chr<<":"<<extendedUTRBound.k1<<"-"<<extendedUTRBound.k2<<op.excelHyperLinkSuffix<<"\t";
		(*os)<<extendedDen.k1<<"\t";
		(*os)<<extendedDen.k2<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<coreDen.k1<<"\t";
		(*os)<<coreDen.k2<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"1"<<"\t";
		(*os)<<"1"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<coreDen.k1<<"\t";
		(*os)<<coreDen.k2<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<"0"<<"\t";
		(*os)<<1<<"\t";
		(*os)<<1<<"\t";
		(*os)<<extendedDen.k1<<"\t";
		(*os)<<extendedDen.k2<<"\t";
		(*os)<<coreDen.k1<<"\t";
		(*os)<<coreDen.k2<<"\t";
		(*os)<<(coreDen.k1+extendedDen.k1)<<"\t";
		(*os)<<extendedDen.k1<<"\t";
		(*os)<<coreDen.k1<<"\t";

		if(extendedDen.k1+coreDen.k1==0)
			(*os)<<"nan"<<"\t";
		else
			(*os)<<(extendedDen.k1/double(extendedDen.k1+coreDen.k1))<<"\t";

		if(extendedDen.k2==0)
			(*os)<<"nan"<<"\t";
		else
			(*os)<<(double(extendedDen.k1)/extendedDen.k2)<<"\t";

		if(coreDen.k2==0)
			(*os)<<"nan"<<"\t";
		else
			(*os)<<(double(coreDen.k1)/coreDen.k2)<<"\t";

		if(bayesPsi<0)
			(*os)<<"nan"<<endl;
		else
			(*os)<<bayesPsi<<endl;

	}
	if(osSeq && raf)
	{
		(*osSeq)<<"Alt3UTR"<<"\t";
		(*osSeq)<<locusName<<"\t";
		(*osSeq)<<locusName<<":"<<utrGroupIndex<<":"<<yindex<<"\t"; //egstring
		(*osSeq)<<jnxstring<<"\t";
		(*osSeq)<<chr<<"\t";
		(*osSeq)<<strand<<"\t";
		(*osSeq)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"/"<<(coreUTRBound.k1+1)<<"-"<<coreUTRBound.k2<<"\t";
		(*osSeq)<<"1"<<"\t"; //the cobound flag
		(*osSeq)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"/"<<(coreUTRBound.k1+1)<<"-"<<coreUTRBound.k2<<"\t"; //for the bound exonic slot
		(*osSeq)<<(extendedUTRBound.k1+1)<<"-"<<extendedUTRBound.k2<<"\t";//incSpecString
		(*osSeq)<<""<<"\t"; //excSpecString

		//(*osSeq)<<"=HYPERLINK(http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position="<<chr<<":"<<extendedUTRBound.k1<<"-"<<extendedUTRBound.k2<<",UCSC Browse)"<<"\t";
		(*osSeq)<<op.excelHyperLinkPrefix<<chr<<":"<<extendedUTRBound.k1<<"-"<<extendedUTRBound.k2<<op.excelHyperLinkSuffix<<"\t";

		string gseq=raf->get(extendedUTRBound.k1-1,extendedUTRBound.k2);
		if(strand==GffEntry::FORWARD)
		{
			(*osSeq)<<gseq;
		}
		else
		{
			(*osSeq)<<::reverse_complement(gseq);
		}
		//raf->transfer(*osSeq,extendedUTRBound.k1-1,extendedUTRBound.k2);
		(*osSeq)<<"\t";
		//raf->transfer(*osSeq,coreUTRBound.k1-1,coreUTRBound.k2);
	     gseq=raf->get(coreUTRBound.k1-1,coreUTRBound.k2);
		if(strand==GffEntry::FORWARD)
		{
			(*osSeq)<<gseq;
		}
		else
		{
					(*osSeq)<<::reverse_complement(gseq);
		}

		(*osSeq)<<"\t";
		(*osSeq)<<"\t\t\t\t\t";
		(*osSeq)<<endl;

	}

	//cerr<<"o8"<<endl;

}

//splidar_calAEPPerChromGeneric(Splidar_OpFlag op,int minDistanceA,int maxDistanceA,int minDistanceB,int maxDistanceB,string chr,RandomAccessFile*fChrSeq,const vector<string>& fSelexa,const vector<string>& fJnxSelexa,int readMapFormat,string snipFile,string seqoutfile,string transcriptome, string jnxinfoFile,int minSpan,int readLength,bool useUniquelyMappablePosition,string uniquelyMappableChrRef,string uniquelyMappableCompactPositionFilePrefix,string GffFileName,string GBlockFile,string lociFile)

void splidar_expressA3UTRPerChromGeneric(Splidar_OpFlag op,string chr,const vector<string>& fSelexa,int readMapFormat,ostream* os,ostream* osSeq,RandomAccessFile *fRAF,string transcriptome,int readLength)
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

				if(bound.first!=bound.second) //check not a singleton block
				{
					//two contaguous blocks qualify a UTR group

					bound.second++;
					UTRgroups.push_back(bound);
				}

				bound.first=i1;
				bound.second=bound.first;

			}

			prevBlock=curBlock;
		}

		//cerr<<"c"<<endl;

		//for the left over groups
		if(bound.first!=bound.second)
		{
			//two contaguous blocks qualify a UTR group

			bound.second++;
			UTRgroups.push_back(bound);
		}

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

			if(strand==GffEntry::FORWARD)
			{
				I1 gf;
				I1 x;
				I1 y;
				I1 ef;
				I1 end;
				//cerr<<"now the forward"<<endl;
				gf=w->first;//V.begin();
				end=w->second;//V.end();

				//cerr<<"forward"<<endl;
				//cerr<<"gf:"<<(gf-firstBlock)<<"end:"<<(end-firstBlock)<<endl;

				//// for this particular UTR group of blocks
				////     /first block
				////   [ gf |  |    |   ] end
				////          x--------
				////         ef=x
				////             y----------
				////      c  x/e  y            (yindex=2)
				////      c  x/e  e   y        (yindex=3)
				////      c  x/e  e   e    y   (yindex=4)
				///       c   c  x/e  y        (yindex=3)
				///       c   c  x/e  e    y   (yindex=4)
				///       c   c  c    x/e  y   (yindex=4)

				for(x=gf+1;x!=end;x++)
				{
					ef=x;
					for(y=ef+1;y!=end+1;y++)
					{
						I2 core(gf,x);
						I2 extended(ef,y);

						int yindex=y-firstBlock;

					//	printVR(firstBlock,gf,x,false);
					//	//cerr<<"| ";
					//	printVR(firstBlock,ef,y);
						//printVR(firstBlock,gf,x,false);
						//cerr<<"| ";
						//printVR(firstBlock,ef,y);
						//printVR(gf,x,false);
						//cerr<<"| ";
						//printVR(ef,y);
						splidar_expressA3UTRPerChrom_output(op,os, osSeq, fRAF, geneName, chr, strand,core, extended,utrGroupIndex,yindex,readLength);
					}
				}

			}
			else if(strand==GffEntry::REVERSE)
			{
				RI1 rgf;
				RI1 rx;
				RI1 ry;
				RI1 ref;
				RI1 rend;

				rgf=RI1(w->second);//V.rbegin();
				rend=RI1(w->first);//V.rend();

				//cerr<<"reverse!"<<endl;
				//cerr<<"gf:"<<(I1(rgf.base())-firstBlock)<<"end:"<<(I1(rend.base())-firstBlock)<<endl;

				printVR2(w->first,w->second);
				printVR2(I1(rend.base()),I1(rgf.base()));

				for(rx=rgf+1;rx!=rend;rx++)
				{
					//cerr<<"rx:"<<(I1(rx.base())-firstBlock)<<endl;
					ref=rx;

					for(ry=ref+1;ry!=rend+1;ry++)
					{
						//cerr<<"rx="<<printVR3(rx)<<endl;

						I1 gf(rx.base());
						//cerr<<"gf="<<printVR3(gf)<<endl;

						I1 x(rgf.base());
						I1 ef(ry.base());
						I1 y(ref.base());

						//cerr<<"core="<<printVR3(gf)<<" to "<<printVR3(x-1)<<endl;
						//cerr<<"extended="<<printVR3(ef)<<" to "<<printVR3(y-1)<<endl;
						I2 core(gf,x);
						I2 extended(ef,y);

						int yindex=ef-firstBlock;

						//printVR(firstBlock,gf,x,false);
						//cerr<<"| ";
						//printVR(firstBlock,ef,y);

						splidar_expressA3UTRPerChrom_output(op,os, osSeq, fRAF,geneName, chr, strand,core, extended,utrGroupIndex,yindex,readLength);


						//printVR(rx.base(),rgf.base(),false);
						//cerr<<"| ";
						//printVR(ry.base(),ref.base());
					}
				}

			}
			else
			{
				cerr<<"serious error!: strand is neither forward nor reverse!!";
			}
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


void splidar_expressA3UTRGeneric
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

	GffEntry::Loader loaderUMP(NULL);
	loaderUMP.loadUniquelyMappablePosition(readLength,uniquelyMappableChrRef,uniquelyMappableCompactPositionFilePrefix,COMPACTPOSITIONSUFFIX,true,false);



	for(unsigned int i=0;i<chromLabels.size();i++)
	{



		string chr=chromLabels[i];

		cerr<<"cal 3UTR "<<i<<":"<<chr<<endl;
		RandomAccessFile *fRAF=NULL;
		ofstream *foutSeq=NULL;
		ofstream *fout=NULL;

	    map<string,vector<string> >::const_iterator genomeReadFilesPerChromI=genomeReadFiles->find(chr);

		if(genomeReadFilesPerChromI==genomeReadFiles->end() || genomeReadFilesPerChromI->second.size()==0)
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
		splidar_expressA3UTRPerChromGeneric(op,chr,genomeReadFilesPerChromI->second,readMapFormat,fout,foutSeq,fRAF,transcriptome,readLength);

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

#endif
