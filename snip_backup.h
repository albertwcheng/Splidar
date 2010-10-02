/*
 * snip_backup.h
 *
 *  Created on: Jan 22, 2010
 *      Author: albertwcheng
 */

#ifndef SNIP_BACKUP_H_
#define SNIP_BACKUP_H_

void getConstitutiveExonBed(string chr,double ratioRep,bool ambigCheck)
{
//	ofstream snip(snipFile.c_str());

	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}


	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}


	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;



	multimap<GffEntry::GBlockPtr,int> blockCountMap;

	/* algorithm get gblock of freq/totalTranscripts > ratioRep:
	foreach locus l
		foreach transcript t in l
			foreach exon e in t
				foreach block b in e
					add b to blockCountMap <block,count>
					and increment count by 1
				rof
			rof
		rof
	rof
	*/

	int lc=0;
	int tc=chrMapLocus->size();

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		vector<GffEntry*>& transcripts=locus->transcripts;
		blockCountMap.clear();

		int rtranscriptCount=transcripts.size();

		////new
		int maxExonCount=0;



		//multimap<int,GffEntry*> orderedTranscripts;

		int *ecounts=new int[rtranscriptCount];
		//set<int> lengthsECount;

		int ei=0;

		for(vector<GffEntry*>::iterator t=transcripts.begin();t!=transcripts.end();t++)
		{
			int thisExonCount=(*t)->exonCount;
			if(thisExonCount>maxExonCount)
				maxExonCount=thisExonCount;

			GffEntry* transcript=(*t);
			ecounts[ei++]=transcript->exonCount;

		}

		int medianIndex=int(ceil((double)rtranscriptCount/2)-1);
		nth_element(ecounts,ecounts+medianIndex,ecounts+rtranscriptCount);


		int LowerThresholdExonCount=ecounts[medianIndex];//MAX(maxExonCount-3,1);

		//if(maxExonCount>1 && LowerThresholdExonCount==1)
			//LowerThresholdExonCount=2;


		delete[] ecounts;

		int transcriptCount=0;

		////new
		cerr<<"processing locus "<<(*locus->names.begin())<<endl;
		cerr<<"lower threshold exoncount="<<LowerThresholdExonCount<<endl;
		for(vector<GffEntry*>::iterator t=transcripts.begin();t!=transcripts.end();t++)
		{



			GffEntry* transcript=(*t);


			//new!! add to the expression code??
			if(transcript->chrom!=locus->chr)
				continue;


			////new
			if(transcript->exonCount<LowerThresholdExonCount)
				continue;

			transcriptCount++;
			////new
			cerr<<"transcript "<<transcript->name<<" of "<<transcriptCount<<" with "<<transcript->exonCount<<" exons"<<endl;

			for(int e=0;e<transcript->exonCount;e++)
			{

				GffEntry::ExonPtr exon=transcript->exons[e];
				if(!exon->blocks)
				{
					cerr<<">>>>>>>>>>>>>>>>>>>no block for exon"<<endl;
					continue;
				}
				for(vector<GffEntry::GBlockPtr>::iterator b=exon->blocks->begin();b!=exon->blocks->end();b++)
				{
					GffEntry::GBlockPtr block=(*b);
					cerr<<"block:"<<block->getStart1()<<"-"<<block->getEnd1()<<"\t";
					if(ambigCheck && block->getNaiveAmbiguity())
					{
						cerr<<"ambiguous"<<endl;
						continue;
					}

					map<GffEntry::GBlockPtr,int>::iterator i=blockCountMap.find(block);
					if(i==blockCountMap.end())
					{
						cerr<<"inserted"<<endl;
						blockCountMap.insert(map<GffEntry::GBlockPtr,int>::value_type(block,1));
					}
					else
					{
						cerr<<"incremented"<<endl;
						(*i).second++;
					}
				}//rof b
			}//rof e
		}//rof t

		//now go through all blocks in the BlockCountMap
		//find contaguous pieces of blocks with count/transcriptTotal>ratioRep
		//get number of possible positions for each piece and the associated read counts
		//add to the records the pos and the reads

		map<GffEntry::GBlockPtr,int>::iterator lMarker;

		string _chrom=locus->chr;
		char _strand=locus->strand;
		string _name;
		int _start=INT_MAX;
		int _end=INT_MIN;
		vector<int> _starts;
		vector<int> _sizes;
		//vector<int> _scores;

		/*for(set<string>::iterator i=locus->names.begin();i!=locus->names.end();i++)
		{
			snip<<(*i)<<",";
		}*/ //changed to below

		set<string>::iterator nameI=locus->names.begin();
		if(nameI!=locus->names.end())
		{
			_name=(*nameI);
		}

		//snip<<"\t";



		//output basic info for block

		//find the first eligible block

		for(lMarker=blockCountMap.begin();lMarker!=blockCountMap.end();lMarker++)
		{
			int count=(*lMarker).second;
			GffEntry::GBlockPtr pbl=(*lMarker).first;
			cerr<<"block "<<pbl->getStart1()<<"-"<<pbl->getEnd1()<<"\t"<<"representing "<<count<<"/"<<transcriptCount<<" of transcript ";
			if((double)count/transcriptCount>=ratioRep)
			{
				cerr<<" accepted"<<endl;
				break;
			}else
			{
				cerr<<" ignored: low representative"<<endl;
			}

		}

		//lMarker reached end
		if(lMarker==blockCountMap.end())
		{
			cerr<<"No Blocks are good for locus "<<(*locus->names.begin())<<endl;
			//snip<<"NBG\t;";
			//snip<<endl;
			continue;
		}

		_start=lMarker->first->start0;
		_end=_start;

		//an imaginative end
		map<GffEntry::GBlockPtr,int>::iterator gi;
		gi=lMarker;

		int prevEnd=(*gi).first->start0;

		//int readcount=0;
		//int pos=0;

		string blockIResult;
		//int nBlockCounted=0;


		map<GffEntry::GBlockPtr,int>::iterator prevGI;

	//	int bn=0;
		//bool prevIsMinor=false;

		for(;gi!=blockCountMap.end();gi++)
		{
			//bn++;
//			cerr<<"a"<<endl;

			int count=(*gi).second;
//			cerr<<"b"<<endl;
			GffEntry::GBlockPtr pbl=(*gi).first;
//			cerr<<"c"<<endl;
			//moved  bool prevIsMinor=false;

			cerr<<"block:"<<pbl->getStart1()<<"-"<<pbl->getEnd1()<<"\t"<<"representing "<<count<<"/"<<transcriptCount<<" of transcript ";

			bool isMinor=((double)count/transcriptCount<ratioRep);




			if(!isMinor)
			{
				cerr<<" accepted"<<endl;
				_starts.push_back(pbl->getStart0()-_start);
				_sizes.push_back(pbl->getLength());
				//_scores.push_back(((double)count/transcriptCount)*1000);
				_end=pbl->getEnd1();
			}
			else
			{
				cerr<<" ignored: low representative"<<endl;
			}



			prevGI=gi;

			prevEnd=pbl->end1;
//			cerr<<"e"<<endl;
			if(isMinor)
				prevGI++;
//			cerr<<"f"<<endl;
		}

		if(_end==_start)
			continue;

		cout<<_chrom<<"\t"<<_start<<"\t"<<_end<<"\t"<<_name<<"\t"<<0<<"\t"<<_strand<<"\t"<<_start<<"\t"<<_end<<"\t"<<"0,0,0"<<"\t"<<_starts.size()<<"\t";
		for(vector<int>::iterator kk=_sizes.begin();kk!=_sizes.end();kk++)
		{
			cout<<StringUtil::str(*kk)<<",";
		}
		cout<<"\t";

		for(vector<int>::iterator kk=_starts.begin();kk!=_starts.end();kk++)
		{
			cout<<StringUtil::str(*kk)<<",";
		}
		cout<<endl;


	}//rof l



}



void calculateSEn(Splidar_OpFlag op,int minDistance,int maxDistance,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{


	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());
		//if(locusName!="CD44")
		//	continue;


		//cerr<<"CD44 found!"<<endl;
		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously


		//cerr<<" printing exon group: "<<locusName<<endl;
		::printExonTree(cerr,locus);

		//now calculate SEn
		set<NESuperGroupPtr> rightSuperGroups;

		locus->root->getSuperGroups(NULL,&rightSuperGroups);

		//cerr<<"done get supergroups"<<endl;


		for(set<NESuperGroupPtr>::iterator i=rightSuperGroups.begin();i!=rightSuperGroups.end();i++)
		{
				//do something with range;

				NESuperGroup* sgroup=*i;
				cerr<<"passing right super group"<<sgroup->getID()<<endl;
				SEn_SplidarGraph sgraph(op,minDistance,maxDistance,fout,fSeqOut,rafIn,sgroup,locus,locusName,true,readLength);
				sgraph.enterGraphLoop();
		}



		//cerr<<"done"<<endl;
		//cerr<<"****"<<endl;

		break;

	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}



void calculateSE2(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{


	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously


		//cerr<<" printing exon group: "<<locusName<<endl;
		//::printExonTree(cerr,locus);

		//now calculate MXE
		NExonGroup::IntExonPtrMMap dBoundMap;

	 	//[ (]) => Exon

		locus->root->getBoundariesMap(NULL,&dBoundMap);
		//cerr<<"done get bound"<<endl;

		NExonGroup::IntExonPtrMMapDI range;
		NExonGroup::IntExonPtrMMapI i=dBoundMap.begin();


		if(i==dBoundMap.end())
		{
			cerr<<"Strange Error: i==dBounMap.end() at the very beginning"<<endl;
			continue;
		}

		int prevB=i->first;

		range.first=i;
		range.second=(++i);
		//cerr<<"now do the job"<<endl;

		while(i!=dBoundMap.end())
		{
			if(i->first!=prevB)
			{
				//do something with range;
				SE_SplidarGraph sgraph(op,fout,fSeqOut,rafIn,range,locus,locusName,true,readLength);
				sgraph.enterGraphLoop();
				prevB=i->first;
				range.first=i;
			}

			range.second=(++i);
		}



		if(range.first!=dBoundMap.end())
		{
			//do something with range (remaining) here
			SE_SplidarGraph sgraph(op,fout,fSeqOut,rafIn,range,locus,locusName,true,readLength);
			sgraph.enterGraphLoop();
		}

		//cerr<<"done"<<endl;
		//cerr<<"****"<<endl;



	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}


/*
void _calculateSE2(string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream snip(snipFile.c_str());


	ofstream *fSeqOut=NULL;


	if(seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}


	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//cerr<<"****"<<endl;

		pair<int,NExonGroup::NExonGroupPtr> exongNext;
		NExonGroup::NExonGroupIterator i=locus->root->getExonGroupIterator(true);
		while(NExonGroup::NExonGroupIterator::isValidItem(exongNext=i.nextItem()))
		{
			//cerr<<"prcessing exonGroup "<<exongNext.second->sid<<" with coord=" <<exongNext.second->getBound()<<endl;
			SE_SplidarGraph sgraph(&fout,fSeqOut,rafIn,exongNext.second,locus,locusName,true,39);
			sgraph.enterGraphLoop();
		}

	}

	fout.close();
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}*/

/*void _calculateMXE(string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream snip(snipFile.c_str());
	ofstream *fSeqOut=NULL;

	if(seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

		}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//cerr<<"****"<<endl;

		pair<int,NExonGroup::NExonGroupPtr> exongNext;
		NExonGroup::NExonGroupIterator i=locus->root->getExonGroupIterator(true);
		while(NExonGroup::NExonGroupIterator::isValidItem(exongNext=i.nextItem()))
		{
			//cerr<<"prcessing exonGroup "<<exongNext.second->sid<<" with coord=" <<exongNext.second->getBound()<<endl;
			MXE_SplidarGraph sgraph(&fout,fSeqOut,rafIn,exongNext.second,locus,locusName,true,39);
			sgraph.enterGraphLoop();
		}

	}



	fout.close();
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}
*/
//ici;
void calculateMXE(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{

	ofstream *fout=NULL;

	ofstream *fSeqOut=NULL;

	if(op.getCount && snipFile.length()>0)
	{
		fout=new ofstream(snipFile.c_str());
	}


	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount&&!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//now calculate MXE
		NExonGroup::IntExonPtrMMap dBoundMap;

		locus->root->getBoundariesMap(NULL,&dBoundMap);
		cerr<<"done get bound"<<endl;

		NExonGroup::IntExonPtrMMapDI range;
		NExonGroup::IntExonPtrMMapI i=dBoundMap.begin();


		if(i==dBoundMap.end())
		{
			cerr<<"Strange Error: i==dBounMap.end() at the very beginning"<<endl;
			continue;
		}

		int prevB=i->first;

		range.first=i;
		range.second=(++i);
		cerr<<"now do the job"<<endl;

		while(i!=dBoundMap.end())
		{
			if(i->first!=prevB)
			{
				//do something with range;
				MXE_SplidarGraph sgraph(op,fout,fSeqOut,rafIn,range,locus,locusName,true,readLength); //
				sgraph.enterGraphLoop();
				prevB=i->first;
				range.first=i;
			}

			range.second=(++i);
		}



		if(range.first!=dBoundMap.end())
		{
			//do something with range (remaining) here
			MXE_SplidarGraph sgraph(op,fout,fSeqOut,rafIn,range,locus,locusName,true,readLength);
			sgraph.enterGraphLoop();
		}

	//	cerr<<"done"<<endl;
		//cerr<<"****"<<endl;



	}


	if(fout)
	{


	fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}

void calculateRI(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream *fout=NULL;
	ofstream *fSeqOut=NULL;

	if(op.getCount && snipFile.length()>0)
	{
		fout=new ofstream(snipFile.c_str());
	}


	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//cerr<<"****"<<endl;

		RI_Splidar rigraph(op,fout,fSeqOut,rafIn,locus->root,locus,locusName,readLength);
		//sgraph.enterGraphLoop();


	}


	if(fout)
	{
		fout->close();
		delete fout;

	}

	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}





void calculateATE(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream *fout=NULL;

	ofstream *fSeqOut=NULL;

	if(op.getCount && snipFile.length()>0)
	{
		fout=new ofstream(snipFile.c_str());
	}


	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}




	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//cerr<<"****"<<endl;

		ATE_Splidar ategraph(op,fout,fSeqOut,rafIn,locus,locusName,readLength);
		//sgraph.enterGraphLoop();


	}


	if(fout)
	{


	fout->close();
	delete fout;
	}


	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}


void calculateA53SS(Splidar_OpFlag op,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream *fout=NULL;

	ofstream *fSeqOut=NULL;

	if(op.getCount && snipFile.length()>0)
	{
		fout=new ofstream(snipFile.c_str());
	}


	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously

		//cerr<<" printing exon group"<<endl;
		//::printExonTree(cerr,locus);

		//cerr<<"****"<<endl;

		A53SS_Splidar a53sssplidar(op,fout,fSeqOut,rafIn,locus,locusName,readLength);
		//sgraph.enterGraphLoop();


	}



	if(fout)
	{


	fout->close();
	delete fout;
	}

	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}
}




void calculateSE(string snipFile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{
	ofstream snip(snipFile.c_str());

	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(!GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;

	string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
	string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

	int lc=0;
	int tevent=0;
	int tc=chrMapLocus->size();

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;


		//vector<SExonTraversalThread*> memBin;


		for(int ei=0;ei<locus->exonCount;ei++)
		{
			map<GffEntry::ExonPtr,multimap<int,SExonTraversalThread>*  > d2path;
			GffEntry::ExonPtr exon=locus->exons[ei];

			//start from this exons, go distance 2;
			queue<SExonTraversalThread> q;
			q.push(SExonTraversalThread(exon,0));

			while(!q.empty())
			{
				SExonTraversalThread& ut=q.front();

				GffEntry::ExonPtr u=ut.u;
				ut.exons.push_back(u);

				if(ut.distance>0)
				{
					multimap<int,SExonTraversalThread>* mmap;
				map<GffEntry::ExonPtr,multimap<int,SExonTraversalThread>*  >::iterator i=d2path.find(u);
					if(i==d2path.end())
					{
						mmap=new multimap<int,SExonTraversalThread>;
						d2path.insert(map<GffEntry::ExonPtr,multimap<int,SExonTraversalThread>*  >::value_type(u,mmap));

					}else
						mmap=(*i).second;

					mmap->insert(multimap<int,SExonTraversalThread>::value_type(ut.distance,ut));
				}

				if(ut.distance<2 && u->outJnxs)
				{
					for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator i=u->outJnxs->begin();i!=u->outJnxs->end();i++)
					{
						GffEntry::ExonPtr v=(*i).first;
						GffEntry::Jnx* jnx=(*i).second;
						SExonTraversalThread nt(ut);
						nt.u=v;
						nt.distance++;
						nt.jnxs.push_back(jnx);
						q.push(nt);
					}

				}



				q.pop();
			}

			//for Casette Exon:

			for(map<GffEntry::ExonPtr,multimap<int,SExonTraversalThread>*  >::iterator i=d2path.begin();i!=d2path.end();i++)
			{
				GffEntry::ExonPtr v=(*i).first;
				multimap<int,SExonTraversalThread>* mmap=(*i).second;

				ETDIterator ones=mmap->equal_range(1);
				ETDIterator twos=mmap->equal_range(2);

				for(ETSIterator oi=ones.first;oi!=ones.second;oi++)
				{
					for(ETSIterator ti=twos.first;ti!=twos.second;ti++)
					{
						//for twos: jnxs[0]+exon[1]+jnx[1];
						//for total: twos+jnx[0]+exon[0]+exon[1];

						SExonTraversalThread &ostt=(*oi).second;
						SExonTraversalThread &tstt=(*ti).second;

						GffEntry::Jnx* tj0=tstt.jnxs[0];
						GffEntry::ExonPtr te0=tstt.exons[0];
						GffEntry::ExonPtr te1=tstt.exons[1];
						GffEntry::ExonPtr te2=tstt.exons[2];
						GffEntry::Jnx* tj1=tstt.jnxs[1];

						GffEntry::Jnx* oj0=ostt.jnxs[0];
						GffEntry::ExonPtr oe0=ostt.exons[0];
						GffEntry::ExonPtr oe1=ostt.exons[1];


						KeyPair<int,int> tj0stat=tj0->getDensity(readLength,useUniquelyMappablePosition);
						KeyPair<int,int> te1stat=te1->getDensity(readLength,useUniquelyMappablePosition);
						KeyPair<int,int> tj1stat=tj1->getDensity(readLength,useUniquelyMappablePosition);
						KeyPair<int,int> oj0stat=oj0->getDensity(readLength,useUniquelyMappablePosition);
						KeyPair<int,int> oe0stat=oe0->getDensity(readLength,useUniquelyMappablePosition);
						KeyPair<int,int> oe1stat=oe1->getDensity(readLength,useUniquelyMappablePosition);

						double twoReads=tj0stat.k1+te1stat.k1+tj1stat.k1;
						double twoPos=tj0stat.k2+te1stat.k2+tj1stat.k2;

						//double totalReads=twoReads+oj0stat.k1+oe0stat.k1+oe1stat.k1;
						//double totalPos=twoPos+oj0stat.k2+oe0stat.k2+oe1stat.k2;

						double oneReads=oj0stat.k1;//+oe0stat.k1+oe1stat.k1;
						double onePos=oj0stat.k2;//+oe0stat.k2+oe1stat.k2;

						double totalReads=twoReads+oneReads+oe0stat.k1+oe1stat.k1;
						double totalPos=twoPos+oneReads+oe0stat.k2+oe1stat.k2;

						if(twoPos==0 || onePos==0)
						{
							cerr<<"Event cannot be seen, no uniquely mappable position for either isoform:";
							cerr<<oe0->chr<<":"<<oe0->getStart1()<<"-"<<oe0->getEnd1()<<"\t";
							cerr<<te1->chr<<":"<<te1->getStart1()<<"-"<<te1->getEnd1()<<"\t";
							cerr<<oe1->chr<<":"<<oe1->getStart1()<<"-"<<oe1->getEnd1()<<"\t";
							cerr<<endl;
							continue;
						}

						double twoDensity=twoReads/twoPos;
						double oneDensity=oneReads/onePos;

						tevent++;

						snip<<tevent<<"\t";

						snip<<(*(locus->names.begin()))<<"\t";
						snip<<locus->chr<<"\t";
						//snip<<te0->getStart1()<<"-"<<te0->getEnd1()<<"\t";
						snip<<oe0->getStart1()<<"-"<<oe0->getEnd1()<<"\t";
						snip<<te1->getStart1()<<"-"<<te1->getEnd1()<<"\t";
						snip<<oe1->getStart1()<<"-"<<oe1->getEnd1()<<"\t";
						if(rafIn)
						{
							rafIn->transfer(snip,oe0->getStart0(),oe0->getEnd1());
							snip<<"<";
							rafIn->transfer(snip,te1->getStart0(),te1->getEnd1());
							snip<<">";
							rafIn->transfer(snip,oe1->getStart0(),oe1->getEnd1());
							snip<<"\t";
						}
					//	snip<<te2->getStart1()<<"-"<<te2->getEnd1()<<"\t";
						snip<<EXCELHyperLinkPrefix<<locus->chr<<":"<<oe0->getStart1()<<"-"<<oe1->getEnd1()<<EXCELHyperLinkSuffix<<"\t";

						snip<<twoReads<<"\t"<<twoPos<<"\t"<<oneReads<<"\t"<<onePos<<"\t"<<(oneReads+twoReads)<<"\t";
						snip<<twoDensity<<"\t"<<oneDensity<<"\t";


						if( twoPos==0 || totalPos==0)
						{
							//snip<<"N/A";
							//snip<<"\t";
							snip<<"N/A";
						}
						else if(totalReads==0)
						{
							//snip<<"N/E";
							//snip<<"\t";
							snip<<"N/E";
						}
						else
						{

							//double EIL=(twoReads/twoPos)/(totalReads/totalPos);
							//snip<<"\t";
							double EIL=twoDensity/(oneDensity+twoDensity);
							snip<<EIL;
						}

						snip<<"\t";

						snip<<int(twoReads)<<"\t"; //NI
						snip<<int(totalReads)<<"\t"; //NE+

						snip<<tj0stat.k1<<"\t";
						snip<<tj0stat.k2<<"\t";
						snip<<te1stat.k1<<"\t";
						snip<<te1stat.k2<<"\t";
						snip<<tj1stat.k1<<"\t";
						snip<<tj1stat.k2<<"\t";
						snip<<oj0stat.k1<<"\t";
						snip<<oj0stat.k2<<"\t";
						snip<<oe0stat.k1<<"\t";
						snip<<oe0stat.k2<<"\t";
						snip<<oe1stat.k1<<"\t";
						snip<<oe1stat.k2;


						snip<<endl;

					}
				}



				delete mmap;
			}
		}






	}//rof l




	///////here
	snip.close();

}

/*void calculateSEn(Splidar_OpFlag op,int minDistance,int maxDistance,string snipFile,string seqfile,string chr,int readLength,bool useUniquelyMappablePosition,RandomAccessFile* rafIn)
{


	ofstream *fout=NULL;
	if(op.getCount && snipFile.length()>0)
	{

		fout=new ofstream(snipFile.c_str());

	}
	ofstream *fSeqOut=NULL;

	if(op.getSequence && seqfile.length()>0)
	{
		fSeqOut=new ofstream(seqfile.c_str());

	}



	if(!GffEntry::isGffLoaded())
	{
		die("Gff Not Loaded");
	}

	if(!GffEntry::isLocusLoaded())
	{
		die("Loci Not Loaded");
	}

	if(op.getCount && !GffEntry::isSelexaMatchesLoaded())
	{
		die("Selexa Matches Generally Not Loaded");
	}

	if(!GffEntry::isGBlocksLoaded())
	{
		die("GBlocks not loaded");
	}
	//now here comes something!!

	//for each Locus

	multimap<GffEntry::ExonPtr,GffEntry::Locus*> * chrMapLocus;

	map<string, multimap<GffEntry::ExonPtr,GffEntry::Locus*> * >::iterator mi=GffEntry::loci.find(chr);

	if(mi==GffEntry::loci.end())
	{
		die("No Locus Defined");
	}

	chrMapLocus=(*mi).second;


	int lc=0;

	int tc=chrMapLocus->size();

	//ofstream fout(snipFile.c_str());

	int jobID=-1;

	for(multimap<GffEntry::ExonPtr,GffEntry::Locus*>::iterator l=chrMapLocus->begin();l!=chrMapLocus->end();l++)
	{
		jobID--;
		lc++;
		if(lc%100==1)
			cerr<<"processing locus "<<lc<<" of "<<tc<<endl;

		GffEntry::Locus* locus=(*l).second;
		string locusName=(*locus->names.begin());

		//cerr<<" assign exon group " <<endl;
		::assignExonGroupPerLocus(locus,jobID,false); //save time: don't validate exon group. already done previously


		//cerr<<" printing exon group: "<<locusName<<endl;
		//::printExonTree(cerr,locus);

		//now calculate MXE
		NExonGroup::IntExonPtrMMap dBoundMap;

		locus->root->getBoundariesMap(NULL,&dBoundMap);
		cerr<<"done get bound"<<endl;

		NExonGroup::IntExonPtrMMapDI range;
		NExonGroup::IntExonPtrMMapI i=dBoundMap.begin();


		if(i==dBoundMap.end())
		{
			cerr<<"Strange Error: i==dBounMap.end() at the very beginning"<<endl;
			continue;
		}

		int prevB=i->first;

		range.first=i;
		range.second=(++i);
		cerr<<"now do the job"<<endl;

		while(i!=dBoundMap.end())
		{
			if(i->first!=prevB)
			{
				//do something with range;
				SEn_SplidarGraph sgraph(op,minDistance,maxDistance,fout,fSeqOut,rafIn,range,locus,locusName,true,39);
				sgraph.enterGraphLoop();
				prevB=i->first;
				range.first=i;
			}

			range.second=(++i);
		}



		if(range.first!=dBoundMap.end())
		{
			//do something with range (remaining) here
			SEn_SplidarGraph sgraph(op,minDistance,maxDistance,fout,fSeqOut,rafIn,range,locus,locusName,true,39);
			sgraph.enterGraphLoop();
		}

		cerr<<"done"<<endl;
		//cerr<<"****"<<endl;



	}


	if(fout)
	{
		fout->close();
	delete fout;
	}
	if(fSeqOut)
	{
		fSeqOut->close();

		delete fSeqOut;

	}



}*/

#endif /* SNIP_BACKUP_H_ */
