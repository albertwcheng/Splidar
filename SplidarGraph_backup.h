/*
 * SplidarGraph_backup.h
 *
 *  Created on: Jan 19, 2010
 *      Author: awcheng
 */

#ifndef SPLIDARGRAPH_BACKUP_H_
#define SPLIDARGRAPH_BACKUP_H_
/*

////**** Actual SplidarGraph class selector:

#define SG_JUNC_HASH 1
#define SG_REF_EG_HASH 2
#define SG_TARGET_EG_HASH 3

//here choose the HASH class

#define SE_HASH	SG_JUNC_HASH //SG_REF_EG_HASH
#define MXE_HASH  SG_JUNC_HASH //SG_REF_EG_HASH
#define AEP_HASH  SG_REF_EG_HASH
#define ATE_HASH  SG_JUNC_HASH//SG_TARGET_EG_HASH
#define A53SS_HASH SG_JUNC_HASH//SG_TARGET_EG_HASH


#if SE_HASH == SG_JUNC_HASH
	typedef SE_SplidarGraph_J SE_SplidarGraph;
#else

#endif

#if MXE_HASH == SG_JUNC_HASH
	typedef MXE_SplidarGraph_J MXE_SplidarGraph;
#else

#endif

#if ATE_HASH == SG_JUNC_HASH
	typedef ATE_Splidar_J ATE_Splidar;
#else
	typedef ATE_Splidar_TEG ATE_Splidar;
#endif

#if A53SS_HASH == SG_JUNC_HASH
	typedef A53SS_Splidar_J A53SS_Splidar;
#else
	typedef A53SS_Splidar_TEG A53SS_Splidar;
#endif
*/



//not checked
class ATE_Splidar_J: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	inline ATE_Splidar_J(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength):GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf), op(_op)
	{
		//for each transcript
		//add to second exon into secondExonSet

		//set { fixedCoord => variantExon => commonExon /*common is just to avoid the same exact set of exons*/ }
		set<ExonPairWJnx > twoMapAG5E; //genomic alternative 5' exon
		set<ExonPairWJnx > twoMapAG3E; //genomic alternative 3' exon
		typedef set<ExonPairWJnx>::iterator ExonPairWJnxI;
		typedef pair<ExonPairWJnxI,ExonPairWJnxI> ExonPairWJnxDI;
		typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
		typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
		typedef pair<sg_i,sg_i> Dsg_i;

		for(vector<GffEntry*>::iterator it=locus->transcripts.begin();it!=locus->transcripts.end();it++)
		{
			GffEntry*transcript=*it;
			int exonCount=transcript->exonCount;
			if(exonCount<2)
				continue;

			GffEntry::ExonPtr secondExon=transcript->exons[1];
			GffEntry::ExonPtr firstExon=transcript->exons[0];

			GffEntry::ExonPtr lastExon=transcript->exons[exonCount-1]; //last
			GffEntry::ExonPtr meiyiExon=transcript->exons[exonCount-2]; //secondLast



			if(firstExon->outJnxs)
			{
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=firstExon->outJnxs->find(secondExon);
				if(ise==firstExon->outJnxs->end())
				{
					cerr<<"strange error: first and second exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					continue;
				}

				GffEntry::Jnx* jnx=ise->second;
				//the common exon is the second exon;
				//the variant exon is the first exon;
				//fixedcoord is the start of the second exon;
				twoMapAG5E.insert(ExonPairWJnx(secondExon->getStart1(),firstExon,secondExon,jnx));

			}else
			{
				cerr<<"strange firstExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
			}




			if(meiyiExon->outJnxs)
			{
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=meiyiExon->outJnxs->find(lastExon);
				if(ise==meiyiExon->outJnxs->end())
				{
					cerr<<"strange error: meiyi and last exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					continue;
				}

				//the common exon is the meiyi exon;
				//the variant exon is the last exon;
				//fixedcoord is the end of meiyi exon;

				GffEntry::Jnx* jnx=ise->second;
				twoMapAG3E.insert(ExonPairWJnx(meiyiExon->getEnd1(),lastExon,meiyiExon,jnx));

			}else
			{
				cerr<<"strange meiyiExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
			}




		}


		//process g5'
		//cerr<<"size of twoMapAG5E="<<twoMapAG5E.size()<<endl;

		set<ExonPairWJnx>::iterator i=twoMapAG5E.begin(); ///
		ExonPairWJnxDI curRange;
		int prevB;

		while(i!=twoMapAG5E.end()) ///
		{
			prevB=i->fixedCoord;
		//	cerr<<i->fixedCoord<<"\t";
			curRange.first=i;
			curRange.second=(++i);

			while(i!=twoMapAG5E.end()) ///
			{
				if(i->fixedCoord!=prevB)
					break;

				curRange.second=(++i);
			}

			//now the curRange contains all the stuff with the same fixed Coord, and the variant Exon are sorted by start



			for(ExonPairWJnxI x=curRange.first;x!=curRange.second;x++)
				{
					ExonPairWJnxI y=x;
					y++;
					for(;y!=curRange.second;y++)
					{
					//	cerr<<"in business:";
						//just need to swap here!
						ExonPairWJnx proximaliev=*y;
						ExonPairWJnx distaliev=*x;

						GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
						KeyPair<int,int> proximalBound=proximalExon->getBound();
						GffEntry::Jnx* proximalJnx=proximaliev.jnx;

						GffEntry::ExonPtr distalExon=distaliev.variantExon;
						KeyPair<int,int> distalBound=distalExon->getBound();
						GffEntry::Jnx* distalJnx=distaliev.jnx;

						if(::areOverlapping(proximalBound,distalBound))
						{
					//		cerr<<"overlapping "<<proximalBound<<" vs "<<distalBound<<endl;
							continue; //the two variable exons overlap, ignore
						}

					//	cerr<<"ok! continue!"<<endl;

						TrafficInfoAS xinInfo;
						TrafficInfoAS yexInfo;

						xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximalExon->getEnd1())+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand;
						yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distalExon->getEnd1())+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand;
						xinInfo.gsid=proximalExon->exonGroup->sid+">"+proximaliev.commonExon->exonGroup->relRootString;
						yexInfo.gsid=distalExon->exonGroup->sid+">"+distaliev.commonExon->exonGroup->relRootString;


						KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());



						xinInfo.bounds.insert(proximalExon->getBound());
						xinInfo.bounds.insert(commonBound);

						yexInfo.bounds.insert(distalExon->getBound());
						yexInfo.bounds.insert(commonBound);

						if(this->op.getCount && this->hasDataOut())
						{


							Dsg_i commonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.commonExon->blocks,commonBound,true);
							KeyPair<int,int> commonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(commonBlocks,readLength);


							addDensityVectors(xinInfo.flankingExonsFlow,commonDen);
							addDensityVectors(yexInfo.flankingExonsFlow,commonDen);
							addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
							addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
							addDensityVectors(xinInfo.middleExonsFlow,proximalExon->getDensity(readLength,true));
							addDensityVectors(yexInfo.middleExonsFlow,distalExon->getDensity(readLength,true));


							xinInfo.JRCheckString=StringUtil::str(xinInfo.jnxsFlow.k1);
							xinInfo.JPCheckString=StringUtil::str(xinInfo.jnxsFlow.k2);
							xinInfo.JRF=xinInfo.jnxsFlow.k1;
							xinInfo.JPF=xinInfo.jnxsFlow.k2;

							yexInfo.JRCheckString=StringUtil::str(yexInfo.jnxsFlow.k1);
							yexInfo.JPCheckString=StringUtil::str(yexInfo.jnxsFlow.k2);
							yexInfo.JRF=yexInfo.jnxsFlow.k1;
							yexInfo.JPF=yexInfo.jnxsFlow.k2;


							addDensityVectors(xinInfo.noFlankingInfo,xinInfo.middleExonsFlow);
							addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
							xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;
						//	xinInfo.JRF=xinInfo.jnxsFlow.k1;
						//	xinInfo.JPF=xinInfo.jnxsFlow.k2;

							addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
							addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow);
							yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;
						//	yexInfo.JRF=yexInfo.jnxsFlow.k1;
						//	yexInfo.JPF=yexInfo.jnxsFlow.k2;

						}


						this->outData((this->strand==GffEntry::FORWARD)?"AFE":"ALE",xinInfo,yexInfo);

					}
				}

		}




		//process g3'
		i=twoMapAG3E.begin(); ///
		//ExonPairWJnxDI curRange;


		while(i!=twoMapAG3E.end()) ///
		{
			prevB=i->fixedCoord;

			curRange.first=i;
			curRange.second=(++i);

			while(i!=twoMapAG3E.end()) ///
			{
				if(i->fixedCoord!=prevB)
					break;

				curRange.second=(++i);
			}

			//now the curRange contains all the stuff with the same fixed Coord, and the variant Exon are sorted by start



			for(ExonPairWJnxI x=curRange.first;x!=curRange.second;x++)
				{
					ExonPairWJnxI y=x;
					y++;
					for(;y!=curRange.second;y++)
					{
						//just need to swap here!
						ExonPairWJnx proximaliev=*x;
						ExonPairWJnx distaliev=*y;

						GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
						KeyPair<int,int> proximalBound=proximalExon->getBound();
						GffEntry::Jnx* proximalJnx=proximaliev.jnx;

						GffEntry::ExonPtr distalExon=distaliev.variantExon;
						KeyPair<int,int> distalBound=distalExon->getBound();
						GffEntry::Jnx* distalJnx=distaliev.jnx;

						if(::areOverlapping(proximalBound,distalBound))
							continue; //the two variable exons overlap, ignore

						TrafficInfoAS xinInfo;
						TrafficInfoAS yexInfo;

						xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximalExon->getStart1())+":"+this->strand;
						yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distalExon->getStart1())+":"+this->strand;
						xinInfo.gsid=proximaliev.commonExon->exonGroup->relRootString+">"+proximalExon->exonGroup->sid; //flip!
						yexInfo.gsid=distaliev.commonExon->exonGroup->relRootString+">"+distalExon->exonGroup->sid;


						KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());

						Dsg_i commonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.commonExon->blocks,commonBound,true);
						KeyPair<int,int> commonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(commonBlocks,readLength);

						xinInfo.bounds.insert(commonBound);
						xinInfo.bounds.insert(proximalExon->getBound());

						yexInfo.bounds.insert(commonBound);
						yexInfo.bounds.insert(distalExon->getBound());


						addDensityVectors(xinInfo.flankingExonsFlow,commonDen);
						addDensityVectors(yexInfo.flankingExonsFlow,commonDen);
						addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
						addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
						addDensityVectors(xinInfo.middleExonsFlow,proximalExon->getDensity(readLength,true));
						addDensityVectors(yexInfo.middleExonsFlow,distalExon->getDensity(readLength,true));


						xinInfo.JRCheckString=StringUtil::str(xinInfo.jnxsFlow.k1);
						xinInfo.JPCheckString=StringUtil::str(xinInfo.jnxsFlow.k2);
						xinInfo.JRF=xinInfo.jnxsFlow.k1;
						xinInfo.JPF=xinInfo.jnxsFlow.k2;

						yexInfo.JRCheckString=StringUtil::str(yexInfo.jnxsFlow.k1);
						yexInfo.JPCheckString=StringUtil::str(yexInfo.jnxsFlow.k2);
						yexInfo.JRF=yexInfo.jnxsFlow.k1;
						yexInfo.JPF=yexInfo.jnxsFlow.k2;



						addDensityVectors(xinInfo.noFlankingInfo,xinInfo.middleExonsFlow);
						addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
						xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;
				//		xinInfo.JRF=xinInfo.jnxsFlow.k1;
				//		xinInfo.JPF=xinInfo.jnxsFlow.k2;

						addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
						addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow);
						yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;
				//		yexInfo.JRF=yexInfo.jnxsFlow.k1;
				//		yexInfo.JPF=yexInfo.jnxsFlow.k2;

						this->outData((this->strand==GffEntry::REVERSE)?"AFE":"ALE",xinInfo,yexInfo);

					}
				}

		}




	}

	inline ATE_Splidar_J(ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength, bool _oldconstructor):GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf)
	{
		//for each transcript
		//add to second exon into secondExonSet

		//secondExon => (jnx, firstExon);
		map<GffEntry::ExonPtr, set<InEdgeAndVertex>* > twoMapAG5E; //genomic alternative 5' exon
		map<GffEntry::ExonPtr, set<InEdgeAndVertex>* > twoMapAG3E; //genomic alternative 3' exon


		for(vector<GffEntry*>::iterator it=locus->transcripts.begin();it!=locus->transcripts.end();it++)
		{
			GffEntry*transcript=*it;
			int exonCount=transcript->exonCount;
			if(exonCount<2)
				continue;

			GffEntry::ExonPtr secondExon=transcript->exons[1];
			GffEntry::ExonPtr firstExon=transcript->exons[0];

			GffEntry::ExonPtr lastExon=transcript->exons[exonCount-1]; //last
			GffEntry::ExonPtr meiyiExon=transcript->exons[exonCount-2]; //secondLast


			//hash the g5' alt exon
			set<InEdgeAndVertex>* ies;
			map<GffEntry::ExonPtr, set<InEdgeAndVertex>* >::iterator i=twoMapAG5E.find(secondExon);
			if(i==twoMapAG5E.end())
			{
				ies=new set<InEdgeAndVertex>;
				twoMapAG5E.insert(map<GffEntry::ExonPtr, set<InEdgeAndVertex>* >::value_type(secondExon,ies));
			}
			else
				ies=i->second;


			if(firstExon->outJnxs)
			{
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=firstExon->outJnxs->find(secondExon);
				if(ise==firstExon->outJnxs->end())
				{
					cerr<<"strange error: first and second exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					continue;
				}

				GffEntry::Jnx* jnx=ise->second;
				ies->insert(InEdgeAndVertex(jnx,firstExon));

			}else
			{
				cerr<<"strange firstExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
			}

			//hash the g3' alt exon
			 i=twoMapAG3E.find(meiyiExon);
			if(i==twoMapAG3E.end())
			{
				ies=new set<InEdgeAndVertex>;
				twoMapAG3E.insert(map<GffEntry::ExonPtr, set<InEdgeAndVertex>* >::value_type(meiyiExon,ies));
			}
			else
				ies=i->second;

			if(meiyiExon->outJnxs)
			{
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=meiyiExon->outJnxs->find(lastExon);
				if(ise==meiyiExon->outJnxs->end())
				{
					cerr<<"strange error: meiyi and last exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					continue;
				}

				GffEntry::Jnx* jnx=ise->second;
				ies->insert(InEdgeAndVertex(jnx,lastExon));

			}else
			{
				cerr<<"strange meiyiExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
			}




		}


		//process g5'
		for(map<GffEntry::ExonPtr, set<InEdgeAndVertex>* >::iterator i=twoMapAG5E.begin();i!=twoMapAG5E.end();i++)
		{
			GffEntry::ExonPtr commonExon=i->first;
			set<InEdgeAndVertex>* ies=i->second;
			if(ies->size()<2) //no alternative, only one
				continue;

			//x is proximal, y is distal
			for(set<InEdgeAndVertex>::reverse_iterator x=ies->rbegin();x!=ies->rend();x++)
			{
				set<InEdgeAndVertex>::reverse_iterator y=x;
				y++;
				for(;y!=ies->rend();y++)
				{
					InEdgeAndVertex proximaliev=*x;
					InEdgeAndVertex distaliev=*y;

					GffEntry::ExonPtr proximalExon=proximaliev.exon;
					KeyPair<int,int> proximalBound=proximalExon->getBound();
					GffEntry::Jnx* proximalJnx=proximaliev.inJnx;

					GffEntry::ExonPtr distalExon=distaliev.exon;
					KeyPair<int,int> distalBound=distalExon->getBound();
					GffEntry::Jnx* distalJnx=distaliev.inJnx;

					if(::areOverlapping(proximalBound,distalBound))
						continue; //the two variable exons overlap, ignore

					TrafficInfoAS xinInfo;
					TrafficInfoAS yexInfo;

					xinInfo.gsid=proximalExon->exonGroup->sid+"_"+commonExon->exonGroup->sid;
					yexInfo.gsid=distalExon->exonGroup->sid+"_"+commonExon->exonGroup->sid;

					xinInfo.bounds.insert(proximalExon->getBound());
					xinInfo.bounds.insert(commonExon->getBound());

					yexInfo.bounds.insert(distalExon->getBound());
					yexInfo.bounds.insert(commonExon->getBound());

					addDensityVectors(xinInfo.flankingExonsFlow,commonExon->getDensity(readLength,true));
					addDensityVectors(yexInfo.flankingExonsFlow,commonExon->getDensity(readLength,true));
					addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
					addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
					addDensityVectors(xinInfo.middleExonsFlow,proximalExon->getDensity(readLength,true));
					addDensityVectors(yexInfo.middleExonsFlow,distalExon->getDensity(readLength,true));

					addDensityVectors(xinInfo.noFlankingInfo,xinInfo.middleExonsFlow);
					addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
					xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;

					addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
					addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow);
					yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;

					this->outData((this->strand==GffEntry::FORWARD)?"AFE":"ALE",xinInfo,yexInfo);

				}
			}

		} ///

		//process g3'
		for(map<GffEntry::ExonPtr, set<InEdgeAndVertex>* >::iterator i=twoMapAG3E.begin();i!=twoMapAG3E.end();i++)
		{
			GffEntry::ExonPtr commonExon=i->first; //meiyi exon
			set<InEdgeAndVertex>* ies=i->second;
			if(ies->size()<2) //no alternative, only one
				continue;

			//x is proximal, y is distal
			for(set<InEdgeAndVertex>::iterator x=ies->begin();x!=ies->end();x++)
			{
				set<InEdgeAndVertex>::iterator y=x;
				y++;
				for(;y!=ies->end();y++)
				{
					InEdgeAndVertex proximaliev=*x;
					InEdgeAndVertex distaliev=*y;

					GffEntry::ExonPtr proximalExon=proximaliev.exon;
					KeyPair<int,int> proximalBound=proximalExon->getBound();
					GffEntry::Jnx* proximalJnx=proximaliev.inJnx;

					GffEntry::ExonPtr distalExon=distaliev.exon;
					KeyPair<int,int> distalBound=distalExon->getBound();
					GffEntry::Jnx* distalJnx=distaliev.inJnx;

					if(::areOverlapping(proximalBound,distalBound))
						continue; //the two variable exons overlap, ignore

					TrafficInfoAS xinInfo;
					TrafficInfoAS yexInfo;

					xinInfo.gsid=proximalExon->exonGroup->sid+"_"+commonExon->exonGroup->sid;
					yexInfo.gsid=distalExon->exonGroup->sid+"_"+commonExon->exonGroup->sid;

					xinInfo.bounds.insert(proximalExon->getBound());
					xinInfo.bounds.insert(commonExon->getBound());

					yexInfo.bounds.insert(distalExon->getBound());
					yexInfo.bounds.insert(commonExon->getBound());

					addDensityVectors(xinInfo.flankingExonsFlow,commonExon->getDensity(readLength,true));
					addDensityVectors(yexInfo.flankingExonsFlow,commonExon->getDensity(readLength,true));
					addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
					addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
					addDensityVectors(xinInfo.middleExonsFlow,proximalExon->getDensity(readLength,true));
					addDensityVectors(yexInfo.middleExonsFlow,distalExon->getDensity(readLength,true));

					addDensityVectors(xinInfo.noFlankingInfo,xinInfo.middleExonsFlow);
					addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
					xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;

					addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
					addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow);
					yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;

					this->outData((this->strand==GffEntry::REVERSE)?"AFE":"ALE",xinInfo,yexInfo);

				}
			}

		} ///



	}
};


class MXE_SplidarGraph_J : public SplidarGraph {
public:
	inline MXE_SplidarGraph_J(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NExonGroup::NExonGroupPtr _exongroup, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength) :
		SplidarGraph(_op,_fout, _foutSeq, _raf, _exongroup, _locus, _locusName, _checkLocusName, 2, _readLength) {

	}

	inline MXE_SplidarGraph_J(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NExonGroup::IntExonPtrMMapDI _range, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength) :
		SplidarGraph(_op,_fout, _foutSeq, _raf, _range, _locus, _locusName, _checkLocusName, 2, _readLength) {

	}

	void endGraph() {
		//cerr<<"end graph SE"<<endl;

		string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		//dnCoord->thread
		map<int,vector<ExonPostItem*>* > twoThread;
		//cerr<<"a";

		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			//just want the twos;

			ExonPostRow::I epri=epr->find(2);


			if(epri==epr->end())
			{
				//no twos, ignore;
				continue;
			}

			//cerr<<"b"<<endl;
			//now register
			int dnCoord=terminalExon->start;

			map<int,vector<ExonPostItem*>* >::iterator miv=twoThread.find(dnCoord);

			vector<ExonPostItem*>* epiV;
			if(miv==twoThread.end())
			{
				epiV=new vector<ExonPostItem*>;
				twoThread.insert(map<int,vector<ExonPostItem*>* >::value_type(dnCoord,epiV));
			}
			else
				epiV=miv->second;

			//cerr<<"c"<<endl;
			ExonPost* ep=epri->second;

			for(ExonPost::I epi=ep->begin();epi!=ep->end();epi++)
			{
				epiV->push_back(&(*epi));
			}

			//cerr<<"d"<<endl;
		}

		//cerr<<"e"<<endl;
			//here comes the critical step!
			//for each terminal coordinate

			for (map<int,vector<ExonPostItem*>* >::iterator i1=twoThread.begin();i1!=twoThread.end();i1++) {

				//int dnCoord=i1->first;  //remove 6/6/2009
				vector<ExonPostItem*>* s2=i1->second;

				if(s2->size()<2) continue; //just one path to such dnCoord. ignore

				//cerr<<"f"<<endl;


				for (vector<ExonPostItem*>::iterator x3=s2->begin(); x3!=s2->end(); x3++) {
					vector<ExonPostItem*>::iterator y3=x3;
					y3++;
					for (; y3!=s2->end(); y3++) {
						//cerr<<"0g"<<endl;
						ExonPostItem* one_epi=*x3;
						ExonPostItem* two_epi=*y3;
						//cerr<<"g"<<endl;



						AESPathThread* one_thr=one_epi->exonThread;
						AESPathThread* two_thr=two_epi->exonThread;
						AESPathThread::backtrack_iterator one_right=
								one_thr->rbegin(one_epi->vertexPointer);
						AESPathThread::backtrack_iterator one_end1=
								one_thr->rend();
						AESPathThread::backtrack_iterator two_right=
								two_thr->rbegin(two_epi->vertexPointer);
						AESPathThread::backtrack_iterator two_end1=
								two_thr->rend();

						AESPathThread::backtrack_iterator one_center=one_right;
						AESPathThread::backtrack_iterator two_center=two_right;
						one_center++;
						two_center++;

						//cerr<<"h"<<endl;


						AESPathThread::backtrack_iterator one_left=one_center;
						one_left++;
						AESPathThread::backtrack_iterator two_left=two_center;
						two_left++;

						//cerr<<"i"<<endl;
						//if center exons overlap, ignored;
						if(::areOverlapping(one_center->exon->getBound(),two_center->exon->getBound()))
							continue;
						//cerr<<"i1"<<endl;
						//inclusion form is the g

						string exonCoordPathString;
						string onegsidString;
						string twogsidString;


						TrafficInfoAS one_info1;
						TrafficInfoAS two_info1;

						typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
						typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
						typedef pair<sg_i,sg_i> Dsg_i;


						string chr=one_center->exon->chr;

						//shouildn't use relRootString anymore

						one_info1.gsid=one_left->exon->exonGroup->relRootString+">"+one_center->exon->exonGroup->sid+">"+one_right->exon->exonGroup->relRootString;
						two_info1.gsid=two_left->exon->exonGroup->relRootString+">"+two_center->exon->exonGroup->sid+">"+two_right->exon->exonGroup->relRootString;




						KeyPair<int,int> leftOverlapBound=::overlapBound(one_left->exon->getBound(),two_left->exon->getBound());
						KeyPair<int,int> rightOverlapBound=::overlapBound(one_right->exon->getBound(),two_right->exon->getBound());

						one_info1.jnxstring=chr+":"+StringUtil::str(leftOverlapBound.k2+1)+this->strand+">"+chr+":"+StringUtil::str(one_center->exon->getStart1())+this->strand+","+chr+":"+StringUtil::str(one_center->exon->getEnd1())+this->strand+">"+chr+":"+StringUtil::str(rightOverlapBound.k1)+this->strand;
						two_info1.jnxstring=chr+":"+StringUtil::str(leftOverlapBound.k2+1)+this->strand+">"+chr+":"+StringUtil::str(two_center->exon->getStart1())+this->strand+","+chr+":"+StringUtil::str(two_center->exon->getEnd1())+this->strand+">"+chr+":"+StringUtil::str(rightOverlapBound.k1)+this->strand;



						one_info1.bounds.insert(leftOverlapBound);
						one_info1.bounds.insert(one_center->exon->getBound());
						one_info1.bounds.insert(rightOverlapBound);

						two_info1.bounds.insert(leftOverlapBound);
						two_info1.bounds.insert(two_center->exon->getBound());
						two_info1.bounds.insert(rightOverlapBound);



						if(this->op.getCount && this->hasDataOut())
						{
						//cerr<<"k"<<endl;
						Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_left->exon->blocks,leftOverlapBound,true);
						Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_right->exon->blocks,rightOverlapBound,true);
						//cerr<<"l"<<endl;

						//cerr<<"m"<<endl;
						KeyPair<int,int> leftDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(leftBlocks,readLength);
						//cerr<<"m2"<<endl;
						KeyPair<int,int> one_middleDen=one_center->exon->getDensity(readLength,true);
						KeyPair<int,int> two_middleDen=two_center->exon->getDensity(readLength,true);
						//cerr<<"m3"<<endl;
						KeyPair<int,int> rightDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(rightBlocks,readLength);
						//cerr<<"n"<<endl;
						//////here-->
							//one_info1;
							//right flanking exon, and jnx before.
							{

								KeyPair<int,int> jnxDen=one_right->inJnx->getDensity(readLength,true);

								one_info1.JRCheckString=StringUtil::str(jnxDen.k1);
								one_info1.JPCheckString=StringUtil::str(jnxDen.k2);
								one_info1.JRF*=jnxDen.k1;
								one_info1.JPF*=jnxDen.k2;

								addDensityVectors(one_info1.jnxsFlow,jnxDen);
								addDensityVectors(one_info1.noFlankingInfo,jnxDen);
								addDensityVectors(one_info1.withFlankingInfo,jnxDen);


								addDensityVectors(one_info1.flankingExonsFlow,rightDen);

								//cerr<<"gt"<<"<a1"<<endl;
							}
							 //left flanking exon (no jnx before)
							{
								//cerr<<"gt"<<"a3"<<endl;

								addDensityVectors(one_info1.flankingExonsFlow,leftDen);

								//cerr<<"gt"<<"<a3"<<endl;
							}
							 //middle exon and jnx to the left
							{
								//cerr<<"gt"<<"a2"<<endl;

								KeyPair<int,int> jnxDen=one_center->inJnx->getDensity(readLength,true);

								one_info1.JRCheckString=StringUtil::str(jnxDen.k1)+","+one_info1.JRCheckString;
								one_info1.JPCheckString=StringUtil::str(jnxDen.k2)+","+one_info1.JPCheckString;
								one_info1.JRF*=jnxDen.k1;
								one_info1.JPF*=jnxDen.k2;

								addDensityVectors(one_info1.noFlankingInfo,jnxDen);
								addDensityVectors(one_info1.noFlankingInfo,one_middleDen);
								addDensityVectors(one_info1.withFlankingInfo,jnxDen);
								addDensityVectors(one_info1.withFlankingInfo,one_middleDen);

								addDensityVectors(one_info1.middleExonsFlow,one_middleDen);
								addDensityVectors(one_info1.jnxsFlow,jnxDen);

								//cerr<<"gt"<<"<a2"<<endl;
							}


							//two_info1;
							//right flanking exon, and jnx before.
							{

								KeyPair<int,int> jnxDen=two_right->inJnx->getDensity(readLength,true);

								two_info1.JRCheckString=StringUtil::str(jnxDen.k1);
								two_info1.JPCheckString=StringUtil::str(jnxDen.k2);
								two_info1.JRF*=jnxDen.k1;
								two_info1.JPF*=jnxDen.k2;

								addDensityVectors(two_info1.jnxsFlow,jnxDen);
								addDensityVectors(two_info1.noFlankingInfo,jnxDen);
								addDensityVectors(two_info1.withFlankingInfo,jnxDen);

								addDensityVectors(two_info1.flankingExonsFlow,rightDen);

								//cerr<<"gt"<<"<a1"<<endl;
							}
							 //left flanking exon (no jnx before)
							{
								//cerr<<"gt"<<"a3"<<endl;

								addDensityVectors(two_info1.flankingExonsFlow,leftDen);

								//cerr<<"gt"<<"<a3"<<endl;
							}
							 //middle exon and jnx to the left
							{
								//cerr<<"gt"<<"a2"<<endl;

								KeyPair<int,int> jnxDen=two_center->inJnx->getDensity(readLength,true);

								two_info1.JRCheckString=StringUtil::str(jnxDen.k1)+","+two_info1.JRCheckString;
								two_info1.JPCheckString=StringUtil::str(jnxDen.k2)+","+two_info1.JPCheckString;
								two_info1.JRF*=jnxDen.k1;
								two_info1.JPF*=jnxDen.k2;

								addDensityVectors(two_info1.noFlankingInfo,jnxDen);
								addDensityVectors(two_info1.noFlankingInfo,two_middleDen);
								addDensityVectors(two_info1.withFlankingInfo,jnxDen);
								addDensityVectors(two_info1.withFlankingInfo,two_middleDen);

								addDensityVectors(two_info1.middleExonsFlow,two_middleDen);
								addDensityVectors(two_info1.jnxsFlow,jnxDen);

								//cerr<<"gt"<<"<a2"<<endl;
							}


						}



						//cerr<<"0j"<<endl;




						//cerr<<"j"<<endl;


						//cerr<<"k"<<endl;

						TrafficInfoAS* one_info;
						TrafficInfoAS* two_info;
						if(one_center->exon->getBound()<two_center->exon->getBound())
						{
							//one is genomic 5'
							//two is genomic 3'
							if(strand==GffEntry::FORWARD)
							{
								//two need to be transcript 5' so swap
								two_info=&one_info1;
								one_info=&two_info1;
							}
							else
							{
								one_info=&one_info1;
								two_info=&two_info1;
							}
						}
						else
						{
							//two is genomic 5'
							//one is genomic 3'
							if(strand==GffEntry::REVERSE)
							{
								//two need to be transcript 5' so swap
								two_info=&one_info1;
								one_info=&two_info1;
							}
							else
							{
								one_info=&one_info1;
								two_info=&two_info1;
							}
						}

						//cerr<<"l"<<endl;

						this->outData("MXE",*two_info,*one_info);

						//cerr<<"m"<<endl;
					}
				}

			}

		//now free memory from the registry
			//cerr<<"n"<<endl;
			for(map<int,vector<ExonPostItem*>* >::iterator i=twoThread.begin();i!=twoThread.end();i++)
			{
				delete i->second;
			}

			//cerr<<"o"<<endl;


	}


	void _endGraph() {
		//cerr<<"end graph SE"<<endl;

		string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			SrcExonPostStruct * srcS=epr->getSrcExonStructure();

			//now here comes the easier but critical step!

			for (SrcExonPostStruct::I1 i1=srcS->begin(); i1!=srcS->end(); i1++) {

				GffEntry::ExonPtr srcExon=i1->first;
				SrcExonPostStruct::S2* s2=i1->second;


				SrcExonPostStruct::I2 twoAway=s2->find(2);

				if (twoAway==s2->end()) {
					continue;
				}

				SrcExonPostStruct::S3* twos=twoAway->second;

				for (SrcExonPostStruct::I3 x3=twos->begin(); x3!=twos->end(); x3++) {
					SrcExonPostStruct::I3 y3=x3;
					y3++;
					for (; y3!=twos->end(); y3++) {
						ExonPostItem* one_epi=*x3;
						ExonPostItem* two_epi=*y3;




						AESPathThread* one_thr=one_epi->exonThread;
						AESPathThread* two_thr=two_epi->exonThread;
						AESPathThread::backtrack_iterator one_beg=
								one_thr->rbegin(one_epi->vertexPointer);
						AESPathThread::backtrack_iterator one_end1=
								one_thr->rend();
						AESPathThread::backtrack_iterator two_beg=
								two_thr->rbegin(two_epi->vertexPointer);
						AESPathThread::backtrack_iterator two_end1=
								two_thr->rend();

						AESPathThread::backtrack_iterator one_center=one_beg;
						AESPathThread::backtrack_iterator two_center=two_beg;
						one_center++;
						two_center++;

						//if center exons overlap, ignored;
						if(::areOverlapping(one_center->exon->getBound(),two_center->exon->getBound()))
							continue;

						//inclusion form is the g

						string exonCoordPathString;
						string onegsidString;
						string twogsidString;

						TrafficInfoAS one_info1=getTrafficAS(one_beg,
								one_end1,false);
						TrafficInfoAS two_info1=getTrafficAS(two_beg,
								two_end1,false); //do not include common regions for NE+

						TrafficInfoAS* one_info;
						TrafficInfoAS* two_info;
						if(one_center->exon->getBound()<two_center->exon->getBound())
						{
							//one is genomic 5'
							//two is genomic 3'
							if(strand==GffEntry::FORWARD)
							{
								//two need to be transcript 5' so swap
								two_info=&one_info1;
								one_info=&two_info1;
							}
							else
							{
								one_info=&one_info1;
								two_info=&two_info1;
							}
						}
						else
						{
							//two is genomic 5'
							//one is genomic 3'
							if(strand==GffEntry::REVERSE)
							{
								//two need to be transcript 5' so swap
								two_info=&one_info1;
								one_info=&two_info1;
							}
							else
							{
								one_info=&one_info1;
								two_info=&two_info1;
							}
						}



						this->outData("MXE",*two_info,*one_info);

					}
				}

			}

			delete srcS;

		}
	}

};


/*backup 1/6/2010
class RI_Splidar: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	inline RI_Splidar(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::NExonGroupPtr _exongroup,GffEntry::Locus* _locus,string _locusName,int _readLength):GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf), op(_op)
	{
		NExonGroup::NExonGroupIterator ngi=_exongroup->getExonGroupIterator(false,TRAVEL_DOWNTO_LEAVES);
		pair<int,NExonGroup::NExonGroupPtr> exongNext;
		while(NExonGroup::NExonGroupIterator::isValidItem(exongNext=ngi.nextItem()))
		{
			NExonGroup::NExonGroupPtr curExonGroup=exongNext.second;
			if(curExonGroup->hasChildren() && curExonGroup->levelExons.size()>0) //exclude root. have both itself and children
			{



				if(curExonGroup->children->size()<2)
				{
					cerr<<"Strange Error: has only 1 child?"<<endl;
				}

				//
				//  [       ]
				// [  ]  [  ]





				NExonGroup::ChildI leftI=curExonGroup->children->begin();
				NExonGroup::ChildI rightI=leftI;
				rightI++;
				NExonGroup::ChildI rightBoundIE=curExonGroup->children->rbegin().base();

				//exon groups
				//[ curExonGroup                               ]
				//[leftI] [rightI] ...... [   ] [ rightBoundIE ] end
				// >>>>>>                 [leftI] [rightI]

				for(;rightI!=rightBoundIE;leftI++,rightI++)
				{
					NExonGroup::NExonGroupPtr leftExonGroup=*leftI;
					NExonGroup::NExonGroupPtr rightExonGroup=*rightI;


					//each pair of exons between the two exon groups
					for(NExonGroup::ExonI lei=leftExonGroup->levelExons.begin();lei!=leftExonGroup->levelExons.end();lei++)
					{
						GffEntry::ExonPtr leftExon=*lei;

						if(!leftExon->outJnxs)
							continue;

						for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ji=leftExon->outJnxs->begin();ji!=leftExon->outJnxs->end();ji++)
						{
							string chr;

							GffEntry::ExonPtr rightExon=ji->first;
							GffEntry::Jnx* jnx=ji->second;

							chr=rightExon->chr;

							if(rightExon->exonGroup!=rightExonGroup) //only right exon of the right Exon group.
								continue;

							//////here

							KeyPair<int,int> parentConcatBound(INT_MAX,INT_MIN);

							//cerr<<"parentBound:"<<parentBound<<endl;

							KeyPair<int,int> leftExonBound=leftExon->getBound();
							KeyPair<int,int> rightExonBound=rightExon->getBound();

							cerr<<"leftExon Bound:"<<leftExonBound<<endl;
							cerr<<"rightExon Bound:"<<rightExonBound<<endl;

							set<GffEntry::GBlockPtr> parentBlocksConcat;

							for(NExonGroup::ExonI pei=curExonGroup->levelExons.begin();pei!=curExonGroup->levelExons.end();pei++)
							{
								GffEntry::ExonPtr curParentExon=*pei;
								KeyPair<int,int> curParentBound=curParentExon->getBound();
								if(areOverlapping(leftExonBound,curParentBound) && areOverlapping(rightExonBound,curParentBound)) //cur parent overlaps with the current child pairs
								{
									//push parents all blocks into the parentBlocksConcat;
									parentBlocksConcat.insert(curParentExon->blocks->begin(),curParentExon->blocks->end());
									::expandBoundMutFirst(parentConcatBound,curParentBound);

								}
							}

							if(parentBlocksConcat.size()<1)
							{
								cerr<<"No parent exons overlap with both left and right exons"<<endl;
								continue;
							}

							KeyPair<int,int> leftOverlapBound=::overlapBound(parentConcatBound,leftExonBound);
							KeyPair<int,int> rightOverlapBound=::overlapBound(parentConcatBound,rightExonBound);
							KeyPair<int,int> middleBound(leftOverlapBound.k2,rightOverlapBound.k1);



							::TrafficInfoAS inInfo;
							::TrafficInfoAS exInfo;

							inInfo.bounds.insert(KeyPair<int,int>(leftOverlapBound.k1,rightOverlapBound.k2));
							exInfo.bounds.insert(leftOverlapBound);
							exInfo.bounds.insert(rightOverlapBound);



							inInfo.jnxstring="";
						//	exInfo.jnxstring=StringUtil::str(parentConcatBound.k1)+"-"+StringUtil::str(parentConcatBound.k2);
							exInfo.jnxstring=chr+":"+StringUtil::str(leftOverlapBound.k2)+this->strand+">"+chr+":"+StringUtil::str(rightOverlapBound.k1+1)+this->strand;

						//	inInfo.gsid=StringUtil::str(parentConcatBound.k1)+"-"+StringUtil::str(parentConcatBound.k2);
						//	exInfo.gsid=leftExonGroup->sid+"_"+rightExonGroup->sid;
							inInfo.gsid="";
							exInfo.gsid=leftExonGroup->sid+">"+rightExonGroup->sid;

							//cerr<<"1:"<<leftParentExon->getBound()<<"\t"<<leftOverlapBound<<endl;
							//cerr<<"2:"<<leftParentExon->getBound()<<"\t"<<middleBound<<endl;
							//cerr<<"3:"<<rightParentExon->getBound()<<"\t"<<rightOverlapBound<<endl;
							typedef set<GffEntry::GBlockPtr>::const_iterator sg_i;
							typedef set<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
							typedef pair<sg_i,sg_i> Dsg_i;


							if(this->op.getCount && this->hasDataOut())
							{


								Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocksConcat,leftOverlapBound,true);
								Dsg_i middleBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocksConcat,middleBound,true);
								Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocksConcat,rightOverlapBound,true);





								//cerr<<"a"<<endl;
								KeyPair<int,int> R0Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(leftBlocks,readLength);
								//cerr<<"a2"<<endl;
								KeyPair<int,int> R1Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(middleBlocks,readLength); //both left and right parent exons should span through the middle
								//cerr<<"a3"<<endl;
								KeyPair<int,int> R2Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(rightBlocks,readLength);
								//cerr<<"b"<<endl;
								//////here-->

								KeyPair<int,int> jnxDen=jnx->getDensity(readLength,true);

								exInfo.JRCheckString=StringUtil::str(R1Den.k1);
								exInfo.JPCheckString=StringUtil::str(R1Den.k2);
								exInfo.JPF=R1Den.k2;
								exInfo.JRF=R1Den.k1;

								inInfo.JRCheckString=StringUtil::str(jnxDen.k1);
								inInfo.JPCheckString=StringUtil::str(jnxDen.k2);
								inInfo.JRF=jnxDen.k1;
								inInfo.JPF=jnxDen.k2;



								addDensityVectors(inInfo.flankingExonsFlow,R0Den);
								addDensityVectors(exInfo.flankingExonsFlow,R0Den);
								addDensityVectors(inInfo.flankingExonsFlow,R2Den);
								addDensityVectors(exInfo.flankingExonsFlow,R2Den);

								addDensityVectors(inInfo.middleExonsFlow,R1Den);
								addDensityVectors(exInfo.jnxsFlow,jnxDen);

								addDensityVectors(inInfo.noFlankingInfo,inInfo.middleExonsFlow);
								addDensityVectors(inInfo.withFlankingInfo,inInfo.middleExonsFlow);
								addDensityVectors(inInfo.withFlankingInfo,inInfo.flankingExonsFlow);

								addDensityVectors(exInfo.noFlankingInfo,exInfo.jnxsFlow);
								addDensityVectors(exInfo.withFlankingInfo,exInfo.jnxsFlow);
								addDensityVectors(exInfo.withFlankingInfo,exInfo.flankingExonsFlow);

							}


							this->outData("RI",inInfo,exInfo);








						}
					}
				}
			}
		}
	}
};
*/

class SE_SplidarGraph_J : public SplidarGraph {
public:
	inline SE_SplidarGraph_J(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NExonGroup::NExonGroupPtr _exongroup, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength) :
		SplidarGraph(_op,_fout, _foutSeq, _raf, _exongroup, _locus, _locusName, _checkLocusName, 2, _readLength) {

	}

	inline SE_SplidarGraph_J(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NExonGroup::IntExonPtrMMapDI _range, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength) :
		SplidarGraph(_op,_fout, _foutSeq, _raf, _range, _locus, _locusName, _checkLocusName, 2, _readLength) {

	}

	void endGraph() {
		//cerr<<"end graph SE"<<endl;

		string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		//dnCoord->thread
		map<int,multimap<int,ExonPostItem*>* > oneTwoThread;
		//cerr<<"a";

		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			//just want the twos;

			ExonPostRow::I epri2=epr->find(2);
			ExonPostRow::I epri1=epr->find(1);

			if(epri1==epr->end() && epri2==epr->end())
			{
				//no ones no twos, ignore; //can come other round, so use and instead of or
				continue;
			}

			//cerr<<"b"<<endl;
			//now register
			int dnCoord=terminalExon->getStart1();

			map<int,multimap<int,ExonPostItem*>* >::iterator miv=oneTwoThread.find(dnCoord);

			multimap<int,ExonPostItem*>* epiV;
			if(miv==oneTwoThread.end())
			{
				epiV=new multimap<int,ExonPostItem*>;
				oneTwoThread.insert(map<int,multimap<int,ExonPostItem*>* >::value_type(dnCoord,epiV));
			}
			else
				epiV=miv->second;

			//cerr<<"c"<<endl;

			if(epri1!=epr->end())
			{
				ExonPost* ep1=epri1->second;
				for(ExonPost::I epi=ep1->begin();epi!=ep1->end();epi++)
				{
					epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}

			if(epri2!=epr->end())
			{
				ExonPost* ep2=epri2->second;
				for(ExonPost::I epi=ep2->begin();epi!=ep2->end();epi++)
				{
				epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}

			//cerr<<"d"<<endl;
		}

		//cerr<<"e"<<endl;
			//here comes the critical step!
			//for each terminal coordinate

			for (map<int,multimap<int,ExonPostItem*>* >::iterator i1=oneTwoThread.begin();i1!=oneTwoThread.end();i1++) {

				int dnCoord=i1->first;
				multimap<int,ExonPostItem*>* s2=i1->second;

				//all twos versus one!
				//cerr<<"f"<<endl;
				typedef multimap<int,ExonPostItem*>::iterator II;
				typedef pair<II,II> DII;

				DII diiTwos=s2->equal_range(2);
				DII diiOnes=s2->equal_range(1);

				for (II x3=diiOnes.first; x3!=diiOnes.second; x3++) {
					for (II y3=diiTwos.first; y3
							!=diiTwos.second; y3++) {

						//cerr<<"g"<<endl;
						ExonPostItem* one_epi=x3->second;
						ExonPostItem* two_epi=y3->second;

						AESPathThread* one_thr=one_epi->exonThread;
						AESPathThread* two_thr=two_epi->exonThread;

						AESPathThread::backtrack_iterator one_right=
								one_thr->rbegin(one_epi->vertexPointer);

						AESPathThread::backtrack_iterator one_left(one_right);
						one_left++;

						AESPathThread::backtrack_iterator two_right=
								two_thr->rbegin(two_epi->vertexPointer);

						AESPathThread::backtrack_iterator two_center(two_right);
						two_center++;

						AESPathThread::backtrack_iterator two_left(two_center);
						two_left++;

						//cerr<<"h"<<endl;
						//now get the bound of the common area
						KeyPair<int,int> leftOverlapBound=::overlapBound(one_left->exon->getBound(),two_left->exon->getBound());
						KeyPair<int,int> rightOverlapBound=::overlapBound(one_right->exon->getBound(),two_right->exon->getBound());

						//cerr<<"i"<<endl;
						::TrafficInfoAS inInfo;
						::TrafficInfoAS exInfo;

						int uCoord=one_left->exon->getEnd1();

						//cerr<<"j"<<endl;
						inInfo.bounds.insert(leftOverlapBound);
						inInfo.bounds.insert(two_center->exon->getBound());
						inInfo.bounds.insert(rightOverlapBound);
						exInfo.bounds.insert(leftOverlapBound);
						exInfo.bounds.insert(rightOverlapBound);

						string chr=two_center->exon->chr;

						inInfo.jnxstring=chr+":"+StringUtil::str(uCoord)+this->strand+">"+StringUtil::str(two_center->exon->getStart1())+","+chr+":"+StringUtil::str(two_center->exon->getEnd1())+strand+">"+chr+":"+StringUtil::str(dnCoord)+this->strand;
						exInfo.jnxstring=chr+":"+StringUtil::str(uCoord)+this->strand+">"+chr+":"+StringUtil::str(dnCoord)+this->strand;


						//shouldn't use relRootString anymore
						exInfo.gsid=one_left->exon->exonGroup->relRootString+">"+one_right->exon->exonGroup->relRootString;
						inInfo.gsid=one_left->exon->exonGroup->relRootString+">"+two_center->exon->exonGroup->sid+">"+one_right->exon->exonGroup->relRootString;

						if(one_left->exon->exonGroup->relRootString=="")
							cerr<<"error: relRootString==empty"<<endl;
						if(one_right->exon->exonGroup->relRootString=="")
							cerr<<"error: relRootString==empty"<<endl;

						//cerr<<"1:"<<leftOverlapBound<<endl;
						//cerr<<"3:"<<rightOverlapBound<<endl;

						typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
						typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
						typedef pair<sg_i,sg_i> Dsg_i;

						//KeyPair<int,int> leftOverlapBound=::overlapBound(one_left->exon->getBound(),two_left->exon->getBound());
						//KeyPair<int,int> rightOverlapBound=::overlapBound(one_right->exon->getBound(),two_right->exon->getBound());

						if(this->op.getCount && this->hasDataOut())
						{
						//cerr<<"k"<<endl;
						Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_left->exon->blocks,leftOverlapBound,true);
						Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_right->exon->blocks,rightOverlapBound,true);
						//cerr<<"l"<<endl;

						//cerr<<"m"<<endl;
						KeyPair<int,int> leftDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(leftBlocks,readLength);
						//cerr<<"m2"<<endl;
						KeyPair<int,int> middleDen=two_center->exon->getDensity(readLength,true);
						//cerr<<"m3"<<endl;
						KeyPair<int,int> rightDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(rightBlocks,readLength);
						//cerr<<"n"<<endl;
						//////here-->


//						string exonCoordPathString;
//						string onegsidString;
//						string twogsidString;



						     //right flanking exon, and jnx before.
							{
								////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"o"<<endl;
								KeyPair<int,int> jnxDen=one_right->inJnx->getDensity(readLength,true);

								addDensityVectors(exInfo.jnxsFlow,jnxDen);
								addDensityVectors(exInfo.noFlankingInfo,jnxDen);
								addDensityVectors(exInfo.withFlankingInfo,jnxDen);

								addDensityVectors(exInfo.flankingExonsFlow,rightDen);

								//if(includeFlankExon)
								addDensityVectors(exInfo.withFlankingInfo,rightDen);

								exInfo.JPF*=jnxDen.k2;
								exInfo.JRF*=jnxDen.k1;
								exInfo.JRCheckString=StringUtil::str(jnxDen.k1);
								exInfo.JPCheckString=StringUtil::str(jnxDen.k2);

								addDensityVectors(exInfo.flankingExonsFlow,leftDen);

								//if(includeFlankExon)
								addDensityVectors(exInfo.withFlankingInfo,leftDen);
								//cerr<<"p"<<endl;
								////cerr<<"gt"<<"<a1"<<endl;
							}
							{
								//center stuff
								//cerr<<"s"<<endl;
								KeyPair<int,int> jnxDen=two_center->inJnx->getDensity(readLength,true);
								inInfo.JPF*=jnxDen.k2;
								inInfo.JRF*=jnxDen.k1;

								inInfo.JPCheckString=StringUtil::str(jnxDen.k2);
								inInfo.JRCheckString=StringUtil::str(jnxDen.k1);

								addDensityVectors(inInfo.noFlankingInfo,jnxDen);
								addDensityVectors(inInfo.noFlankingInfo,middleDen);
								addDensityVectors(inInfo.withFlankingInfo,jnxDen);
								addDensityVectors(inInfo.withFlankingInfo,middleDen);

								addDensityVectors(inInfo.middleExonsFlow,middleDen);
								addDensityVectors(inInfo.jnxsFlow,jnxDen);
								////cerr<<"gt"<<"<a2"<<endl;
								//cerr<<"t"<<endl;
							}
							{
								////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"q"<<endl;
								KeyPair<int,int> jnxDen=two_right->inJnx->getDensity(readLength,true);

								addDensityVectors(inInfo.jnxsFlow,jnxDen);
								addDensityVectors(inInfo.noFlankingInfo,jnxDen);
								addDensityVectors(inInfo.withFlankingInfo,jnxDen);

								addDensityVectors(inInfo.flankingExonsFlow,rightDen);

								//if(includeFlankExon)
								addDensityVectors(inInfo.withFlankingInfo,rightDen);

								inInfo.JPF*=jnxDen.k2;
								inInfo.JRF*=jnxDen.k1;

								inInfo.JPCheckString+=string(",")+StringUtil::str(jnxDen.k2);
								inInfo.JRCheckString+=string(",")+StringUtil::str(jnxDen.k1);

								addDensityVectors(inInfo.flankingExonsFlow,leftDen);

								//if(includeFlankExon)
								addDensityVectors(inInfo.withFlankingInfo,leftDen);
								//cerr<<"r"<<endl;
							}
						}

						//cerr<<"u"<<endl;
						this->outData("SE",inInfo,exInfo);
						//cerr<<"v"<<endl;
					}
				}

			}

		//now free memory from the registry
			////cerr<<"n"<<endl;
			for(map<int,multimap<int,ExonPostItem*>* >::iterator i=oneTwoThread.begin();i!=oneTwoThread.end();i++)
			{
				delete i->second;
			}

			////cerr<<"o"<<endl;


	}


	void _endGraph() {
		//cerr<<"end graph SE"<<endl;

		string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			SrcExonPostStruct * srcS=epr->getSrcExonStructure();

			//now here comes the easier but critical step!

			for (SrcExonPostStruct::I1 i1=srcS->begin(); i1!=srcS->end(); i1++) {

				GffEntry::ExonPtr srcExon=i1->first;
				SrcExonPostStruct::S2* s2=i1->second;

				SrcExonPostStruct::I2 oneAway=s2->find(1);
				SrcExonPostStruct::I2 twoAway=s2->find(2);

				if (oneAway==s2->end() || twoAway==s2->end()) {
					continue;
				}

				SrcExonPostStruct::S3* ones=oneAway->second;
				SrcExonPostStruct::S3* twos=twoAway->second;

				for (SrcExonPostStruct::I3 x3=ones->begin(); x3!=ones->end(); x3++) {
					for (SrcExonPostStruct::I3 y3=twos->begin(); y3
							!=twos->end(); y3++) {
						ExonPostItem* one_epi=*x3;
						ExonPostItem* two_epi=*y3;
						AESPathThread* one_thr=one_epi->exonThread;
						AESPathThread* two_thr=two_epi->exonThread;
						AESPathThread::backtrack_iterator one_beg=
								one_thr->rbegin(one_epi->vertexPointer);
						AESPathThread::backtrack_iterator one_end1=
								one_thr->rend();
						AESPathThread::backtrack_iterator two_beg=
								two_thr->rbegin(two_epi->vertexPointer);
						AESPathThread::backtrack_iterator two_end1=
								two_thr->rend();

						string exonCoordPathString;
						string onegsidString;
						string twogsidString;

						TrafficInfoAS one_info=getTrafficAS(one_beg,
								one_end1);
						TrafficInfoAS two_info=getTrafficAS(two_beg,
								two_end1);


						this->outData("SE",two_info,one_info);

					}
				}

			}

			delete srcS;

		}
	}

};

class SEn_SplidarGraph : public SplidarGraph {

private:
	int minDistance;
	int maxDistance;
	NESuperGroupPtr startGroup;
public:


	inline SEn_SplidarGraph(Splidar_OpFlag _op,int _minDistance,int _maxDistance,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NESuperGroupPtr  _sgroup, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength) :
				startGroup(_sgroup),SplidarGraph(_op,_fout, _foutSeq, _raf, _sgroup, _locus, _locusName, _checkLocusName, _maxDistance, _readLength), minDistance(_minDistance),maxDistance(_maxDistance) {
		SG_DEBUG2("processing startGroup="<<_sgroup->getID()<<endl);
	}

	void endGraph() {
		SG_DEBUG2("end graph SEn:"<<startGroup->getID()<<endl);

		//*******CHANGE PENDING CONSIDERATION allow setting genome for the hyperlink


		//string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		//string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		//dn LeftSuperGroup ->thread
		map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* > oneNThread; //dsExon->SuperGroup => [ length-of-SEn-inclusionPath => ExonPostItem (a record that points to the exon thread and the location freeze of that thread) ]
		//cerr<<"a";

		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			//just want the twos;

			//ExonPostRow::I epri2=epr->find(2);
			ExonPostRow::I epriMaxEnd=epr->upper_bound(maxDistance);
			ExonPostRow::I epriMinStart=epr->lower_bound(minDistance);

			SG_DEBUG2("from "<<startGroup->getID()<<" reaching:"<< terminalExon->egsid<<" of "<<terminalExon->exonGroup->sid<<endl);




			ExonPostRow::I epri1=epr->find(1);

			if(epri1==epr->end() && epriMinStart==epr->end())
			{
				//no ones and no min-skipped-count, ignore; //???two can come other round
				continue;
			}

			////cerr<<"b"<<endl;
			//now register

			//int dnCoord=terminalExon->getStart1(); //remove 6/6/2009


			NESuperGroupPtr  dnSGroup=terminalExon->exonGroup->leftSuperGroup;

			SG_DEBUG2("jumping "<<startGroup->getID()<<" to "<<dnSGroup->getID()<<endl);

			map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator miv=oneNThread.find(dnSGroup);

			multimap<int,ExonPostItem*>* epiV;
			if(miv==oneNThread.end())
			{
				epiV=new multimap<int,ExonPostItem*>;
				oneNThread.insert(map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::value_type(dnSGroup,epiV));
			}
			else
				epiV=miv->second;

			////cerr<<"c"<<endl;


			//add all the one-degree items
			if(epri1!=epr->end())
			{
				ExonPost* ep1=epri1->second;
				for(ExonPost::I epi=ep1->begin();epi!=ep1->end();epi++)
				{

					epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}


			//add all the N-degree items
			for(ExonPostRow::I eprin=epriMinStart;eprin!=epriMaxEnd;eprin++)
			{
				ExonPost* epn=eprin->second;
				for(ExonPost::I epi=epn->begin();epi!=epn->end();epi++)
				{
					//cerr<<"adding time "<<epi->time<<endl;
					epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}

			////cerr<<"d"<<endl;
		}

		////cerr<<"e"<<endl;
			//here comes the critical step!
			//for each terminal coordinate

			for (map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator i1=oneNThread.begin();i1!=oneNThread.end();i1++) {

				//int dnCoord=i1->first;
				multimap<int,ExonPostItem*>* s2=i1->second; //s2 is the record for all threads ending in current i1->first,aka dnExonGroup

				//all twos versus one!
				////cerr<<"f"<<endl;
				typedef multimap<int,ExonPostItem*>::iterator II;
				typedef pair<II,II> DII;

				//DII diiTwos=s2->equal_range(2);
				DII diiOnes=s2->equal_range(1);
				DII diiNs;
				diiNs.first=diiOnes.second;
				diiNs.second=s2->end();


				for (II x3=diiOnes.first; x3!=diiOnes.second; x3++) {
					for (II y3=diiNs.first; y3
							!=diiNs.second; y3++) {

						////cerr<<"g"<<endl;



						ExonPostItem* one_epi=x3->second;
						ExonPostItem* N_epi=y3->second;

						AESPathThread* one_thr=one_epi->exonThread;
						AESPathThread* N_thr=N_epi->exonThread;

						AESPathThread::backtrack_iterator one_right=
								one_thr->rbegin(one_epi->vertexPointer);

						AESPathThread::backtrack_iterator one_left(one_right);
						one_left++;

						AESPathThread::backtrack_iterator N_right=N_thr->rbegin(N_epi->vertexPointer);

						AESPathThread::backtrack_iterator N_i(N_right);
						AESPathThread::backtrack_iterator N_end=N_thr->rend();

						AESPathThread::backtrack_iterator N_left(N_right); //gonna change anyway


						deque<InEdgeAndVertex> ievs_N;

						::TrafficInfoAS inInfo;
						::TrafficInfoAS exInfo;
						//start from left of N_right to left end


						string inMiddleJnxString;
						string inMiddleEgString;

						string chr=one_left->exon->chr;

						//setup inclusion bounds and exon strings
						//N_i++ initially to skip including the rightmost exon
						//N_i!=N_end as terminal condition to avoid including the leftmost exon
						//from left to right
						//iev N_i is iev contains edge to dst node and the dst node

						for(N_i++;N_i!=N_end;N_i++)
						{
							if(N_i.hasAtLeastOneMore())
							{
								//not the last one
								const InEdgeAndVertex& iev=*N_i;
								inMiddleJnxString=chr+":"+StringUtil::str(iev.exon->getStart1())+this->strand+","+chr+":"+StringUtil::str(iev.exon->getEnd1())+":"+this->strand+">"+inMiddleJnxString;
								inMiddleEgString=iev.exon->exonGroup->sid+">"+inMiddleEgString;
								inInfo.bounds.insert(iev.exon->getBound());
								ievs_N.push_front(iev);
							}
							else
							{
								//the last one
								N_left=N_i;
							}
						}


						if(!areOverlapping(one_left->exon->getBound(),N_left->exon->getBound()))
							continue;

						//*****CHANGE PENDING CONSIDERATION: Require one_left->exon and N_left->exon to be in the same rightSuperGroup*****
						//
						 if (one_left->exon->exonGroup->rightSuperGroup != N_left->exon->exonGroup->rightSuperGroup)
							 continue;
						//
						//////




						if(!areOverlapping(one_right->exon->getBound(),N_right->exon->getBound()))
							continue;




						//cerr<<"h"<<endl;

						//now get the bound of the common area
						KeyPair<int,int> leftOverlapBound=::overlapBound(one_left->exon->getBound(),N_left->exon->getBound());
						KeyPair<int,int> rightOverlapBound=::overlapBound(one_right->exon->getBound(),N_right->exon->getBound());

						//cerr<<"i"<<endl;

						SG_DEBUG2("SEn event:"<<one_left->exon->egsid<<":"<<one_right->exon->egsid<<endl);


						//cerr<<"j"<<endl;
						inInfo.bounds.insert(leftOverlapBound);
						//inInfo.bounds.insert(two_center->exon->getBound());
						inInfo.bounds.insert(rightOverlapBound);
						exInfo.bounds.insert(leftOverlapBound);
						exInfo.bounds.insert(rightOverlapBound);

						int uCoord=leftOverlapBound.k2;
						int dnCoord=rightOverlapBound.k1+1; //convert to 1-base coordinate

						//finish up the jnxstring annotation
						string leftjnxstring=chr+":"+StringUtil::str(uCoord)+this->strand;
						string rightjnxstring=chr+":"+StringUtil::str(dnCoord)+this->strand;
						//inInfo.jnxstring=chr+":"+StringUtil::str(uCoord)+this->strand+">"+StringUtil::str(two_center->exon->getStart1())+","+chr+":"+StringUtil::str(two_center->exon->getEnd1())+strand+">"+chr+":"+StringUtil::str(dnCoord)+this->strand;
						exInfo.jnxstring=leftjnxstring+">"+rightjnxstring;
						inInfo.jnxstring=leftjnxstring+">"+inMiddleJnxString+rightjnxstring;


						//finish up the exongroupstring annotation
						string leftegstring=one_left->exon->exonGroup->rightSuperGroup->getID();
						string rightegstring=one_right->exon->exonGroup->leftSuperGroup->getID();
						exInfo.gsid=leftegstring+">"+rightegstring;
						inInfo.gsid=leftegstring+">"+inMiddleEgString+rightegstring;

						//inInfo.gsid=one_left->exon->exonGroup->relRootString+">"+two_center->exon->exonGroup->sid+">"+one_right->exon->exonGroup->relRootString;


						///**** CHANGE PENDING: REMOVAL: some root exongroup can have level exons which have relRootString==""
						//if(one_left->exon->exonGroup->relRootString=="")
						//	cerr<<"error: relRootString==empty"<<endl;
						//if(one_right->exon->exonGroup->relRootString=="")
						//	cerr<<"error: relRootString==empty"<<endl;
						//
						////*****




						//cerr<<"1:"<<leftOverlapBound<<endl;
						//cerr<<"3:"<<rightOverlapBound<<endl;

						typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
						typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
						typedef pair<sg_i,sg_i> Dsg_i;

						if(this->op.getCount && this->hasDataOut())
						{
						//cerr<<"k"<<endl;
						Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_left->exon->blocks,leftOverlapBound,true);
						Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*one_right->exon->blocks,rightOverlapBound,true);
						//cerr<<"l"<<endl;

						//cerr<<"m"<<endl;
						KeyPair<int,int> leftDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(leftBlocks,readLength);
						//cerr<<"m2"<<endl;

						//cerr<<"m3"<<endl;
						KeyPair<int,int> rightDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(rightBlocks,readLength);
						//cerr<<"n"<<endl;
						//////here-->


//						string exonCoordPathString;
//						string onegsidString;
//						string twogsidString;
						SG_DEBUG2(one_left->exon->egsid<<":left Den:"<<leftDen<<endl);
						SG_DEBUG2(one_right->exon->egsid<<":right Den:"<<rightDen<<endl);


						     //right flanking exon, and jnx before.
							{
								//////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"o"<<endl;
								KeyPair<int,int> jnxDen=one_right->inJnx->getDensity(readLength,true);  //jxnDen= [ ]V*[ ]

								addDensityVectors(exInfo.jnxsFlow,jnxDen);
								addDensityVectors(exInfo.noFlankingInfo,jnxDen);
								addDensityVectors(exInfo.withFlankingInfo,jnxDen);

								addDensityVectors(exInfo.flankingExonsFlow,rightDen);  //rightDen= [ ]V[*] = []^[]^[]^[*]

								//if(includeFlankExon)
								addDensityVectors(exInfo.withFlankingInfo,rightDen);

								exInfo.JPF*=jnxDen.k2;
								exInfo.JRF*=jnxDen.k1;
								exInfo.JRCheckString=StringUtil::str(jnxDen.k1);
								exInfo.JPCheckString=StringUtil::str(jnxDen.k2);

								addDensityVectors(exInfo.flankingExonsFlow,leftDen); //leftDen= [*]V[ ] = [*]^[]^[]^[]

								//if(includeFlankExon)
								addDensityVectors(exInfo.withFlankingInfo,leftDen);
								//cerr<<"p"<<endl;
								////cerr<<"gt"<<"<a1"<<endl;
							}
							{
								//center stuff
								//cerr<<"s"<<endl;


								//for each of the middle exons and jnx
								for(deque<InEdgeAndVertex>::iterator i=ievs_N.begin();i!=ievs_N.end();i++)
								{

									////*** PENDING CONSIDERING: REMOVAL
									//if(!i->inJnx)
									//{
									//	cerr<<"you are shitting me with in Jnx == NULL"<<endl;
									//}
									/////******



									KeyPair<int,int> middleDen=i->exon->getDensity(readLength,true); //middleDen= [ ](^[ ]^[*])^[ ]
									KeyPair<int,int> jnxDen=i->inJnx->getDensity(readLength,true); //jnxDen= [ ](^[ ]^*[ ])^[ ]

									SG_DEBUG2(i->exon->egsid<<":"<<middleDen<<"\t");

									inInfo.JPF*=jnxDen.k2;
									inInfo.JRF*=jnxDen.k1;

									inInfo.JPCheckString+=StringUtil::str(jnxDen.k2)+",";
									inInfo.JRCheckString+=StringUtil::str(jnxDen.k1)+",";

									addDensityVectors(inInfo.noFlankingInfo,jnxDen);
									addDensityVectors(inInfo.noFlankingInfo,middleDen);
									addDensityVectors(inInfo.withFlankingInfo,jnxDen);
									addDensityVectors(inInfo.withFlankingInfo,middleDen);

									addDensityVectors(inInfo.middleExonsFlow,middleDen);
									addDensityVectors(inInfo.jnxsFlow,jnxDen);
								}
								////cerr<<"gt"<<"<a2"<<endl;
								//cerr<<"t"<<endl;
							}

							SG_DEBUG2(endl);


							{
								////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"q"<<endl;
								KeyPair<int,int> jnxDen=N_right->inJnx->getDensity(readLength,true); //jnxDen= [ ]^[ ]^[ ]^*[ ]

								addDensityVectors(inInfo.jnxsFlow,jnxDen);
								addDensityVectors(inInfo.noFlankingInfo,jnxDen);
								addDensityVectors(inInfo.withFlankingInfo,jnxDen);

								addDensityVectors(inInfo.flankingExonsFlow,rightDen);

								//if(includeFlankExon)
								addDensityVectors(inInfo.withFlankingInfo,rightDen);

								inInfo.JPF*=jnxDen.k2;
								inInfo.JRF*=jnxDen.k1;

								inInfo.JPCheckString+=StringUtil::str(jnxDen.k2);
								inInfo.JRCheckString+=StringUtil::str(jnxDen.k1);

								addDensityVectors(inInfo.flankingExonsFlow,leftDen);

								//if(includeFlankExon)
								addDensityVectors(inInfo.withFlankingInfo,leftDen);
								//cerr<<"r"<<endl;
							}
						}

						//cerr<<"u"<<endl;
						this->outData("SE",inInfo,exInfo);
						//cerr<<"v"<<endl;
					}
				}

			}

		//now free memory from the registry
			//cerr<<"n"<<endl;
			for(map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator i=oneNThread.begin();i!=oneNThread.end();i++)
			{
				delete i->second;
			}

			//cerr<<"o"<<endl;


	}




};


#endif /* SPLIDARGRAPH_BACKUP_H_ */
