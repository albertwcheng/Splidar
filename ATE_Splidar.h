/*
 * ATE_SplidarGraph.h
 *
 *  Created on: Jan 19, 2010
 *      Author: awcheng
 */

#ifndef ATE_SPLIDARGRAPH_H_
#define ATE_SPLIDARGRAPH_H_


class ATE_Splidar_REGJ: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	inline ATE_Splidar_REGJ(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry):op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{
		//for each transcript
		//add to second exon into secondExonSet

		//map anchor left/right supergroup :> set { fixedCoord => variantExon => commonExon /*common is just to avoid the same exact set of exons*/ }

		//typedef map<NESuperGroupPtr,set<ExonPairWJnx>* > ATEHashMap;
		class ATEHashMap: public map<NESuperGroupPtr,set<ExonPairWJnx>* >
		{
		public:
			~ATEHashMap()
			{
				for(iterator i=this->begin();i!=this->end();i++)
				{
					delete i->second;
				}
			}
		};

		typedef ATEHashMap::iterator ATEHashMapI;
		typedef ATEHashMap::value_type ATEHashMapV;

		ATEHashMap twoMapAG5E; //genomic alternative 5' exon
		ATEHashMap twoMapAG3E; //genomic alternative 3' exon
		typedef set<ExonPairWJnx>::iterator ExonPairWJnxI;
		//typedef pair<ExonPairWJnxI,ExonPairWJnxI> ExonPairWJnxDI;
		typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
		typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
		typedef pair<sg_i,sg_i> Dsg_i;

		for(vector<GffEntry*>::iterator it=locus->transcripts.begin();it!=locus->transcripts.end();it++)
		{
			GffEntry*transcript=*it;
			int exonCount=transcript->exonCount;

			//if single-exon transcript, ignore

			if(exonCount<2)
				continue;

			GffEntry::ExonPtr secondExon=transcript->exons[1];
			GffEntry::ExonPtr firstExon=transcript->exons[0];

			GffEntry::ExonPtr lastExon=transcript->exons[exonCount-1]; //last
			GffEntry::ExonPtr meiyiExon=transcript->exons[exonCount-2]; //secondLast

			//no need to check usageFlag or locus name consistency between exons, because it is by transcript

			if(firstExon->outJnxs)
			{
				//outJnx is [ dstExon :> jnx ]

				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=firstExon->outJnxs->find(secondExon);
				if(ise==firstExon->outJnxs->end())
				{
					cerr<<"strange error: first and second exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					die("");
					continue;
				}

				GffEntry::Jnx* jnx=ise->second;
				//the common exon is the second exon;
				//the variant exon is the first exon;

				set<ExonPairWJnx> *epwjset;

				ATEHashMapI ATEHashMapi=twoMapAG5E.find(secondExon->exonGroup->leftSuperGroup);

				if(ATEHashMapi==twoMapAG5E.end())
				{
					epwjset=new set<ExonPairWJnx>;
					twoMapAG5E.insert(ATEHashMapV(secondExon->exonGroup->leftSuperGroup,epwjset));
				}
				else
					epwjset=ATEHashMapi->second;

				epwjset->insert(ExonPairWJnx(/*dummy*/0,firstExon,secondExon,jnx));

			}else
			{
				cerr<<"strange firstExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
				die("");
			}




			if(meiyiExon->outJnxs)
			{
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ise=meiyiExon->outJnxs->find(lastExon);
				if(ise==meiyiExon->outJnxs->end())
				{
					cerr<<"strange error: meiyi and last exon are not connected by a jnx:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
					die("");
					continue;
				}

				//the common exon is the meiyi exon;
				//the variant exon is the last exon;
				GffEntry::Jnx* jnx=ise->second;

				set<ExonPairWJnx> *epwjset;

				ATEHashMapI ATEHashMapi=twoMapAG3E.find(meiyiExon->exonGroup->rightSuperGroup);

				if(ATEHashMapi==twoMapAG3E.end())
				{
					epwjset=new set<ExonPairWJnx>;
					twoMapAG3E.insert(ATEHashMapV(meiyiExon->exonGroup->rightSuperGroup,epwjset));
				}
				else
					epwjset=ATEHashMapi->second;

				epwjset->insert(ExonPairWJnx(/*dummy*/0,lastExon,meiyiExon,jnx));

			}else
			{
				cerr<<"strange meiyiExon.outjnx == Null:"<<locusName<<":"<<firstExon->exonID<<" to "<<secondExon->exonID<<" with bound "<<firstExon->getBound()<<" to "<<secondExon->getBound()<<endl;
				die("");
			}




		}

		////HERE

		//process g5'
		//cerr<<"size of twoMapAG5E="<<twoMapAG5E.size()<<endl;

		//set<ExonPairWJnx>::iterator i=twoMapAG5E.begin(); ///
		//ExonPairWJnxDI curRange;
		//int prevB;

		//while(i!=twoMapAG5E.end()) ///

		for(ATEHashMapI i=twoMapAG5E.begin();i!=twoMapAG5E.end();i++)
		{
			//prevB=i->fixedCoord;
		//	cerr<<i->fixedCoord<<"\t";
			//curRange.first=i;
			//curRange.second=(++i);

			NESuperGroupPtr sgroup=i->first;
			set<ExonPairWJnx> * epwjset=i->second;

			SG_DEBUG3("processing 5' ATE with 3' exon super group "<<sgroup->getID());



			//sorted by variable exon

			for(ExonPairWJnxI x=epwjset->begin();x!=epwjset->end();x++)
				{
					ExonPairWJnxI y=x;
					y++;
					for(;y!=epwjset->end();y++)
					{
					//	cerr<<"in business:";
						//just need to swap here!


						//  P    D
						//[ x ]------[]
						//       [y]-[]

						ExonPairWJnx proximaliev=*y;
						ExonPairWJnx distaliev=*x;

						GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
						KeyPair<int,int> proximalBound=proximalExon->getBound();
						GffEntry::Jnx* proximalJnx=proximaliev.jnx;

						GffEntry::ExonPtr distalExon=distaliev.variantExon;
						KeyPair<int,int> distalBound=distalExon->getBound();
						GffEntry::Jnx* distalJnx=distaliev.jnx;


						if(op.commonFlankingCobound)
						{
							if(proximaliev.commonExon->getBound().k1 != distaliev.commonExon->getBound().k1)
								continue;
						}

						if(::areOverlapping(proximalBound,distalBound))
						{
					//		cerr<<"overlapping "<<proximalBound<<" vs "<<distalBound<<endl;
							continue; //the two variable exons overlap, ignore
							//this can happen,e.g.,
							// [    ]-----[]
							//   [   ]----[]
						}

					//	cerr<<"ok! continue!"<<endl;

						TrafficInfoAS xinInfo;
						TrafficInfoAS yexInfo;


						//jnxstring
						//( ]-----[      )
						//    ( ]---[   )

						xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximalExon->getEnd1())+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximaliev.commonExon->getStart1())+":"+this->strand;
						yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distalExon->getEnd1())+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distaliev.commonExon->getStart1())+":"+this->strand;

						//the commonExon's leftsupergroups are the same anyawy
						string xinCommonEG;
						string yexCommonEG;

						if(op.egStringOutputCommonFlankingAsCoord)
						{
							xinCommonEG=proximaliev.commonExon->chr+":"+StringUtil::str(proximaliev.commonExon->getBound().k1+1)+":"+this->strand;
							yexCommonEG=distaliev.commonExon->chr+":"+StringUtil::str(distaliev.commonExon->getBound().k1+1)+":"+this->strand;
						}
						else
						{
							xinCommonEG=proximaliev.commonExon->exonGroup->leftSuperGroup->getID();
							yexCommonEG=distaliev.commonExon->exonGroup->leftSuperGroup->getID();
						}
							xinInfo.gsid=proximalExon->exonGroup->rightSuperGroup->getID()+">"+xinCommonEG;
							yexInfo.gsid=distalExon->exonGroup->rightSuperGroup->getID()+">"+yexCommonEG;

						//xinInfo.gsid=proximalExon->exonGroup->sid+">"+proximaliev.commonExon->exonGroup->relRootString;
						//yexInfo.gsid=distalExon->exonGroup->sid+">"+distaliev.commonExon->exonGroup->relRootString;

						KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());



						xinInfo.bounds.insert(proximalExon->getBound());
						xinInfo.bounds.insert(commonBound);

						yexInfo.bounds.insert(distalExon->getBound());
						yexInfo.bounds.insert(commonBound);


						if(op.avoidIsoBound)
						{

							string isoboundKeyL=xinInfo.getCoordPathString(NULL);
							string isoboundKeyR=yexInfo.getCoordPathString(NULL);
							string isoboundKey=(isoboundKeyL<isoboundKeyR)?(isoboundKeyL+"/"+isoboundKeyR):(isoboundKeyR+"/"+isoboundKeyL);

							if(isoBoundRegistry.find(isoboundKey)!=isoBoundRegistry.end())
							{
								continue;
							}
							else
							{

								isoBoundRegistry.insert(isoboundKey);
							}



						}

						SplidarOutputExtraData extraDat;
						extraDat.commonExonIsCobound=true;

						//              CC
						//( ]-----------[  )
						//    (*]----[xx  )
						KeyPair<int,int> xinCommonExonSpecBound(proximaliev.commonExon->getBound().k1,commonBound.k1+1);
						xinInfo.specBounds.insert(proximalExon->getBound());
						if(len01(xinCommonExonSpecBound)>0){
							xinInfo.specBounds.insert(xinCommonExonSpecBound);
							extraDat.commonExonIsCobound=false;
						}
						//          CCCC
						//(*]-----[y     )
						//    ( ]---[   )
						KeyPair<int,int> yexCommonExonSpecBound(distaliev.commonExon->getBound().k1,commonBound.k1+1);
						yexInfo.specBounds.insert(distalExon->getBound());
						if(len01(yexCommonExonSpecBound)>0){
							yexInfo.specBounds.insert(yexCommonExonSpecBound);
							extraDat.commonExonIsCobound=false;
						}


						xinInfo.boundsExonic.insert(proximalExon->getBound());
						xinInfo.boundsExonic.insert(proximaliev.commonExon->getBound());

						yexInfo.boundsExonic.insert(distalExon->getBound());
						yexInfo.boundsExonic.insert(distaliev.commonExon->getBound());



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


						this->outData((this->strand==GffEntry::FORWARD)?"AFE":"ALE",xinInfo,yexInfo,extraDat);

					}
				}

		}




		//process g3'




		for(ATEHashMapI i=twoMapAG3E.begin();i!=twoMapAG3E.end();i++)
		{


			NESuperGroupPtr sgroup=i->first;
			set<ExonPairWJnx> * epwjset=i->second;

			SG_DEBUG3("processing 3' ATE with 5' exon super group "<<sgroup->getID());




			for(ExonPairWJnxI x=epwjset->begin();x!=epwjset->end();x++)
				{
					ExonPairWJnxI y=x;
					y++;
					for(;y!=epwjset->end();y++)
					{
						//just need to swap here!

						// [  ]---- [x ]
						//  [ ]----------[y ]
						//


						ExonPairWJnx proximaliev=*x;
						ExonPairWJnx distaliev=*y;

						GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
						KeyPair<int,int> proximalBound=proximalExon->getBound();
						GffEntry::Jnx* proximalJnx=proximaliev.jnx;

						GffEntry::ExonPtr distalExon=distaliev.variantExon;
						KeyPair<int,int> distalBound=distalExon->getBound();
						GffEntry::Jnx* distalJnx=distaliev.jnx;

						if(op.commonFlankingCobound)
						{
							if(proximaliev.commonExon->getBound().k2 != distaliev.commonExon->getBound().k2)
								continue;
						}

						if(::areOverlapping(proximalBound,distalBound))
						{
							continue; //the two variable exons overlap, ignore
						}

						TrafficInfoAS xinInfo;
						TrafficInfoAS yexInfo;

						//jnxstring
						//   ]-----[  )
						//  ]---------------[ )
						//


						xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximaliev.commonExon->getEnd1())+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximalExon->getStart1())+":"+this->strand;
						yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distaliev.commonExon->getEnd1())+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distalExon->getStart1())+":"+this->strand;

						//egstring
						//  RSG-----LSG
						//  RSG------------LSG
						//
						// the common exon's right supergroup are the same anyway! (by Hashing)

						string xinCommonEG;
						string yexCommonEG;

						if(op.egStringOutputCommonFlankingAsCoord)
						{
							xinCommonEG=proximaliev.commonExon->chr+":"+StringUtil::str(proximaliev.commonExon->getBound().k2)+":"+this->strand;
							yexCommonEG=distaliev.commonExon->chr+":"+StringUtil::str(distaliev.commonExon->getBound().k2)+":"+this->strand;
						}
						else
						{
							xinCommonEG=proximaliev.commonExon->exonGroup->leftSuperGroup->getID();
							yexCommonEG=distaliev.commonExon->exonGroup->leftSuperGroup->getID();
						}


						xinInfo.gsid=xinCommonEG+">"+proximalExon->exonGroup->leftSuperGroup->getID(); //flip!
						yexInfo.gsid=yexCommonEG+">"+distalExon->exonGroup->leftSuperGroup->getID();


						KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());

						Dsg_i commonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.commonExon->blocks,commonBound,true);
						KeyPair<int,int> commonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(commonBlocks,readLength);

						xinInfo.bounds.insert(commonBound);
						xinInfo.bounds.insert(proximalExon->getBound());

						yexInfo.bounds.insert(commonBound);
						yexInfo.bounds.insert(distalExon->getBound());

						if(op.avoidIsoBound)
						{

							string isoboundKeyL=xinInfo.getCoordPathString(NULL);
							string isoboundKeyR=yexInfo.getCoordPathString(NULL);
							string isoboundKey=(isoboundKeyL<isoboundKeyR)?(isoboundKeyL+"/"+isoboundKeyR):(isoboundKeyR+"/"+isoboundKeyL);

							if(isoBoundRegistry.find(isoboundKey)!=isoBoundRegistry.end())
							{
								continue;
							}
							else
							{

								isoBoundRegistry.insert(isoboundKey);
							}



						}

						SplidarOutputExtraData extraDat;
						extraDat.commonExonIsCobound=true;

						xinInfo.boundsExonic.insert(proximaliev.commonExon->getBound());
						xinInfo.boundsExonic.insert(proximalExon->getBound());

						yexInfo.boundsExonic.insert(distaliev.commonExon->getBound());
						yexInfo.boundsExonic.insert(distalExon->getBound());

						//CCC
						//   x]-----[***)
						//  ]---------------[ )
						//
						KeyPair<int,int> xinCommonExonSpecBound(commonBound.k2,proximaliev.commonExon->getBound().k2);
						xinInfo.specBounds.insert(proximalExon->getBound());
						if(len01(xinCommonExonSpecBound)>0){
							xinInfo.specBounds.insert(xinCommonExonSpecBound);
							extraDat.commonExonIsCobound=false;
						}
						//CCC
						//  ]-------[ )
						//   yyy]---------------[**)
						//
						KeyPair<int,int> yexCommonExonSpecBound(commonBound.k2,distaliev.commonExon->getBound().k2);
						yexInfo.specBounds.insert(distalExon->getBound());
						if(len01(yexCommonExonSpecBound)>0){
							yexInfo.specBounds.insert(yexCommonExonSpecBound);
							extraDat.commonExonIsCobound=false;
						}

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

						this->outData((this->strand==GffEntry::REVERSE)?"AFE":"ALE",xinInfo,yexInfo,extraDat);

					}
				}

		}




	}


};
#endif /* ATE_SPLIDARGRAPH_H_ */
