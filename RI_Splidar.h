/*
 * RI_Splidar.h
 *
 *  Created on: Jan 19, 2010
 *      Author: awcheng
 */

#ifndef RI_SPLIDAR_H_
#define RI_SPLIDAR_H_




class RI_Splidar_REG: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	inline RI_Splidar_REG(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry,bool parentsNeedToBeSplicedEitherSide,bool parentsNeedToBeSplicedOnBothSides, bool parentsNeedToBeSplicedBothSidesInSameTranscript,bool retainedIntronsHaveToBeFreeOfExons,bool useUnambiguousRegionsOfRetainedIntron): op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{

		//^[       ]^
		//  [ ]^[]

		string chr=_locus->chr;

		//use the visited flag of exon to remember whether it is a good parent

		//now reset visited flag to false, everyone is not a good parent


		for(set<GffEntry::ExonPtr>::iterator exonI=_locus->exonSet.begin();exonI!=_locus->exonSet.end();exonI++)
		{

			GffEntry::ExonPtr exon=*exonI;


			if(parentsNeedToBeSplicedBothSidesInSameTranscript)
			{
							exon->visited=false;
			}
			else if(parentsNeedToBeSplicedOnBothSides)
			{
				if(exon->hasBothInAndOut())
				{
					exon->visited=true;
				}
				else
				{
					exon->visited=false;
				}
			}else if(parentsNeedToBeSplicedEitherSide)
			{
				if(exon->hasInOrOut())
				{
					exon->visited=true;
				}
				else
				{
					exon->visited=false;
				}
			}
			else
			{
				exon->visited=true;
			}

		}



		//now for each transcript:
		if(parentsNeedToBeSplicedBothSidesInSameTranscript)
		{
			for(vector<GffEntry*>::iterator transcriptI=_locus->transcripts.begin();transcriptI!=locus->transcripts.end();transcriptI++)
			{
				GffEntry* transcript=*transcriptI;


				//[ ]-[ T ]-[T ]-[T ]-[  ]
				for(int rank=1;rank<transcript->exonCount-1;rank++)
				{
					transcript->exons[rank]->visited=true;
				}
			}

		}



		for(set<GffEntry::ExonPtr>::iterator leftExonI=_locus->exonSet.begin();leftExonI!=_locus->exonSet.end();leftExonI++)
		{
			GffEntry::ExonPtr leftExon=*leftExonI;
			if(!leftExon->outJnxs)
				continue;

			for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ievI=leftExon->outJnxs->begin();ievI!=leftExon->outJnxs->end();ievI++)
			{
				GffEntry::ExonPtr rightExon=ievI->first;

				if(leftExon->usageFlag!=rightExon->usageFlag)
				{
					continue;
				}

				GffEntry::Jnx* jnx=ievI->second;


				KeyPair<int,int> leftExonBound=leftExon->getBound();
				KeyPair<int,int> rightExonBound=rightExon->getBound();

				KeyPair<int,int> middleBound(leftExonBound.k2,rightExonBound.k1);

				//retainedIntronsHaveToBeFreeOfExons: discard event if
				//[              ]
				//[ ]----------[]
				//     [  ]
				//


				//useUmabiguousRegionsOfRetainedIntron (use only marked region for inclusion reads and pos
				//
				//[    xxxxxx    ]    }inc
				//[ ]---------[  ]   }exc
				// [  ]--------------
				//           [   ]---




				BoundSet middleBoundSet(middleBound);
				set<GffEntry::ExonPtr> parentExons;

				if(retainedIntronsHaveToBeFreeOfExons || useUnambiguousRegionsOfRetainedIntron)
				{
					//KeyPair<int,int> newMiddleBound(middleBound);
					bool violateRetainedIntronsHaveToBeFreeOfExons=false;

					for(set<GffEntry::ExonPtr>::iterator internalExonI=_locus->exonSet.begin();internalExonI!=_locus->exonSet.end();internalExonI++)
					{
						GffEntry::ExonPtr internalExon=*internalExonI;
						if( internalExon==leftExon || internalExon==rightExon)
							continue;

						if(areOverlapping(leftExon->getBound(),internalExon->getBound()) && areOverlapping(rightExon->getBound(),internalExon->getBound()))
						{
							//a potential parent exon, immune;
							/*if(parentsNeedToBeSpliced)
							{
								//whether parent is spliced? if not, and need to check, ignore
								if(internalExon->hasBothInAndOut()) //chris says need both end spliced... so instead of internalExon->hasInOrOut()
								{

									parentExons.insert(internalExon);
								}
							}*/


							if(internalExon->visited) //this visited flag is used to remember whether it is a good candidate of parent by the above criteria
							{
								parentExons.insert(internalExon);
							}


							continue;
						}


						if(internalExon->getStart1()<=leftExon->getEnd1() && internalExon->getEnd1()>leftExon->getEnd1())
						{
							//  [ ]------[ ]
							//   [iE  ]
							// [ iE  ]
							if(useUnambiguousRegionsOfRetainedIntron)
							{
								middleBoundSet.subtractInPlaceBy(BoundSet(internalExon->getBound()));
							}


						}
						else if(internalExon->getStart1()>leftExon->getEnd1() && internalExon->getEnd1()<rightExon->getStart1())
						{
							// [ ]-----[]
							//     []
							//
							//
							if(retainedIntronsHaveToBeFreeOfExons)
							{
								violateRetainedIntronsHaveToBeFreeOfExons=true;
								break;
							}

							if(useUnambiguousRegionsOfRetainedIntron)
							{
								middleBoundSet.subtractInPlaceBy(BoundSet(internalExon->getBound()));
							}

						}else if(internalExon->getStart1()<rightExon->getStart1() && internalExon->getStart1()>=rightExon->getStart1())
						{

							// [ ]-----[    ]
							//       [  iE   ]
							//        [ iE]
							if(useUnambiguousRegionsOfRetainedIntron)
							{
								middleBoundSet.subtractInPlaceBy(BoundSet(internalExon->getBound()));
							}

						}else
						{
							//not an internal exon at all!~
						}






					}

					//the bound masks essentially mask the whole middle bound~!~! ignore
					if(middleBoundSet.isEmpty())
					{
						cerr<<"middle bound set is empty"<<this->locus->chr<<":"<<(middleBound.k1+1)<<"-"<<middleBound.k2<<endl;
						continue;
					}

					if(violateRetainedIntronsHaveToBeFreeOfExons)
						continue;


				}



				//now we have a pair of exon, find the parent now
				for(set<GffEntry::ExonPtr>::iterator parentExonI=parentExons.begin();parentExonI!=parentExons.end();parentExonI++)
				{


					GffEntry::ExonPtr parentExon=*parentExonI;



					//now these are good set of parent exon and left,right exons and the jnx connecting left and right exons


					///[              ]
					///  [  ]-----[  ]

					///[              ]     parentBound
					///  [  ]               leftOverlapBound
					///           [  ]      rightOverlapBound
					///     [     ]         middleBound (the retained intron)
					///  [           ]      eventBound


					KeyPair<int,int> parentBound=parentExon->getBound();

					KeyPair<int,int> leftOverlapBound=::overlapBound(parentBound,leftExonBound);
					KeyPair<int,int> rightOverlapBound=::overlapBound(parentBound,rightExonBound);

					KeyPair<int,int> eventBound(leftOverlapBound.k1,rightOverlapBound.k2);






					::TrafficInfoAS inInfo;
					::TrafficInfoAS exInfo;

					inInfo.bounds.insert(eventBound);
					exInfo.bounds.insert(leftOverlapBound);
					exInfo.bounds.insert(rightOverlapBound);


					if(op.avoidIsoBound)
					{

						string isoboundKeyL=inInfo.getCoordPathString(NULL);
						string isoboundKeyR=exInfo.getCoordPathString(NULL);
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

					inInfo.boundsExonic.insert(parentBound);

					exInfo.boundsExonic.insert(leftExonBound);
					exInfo.boundsExonic.insert(rightExonBound);


					inInfo.specBounds=middleBoundSet._set;

					SplidarOutputExtraData extraDat;
					extraDat.commonExonIsCobound=true;


					inInfo.jnxstring="";
					exInfo.jnxstring=chr+":"+StringUtil::str(leftOverlapBound.k2)+this->strand+">"+chr+":"+StringUtil::str(rightOverlapBound.k1+1)+this->strand;

					inInfo.gsid="";
					exInfo.gsid=leftExon->exonGroup->rightSuperGroup->getID()+">"+rightExon->exonGroup->leftSuperGroup->getID();

					//cerr<<"1:"<<leftParentExon->getBound()<<"\t"<<leftOverlapBound<<endl;
					//cerr<<"2:"<<leftParentExon->getBound()<<"\t"<<middleBound<<endl;
					//cerr<<"3:"<<rightParentExon->getBound()<<"\t"<<rightOverlapBound<<endl;
					typedef set<GffEntry::GBlockPtr>::const_iterator sg_i;
					typedef set<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
					typedef pair<sg_i,sg_i> Dsg_i;


					set<GffEntry::GBlockPtr> parentBlocks;
					for(vector<GffEntry::GBlockPtr>::iterator pi=parentExon->blocks->begin();pi!=parentExon->blocks->end();pi++)
						{
						  parentBlocks.insert(*pi);

						}

					if(this->op.getCount && this->hasDataOut())
					{


						Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocks,leftOverlapBound,true);

						Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocks,rightOverlapBound,true);





						//cerr<<"a"<<endl;
						KeyPair<int,int> R0Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(leftBlocks,readLength);
						//cerr<<"a2"<<endl;


						KeyPair<int,int> R1Den(0,0);

						for(BoundSet::iterator bi=middleBoundSet.begin();bi!=middleBoundSet.end();bi++)
						{
							Dsg_i middleBlocks=GffEntry::Exon::getBlocksInRangeGeneric<set<GffEntry::GBlockPtr>, sg_i>(parentBlocks,*bi,true);

							KeyPair<int,int> newR1Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(middleBlocks,readLength); //both left and right parent exons should span through the middle
							::addDensityVectors(R1Den,newR1Den);

						}
						//cerr<<"a3"<<endl;
						KeyPair<int,int> R2Den=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(rightBlocks,readLength);
						//cerr<<"b"<<endl;
						//////here-->


						//[R0|   R1   |R2]
						//[  ]--------[  ]

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


					this->outData("RI",inInfo,exInfo,extraDat);




				}


			}
		}


	}
};



#endif /* RI_SPLIDAR_H_ */
