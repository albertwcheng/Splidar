/*
 * RI_Splidar_backup.h
 *
 *  Created on: Feb 4, 2010
 *      Author: awcheng
 */

#ifndef RI_SPLIDAR_BACKUP_H_
#define RI_SPLIDAR_BACKUP_H_
class RI_Splidar: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	inline RI_Splidar(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NExonGroup::NExonGroupPtr _exongroup,GffEntry::Locus* _locus,string _locusName,int _readLength,bool parentsNeedToBeSpliced=true): op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{

		//^[       ]^
		//  [ ]^[]




/*		NExonGroup::NExonGroupIterator ngi=_exongroup->getExonGroupIterator(false,TRAVEL_DOWNTO_LEAVES);
		pair<int,NExonGroup::NExonGroupPtr> exongNext;
		while(NExonGroup::NExonGroupIterator::isValidItem(exongNext=ngi.nextItem()))
		{
			NExonGroup::NExonGroupPtr curExonGroup=exongNext.second;
			if(curExonGroup->hasChildren() && curExonGroup->levelExons.size()>0) //exclude root. have both itself and children
			{

				//curExonGroup is now the "parent"

				if(curExonGroup->children->size()<2)
				{
					//SG_DEBUG1
					die("Strange Error: has only 1 child?");
				}

				//
				//  [       ] curExonGroup
				// [  ]  [  ]



				bool noExonsHaveAJnx=true;
				if(parentsNeedToBeSpliced)
				{
					//check ignore this exon group if it has no junction coming in or going out, potential unspliced transcript contamination
					//
					//
					//

					for(NExonGroup::ExonI i=curExonGroup->levelExons.begin();i!=curExonGroup->levelExons.end();i++)
					{
						GffEntry::ExonPtr e=*i;
						if(e->hasInOrOut())
						{
							noExonsHaveAJnx=false;
						}
					}

					if(noExonsHaveAJnx)
						continue;
				}

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


					// [        @         ]
					// [  *  ]---------[*   ]
					// [ ] [* ]-------[*  ]-[  ]
					// [][][][ *]----[*][][]
					//from leftExonGroup, travel all the right most branch and get all the *left group exons*
					//from rightExonGroup, travel all the left most branch and get all the *right group exons*
					//all pairs of left and right exons, see whether anyone are connected.
					//for exon u in *left group exons*
					// for jnx j,dest v in u->jnx
					//  if v is in *right group exons*
					//    work on u,v

					//now find all left group exons and right group exons
					set<GffEntry::ExonPtr> leftGroupExons;
					set<NExonGroup::NExonGroupPtr > rightGroupExonGroups;

					//left group exons
					NExonGroup::NExonGroupPtr egTraveler;

					egTraveler=leftExonGroup;
					while(egTraveler)
					{

						for(NExonGroup::ExonI lei=egTraveler->levelExons.begin();lei!=egTraveler->levelExons.end();lei++)
						{
							leftGroupExons.insert(*lei);
						}

						//down one level getting rightmost child exongroup


						egTraveler=egTraveler->getRightMostChild();
					}

					//also find the rightSuperGroups with depth < current leftExonGroup
					NESuperGroupPtr rightSuperGroupOfLeftExonGroup=leftExonGroup->rightSuperGroup;
					for(NESuperGroup::iterator sgi=rightSuperGroupOfLeftExonGroup->begin();
						sgi!=rightSuperGroupOfLeftExonGroup->end();
						sgi++)
					{
						NExonGroup* eg=*sgi;

						if(eg->depth>=leftExonGroup->depth)
							break;

						for(NExonGroup::ExonI lei=eg->levelExons.begin();lei!=eg->levelExons.end();lei++)
						{
							leftGroupExons.insert(*lei);
						}
					}


					egTraveler=rightExonGroup;
					while(egTraveler)
					{

						rightGroupExonGroups.insert(egTraveler);

						//down one level getting rightmost child exongroup


						egTraveler=egTraveler->getLeftMostChild();
					}

					NESuperGroupPtr leftSuperGroupOfRightExonGroup=rightExonGroup->leftSuperGroup;
					for(NESuperGroup::iterator sgi=leftSuperGroupOfRightExonGroup->begin();
						sgi!=leftSuperGroupOfRightExonGroup->end();
						sgi++)
					{
						NExonGroup* eg=*sgi;
						if(eg->depth>=rightExonGroup->depth)
							break;

						rightGroupExonGroups.insert(eg);
					}


					//each pair of exons between the two exon groups
					for(set<GffEntry::ExonPtr>::iterator lei=leftGroupExons.begin();lei!=leftGroupExons.end();lei++)
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

							if(!key_exists(rightGroupExonGroups,rightExon->exonGroup)) //only the right exons from the registered right exon groups
								continue;

							//////here

							KeyPair<int,int> parentConcatBound(INT_MAX,INT_MIN);

							//cerr<<"parentBound:"<<parentBound<<endl;

							KeyPair<int,int> leftExonBound=leftExon->getBound();
							KeyPair<int,int> rightExonBound=rightExon->getBound();

							SG_DEBUG2("leftExon Bound:"<<leftExonBound<<endl);
							SG_DEBUG2("rightExon Bound:"<<rightExonBound<<endl);

							set<GffEntry::GBlockPtr> parentBlocksConcat;

							for(NExonGroup::ExonI pei=curExonGroup->levelExons.begin();pei!=curExonGroup->levelExons.end();pei++)
							{
								GffEntry::ExonPtr curParentExon=*pei;
								if(parentsNeedToBeSpliced && !curParentExon->hasInOrOut()) //again check this condition to include spliced transcript exon only for parents.
								{
									continue;
								}
								KeyPair<int,int> curParentBound=curParentExon->getBound();
								if(areOverlapping(leftExonBound,curParentBound) && areOverlapping(rightExonBound,curParentBound)) //cur parent overlaps with the current child pairs
								{
									//
									//[   ]
									//[] []
									//push parents all blocks into the parentBlocksConcat;
									parentBlocksConcat.insert(curParentExon->blocks->begin(),curParentExon->blocks->end());

									//expand the bound of parentConcatBound by curParentBound in place (to parentConcatBound)
									::expandBoundMutFirst(parentConcatBound,curParentBound);

								}
							}

							if(parentBlocksConcat.size()<1)
							{


								//condition like
								//  [       ]
								// [ ]     [   ]
								//[]----------[]

								SG_DEBUG1("No parent exons overlap with both left and right exons"<<endl);
								continue;


							}


							///   [          ]
							///  [         ]
							///[            ]
							/// [   ]-----[   ]

							///[             ]     parentConcatBound
							///  [  ]              leftOverlapBound
							///           [  ]     rightOverlapBound
							///     [     ]        middleBound (the retained intron)
							///  [           ]     eventBound

							KeyPair<int,int> leftOverlapBound=::overlapBound(parentConcatBound,leftExonBound);
							KeyPair<int,int> rightOverlapBound=::overlapBound(parentConcatBound,rightExonBound);
							KeyPair<int,int> middleBound(leftOverlapBound.k2,rightOverlapBound.k1);
							KeyPair<int,int> eventBound(leftOverlapBound.k1,rightOverlapBound.k2);



							::TrafficInfoAS inInfo;
							::TrafficInfoAS exInfo;

							inInfo.bounds.insert(eventBound);
							exInfo.bounds.insert(leftOverlapBound);
							exInfo.bounds.insert(rightOverlapBound);

							//inInfo.boundsActual=inInfo.bounds;
							//exInfo.boundsActual.insert(leftExonBound);
							//exInfo.boundsActual.insert(rightExonBound);
							inInfo.specBounds.insert(middleBound);

							SplidarOutputExtraData extraDat;
							extraDat.commonExonIsCobound=true;


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


							this->outData("RI",inInfo,exInfo,extraDat);








						}
					}
				}
			}
		}*/
	}

#endif /* RI_SPLIDAR_BACKUP_H_ */
