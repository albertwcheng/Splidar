/*
 * A53SS_SplidarGraph.h
 *
 *  Created on: Jan 19, 2010
 *      Author: awcheng
 */

#ifndef A53SS_SPLIDARGRAPH_H_
#define A53SS_SPLIDARGRAPH_H_

#define A53SS_LEFT true
#define A53SS_RIGHT false

class A53SSHashStruct
{
public:
	int fixedCoord;
	int variantCoord;
	GffEntry::ExonPtr variantExon;
	GffEntry::ExonPtr commonExon;
	GffEntry::JnxPtr jnx;
	inline A53SSHashStruct(int _fixedCoord,int _variantCoord,GffEntry::ExonPtr _variantExon,GffEntry::ExonPtr _commonExon,GffEntry::JnxPtr _jnx)
		:fixedCoord(_fixedCoord),variantCoord(_variantCoord),variantExon(_variantExon),commonExon(_commonExon),jnx(_jnx){}

	inline bool operator < (const A53SSHashStruct& right) const
	{
		if(fixedCoord==right.fixedCoord)
		{
			if(variantCoord==right.variantCoord)
			{
				if(variantExon==right.variantExon)
				{
					return commonExon<right.commonExon;
				}
				else return variantExon<right.variantExon;
			}
			else
				return variantCoord<right.variantCoord;
		}
		else
			return fixedCoord<right.fixedCoord;
	}
	inline bool operator > (const A53SSHashStruct& right) const
	{
		if(fixedCoord==right.fixedCoord)
		{
			if(variantCoord==right.variantCoord)
			{
				if(variantExon==right.variantExon)
				{
					return commonExon>right.commonExon;
				}
				else return variantExon>right.variantExon;
			}
			else
				return variantCoord>right.variantCoord;
		}
		else
			return fixedCoord>right.fixedCoord;
	}

	inline bool operator == ( const A53SSHashStruct& right) const
	{
		return fixedCoord==right.fixedCoord && variantCoord==right.variantCoord && variantExon==right.variantExon && commonExon==right.commonExon;
	}

	inline bool operator !=(const A53SSHashStruct& right) const
	{
		return !(*this==right);
	}

	inline bool operator <= ( const A53SSHashStruct& right) const
	{
		return !(*this>right);
	}

	inline bool operator >= (const A53SSHashStruct& right) const
	{
		return !(*this<right);
	}
};




class A53SS_Splidar_REGJ: public GenericSplidarOutputFormat
{
private:
	Splidar_OpFlag op;
public:

	/*inline void countNum(const char* str,int& L,int& M,int &R)
	{

		char c;L=0;M=0;R=0;
		while((c=*(str++))!='\0')
		{
			switch(c)
			{
			case 'L':
				L++;
				break;
			case 'M':
				M++;
				break;
			case 'R':
				R++;
				break;

			}
		}
	}*/



	/*inline bool theyAreTrueExonVariantOn(NExonGroup::NExonGroupPtr a,NExonGroup::NExonGroupPtr b,bool left) //just a switch: left => true; right => false;
	{
		const string& asid=a->sid;
		const string& bsid=b->sid;
		int alen=asid.length();
		int blen=bsid.length();
		const char* asidc=asid.c_str();
		const char* bsidc=bsid.c_str();
		const char *prefix;
		if(alen<blen)
		{	//a is probably a prefix
			prefix=StringUtil::isPrefix(asidc,bsidc);
		}else
		{
			prefix=StringUtil::isPrefix(bsidc,asidc);
		}

		if(!prefix)
		{
			cerr<<"strange: neither is prefix of another:"<<asidc<<" vs "<<bsidc<<endl;
			return false;
		}

		int L;
		int M;
		int R;

		this->countNum(prefix,L,M,R);

		if(left)
		{
			//focus on left
			//--[ ]
			//---[]

			return (R==0 && M==0);
		}
		else
		{
			//focus on right
			//[ ]--
			//[]---

			return (L==0 && M==0);
		}


	}*/

	inline bool theyAreTrueExonVariantOn(GffEntry::ExonPtr a,GffEntry::ExonPtr b,bool left) //just a switch: left => true; right => false;
	{
		//return this->theyAreTrueExonVariantOn(a->exonGroup,b->exonGroup,left);
		if(left)
		{
			return a->exonGroup->leftSuperGroup==b->exonGroup->leftSuperGroup;
		}
		else
		{
			return a->exonGroup->rightSuperGroup==b->exonGroup->rightSuperGroup;
		}
	}

	inline A53SS_Splidar_REGJ(Splidar_OpFlag _op,ofstream* _fout,ofstream* _foutSeq,RandomAccessFile* _raf,NESuperGroupPtr  _sgroup,GffEntry::Locus* _locus,string _locusName,int _readLength,set<string>& isoBoundRegistry):op(_op),GenericSplidarOutputFormat( _readLength, _locus, _locusName,  _fout, _foutSeq,_raf, _op.excelHyperLinkPrefix,_op.excelHyperLinkSuffix,_op.seqGetI5Len,_op.seqGetI3Len,_op.seqGetE5Len,_op.seqGetE3Len)
	{
		//for each transcript
		//add to second exon into secondExonSet



		//set { fixedCoord => variantCoord => variantExon => commonExon /*common is just to avoid the same exact set of exons*/ }
		//set<A53SSHashStruct > twoMapAG5E; //genomic alternative 5' exon
		//set<A53SSHashStruct > twoMapAG3E; //genomic alternative 3' exon
		typedef set<A53SSHashStruct> A53SSHash;
		typedef pair<A53SSHash*,A53SSHash*> A53SSHashPair;
		typedef map<NESuperGroupPtr,  A53SSHashPair> A53SSBigHash;
		typedef A53SSBigHash::iterator A53SSBigHashI;

		typedef A53SSHash::iterator A53SSHashStructI;
		typedef pair<A53SSHashStructI,A53SSHashStructI> A53SSHashStructDI;
		typedef vector<GffEntry::GBlockPtr>::iterator sg_i;
		typedef vector<GffEntry::GBlockPtr>::reverse_iterator sg_ri;
		typedef pair<sg_i,sg_i> Dsg_i;

		A53SSBigHash bigHash;

		//hashing phase: goto each exon in the RightSuperGroup of left exon, go to the right exon by traversing an edge and then hashing the by the leftsupergroup of the right exon

		//NExonGroup::ExonIterator exoni=locus->root->getAllLevelsExonIterator();
		GffEntry::ExonPtr leftExon;

		for(NESuperGroup::iterator exonGroupI=_sgroup->begin();exonGroupI!=_sgroup->end();exonGroupI++)
		{
			NExonGroup* exonGroup=*exonGroupI;





			for(NExonGroup::ExonI exonI=exonGroup->levelExons.begin();exonI!=exonGroup->levelExons.end();exonI++)
			{
				leftExon=*exonI;

				if(!leftExon->outJnxs) //this exon has no out edge
					continue;


				for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator i=leftExon->outJnxs->begin();i!=leftExon->outJnxs->end();i++)
				{
					GffEntry::ExonPtr rightExon=i->first;

					A53SSHash* twoMapAG5E;
					A53SSHash* twoMapAG3E;


					//cerr<<"usageflag->usageflag "<<leftExon->usageFlag<<"->"<<rightExon->usageFlag<<endl;

					if(rightExon->usageFlag!=leftExon->usageFlag)
						continue; //avoid interlocus links

					NESuperGroupPtr thisRightExonLeftSuperGroup=rightExon->exonGroup->leftSuperGroup;

					A53SSBigHashI bigHashI=bigHash.find(thisRightExonLeftSuperGroup);

					if(bigHashI==bigHash.end()) //
					{
						twoMapAG5E=new A53SSHash;
						twoMapAG3E=new A53SSHash;
						bigHash.insert(A53SSBigHash::value_type(thisRightExonLeftSuperGroup,A53SSHashPair(twoMapAG5E,twoMapAG3E)));
					}else
					{
						twoMapAG5E=bigHashI->second.first;
						twoMapAG3E=bigHashI->second.second;
					}

					GffEntry::JnxPtr jnx=i->second;

					//A53SSHashStruct(fixedCoord,VarCoord,VarExon,AnchorExon)

					twoMapAG5E->insert(A53SSHashStruct(rightExon->getStart1(),leftExon->getEnd1(),leftExon,rightExon,jnx));
					//( V v------f A )
					//(  V  v----f A )


					twoMapAG3E->insert(A53SSHashStruct(leftExon->getEnd1(),rightExon->getStart1(),rightExon,leftExon,jnx));
					//( A f-------v V )
					//( A f-----v   V )
				}



			}

		}

		//hash complete, now traverse hash




		//process g5'
		//cerr<<"size of twoMapAG5E="<<twoMapAG5E.size()<<endl;

		for(A53SSBigHashI bigHashI=bigHash.begin();bigHashI!=bigHash.end();bigHashI++)
		{

			NESuperGroupPtr thisRightExonLeftSuperGroup=bigHashI->first;
			A53SSHashPair hashPair=bigHashI->second;

			A53SSHash& twoMapAG5E=*hashPair.first;
			A53SSHash& twoMapAG3E=*hashPair.second;

			A53SSHashStructI i=twoMapAG5E.begin(); ///
			A53SSHashStructDI curRange;
			int prevB;

			while(i!=twoMapAG5E.end()) ///
			{
				prevB=i->fixedCoord;
			//	cerr<<i->fixedCoord<<"\t";
				curRange.first=i;
				curRange.second=(++i);

				if(op.commonFlankingCobound)
				{

					while(i!=twoMapAG5E.end()) ///
					{
						if(i->fixedCoord!=prevB)
							break;

						curRange.second=(++i);
					}

				}else
				{
					//don't care about the fixedCoord, just go to the end
					curRange.second=twoMapAG5E.end();
				}

				//now the curRange contains all the stuff with the same fixed Coord, and the variant Exon are sorted by coordinate of left->end,

				//cerr<<"a"<<endl;

				//cannot assume fixed coord the same

				for(A53SSHashStructI x=curRange.first;x!=curRange.second;x++)
					{
						A53SSHashStructI y=x;
						y++;
						for(;y!=curRange.second;y++)
						{
							//proximal is the one close to anchor exon , this having a larger variant coord (right end)

							// [ P   ]---[]
							// [ D ]-----[]
							A53SSHashStruct proximaliev=*y;
							A53SSHashStruct distaliev=*x;

							//cerr<<"b1"<<endl;
							GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
							KeyPair<int,int> proximalBound=proximalExon->getBound();
							GffEntry::Jnx* proximalJnx=proximaliev.jnx;

							//cerr<<"b2"<<endl;
							GffEntry::ExonPtr distalExon=distaliev.variantExon;
							KeyPair<int,int> distalBound=distalExon->getBound();
							GffEntry::Jnx* distalJnx=distaliev.jnx;

							//cerr<<"b3"<<endl;

							if(!::areOverlapping(proximalBound,distalBound) || proximaliev.variantCoord==distaliev.variantCoord)
							{
								continue; //the two variable exons do not overlap or they have the same variant coord, ignore
							}

							//[ ]--[]
							//[]---[]


							if(!theyAreTrueExonVariantOn(proximalExon,distalExon,A53SS_RIGHT))
							{

								//require
								/*
								 * [   x   ]-----[]
								 * [ ] [x]-------[]
								 *
								 */

								continue;
							}

						//	//cerr<<"ok! continue!"<<endl;

							TrafficInfoAS xinInfo;
							TrafficInfoAS yexInfo;

							//cerr<<"b4"<<endl;


							//jnx string
							//    ]-----[
							//      ]---[
							//   or
							//    ]------["fixed"
							//      ]---["fixed"

							//these fixed coord can be different...
							xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximalExon->getEnd1())+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand;
							yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distalExon->getEnd1())+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand;

							//
							//  RSG ------LSG
							//  RSG   ---LSG
							// or op.egStringOutputCommonFlankingAsCoord==TRUE
							//  RSG -----[
							//  RSG   ---[

							string ifixedSideString;
							string efixedSideString;




							if(op.egStringOutputCommonFlankingAsCoord)
							{
								ifixedSideString=proximaliev.commonExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand;
								efixedSideString=distaliev.commonExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand;
							}
							else
							{
								ifixedSideString=proximaliev.commonExon->exonGroup->leftSuperGroup->getID();
								efixedSideString=distaliev.commonExon->exonGroup->leftSuperGroup->getID();
							}
							xinInfo.gsid=proximalExon->exonGroup->rightSuperGroup->getID()+">"+ifixedSideString;
							yexInfo.gsid=distalExon->exonGroup->rightSuperGroup->getID()+">"+efixedSideString;

							//xinInfo.gsid=proximalExon->exonGroup->relRootString+">"+proximaliev.commonExon->exonGroup->relRootString;
							//yexInfo.gsid=distalExon->exonGroup->relRootString+">"+distaliev.commonExon->exonGroup->relRootString;

							//cerr<<"b5"<<endl;
							KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());

							KeyPair<int,int> inCommonBound=KeyPair<int,int>(proximaliev.commonExon->getBound().k1,commonBound.k2);
							KeyPair<int,int> exCommonBound=KeyPair<int,int>(distaliev.commonExon->getBound().k1,commonBound.k2);
							KeyPair<int,int> inCommonSpecBound=KeyPair<int,int>(proximaliev.commonExon->getBound().k1,commonBound.k1);
							KeyPair<int,int> exCommonSpecBound=KeyPair<int,int>(distaliev.commonExon->getBound().k1,commonBound.k1);


							//cerr<<"b6"<<endl;

							//c=variant common bound, x=variant extra bound, L= variant larger bound, C=Common bound

							//  variant       anchor
							//    LLLLLL
							// (  cccxxx]-------[CCC]
							//   (ccc]--------[  CCC    ]

							KeyPair<int,int> variantCommonBound=::overlapBound(distalBound,proximalBound);
							KeyPair<int,int> variantExtraBound=KeyPair<int,int>(distalBound.k2,proximalBound.k2);

							if(!isValid(variantExtraBound))
							{
								cerr<<"strange error: variantExtraBound invalid: "<<variantExtraBound<<endl;
								die("");
								continue;
							}

							//cerr<<"b7"<<endl;

							KeyPair<int,int> variantLargerBound=KeyPair<int,int>(variantCommonBound.k1,variantExtraBound.k2);

							xinInfo.bounds.insert(variantLargerBound);
							xinInfo.bounds.insert(inCommonBound);

							yexInfo.bounds.insert(variantCommonBound);
							yexInfo.bounds.insert(exCommonBound);

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



							xinInfo.specBounds.insert(variantExtraBound);

							SplidarOutputExtraData extraDat;
							extraDat.commonExonIsCobound=true;

							if(len01(inCommonSpecBound)>0){
								xinInfo.specBounds.insert(inCommonSpecBound);
								extraDat.commonExonIsCobound=false;
							}
							if(len01(exCommonSpecBound)>0){
								yexInfo.specBounds.insert(exCommonSpecBound);
								extraDat.commonExonIsCobound=false;
							}

							//inclusion is proximal




							xinInfo.boundsExonic.insert(proximaliev.variantExon->getBound());
							xinInfo.boundsExonic.insert(proximaliev.commonExon->getBound());

							yexInfo.boundsExonic.insert(distaliev.variantExon->getBound());
							yexInfo.boundsExonic.insert(distaliev.commonExon->getBound());


							//cerr<<"b8"<<endl;
							if(this->op.getCount && this->hasDataOut())
							{
								Dsg_i commonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.commonExon->blocks,commonBound,true);
								KeyPair<int,int> commonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(commonBlocks,readLength);

								//cerr<<"b9"<<endl;
								Dsg_i variantCommonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantCommonBound,true);
								KeyPair<int,int> variantCommonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantCommonBlocks,readLength);
								//cerr<<"b10"<<endl;
								Dsg_i variantLargerBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantLargerBound,true);
								KeyPair<int,int> variantLargerDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantLargerBlocks,readLength);
								//cerr<<"b11"<<endl;
								Dsg_i variantExtraBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantExtraBound,true);
								KeyPair<int,int> variantExtraDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantExtraBlocks,readLength);
								//cerr<<"b12"<<endl;

								addDensityVectors(xinInfo.flankingExonsFlow,commonDen);
								addDensityVectors(yexInfo.flankingExonsFlow,commonDen);
								addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
								addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
								addDensityVectors(xinInfo.middleExonsFlow,variantLargerDen);
								addDensityVectors(yexInfo.middleExonsFlow,variantCommonDen);

								//cerr<<"b13"<<endl;
								xinInfo.JRCheckString=StringUtil::str(xinInfo.jnxsFlow.k1);
								xinInfo.JPCheckString=StringUtil::str(xinInfo.jnxsFlow.k2);
								xinInfo.JRF=xinInfo.jnxsFlow.k1;
								xinInfo.JPF=xinInfo.jnxsFlow.k2;

								yexInfo.JRCheckString=StringUtil::str(yexInfo.jnxsFlow.k1);
								yexInfo.JPCheckString=StringUtil::str(yexInfo.jnxsFlow.k2);
								yexInfo.JRF=yexInfo.jnxsFlow.k1;
								yexInfo.JPF=yexInfo.jnxsFlow.k2;


								addDensityVectors(xinInfo.noFlankingInfo,variantExtraDen); //just the extra area
								addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
								xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;
								addDensityVectors(xinInfo.withFlankingInfo,variantCommonDen);


								//addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
								addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow); //just the jnx!
								yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;
								addDensityVectors(yexInfo.withFlankingInfo,variantCommonDen); //add the common density! for NE+
							}

							this->outData((this->strand==GffEntry::FORWARD)?"A5SS":"A3SS",xinInfo,yexInfo,extraDat);

						}
					}

			}


			//cerr<<"b<"<<endl;

			//process g3'
			i=twoMapAG3E.begin(); ///
			//ExonPairWJnxDI curRange;


			while(i!=twoMapAG3E.end()) ///
			{
				prevB=i->fixedCoord;

				curRange.first=i;
				curRange.second=(++i);

				if(op.commonFlankingCobound)
				{

					while(i!=twoMapAG3E.end()) ///
					{
						if(i->fixedCoord!=prevB)
							break;

						curRange.second=(++i);
					}

				}else
				{
					//don't care about the fixed coord, just go to end
					curRange.second=twoMapAG3E.end();
				}
				//now the curRange contains all the stuff with the same fixed Coord, and the variant Exon are sorted by coordinate of right->start



				for(A53SSHashStructI x=curRange.first;x!=curRange.second;x++)
					{
						A53SSHashStructI y=x;
						y++;
						for(;y!=curRange.second;y++)
						{
							//just need to swap here!
							//cerr<<"c1"<<endl;
							A53SSHashStruct proximaliev=*x;
							A53SSHashStruct distaliev=*y;
							//cerr<<"c2"<<endl;
							GffEntry::ExonPtr proximalExon=proximaliev.variantExon;
							KeyPair<int,int> proximalBound=proximalExon->getBound();
							GffEntry::Jnx* proximalJnx=proximaliev.jnx;
							//cerr<<"c3"<<endl;
							GffEntry::ExonPtr distalExon=distaliev.variantExon;
							KeyPair<int,int> distalBound=distalExon->getBound();
							GffEntry::Jnx* distalJnx=distaliev.jnx;
							//cerr<<"c4"<<endl;
							if(!::areOverlapping(proximalBound,distalBound) || proximaliev.variantCoord==distaliev.variantCoord)
							{
								continue; //the two variable exons do not overlap or they have the same variant coord, ignore
							}

							if(!theyAreTrueExonVariantOn(proximalExon,distalExon,A53SS_LEFT))
							{

								//require
								/*
								 * []--------[   x   ]
								 * []----------[x] [ ]
								 *
								 */

								continue;
							}
						//	//cerr<<"ok! continue!"<<endl;

							//jnx string
							//    ]-----[
							//    ]---[


							TrafficInfoAS xinInfo;
							TrafficInfoAS yexInfo;
							//cerr<<"c5"<<endl;
							xinInfo.jnxstring=proximalExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand+">"+proximalExon->chr+":"+StringUtil::str(proximalExon->getStart1())+":"+this->strand;
							//cerr<<"c5.1"<<endl;
							yexInfo.jnxstring=distalExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand+">"+distalExon->chr+":"+StringUtil::str(distalExon->getStart1())+":"+this->strand;
							//cerr<<"c5.2"<<endl;
							//cerr<<proximaliev.commonExon<<endl;
							//cerr<<proximaliev.commonExon->getBound()<<endl;
							//cerr<<proximaliev.commonExon->exonGroup<<endl;
							//cerr<<proximaliev.commonExon->exonGroup->relRootString<<endl;

							//eg string
							// ]------ LSG
							// ]----   LSG
							// or op.egStringOutputCommonFlankingAsCoord==TRUE
							// RSG]------ LSG
							// RSG]---- LSG

							string ifixedSideString;
							string efixedSideString;

							if(op.egStringOutputCommonFlankingAsCoord)
							{
								ifixedSideString=proximaliev.commonExon->chr+":"+StringUtil::str(proximaliev.fixedCoord)+":"+this->strand;
								efixedSideString=distaliev.commonExon->chr+":"+StringUtil::str(distaliev.fixedCoord)+":"+this->strand;
							}
							else
							{
								ifixedSideString=proximaliev.commonExon->exonGroup->rightSuperGroup->getID();
								efixedSideString=distaliev.commonExon->exonGroup->rightSuperGroup->getID();
							}


							xinInfo.gsid=ifixedSideString+">"+proximalExon->exonGroup->leftSuperGroup->getID();
							yexInfo.gsid=efixedSideString+">"+distalExon->exonGroup->leftSuperGroup->getID();

							//xinInfo.gsid=proximaliev.commonExon->exonGroup->relRootString+">"+proximalExon->exonGroup->relRootString; //flip!
							//cerr<<"c5.3"<<endl;
							//yexInfo.gsid=distaliev.commonExon->exonGroup->relRootString+">"+distalExon->exonGroup->relRootString;
							//cerr<<"c6"<<endl;


							KeyPair<int,int> commonBound=::overlapBound(proximaliev.commonExon->getBound(),distaliev.commonExon->getBound());
							//cerr<<"c7"<<endl;

							KeyPair<int,int> inCommonBound=KeyPair<int,int>(commonBound.k1,proximaliev.commonExon->getBound().k2);
							KeyPair<int,int> exCommonBound=KeyPair<int,int>(commonBound.k1,distaliev.commonExon->getBound().k2);
							KeyPair<int,int> inCommonSpecBound=KeyPair<int,int>(commonBound.k2,proximaliev.commonExon->getBound().k2);
							KeyPair<int,int> exCommonSpecBound=KeyPair<int,int>(commonBound.k2,distaliev.commonExon->getBound().k2);



							KeyPair<int,int> variantCommonBound=::overlapBound(distalBound,proximalBound);
							KeyPair<int,int> variantExtraBound=KeyPair<int,int>(proximalBound.k1,distalBound.k1);


							//c=variant common bound, x=variant extra bound, L= variant larger bound, C=Common bound

							//  anchor     variant
							//             LLLLLL
							// (  CCC]-------[ccc)
							//   (CCC]-----[xxccc   )

							if(!isValid(variantExtraBound))
							{
								cerr<<"strange error: variantExtraBound invalid: "<<variantExtraBound<<endl;
								die("");
								continue;
							}

							KeyPair<int,int> variantLargerBound=KeyPair<int,int>(variantExtraBound.k1,variantCommonBound.k2);

							xinInfo.bounds.insert(inCommonBound);
							xinInfo.bounds.insert(variantLargerBound);

							yexInfo.bounds.insert(exCommonBound);
							yexInfo.bounds.insert(variantCommonBound);

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


							//inclusion is proximal
							xinInfo.boundsExonic.insert(proximaliev.commonExon->getBound());
							xinInfo.boundsExonic.insert(proximaliev.variantExon->getBound());

							yexInfo.boundsExonic.insert(distaliev.commonExon->getBound());
							yexInfo.boundsExonic.insert(distaliev.variantExon->getBound());

							xinInfo.specBounds.insert(variantExtraBound);

							SplidarOutputExtraData extraDat;
							extraDat.commonExonIsCobound=true;

							if(len01(inCommonSpecBound)>0){
								xinInfo.specBounds.insert(inCommonSpecBound);

								extraDat.commonExonIsCobound=false;
							}
							if(len01(exCommonSpecBound)>0){
								yexInfo.specBounds.insert(exCommonSpecBound);
								extraDat.commonExonIsCobound=false;
							}

							//cerr<<"c8"<<endl;
							Dsg_i commonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.commonExon->blocks,commonBound,true);
							KeyPair<int,int> commonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(commonBlocks,readLength);
							//cerr<<"c9"<<endl;
							Dsg_i variantCommonBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantCommonBound,true);
							KeyPair<int,int> variantCommonDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantCommonBlocks,readLength);
							//cerr<<"c10"<<endl;
							Dsg_i variantLargerBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantLargerBound,true);
							KeyPair<int,int> variantLargerDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantLargerBlocks,readLength);
							//cerr<<"c11"<<endl;
							Dsg_i variantExtraBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*proximaliev.variantExon->blocks,variantExtraBound,true);
							KeyPair<int,int> variantExtraDen=GffEntry::GBlock::getDensityOfContigBlocks<sg_i,sg_ri>(variantExtraBlocks,readLength);

							//cerr<<"c12"<<endl;
							addDensityVectors(xinInfo.flankingExonsFlow,commonDen);
							addDensityVectors(yexInfo.flankingExonsFlow,commonDen);
							addDensityVectors(xinInfo.jnxsFlow,proximalJnx->getDensity(readLength,true));
							addDensityVectors(yexInfo.jnxsFlow,distalJnx->getDensity(readLength,true));
							addDensityVectors(xinInfo.middleExonsFlow,variantLargerDen);
							addDensityVectors(yexInfo.middleExonsFlow,variantCommonDen);
							//cerr<<"c3"<<endl;

							xinInfo.JRCheckString=StringUtil::str(xinInfo.jnxsFlow.k1);
							xinInfo.JPCheckString=StringUtil::str(xinInfo.jnxsFlow.k2);
							xinInfo.JRF=xinInfo.jnxsFlow.k1;
							xinInfo.JPF=xinInfo.jnxsFlow.k2;

							yexInfo.JRCheckString=StringUtil::str(yexInfo.jnxsFlow.k1);
							yexInfo.JPCheckString=StringUtil::str(yexInfo.jnxsFlow.k2);
							yexInfo.JRF=yexInfo.jnxsFlow.k1;
							yexInfo.JPF=yexInfo.jnxsFlow.k2;


							addDensityVectors(xinInfo.noFlankingInfo,variantExtraDen); //just the extra area
							addDensityVectors(xinInfo.noFlankingInfo,xinInfo.jnxsFlow);
							xinInfo.withFlankingInfo=xinInfo.noFlankingInfo;
							addDensityVectors(xinInfo.withFlankingInfo,variantCommonDen);


							//addDensityVectors(yexInfo.noFlankingInfo,yexInfo.middleExonsFlow);
							addDensityVectors(yexInfo.noFlankingInfo,yexInfo.jnxsFlow); //just the jnx!
							yexInfo.withFlankingInfo=yexInfo.noFlankingInfo;
							addDensityVectors(yexInfo.withFlankingInfo,variantCommonDen); //add the common density! for NE+

							this->outData((this->strand==GffEntry::REVERSE)?"A5SS":"A3SS",xinInfo,yexInfo,extraDat);

						}



					}



			}
		}
		//now free bigHash
		for(A53SSBigHashI bi=bigHash.begin();bi!=bigHash.end();bi++)
		{
			A53SSHashPair hashPair=bi->second;
			delete hashPair.first;
			delete hashPair.second;
		}




	}
};


#endif /* A53SS_SPLIDARGRAPH_H_ */
