/*
 * AEP_SplidarGraph.h
 *
 *  Created on: Jan 19, 2010
 *      Author: awcheng
 */

#ifndef AEP_SPLIDARGRAPH_H_
#define AEP_SPLIDARGRAPH_H_

///RECHECK LI DOE
class AEP_SplidarGraph_REGJ : public SplidarGraph {

private:
	KeyPair<int,int> distanceA;
	KeyPair<int,int> distanceB;
	NESuperGroupPtr startGroup;
	string eventTypeName;
	set<string>& isoBoundRegistry;
public:
   //string _EXCELHyperLinkPrefix,


	inline AEP_SplidarGraph_REGJ(Splidar_OpFlag _op,int _minDistanceA,int _maxDistanceA,int _minDistanceB,int _maxDistanceB,ofstream* _fout,ofstream* _foutSeq, RandomAccessFile* _raf, NESuperGroupPtr  _sgroup, GffEntry::Locus* _locus,
			string _locusName, bool _checkLocusName, int _readLength, set<string>& _isoBoundRegistry,string _eventTypeName) :
				startGroup(_sgroup),SplidarGraph(_op,_fout, _foutSeq, _raf, _sgroup, _locus, _locusName, _checkLocusName, _maxDistanceB, _readLength), distanceA(_minDistanceA,_maxDistanceA), distanceB(_minDistanceB,_maxDistanceB),isoBoundRegistry(_isoBoundRegistry),eventTypeName(_eventTypeName) {
		SG_DEBUG2("processing startGroup="<<_sgroup->getID()<<endl);

		if(_minDistanceA<=_maxDistanceA && _maxDistanceA<=_minDistanceB && _minDistanceB<=_maxDistanceB && _minDistanceB>=2)
		{
			//ok. criteria met

		}else
		{
			SG_DEBUG1("error AEP(M1,M2,N1,N2) M1<=M2<=N1<=N2 and N2>=2 criteria not met"<<endl);
		}
	}

	void endGraph() {
		SG_DEBUG2("end graph AEP:"<<startGroup->getID()<<endl);

		//*******CHANGE PENDING CONSIDERATION allow setting genome for the hyperlink


		//string EXCELHyperLinkPrefix="=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&acembly=full&mgcGenes=hide&genscan=dense&ensGene=hide&xenoRefGene=hide&mrna=hide&refGene=hide&position=";
		//string EXCELHyperLinkSuffix="\",\"UCSC Browse\")";

		//dn LeftSuperGroup ->thread
		map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* > MNThread; //dsExon->SuperGroup => [ length-of-SEn-inclusionPath => ExonPostItem (a record that points to the exon thread and the location freeze of that thread) ]
		//cerr<<"a";


		//ExonPostItem: terminal Exon and the path, so for each of the terminal exon
		//all paths arriving at the same terminal exon
		for (ExonPostMatrix::I i=this->exonPostMatrix.begin(); i
				!=this->exonPostMatrix.end(); i++) {
			ExonPostRow* epr=i->second;
			GffEntry::ExonPtr terminalExon=i->first;

			//just want the twos;

			//ExonPostRow::I epri2=epr->find(2);
			SG_DEBUG2("from "<<startGroup->getID()<<" reaching:"<< terminalExon->egsid<<" of "<<terminalExon->exonGroup->sid<<endl);

			ExonPostRow::I epriMaxEndA=epr->upper_bound(this->distanceA.k2);
			ExonPostRow::I epriMinStartA=epr->lower_bound(this->distanceA.k1);

			ExonPostRow::I epriMaxEndB=epr->upper_bound(this->distanceB.k2);
			ExonPostRow::I epriMinStartB=epr->lower_bound(this->distanceB.k1);




			if(epriMinStartB==epr->end() && epriMinStartA==epr->end())
			{
				//nothing to work on, ignore
				continue;
			}

			////cerr<<"b"<<endl;
			//now register

			//int dnCoord=terminalExon->getStart1(); //remove 6/6/2009


			NESuperGroupPtr  dnSGroup=terminalExon->exonGroup->leftSuperGroup;

			SG_DEBUG2("jumping "<<startGroup->getID()<<" to "<<dnSGroup->getID()<<endl);

			map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator miv=MNThread.find(dnSGroup);

			multimap<int,ExonPostItem*>* epiV;

			//the dnSG group not found in MNThread, add it.
			if(miv==MNThread.end())
			{
				epiV=new multimap<int,ExonPostItem*>;
				MNThread.insert(map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::value_type(dnSGroup,epiV));
			}
			else
				epiV=miv->second;

			////cerr<<"c"<<endl;


			//add all the M(A)-degree items
			for(ExonPostRow::I eprim=epriMinStartA;eprim!=epriMaxEndA;eprim++)
			{
				ExonPost* epm=eprim->second;
				for(ExonPost::I epi=epm->begin();epi!=epm->end();epi++)
				{

					if(epi->time<distanceA.k1 || epi->time>distanceA.k2)
					{
						cerr<<"error exclusion bound has epi->timenot in range "<<epi->time;
						die_exit(" Abort");
					}
					epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}


			//add all the N(B)-degree items
			for(ExonPostRow::I eprin=epriMinStartB;eprin!=epriMaxEndB;eprin++)
			{
				ExonPost* epn=eprin->second;
				for(ExonPost::I epi=epn->begin();epi!=epn->end();epi++)
				{
					if(epi->time<distanceB.k1 || epi->time>distanceB.k2)
					{
						cerr<<"error inclusion bound has epi->timenot in range "<<epi->time;
						die_exit(" Abort");
					}
					//cerr<<"adding time "<<epi->time<<endl;
					epiV->insert(multimap<int,ExonPostItem*>::value_type(epi->time,&(*epi)));
				}
			}

			////cerr<<"d"<<endl;
		}


		//now we have filled MNThread : dsExon->SuperGroup => [ length-of-SEn-inclusionPath => ExonPostItem (a record that points to the exon thread and the location freeze of that thread) ]

		////cerr<<"e"<<endl;
			//here comes the critical step!
			//for each terminal coordinate

			///// LI DOE

			for (map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator i1=MNThread.begin();i1!=MNThread.end();i1++) {

				//int dnCoord=i1->first;
				multimap<int,ExonPostItem*>* s2=i1->second; //s2 is the record for all threads ending in current i1->first,aka dnExonGroup

				//all Ms vs Ns
				////cerr<<"f"<<endl;
				typedef multimap<int,ExonPostItem*>::iterator II;

				typedef pair<II,II> DII;


				//DII diiTwos=s2->equal_range(2);
				DII diiMs;
				diiMs.first=s2->lower_bound(distanceA.k1);
				diiMs.second=s2->upper_bound(distanceA.k2);
				DII diiNs;
				diiNs.first=s2->lower_bound(distanceB.k1);
				diiNs.second=s2->upper_bound(distanceB.k2);


				for (II x3=diiMs.first; x3!=diiMs.second; x3++) {
					for (II y3=diiNs.first; y3
							!=diiNs.second; y3++) {

						////cerr<<"g"<<endl;
						if(x3==y3)
						{
							//essentially the same thread, ignore
							continue;
						}


						int M_dist=x3->first;
						int N_dist=y3->first;

						if(M_dist>distanceA.k2 || M_dist<distanceA.k1)
						{
							cerr<<"M_dist out of range "<<M_dist<<" ";
							die_exit(" Abort");
						}

						if(N_dist>distanceB.k2 || N_dist< distanceB.k1)
						{
							cerr<<"N_dist out of range "<<N_dist<<" ";
							die_exit(" Abort");
						}

						ExonPostItem* M_epi=x3->second;
						ExonPostItem* N_epi=y3->second;

						AESPathThread* M_thr=M_epi->exonThread;
						AESPathThread* N_thr=N_epi->exonThread;

						AESPathThread::backtrack_iterator M_right=M_thr->rbegin(M_epi->vertexPointer);

						AESPathThread::backtrack_iterator M_i(M_right);
						AESPathThread::backtrack_iterator M_end=M_thr->rend();

						AESPathThread::backtrack_iterator M_left(M_right); //gonna change anyway
						//one_left++;

						AESPathThread::backtrack_iterator N_right=N_thr->rbegin(N_epi->vertexPointer);

						AESPathThread::backtrack_iterator N_i(N_right);
						AESPathThread::backtrack_iterator N_end=N_thr->rend();

						AESPathThread::backtrack_iterator N_left(N_right); //gonna change anyway




						if(!areOverlapping(M_right->exon->getBound(),N_right->exon->getBound()))
							continue;

						if(op.commonFlankingCobound) //require coboundary of flanking exon
						{
							//check for condition for right most exon:
							// -----[
							// -----[
							if(M_right->exon->getBound().k1 != N_right->exon->getBound().k1)
								continue;
						}


						deque<InEdgeAndVertex> ievs_N;
						deque<InEdgeAndVertex> ievs_M;

						::TrafficInfoAS inInfo;
						::TrafficInfoAS exInfo;
						//start from left of N_right to left end


						string inMiddleJnxString;
						string inMiddleEgString;

						string exMiddleJnxString;
						string exMiddleEgString;


						//string chr=one_left->exon->chr;
						string chr=M_right->exon->chr;

						//setup inclusion bounds and exon strings
						//N_i++ initially to skip including the rightmost exon
						//N_i!=N_end as terminal condition to avoid including the leftmost exon
						//from left to right
						//iev N_i is iev contains edge to dst node and the dst node

						//assuming M is ex, N is in

						N_i++;




						int realN_dist=0;
						int realM_dist=0;


						AESPathThread::backtrack_iterator N_i_CheckStart(N_i);

						GffEntry::ExonPtr M_leftSecondExon=NULL;
						GffEntry::ExonPtr N_leftSecondExon=NULL;

						for(;N_i!=N_end;N_i++)
						{
							realN_dist++;
							if(N_i.hasAtLeastOneMore())
							{
								//not the last one
								const InEdgeAndVertex& iev=*N_i;
								inMiddleJnxString=chr+":"+StringUtil::str(iev.exon->getStart1())+":"+this->strand+","+chr+":"+StringUtil::str(iev.exon->getEnd1())+":"+this->strand+">"+inMiddleJnxString;
								inMiddleEgString=iev.exon->exonGroup->sid+">"+inMiddleEgString;
								inInfo.bounds.insert(iev.exon->getBound());
								ievs_N.push_front(iev);


								N_leftSecondExon=iev.exon;
							}
							else
							{
								//the last one
								N_left=N_i;
							}
						}


						if(realN_dist!=N_dist)
						{
							cerr<<"N_dists inconsistent "<<realN_dist<<" " <<N_dist<<" ";
							die_exit("abort");
						}

						AESPathThread::backtrack_iterator N_i_CheckEnd(N_left); //no need to check the N_left


						//the same thing apply to M's, now need to check overlaps of any of the M exons with N exons, if yes, exit!

						bool overlapViolated=false;

						for(M_i++;M_i!=M_end;M_i++)
						{
							realM_dist++;

							if(M_i.hasAtLeastOneMore())
							{
								//not the last one
								const InEdgeAndVertex& iev=*M_i;


								//check overlap

								for(AESPathThread::backtrack_iterator N_check=N_i_CheckStart;N_check!=N_i_CheckEnd;N_check++)
								{
									if(areOverlapping(iev.exon->getBound(),N_check->exon->getBound()))
									{
										//one of the exon from the N paths and the M paths overlap!
										overlapViolated=true;
										break;
									}
								}

								exMiddleJnxString=chr+":"+StringUtil::str(iev.exon->getStart1())+":"+this->strand+","+chr+":"+StringUtil::str(iev.exon->getEnd1())+":"+this->strand+">"+exMiddleJnxString;
								exMiddleEgString=iev.exon->exonGroup->sid+">"+exMiddleEgString;
								exInfo.bounds.insert(iev.exon->getBound());
								ievs_M.push_front(iev);

								M_leftSecondExon=iev.exon;

							}
							else
							{
								//the last one
								M_left=M_i;
							}
						}

						if(realM_dist!=M_dist)
						{
							cerr<<"M_dists inconsistent "<<realM_dist<<" " <<M_dist<<" ";
							die_exit("abort");
						}
						/*cerr<<"AAA:"
						<<M_left->exon->getBound().k2<<">"<<(M_right->exon->getBound().k1+1)
						<<"/"
						<<N_left->exon->getBound().k2<<">"<<(N_right->exon->getBound().k1+1)
						<<endl;*/


						if(overlapViolated) //one of the exon from N paths and the M paths overlap! ignore, continue with other possibilities
							continue;


						///Now More Checks!

						if(!areOverlapping(M_left->exon->getBound(),N_left->exon->getBound()))
						{


							//they shoulbe be overlapping!! else a fatal bug
							//cerr<<this->SpliceTravsersalGraph::locusName<<endl;
							cerr<<M_left->exon->chr<<" ";
							cerr<<"M_left"<<M_left->exon->getBound()<<" ";
							cerr<<"N_left"<<N_left->exon->getBound()<<" ";
							cerr<<"special case ignore: M_left and N_left not overlapping"<<endl;

							continue;
						}

						if(op.commonFlankingCobound) //requires cobound of flanking exons
						{
							//check for condition
							// ]-----
							// ]----

							if(M_left->exon->getBound().k2!=N_left->exon->getBound().k2)
								continue;
						}


						//*****CHANGE PENDING CONSIDERATION: Require one_left->exon and N_left->exon to be in the same rightSuperGroup*****
						//
						//they must be anayway, because caller put in starter exons in the same rightSupergroup
						 if (M_left->exon->exonGroup->rightSuperGroup != N_left->exon->exonGroup->rightSuperGroup)
						 {
							 die("fatal error: AEP M_left,N_left not in the same right Super Group");

							 continue;
						 }
						//
						//////









						//cerr<<"h"<<endl;

						//now get the bound of the common area
						KeyPair<int,int> leftOverlapBound=::overlapBound(M_left->exon->getBound(),N_left->exon->getBound());
						KeyPair<int,int> rightOverlapBound=::overlapBound(M_right->exon->getBound(),N_right->exon->getBound());

						//
						//N (inc)  [   ]----[ ]---[      ]
						//M (exc)   [ ]----------[   ]
						//LOB        x
						//ROB                      xx
						//LIB        xx
						//RIB                      xx
						//LEB        x
						//REB                     xxx
						//IB         xx      x     xx
						//EB         x            xxx
						//LISB        x
						//RISB      -------null----------
						//LESB      -------null----------
						//RESB                    x
						//ISB         x      x
						//ESB                     x

						KeyPair<int,int>  leftIncBound=KeyPair<int,int>(leftOverlapBound.k1,N_left->exon->getBound().k2);
						KeyPair<int,int> rightIncBound=KeyPair<int,int>(N_right->exon->getBound().k1,rightOverlapBound.k2);

						KeyPair<int,int>  leftExcBound=KeyPair<int,int>(leftOverlapBound.k1,M_left->exon->getBound().k2);
						KeyPair<int,int> rightExcBound=KeyPair<int,int>(M_right->exon->getBound().k1,rightOverlapBound.k2);


						KeyPair<int,int> leftIncSpecBound=KeyPair<int,int>(leftOverlapBound.k2,N_left->exon->getBound().k2);
						KeyPair<int,int> rightIncSpecBound=KeyPair<int,int>(N_right->exon->getBound().k1,rightOverlapBound.k1);

						KeyPair<int,int> leftExcSpecBound=KeyPair<int,int>(leftOverlapBound.k2,M_left->exon->getBound().k2);
						KeyPair<int,int> rightExcSpecBound=KeyPair<int,int>(M_right->exon->getBound().k1,rightOverlapBound.k1);




						//cerr<<"i"<<endl;



						SG_DEBUG2("SEmn event:"<<M_left->exon->egsid<<":"<<M_right->exon->egsid<<endl);



						//cerr<<"j"<<endl;

						//now make a copy of bounds (central exon bounds) to specific bounds

						inInfo.specBounds=inInfo.bounds;
						exInfo.specBounds=exInfo.bounds;
						inInfo.boundsExonic=inInfo.bounds;
						exInfo.boundsExonic=exInfo.bounds;

						//now diverge, add in the bounds:

						{
						//inInfo.bounds.insert(leftOverlapBound);
						//inInfo.bounds.insert(rightOverlapBound);
						//exInfo.bounds.insert(leftOverlapBound);
						//exInfo.bounds.insert(rightOverlapBound);
							inInfo.bounds.insert(leftIncBound);
							inInfo.bounds.insert(rightIncBound);
							exInfo.bounds.insert(leftExcBound);
							exInfo.bounds.insert(rightExcBound);



						}

						SplidarOutputExtraData extraDat;
						extraDat.commonExonIsCobound=true;

						//now bounds are completed. Check whether isbounds have been outputed already.
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


						{
							if(len01(leftIncSpecBound)>0){
								inInfo.specBounds.insert(leftIncSpecBound);
								extraDat.commonExonIsCobound=false;
							}

							if(len01(rightIncSpecBound)>0){
								inInfo.specBounds.insert(rightIncSpecBound);
								extraDat.commonExonIsCobound=false;
							}
							if(len01(leftExcSpecBound)>0){
								exInfo.specBounds.insert(leftExcSpecBound);
								extraDat.commonExonIsCobound=false;
							}
							if(len01(rightExcSpecBound)>0){
								exInfo.specBounds.insert(rightExcSpecBound);
								extraDat.commonExonIsCobound=false;
							}

						}

						//P



						{


						inInfo.boundsExonic.insert(N_left->exon->getBound());
						inInfo.boundsExonic.insert(N_right->exon->getBound());
						exInfo.boundsExonic.insert(M_left->exon->getBound());
						exInfo.boundsExonic.insert(M_right->exon->getBound());

						}


						//int uCoord=leftOverlapBound.k2;
						//int dnCoord=rightOverlapBound.k1+1; //convert to 1-base coordinate

						//finish up the jnxstring annotation
						//string leftjnxstring=chr+":"+StringUtil::str(uCoord)+this->strand;
						//string rightjnxstring=chr+":"+StringUtil::str(dnCoord)+this->strand;

						string exleftjnxstring=chr+":"+StringUtil::str(M_left->exon->getEnd1())+":"+this->strand;
						string exrightjnxstring=chr+":"+StringUtil::str(M_right->exon->getStart1())+":"+this->strand;

						string inleftjnxstring=chr+":"+StringUtil::str(N_left->exon->getEnd1())+":"+this->strand;
						string inrightjnxstring=chr+":"+StringUtil::str(N_right->exon->getStart1())+":"+this->strand;

						//inInfo.jnxstring=chr+":"+StringUtil::str(uCoord)+this->strand+">"+StringUtil::str(two_center->exon->getStart1())+","+chr+":"+StringUtil::str(two_center->exon->getEnd1())+strand+">"+chr+":"+StringUtil::str(dnCoord)+this->strand;

						//this doesn't have to be the same!


						exInfo.jnxstring=exleftjnxstring+">"+exMiddleJnxString+exrightjnxstring;
						inInfo.jnxstring=inleftjnxstring+">"+inMiddleJnxString+inrightjnxstring;

						string ileftegstring;
						string irightegstring;
						string eleftegstring;
						string erightegstring;
						//finish up the exongroupstring annotation
						if(op.egStringOutputCommonFlankingAsCoord)
						{
							eleftegstring=StringUtil::str(M_left->exon->getBound().k2);
							erightegstring=StringUtil::str(M_right->exon->getBound().k1+1);
							ileftegstring=StringUtil::str(N_left->exon->getBound().k2);
							irightegstring=StringUtil::str(N_right->exon->getBound().k1+1);
						}
						else
						{
							eleftegstring=M_left->exon->exonGroup->rightSuperGroup->getID();
							erightegstring=M_right->exon->exonGroup->leftSuperGroup->getID();
							ileftegstring=N_left->exon->exonGroup->rightSuperGroup->getID();
							irightegstring=N_right->exon->exonGroup->leftSuperGroup->getID();

							if(ileftegstring!=eleftegstring)
							{
								die("Strange Error: ileftegstring!=eleftegstring");

							}
							if(irightegstring!=erightegstring)
							{
								die("Strange Error: irightegstring!=erightegstring");
							}
						}
						exInfo.gsid=eleftegstring+">"+exMiddleEgString+erightegstring;
						inInfo.gsid=ileftegstring+">"+inMiddleEgString+irightegstring;

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

						//the blocks in M_left,M_right contained in the leftOverlapBound and rightOverlapBound are the same as N_left,N_right respectively
						Dsg_i leftBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*M_left->exon->blocks,leftOverlapBound,true);
						Dsg_i rightBlocks=GffEntry::Exon::getBlocksInRangeGeneric<vector<GffEntry::GBlockPtr>, sg_i>(*M_right->exon->blocks,rightOverlapBound,true);
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
						SG_DEBUG2(M_left->exon->egsid<<":left Den:"<<leftDen<<endl);
						SG_DEBUG2(M_right->exon->egsid<<":right Den:"<<rightDen<<endl);



							{ ///INCLUSION CENTER STUFF
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



									KeyPair<int,int> middleDen=i->exon->getDensity(readLength,true); //middleDen= [ ](^[ ]^[*])^[ ]     => traversal direction
									KeyPair<int,int> jnxDen=i->inJnx->getDensity(readLength,true); //jnxDen= [ ](^[ ]^*[ ])^[ ]        => traversal direction

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

							{ ///EXCLUSION CENTER STUFF
								//center stuff
								//cerr<<"s"<<endl;


								//for each of the middle exons and jnx
								for(deque<InEdgeAndVertex>::iterator i=ievs_M.begin();i!=ievs_M.end();i++)
								{

									////*** PENDING CONSIDERING: REMOVAL
									//if(!i->inJnx)
									//{
									//	cerr<<"you are shitting me with in Jnx == NULL"<<endl;
									//}
									/////******



									KeyPair<int,int> middleDen=i->exon->getDensity(readLength,true); //middleDen= [ ](^[ ]^[*])^[ ]     => traversal direction
									KeyPair<int,int> jnxDen=i->inJnx->getDensity(readLength,true); //jnxDen= [ ](^[ ]^*[ ])^[ ]        => traversal direction

									SG_DEBUG2(i->exon->egsid<<":"<<middleDen<<"\t");

									exInfo.JPF*=jnxDen.k2;
									exInfo.JRF*=jnxDen.k1;

									exInfo.JPCheckString+=StringUtil::str(jnxDen.k2)+",";
									exInfo.JRCheckString+=StringUtil::str(jnxDen.k1)+",";

									addDensityVectors(exInfo.noFlankingInfo,jnxDen);
									addDensityVectors(exInfo.noFlankingInfo,middleDen);
									addDensityVectors(exInfo.withFlankingInfo,jnxDen);
									addDensityVectors(exInfo.withFlankingInfo,middleDen);

									addDensityVectors(exInfo.middleExonsFlow,middleDen);
									addDensityVectors(exInfo.jnxsFlow,jnxDen);
								}
								////cerr<<"gt"<<"<a2"<<endl;
								//cerr<<"t"<<endl;
							}

							SG_DEBUG2(endl);

						     //right flanking exon, and jnx before.
							{ //EXCLUSION FORM
								//////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"o"<<endl;
								KeyPair<int,int> jnxDen=M_right->inJnx->getDensity(readLength,true);  //jxnDen= [ ]V*[ ]

								addDensityVectors(exInfo.jnxsFlow,jnxDen);
								addDensityVectors(exInfo.noFlankingInfo,jnxDen);
								addDensityVectors(exInfo.withFlankingInfo,jnxDen);



								addDensityVectors(exInfo.flankingExonsFlow,rightDen);  //rightDen= [ ]V[*] = []^[]^[]^[*]

								//if(includeFlankExon)
								if(M_dist==1)
									addDensityVectors(exInfo.withFlankingInfo,rightDen);

								exInfo.JPF*=jnxDen.k2;
								exInfo.JRF*=jnxDen.k1;

								exInfo.JRCheckString+=StringUtil::str(jnxDen.k1);
								exInfo.JPCheckString+=StringUtil::str(jnxDen.k2);

								addDensityVectors(exInfo.flankingExonsFlow,leftDen); //leftDen= [*]V[ ] = [*]^[]^[]^[]

								//if(includeFlankExon)
								if(M_dist==1)
									addDensityVectors(exInfo.withFlankingInfo,leftDen);
								//cerr<<"p"<<endl;
								////cerr<<"gt"<<"<a1"<<endl;
							}

							{ //INCLUSION FORM

								////cerr<<"gt"<<"a1"<<endl;
								//cerr<<"q"<<endl;
								KeyPair<int,int> jnxDen=N_right->inJnx->getDensity(readLength,true); //jnxDen= [ ]^[ ]^[ ]^*[ ]

								addDensityVectors(inInfo.jnxsFlow,jnxDen);
								addDensityVectors(inInfo.noFlankingInfo,jnxDen);
								addDensityVectors(inInfo.withFlankingInfo,jnxDen);

								addDensityVectors(inInfo.flankingExonsFlow,rightDen);

								//if(includeFlankExon)
								if(M_dist==1)
									addDensityVectors(inInfo.withFlankingInfo,rightDen);

								inInfo.JPF*=jnxDen.k2;
								inInfo.JRF*=jnxDen.k1;

								inInfo.JPCheckString+=StringUtil::str(jnxDen.k2);
								inInfo.JRCheckString+=StringUtil::str(jnxDen.k1);

								addDensityVectors(inInfo.flankingExonsFlow,leftDen);

								//if(includeFlankExon)
								if(M_dist==1)
									addDensityVectors(inInfo.withFlankingInfo,leftDen);
								//cerr<<"r"<<endl;
							}
						}

						//cerr<<"u"<<endl;


						//now check!

						::TrafficInfoAS *inInfoActual=&inInfo;
						::TrafficInfoAS *exInfoActual=&exInfo;


						//inInfo now is N
						//exInfo now is M

						if(M_dist==N_dist) //MX Epaths
						{
							if(M_leftSecondExon->getBound()<N_leftSecondExon->getBound())
							{
								//M is genomic 5'
								//N is genomic 3'
								if(strand==GffEntry::FORWARD)
								{
									//INCLUSION IS SUPPOSED TO BE TRANSCRIPT 5'
									//now M IS 5' than N => swap such that M~inInfo, N~exInfo
									//two need to be transcript 5' so swap
									inInfoActual=&exInfo;
									exInfoActual=&inInfo;
								}
								else
								{
									//NO SWAPPING
									///one_info=&one_info1;
									///two_info=&two_info1;
								}
							}
							else
							{
								//M is genomic 3'
								//N is genomic 5'
								if(strand==GffEntry::REVERSE)
								{

									//INCLUSION IS SUPPOSED TO BE TRANSCRIPT 5'
									//now M is transcript 5' so swap such that M~inInfo, N~exInfo
									inInfoActual=&exInfo;
									exInfoActual=&inInfo;
								}
								else
								{
									//one_info=&one_info1;
									//two_info=&two_info1;
								}
							}

						}


						this->outData(this->eventTypeName,*inInfoActual,*exInfoActual,extraDat);
						//cerr<<"v"<<endl;
					}
				}

			}

		//now free memory from the registry
			//cerr<<"n"<<endl;
			for(map<NESuperGroupPtr ,multimap<int,ExonPostItem*>* >::iterator i=MNThread.begin();i!=MNThread.end();i++)
			{
				delete i->second;
			}

			//cerr<<"o"<<endl;


	}




};

#endif /* AEP_SPLIDARGRAPH_H_ */
