/*
 * SpliceMMGraph.h
 *
 *  Created on: Dec 1, 2009
 *      Author: awcheng
 */

#ifndef SPLICEMMGRAPH_H_
#define SPLICEMMGRAPH_H_


#include "snipIncludes.h"





#define DENOBJ_JUNCTION 1
#define DENOBJ_GBLOCK 2
#define DENOBJ_TERMINAL 3
#define DENOBJ_COLLAPSED 4

class DenObj;
class DenGraph;

typedef SmartPtr<DenObj> DenObjPtr;

class DenObj
{

public:
	int id;

	inline bool isVisited(int stamp) const
	{
		return this->lastVisitedKey==stamp;
	}
	inline bool attemptVisit(int stamp)
	{
		if ( isVisited(stamp) ){
			return false;
		}

		this->lastVisitedKey=stamp;
		return true;

	}

	int lastVisitedKey;
	GffEntry::GBlockPtr pBlock;


	GffEntry::JnxTagSetPtr pSynJnx;

	static const int Junction;
	static const int Gblock;
	static const int Terminal;
	static const int Collapsed;

	set<DenObjPtr> *rights;
	set<DenObjPtr> *lefts;
	set<DenObjPtr> *tmpLefts;
	set<DenObjPtr> *tmpRights;
	set<DenObjPtr> *collapsedRights; //point to the appropriate rights after collapsed.


	DenObjPtr oneTrack;

	inline bool hasOneTrack() const
	{
		return oneTrack;
	}


	//constructor for a terminal type Den Obj
	inline DenObj():lastVisitedKey(-1),pBlock(NULL),pSynJnx(NULL),lefts(new set<DenObjPtr>),rights(new set<DenObjPtr>),tmpLefts(NULL),tmpRights(NULL),collapsedRights(NULL){}
	inline DenObj(GffEntry::GBlockPtr _pBlock):lastVisitedKey(-1),pBlock(_pBlock),pSynJnx(NULL),lefts(new set<DenObjPtr>),rights(new set<DenObjPtr>),tmpLefts(NULL),tmpRights(NULL),collapsedRights(NULL){}
	inline DenObj(GffEntry::JnxTagSetPtr _pSynJnx): lastVisitedKey(-1),pBlock(NULL),pSynJnx(_pSynJnx),lefts(new set<DenObjPtr>),rights(new set<DenObjPtr>),tmpLefts(NULL),tmpRights(NULL),collapsedRights(NULL){}
	~DenObj()
	{
		delete rights;
		delete lefts;



		if(tmpLefts) //should not exec coz tmpLefts should be transferred
			delete tmpLefts;
		if(tmpRights)
			delete tmpRights;
	}

	string prefix;
	inline void setPrefix(string _prefix)
	{
		prefix=_prefix;
	}

	inline void setID(int _id)
	{
		id=_id;
	}

	inline int getID() const
	{
		return id;
	}

	inline void allocateTmps()
	{
		tmpLefts=new set<DenObjPtr>;
		tmpRights=new set<DenObjPtr>;
	}



	inline void applyTmps()
	{
		if(lefts)
			delete lefts;
		if(rights)
			delete rights;

		lefts=tmpLefts;
		rights=tmpRights;
		tmpLefts=NULL;
		tmpRights=NULL;
	}

	//attempt to add den obj left of or right of the current object, if not existed already, it is added, else, it is destroyed. Returns the Den obj address as desired
    DenObjPtr addLeftEnsured(DenGraph* _parent,DenObjPtr _left);
	DenObjPtr addRightEnsured(DenGraph* _parent,DenObjPtr _right);

	inline void addLeftDirectly(DenObjPtr _left)
	{
		this->lefts->insert(_left);
	}

	inline void addRightDirectly(DenObjPtr _right)
	{
		this->rights->insert(_right);
	}


	inline int getType() const
	{
		if (pBlock)
			return DenObj::Gblock;

		if (pSynJnx)
			return DenObj::Junction;


		return DenObj::Terminal;
	}

	KeyPair<int,int> getBound() const
	{
		switch(this->getType())
		{
		case DENOBJ_JUNCTION:
			return this->pSynJnx->getBound11();

		case DENOBJ_GBLOCK:
			return this->pBlock->getBound11();

		case DENOBJ_TERMINAL:
		default:
			return KeyPair<int,int>(-1,-1);

		}
	}

	inline bool operator > (const DenObj& right) const
	{
		KeyPair<int,int> thisBound=this->getBound();
		KeyPair<int,int> rightBound=right.getBound();
		int thisType=this->getType();
		int rightType=right.getType();

		if(thisBound==rightBound)
		{
			if(thisType==rightType)
			{
				if(thisType==DenObj::Terminal)
				{
					return this>(&right); //Terminal requires absolute same address coz there is only one start and stop. Terminal shouldn't be compared anyway.
				}
				else
				{
					return false;
				}

			}else
			{
				return thisType>rightType;
			}
		}else
		{
			return thisBound>rightBound;
		}

	}

	inline bool operator < (const DenObj& right) const
	{
		KeyPair<int,int> thisBound=this->getBound();
		KeyPair<int,int> rightBound=right.getBound();
		int thisType=this->getType();
		int rightType=right.getType();

		if(thisBound==rightBound)
		{
			if(thisType==rightType)
			{
				if(thisType==DenObj::Terminal)
				{
					return this<(&right); //Terminal requires absolute same address coz there is only one start and stop. Terminal shouldn't be compared anyway.
				}
				else
				{
					return false;
				}

			}else
			{
				return thisType<rightType;
			}
		}else
		{
			return thisBound<rightBound;
		}
	}

	inline bool operator == (const DenObj& right) const
	{
		if(this->getType() != right.getType())
			return false;

		//now the type is the same
		switch(this->getType())
		{
		case DENOBJ_JUNCTION:
			return this->pSynJnx==right.pSynJnx;
		case DENOBJ_TERMINAL:
			return (this)==(&right);
		case DENOBJ_GBLOCK:
			return this->pBlock==right.pBlock;
		default:
			return false;
		}
	}

	inline bool operator != (const DenObj& right) const
	{
		return !(*this==right);
	}
	inline bool operator >=(const DenObj& right) const
	{
		return !(*this<right);
	}
	inline bool operator <=(const DenObj& right) const
	{
		return !(*this>right);
	}

	friend class DenGraph;

	inline static ostream& printBEDHeader(ostream& os,string trackname="SpliceMMGraphDenObjs",string description="SpliceMMGraph Density Objects",string visibility="full",bool useScore=true)
	{
		os<<"track name="<<trackname<<" description=\""<<description<<"\" visibility="<<visibility<<" useScore="<<(useScore?"1":"0")<<endl;
		return os;
	}

	inline string getName(bool simple=false) const
	{
		switch(getType())
		{
			case DENOBJ_JUNCTION:
				return (simple?string(""):prefix)+"J"+StringUtil::str(getID());


			case DENOBJ_GBLOCK:
				return (simple?string(""):prefix)+"B"+StringUtil::str(getID());

				break;
			case DENOBJ_TERMINAL:
				return prefix;
				break;
			default:
				return "";
		}
	}

	inline ostream& printBEDStyle(ostream& os,  string prefix, string chrom, char strand) const
	{
		if (this->getType()==DENOBJ_TERMINAL)
			return os;

		string printName;

		KeyPair<int,int> bound=this->getBound();
		switch(this->getType())
		{
		case DENOBJ_JUNCTION:

			os<<chrom<<" "<<(bound.k1-2)<<" "<<(bound.k2+1)<<" "<<getName()<<" "<<500<<" "<<strand<<" "<<(bound.k1-2)<<" "<<(bound.k2+1)<<" 0,0,0 2 10,10 0,"<<(bound.k2-7-bound.k1)<<endl;
			break;
		case DENOBJ_GBLOCK:

			os<<chrom<<" "<<(bound.k1-1)<<" "<<(bound.k2)<<" "<<getName()<<" "<<1000<<" "<<strand<<" "<<(bound.k1-1)<<" "<<(bound.k2)<<" 0,0,0 1 "<< (bound.k2-bound.k1+1)<<" 0"<<endl;
			break;
		default:
			break;
		}
		return os;
	}





};

inline ostream& operator << (ostream& os,const DenObj& denObj)
{

	KeyPair<int,int> bound=denObj.getBound();

	if(denObj.getType()==DENOBJ_TERMINAL)
		os<<denObj.getName();
	else
		os<<denObj.getName()<< "["<<bound.k1<<"-"<<bound.k2<<"]";


	return os;
}

inline ostream& operator << (ostream& os, DenGraph& graph);

class DenGraph
{
private:
	int traversalKey;
	set <DenObjPtr> denObjPool;
	DenObjPtr start;
	DenObjPtr end;
	bool transferred;
	GffEntry::Locus* locus;

	char gstrand; //original genomic strand of the locus
	char tstrand; //current DenGraph strand wrt transcript strand
	//inline DenGraph(char _gstrand,char _tstrand):gstrand(_gstrand),tstrand(_tstrand),traversalKey(0),start(NULL),end(NULL),transferred(false){}

	typedef KeyPair<DenObjPtr,DenObjPtr> GBlockDenObjPair;
	typedef map<GBlockDenObjPair, GffEntry::JnxTagSetPtr > JnxPool;
	typedef JnxPool::iterator JnxPoolI;

	typedef set<DenObjPtr>::iterator DenObjPoolI;

	JnxPool jnxPool;

public:



	inline int issueNewTraversalKey(){
		return ++traversalKey;
	}
	inline virtual ~DenGraph(){
		if(transferred)
		{
			return;
		}

		for(set<DenObjPtr>::iterator i=denObjPool.begin();i!=denObjPool.end();i++)
		{
			delete *i;
		}

		for(JnxPoolI i=jnxPool.begin();i!=jnxPool.end();i++)
		{
			delete i->second;
		}

		delete start;
		delete end;
	}

	inline DenObjPtr ensureDenObj(DenObjPtr _obj)
	{
		DenObjPtr working=_obj;

		_obj->setID(denObjPool.size()+1);
		_obj->setPrefix(locus->getFirstName()+"_");

		pair<DenGraph::DenObjPoolI,bool> insertResult=this->denObjPool.insert(_obj);
		if (!insertResult.second)
		{
			delete _obj;
			working=*insertResult.first;
		}

		return working;
	}



	inline GffEntry::JnxTagSetPtr registerSynJnx(GffEntry::JnxPtr _jnx,int _jnxMode)
	{
		//get the first block and right block info from jnx by following the exon info
		GffEntry::GBlockPtr leftExonBlock=_jnx->leftExon->getGenomicLastBlock();
		GffEntry::GBlockPtr rightExonBlock=_jnx->rightExon->getGenomicFirstBlock();

		DenObjPtr leftExonBlockDenObj=this->ensureDenObj(new DenObj(leftExonBlock));
		DenObjPtr rightExonBlockDenObj=this->ensureDenObj(new DenObj(rightExonBlock));

		//try to pool out the registered synjnx from jnxPool

		GffEntry::JnxTagSetPtr synJnx;

		GBlockDenObjPair jnxHashKey(leftExonBlockDenObj,rightExonBlockDenObj);
		JnxPoolI synJnxI=this->jnxPool.find(jnxHashKey);

		//add in jnx tag by either union or intersection as defined by jnxmode
		if( synJnxI==this->jnxPool.end())
		{
			synJnx=new GffEntry::JnxTagSet(*_jnx);
			this->jnxPool.insert(JnxPool::value_type(jnxHashKey,synJnx));

		}
		else
		{
			synJnx=synJnxI->second;
			synJnx->addTags(*_jnx,_jnxMode);
		}

		//now connect!
		DenObjPtr pJnxDenObj=this->ensureDenObj(new DenObj(synJnx));

		leftExonBlockDenObj->addRightDirectly(pJnxDenObj);
		pJnxDenObj->addLeftDirectly(leftExonBlockDenObj);
		rightExonBlockDenObj->addLeftDirectly(pJnxDenObj);
		pJnxDenObj->addRightDirectly(rightExonBlockDenObj);

		return synJnx;

	}

	inline GffEntry::Locus* getLocus()
	{
		return locus;
	}

	inline DenGraph(GffEntry::Locus* _locus,int jnxMode):locus(_locus),gstrand(locus->strand),tstrand(locus->strand),traversalKey(0),transferred(false)
	{
		die_exit("DenGraph is not implemented. see SpliceMMGraph.h");

	/*	start=new DenObj(); //create start as a terminal object
		end=new DenObj(); //create end as a terminal object

		start->setID(-1);
		start->setPrefix("START");
		end->setID(-1);
		end->setPrefix("END");

		for(vector<GffEntry*>::iterator i=locus->transcripts.begin();i!=locus->transcripts.end();i++)
		{
				GffEntry* transcript=*i;
				GffEntry::ExonPtr GenomicHeadExon=transcript->getGenomicFirstExon();
				GffEntry::GBlockPtr GenomicHeadBlock=GenomicHeadExon->getGenomicFirstBlock();
				DenObjPtr thisHeadBlockDenObj=start->addRightEnsured(this,new DenObj(GenomicHeadBlock));
				thisHeadBlockDenObj->addLeftDirectly(start);
				GffEntry::ExonPtr GenomicTailExon=transcript->getGenomicLastExon();
				GffEntry::GBlockPtr GenomicTailBlock=GenomicTailExon->getGenomicLastBlock();
				DenObjPtr thisTailBlockDenObj=end->addLeftEnsured(this,new DenObj(GenomicTailBlock));
				thisTailBlockDenObj->addRightDirectly(end);

		}

		//now head and tail blocks are attached to START and END respectively

		NExonGroup::ExonIterator exoni=locus->getExonTree()->getAllLevelsExonIterator();
		GffEntry::ExonPtr curExon;
		while(curExon=exoni.nextItem())
		{
			vector<GffEntry::GBlockPtr>& curBlocks=*curExon->blocks;
			int nblocks=curBlocks.size();
			for(int blocki=0;blocki<nblocks-1;blocki++) //from first block of the exon to the second last block of that exon
			{
				GffEntry::GBlockPtr curBlock=curBlocks[blocki];
				GffEntry::GBlockPtr nextBlock=curBlocks[blocki+1];
				//connect blocks;
				DenObjPtr curBlockDenObj=this->ensureDenObj(new DenObj(curBlock));
				DenObjPtr nextBlockDenObj=curBlockDenObj->addRightEnsured(this,new DenObj(nextBlock));

				nextBlockDenObj->addLeftDirectly(curBlockDenObj);
			}

			//now add junction to record
			if(curExon->hasOutJnxs())
			{
				for(GffEntry::Exon::OutJnxI oji=curExon->outJnxs->begin();oji!=curExon->outJnxs->end();oji++)
				{
						GffEntry::JnxPtr outjnx=oji->second;
						this->registerSynJnx(outjnx,jnxMode);

						//almost forgot, now you need to connect the blocks thru jnx.

				}
			}
		}
*/


	}

	inline ostream& printAllDenObjs(ostream& os)
	{
		for(set <DenObjPtr>::iterator i= denObjPool.begin();i!=denObjPool.end();i++)
		{
			os<<**i<<endl;
		}
		return os;
	}

	inline ostream& printDenObjsinBEDFormat(ostream& os)
	{
		for(set <DenObjPtr>::iterator i= denObjPool.begin();i!=denObjPool.end();i++)
		{
			(*i)->printBEDStyle(os,this->locus->getFirstName(),this->locus->chr,this->locus->strand);
		}
		return os;
	}

	inline ostream& printDenGraphDFSInSIFFormat(ostream& os,bool printCollapsed,string prefixEachLine,string interactionName)
	{



		int kstmp=this->issueNewTraversalKey();

		end->attemptVisit(kstmp); //visit the end node so that it knows to stop.

		typedef KeyPair<DenObjPtr,string> StateMem;

		deque<StateMem> Q;
		Q.push_front(StateMem(start,""));

		while(!Q.empty())
		{
			StateMem stateMemu=Q.front();
			Q.pop_front();

			DenObjPtr u=stateMemu.k1;
			string parentName=stateMemu.k2;


			string thisName=u->getName(true);
			DenObjPtr ori=u;

			if(printCollapsed){
				if(u->hasOneTrack())
				{
					DenObjPtr nu;
					while(nu=u->oneTrack)
					{
						thisName+=nu->getName(true);
						u=nu;
					}
				}
			}

			if(parentName!="")
				os<<prefixEachLine<<parentName<<" "<<interactionName<<" "<<thisName<<endl;



			if(!ori->isVisited(kstmp))
			{
				for(set<DenObjPtr>::iterator v=u->rights->begin();v!=u->rights->end();v++)
				{
					//if(!(*v)->isVisited(kstmp) || *v==end){
						Q.push_front(StateMem(*v,thisName));
					//}

				}
			}

			ori->attemptVisit(kstmp);


		}

		return os;
	}

/*	inline ostream& printDenGraphDFS(ostream& os)
	{



		int kstmp=this->issueNewTraversalKey();

		end->attemptVisit(kstmp); //visit the end node so that it knows to stop.

		typedef KeyPair<DenObjPtr,int> StateMem;

		deque<StateMem> Q;
		Q.push_front(StateMem(start,0));

		while(!Q.empty())
		{
			StateMem stateMemu=Q.front();
			Q.pop_front();

			DenObjPtr u=stateMemu.k1;
			int curLevel=stateMemu.k2;

			if(curLevel<0)
			{
				//revisiting visited node
				curLevel=curLevel*-1;

				for(int x=0;x<curLevel;x++)
				{
					os<<" ";
				}

				os<<"~"<<(*u)<<endl;

				continue;
			}

			u->attemptVisit(kstmp);

			for(int x=0;x<curLevel;x++)
			{
				os<<" ";
			}

			os<<"+"<<(*u);

			if(u->hasOneTrack())
			{
				DenObjPtr nu;
				while(nu=u->oneTrack)
				{
					os<<" "<<(*nu);
					u=nu;
				}
			}

			os<<endl;

			for(set<DenObjPtr>::iterator v=u->rights->begin();v!=u->rights->end();v++)
			{
				if(!(*v)->isVisited(kstmp)){
					Q.push_front(StateMem(*v,curLevel+1));
				}
				else
				{
				Q.push_front(StateMem(*v,-1*(curLevel+1)));
				}
			}




		}

		return os;
	}
*/
	inline ostream& printDenGraphDFS(ostream& os)
	{



		int kstmp=this->issueNewTraversalKey();

		end->attemptVisit(kstmp); //visit the end node so that it knows to stop.

		typedef KeyPair<DenObjPtr,int> StateMem;

		deque<StateMem> Q;
		Q.push_front(StateMem(start,0));

		while(!Q.empty())
		{
			StateMem stateMemu=Q.front();
			Q.pop_front();

			DenObjPtr u=stateMemu.k1;
			int curLevel=stateMemu.k2;

			if(curLevel<0)
			{
				//revisiting visited node
				curLevel=curLevel*-1;

				for(int x=0;x<curLevel;x++)
				{
					os<<" ";
				}

				os<<"~"<<(*u)<<endl;

				if(u->hasOneTrack())
				{
					DenObjPtr nu;
					while(nu=u->oneTrack)
					{
						os<<" "<<(*nu);
						u=nu;
					}
				}

				continue;
			}

			u->attemptVisit(kstmp);

			for(int x=0;x<curLevel;x++)
			{
				os<<" ";
			}

			os<<"+"<<(*u);

			if(u->hasOneTrack())
			{
				DenObjPtr nu;
				while(nu=u->oneTrack)
				{
					os<<" "<<(*nu);
					u=nu;
				}
			}

			os<<endl;

			for(set<DenObjPtr>::iterator v=u->rights->begin();v!=u->rights->end();v++)
			{
				if(!(*v)->isVisited(kstmp)){
					Q.push_front(StateMem(*v,curLevel+1));
				}
				else
				{
				Q.push_front(StateMem(*v,-1*(curLevel+1)));
				}
			}




		}

		return os;
	}

	inline void updateCollapsedRight()
	{
		for(set<DenObjPtr>::iterator i=this->denObjPool.begin();i!=denObjPool.end();i++)
		{
			DenObjPtr ori=*i;
			if(!ori->hasOneTrack())
				continue; //no need to update

			DenObjPtr u=ori;
			DenObjPtr nu;
			while(nu=u->oneTrack)
			{
				u=nu;
			}

			//arriving at a node that has not more oneTrack, i.e., the last node of the track
			//that's what we want. we want the last node's rights
			ori->collapsedRights=u->rights;
		}
	}

	inline void collapseGraphInPlace()
	{
		//do a DFS
		deque<DenObjPtr> Q;


		Q.push_front(start);

		int kstmp=this->issueNewTraversalKey();

		end->attemptVisit(kstmp); //visit the end node so that it knows to stop.



		while(!Q.empty())
		{
			DenObjPtr u=Q.front();
			Q.pop_front();
			u->attemptVisit(kstmp);

			if(u->rights->size()==1)
				u->oneTrack=*u->rights->begin();

			for(set<DenObjPtr>::iterator v=u->rights->begin();v!=u->rights->end();v++)
			{
				if(!(*v)->isVisited(kstmp))
					Q.push_front(*v);
			}

		}
	}

	inline void inverseGraphInPlace()
	{
		DenObjPtr oldStart=start;
		DenObjPtr oldEnd=end;
		end=oldStart;
		start=oldEnd;

		if(this->tstrand==GffEntry::FORWARD)
			this->tstrand=GffEntry::REVERSE;
		else if(this->tstrand==GffEntry::REVERSE)
			this->tstrand=GffEntry::FORWARD;

		oldStart->allocateTmps();
		oldEnd->allocateTmps();

		//needs to invalidate all the one track and collapsedRights

		for(set<DenObjPtr>::iterator i=denObjPool.begin();i!=denObjPool.end();i++)
		{
			DenObjPtr u=*i;
			u->collapsedRights=NULL;
			u->oneTrack=NULL;
		}


		//DFS go leftwards from oldEnd
		int kstmp=this->issueNewTraversalKey();

		//fake visit oldStart such that DFS does not go there
		oldStart->attemptVisit(kstmp);

		deque<DenObjPtr> Q;

		Q.push_front(oldEnd);
		while(!Q.empty())
		{
			DenObjPtr u=Q.front();
			Q.pop_front();

			u->attemptVisit(kstmp);



			for(set<DenObjPtr>::iterator v=u->lefts->begin();v!=u->lefts->end();v++)
			{
				DenObjPtr pv=*v;
				if(!pv->isVisited(kstmp))
				{
					pv->allocateTmps();
					Q.push_front(pv);
				}
				//connect!!
				//pv was left of u, but should be right of u in the future

				pv->tmpLefts->insert(u);
				u->tmpRights->insert(pv);

			}


		}



		//update lefts right from tmp
		for(set<DenObjPtr>::iterator i=this->denObjPool.begin();i!=this->denObjPool.end();i++)
		{
			(*i)->applyTmps();
		}

		start->applyTmps();
		end->applyTmps();
		start->setPrefix("START");
		end->setPrefix("END");

	}
	friend class DenObj;
	friend ostream& operator << (ostream& os, DenGraph& graph);
};

inline ostream& operator << (ostream& os, DenGraph& graph)
{
	os<<"*********"<<endl;
	os<<"genomic (original strand):"<<graph.gstrand<<endl;
	os<<"transcript strand:"<<graph.tstrand<<endl;
	graph.printDenGraphDFS(os);
	os<<"**********"<<endl;
	return os;
}


class SpliceMMGraph: public DenGraph
{
public:
	inline SpliceMMGraph(GffEntry::Locus* _locus,int jnxMode):DenGraph(_locus,jnxMode){}


};

#endif /* SPLICEMMGRAPH_H_ */
