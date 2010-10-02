#ifndef SPLICEGRAPH_H_
#define SPLICEGRAPH_H_

#include "snipIncludes.h"







class InEdgeAndVertex
{
public:
	GffEntry::JnxPtr inJnx;
	GffEntry::ExonPtr exon;
	InEdgeAndVertex(GffEntry::JnxPtr _inJnx,GffEntry::ExonPtr _exon): inJnx(_inJnx),exon(_exon)
	{}
	
	
	inline bool operator < (const InEdgeAndVertex& right) const
	{
		return exon<right.exon;
	}
	
	inline bool operator > (const InEdgeAndVertex& right) const
	{
		return exon>right.exon;
	}
	
	inline bool operator == ( const InEdgeAndVertex& right) const
	{
		return exon==right.exon;
	}
	
	inline bool operator!=(const InEdgeAndVertex& right) const
	{
		return (*this==right);
	}
	
};

class AESPathThreadExonGroup;

class AESPathThread
{
public:
	AESPathThreadExonGroup* parentExonGroupThread;
	vector<InEdgeAndVertex> ievs;
	AESPathThread* root; //the root thread
	AESPathThread* prec; //the preceding thread
	inline AESPathThread(AESPathThreadExonGroup* _parentExonGroupThread,AESPathThread* _prec,AESPathThread* root);
	inline AESPathThread(AESPathThreadExonGroup* _parentExonGroupThread,const InEdgeAndVertex& iev, AESPathThread* _prec,AESPathThread* root);
	
	inline bool hasPrev()
	{
		if(prec==NULL&&ievs.size()==1)
			return false;
		
		return true;
	}
	
	class backtrack_iterator
	{
	public:
		AESPathThread* thr;
		int indexOnThread;
		inline backtrack_iterator(AESPathThread* _thr,int _indexOnThread): thr(_thr),indexOnThread(_indexOnThread)
		{}
		inline void operator++()
		{
			if(indexOnThread>0)
				indexOnThread--;
			else
			{
				thr=thr->prec;
				if(thr)
					indexOnThread=thr->ievs.size()-1;
				else
					indexOnThread=0;
			}
		}
		inline bool operator==(const backtrack_iterator& right)
		{
			return thr==right.thr&& indexOnThread==right.indexOnThread;
		}
		inline bool operator!=(const backtrack_iterator& right)
		{
			return !(*this==right);
		}
		inline void operator++(int)
		{
			++(*this);
		}
		inline const InEdgeAndVertex& operator*()
		{
			return thr->ievs[indexOnThread];
		}
		inline InEdgeAndVertex* operator->()
		{
			return &thr->ievs[indexOnThread];
		}
		inline AESPathThread* getThread()
		{
			return thr;
		}
		inline bool hasAtLeastOneMore()
		{
			//if(thr->prec==NULL&&thr->ievs.size()==1)
			if(thr->prec==NULL && indexOnThread==0)
				return false;
			return true;
		}
	};
	
	
	inline int getCurrentVertexPointer() const
	{
		return ievs.size()-1;
	}
	inline backtrack_iterator rbegin(int indexOnThread=-1)
	{
		return backtrack_iterator(this,(indexOnThread>=0)?indexOnThread:(this->getCurrentVertexPointer()-1));
	}
	
	inline backtrack_iterator rend()
	{
		return backtrack_iterator(NULL,0);
	}
	
	inline GffEntry::ExonPtr getLastExon()
	{
		vector<InEdgeAndVertex>::reverse_iterator i=ievs.rbegin();
		return (i==ievs.rend())?NULL:i->exon;
	}
	
	inline int append(const InEdgeAndVertex& iev)
	{	
		ievs.push_back(iev);
		return ievs.size()-1;
	}


	inline AESPathThread* fork(InEdgeAndVertex iev)
	{
		//fork parent too?
		return new AESPathThread(this->parentExonGroupThread /*transmit it for now, will be updated*/,iev,this,this->root);
	}
	
};


class AESPathThreadExonGroup
{
public:
	AESPathThreadExonGroup* prec;
	vector<NExonGroup::NExonGroupPtr> exongroups;
	typedef pair<int,AESPathThread*> AESPathThreadSnapShot;
	typedef multimap<int,AESPathThreadSnapShot> ExonThreadMap;
	typedef ExonThreadMap::iterator ExonThreadMapI;
	typedef pair<ExonThreadMapI,ExonThreadMapI> DExonThreadMapI;
	ExonThreadMap exonThreads; //holds the exon threads having the same exon group path
	inline AESPathThreadExonGroup(AESPathThreadExonGroup* _prec=NULL): prec(_prec){}
	
	inline bool hasPrev()
	{
		if(prec==NULL&&exongroups.size()==1)
			return false;
		return true;
	}
	
	class backtrack_iterator
	{
	public:
		AESPathThreadExonGroup* thr;
		int indexOnThread;
		inline backtrack_iterator(AESPathThreadExonGroup* _thr,int _indexOnThread): thr(_thr),indexOnThread(_indexOnThread)
		{}
		inline void operator++()
		{
			if(indexOnThread>0)
				indexOnThread--;
			else
			{
				thr=thr->prec;
				if(thr)
					indexOnThread=thr->exongroups.size()-1;
				else
					indexOnThread=0;
			}
		}
		inline bool operator==(const backtrack_iterator & right)
		{
			return thr==right.thr&&indexOnThread==right.indexOnThread;
		}
		inline bool operator!=(const backtrack_iterator & right)
		{
			return !(*this==right);
		}
		inline void operator++(int)
		{
			++(*this);
		}
		inline const NExonGroup::NExonGroupPtr operator*()
		{
			return thr->exongroups[indexOnThread];
		}
		inline const NExonGroup::NExonGroupPtr operator->()
		{
			return thr->exongroups[indexOnThread];
		}
		inline AESPathThreadExonGroup* getThread()
		{
			return thr;
		}
		inline bool hasAtLeastOneMore()
		{
			if(thr->prec==NULL&&thr->exongroups.size()==1)
				return false;
			return true;
		}
	};
	
	inline backtrack_iterator rbegin(int indexOnThread=-1)
	{
		return backtrack_iterator(this,(indexOnThread>=0)?indexOnThread:(exongroups.size()-1));
	}
	
	inline backtrack_iterator rend()
	{
		return backtrack_iterator(NULL,0);
	}
	
	
	
	inline int append(NExonGroup::NExonGroupPtr neg)
	{
		exongroups.push_back(neg);
		return exongroups.size()-1;
	}
	inline void addExonThread(AESPathThread* exonThreadToAdd)
	{
		int indexOnThread=exongroups.size()-1;
		
		exonThreads.insert(ExonThreadMap::value_type(indexOnThread,AESPathThreadSnapShot(exonThreadToAdd->getCurrentVertexPointer(),exonThreadToAdd)));
		exonThreadToAdd->parentExonGroupThread=this;
	}	
	inline AESPathThreadExonGroup(NExonGroup::NExonGroupPtr neg,AESPathThreadExonGroup* _prec=NULL,AESPathThread* founderExonThread=NULL):prec(_prec)
	{
		append(neg);
		if(founderExonThread)
			addExonThread(founderExonThread);
	}
	inline void registerExonThreads(const vector<AESPathThread*>& threadv)
	{
		for(vector<AESPathThread*>::const_iterator i=threadv.begin();i!=threadv.end();i++)
		{
			addExonThread(*i);
		}
	}
	inline AESPathThreadExonGroup(NExonGroup::NExonGroupPtr neg,AESPathThreadExonGroup* _prec,const vector<AESPathThread*>&threadv):prec(_prec)
	{
		
		append(neg);
		registerExonThreads(threadv);
	}
	inline AESPathThreadExonGroup* fork(NExonGroup::NExonGroupPtr neg,vector<AESPathThread*>&threadv)
	{
		return new AESPathThreadExonGroup(neg,this,threadv);
	}
	
	inline AESPathThreadExonGroup* fork(NExonGroup::NExonGroupPtr neg,AESPathThread* founderExonThread)
	{
		return new AESPathThreadExonGroup(neg,this,founderExonThread);
	}
	

	
};

inline AESPathThread::AESPathThread(AESPathThreadExonGroup* _parentExonGroupThread,AESPathThread* _prec,AESPathThread* _root):parentExonGroupThread(_parentExonGroupThread),prec(_prec)
{
	//parentExonGroupThread->addExonThread(this);
	root=(_root)?_root:this;
}

inline AESPathThread::AESPathThread(AESPathThreadExonGroup* _parentExonGroupThread,const InEdgeAndVertex& iev, AESPathThread* _prec,AESPathThread* _root):parentExonGroupThread(_parentExonGroupThread), prec(_prec)
{
	//parentExonGroupThread->addExonThread(this);
	root=(_root)?_root:this;
	append(iev);
}


class ExonPostItem
{
public:
	AESPathThread* exonThread; //to the exonThread triggering the post
	int vertexPointer;  //to the vertex in which the post is at
	int time;
	inline ExonPostItem(AESPathThread* _exonThread,int _vertexPointer,int _time): exonThread(_exonThread),vertexPointer(_vertexPointer),time(_time)
	{
		
	}
};


class ExonGroupPostItem
{
public:
	AESPathThreadExonGroup* exonGroupThread; //to the exonGroup thread triggering the post
	int vertexPointer; //to the vertex in which the post is at
	int time;
	inline ExonGroupPostItem(AESPathThreadExonGroup* _exonGroupThread, int _vertexPointer, int _time): exonGroupThread(_exonGroupThread), vertexPointer(_vertexPointer), time(_time){}
};

class ExonPost: public vector<ExonPostItem> //Exon post in one time at one exon
{
public:
	typedef iterator I;
};

class ExonGroupPost:public vector<ExonGroupPostItem> //ExonGroup Post in one time at one exongroup
{
public:
	typedef iterator I;
};



typedef vector<ExonPostItem*> SrcExonPostStruct_S3;
typedef map<int, SrcExonPostStruct_S3* > SrcExonPostStruct_S2;
typedef map<GffEntry::ExonPtr, SrcExonPostStruct_S2* > SrcExonPostStruct_S1;

class SrcExonPostStruct: public SrcExonPostStruct_S1
{
public:
	
	typedef SrcExonPostStruct_S1 S1;
	typedef SrcExonPostStruct_S2 S2;	
	typedef SrcExonPostStruct_S3 S3;	
	typedef S1::iterator I1;
	typedef S2::iterator I2;
	typedef S3::iterator I3;
	inline ~SrcExonPostStruct()
	{
		for(I1 i1=begin();i1!=end();i1++)
		{
			S2* s2=i1->second;
			for(I2 i2=s2->begin();i2!=s2->end();i2++)
			{
				S3* s3=i2->second;
				delete s3;
			}
			
			delete s2;
		}
	}
	
	inline void add(GffEntry::ExonPtr exon,int time,ExonPostItem* item)
	{
		I1 i1=find(exon);
		S2* s2;
		S3* s3;
		if(i1==end())
		{
			s2=new S2;
			insert(value_type(exon,s2));
		}else
			s2=i1->second;
		
		I2 i2=s2->find(time);
		if(i2==s2->end())
		{
			s3=new S3;
			s2->insert(S2::value_type(time,s3));
		}else
			s3=i2->second;
		
		s3->push_back(item);
	}
};

class ExonPostRow: public map<int,ExonPost*> //time-> Exon Post in one exon
{
public:
	typedef map<int,ExonPost*>::iterator I;
	
	inline void postItem(int time,const ExonPostItem& item)
	{
		
		I i=find(time);
		ExonPost* post;
		if(i==end())
		{
			post=new ExonPost;
			insert(value_type(time,post));
		}else
			post=i->second;
		
		post->push_back(item);
	}
	inline ~ExonPostRow()
	{
		for(I i=begin();i!=end();i++)
		{
			delete i->second;
		}
	}
	
	inline SrcExonPostStruct* getSrcExonStructure()
	{
		SrcExonPostStruct* newStruct=new SrcExonPostStruct;
		
		for(I i=begin();i!=end();i++)
		{
			int time=i->first;
			ExonPost* exonpost=i->second;
			for(ExonPost::I j=exonpost->begin();j!=exonpost->end();j++)
			{
				ExonPostItem *epi=&(*j);
				AESPathThread* exonThread=epi->exonThread;
				AESPathThread* rootThread=exonThread->root;
				GffEntry::ExonPtr srcExon=rootThread->ievs[0].exon;
				newStruct->add(srcExon,time,epi);
			}
		}
		return newStruct;
	}
};

class ExonGroupPostRow: public map<int,ExonGroupPost*> //time-> ExonGroup Post in one exongroup
{
public:
	typedef map<int,ExonGroupPost*>::iterator I;
	inline void postItem(int time,const ExonGroupPostItem& item)
	{
		
		I i=find(time);
		ExonGroupPost* post;
		if(i==end())
		{
			post=new ExonGroupPost;
			insert(value_type(time,post));
		}else
			post=i->second;
		
		post->push_back(item);
	}	
	inline ~ExonGroupPostRow()
	{
		for(I i=begin();i!=end();i++)
		{
			delete i->second;
		}
	}
};

class ExonPostMatrix: public map<GffEntry::ExonPtr, ExonPostRow*>
{
public:
	typedef map<GffEntry::ExonPtr,ExonPostRow*>::iterator I;
	typedef map<GffEntry::ExonPtr,ExonPostRow*>::value_type V;
	typedef pair<ExonPostRow*,bool> IStat;
	
	IStat getRow(GffEntry::ExonPtr exon)
	{
		ExonPostRow* row;
		bool newRow=false;
		I i=find(exon);
		if(i==end())
		{
			newRow=true;
			row=new ExonPostRow;
			insert(V(exon,row));
		}else
			row=i->second;
		
		return IStat(row,newRow);
	}
	inline ~ExonPostMatrix()
	{
		for(I i=begin();i!=end();i++)
		{
			delete i->second;
		}
	}	
};

class ExonGroupPostMatrix: public map<NExonGroup::NExonGroupPtr, ExonGroupPostRow*>
{
public:
	typedef map<NExonGroup::NExonGroupPtr, ExonGroupPostRow*>::iterator I;
	typedef map<NExonGroup::NExonGroupPtr, ExonGroupPostRow*>::value_type V;
	typedef pair<ExonGroupPostRow*,bool> IStat;
	IStat getRow(NExonGroup::NExonGroupPtr exongroup)
	{
		ExonGroupPostRow* row;
		bool newRow=false;
		I i=find(exongroup);
		if(i==end())
		{
			newRow=true;
			row=new ExonGroupPostRow;
			insert(V(exongroup,row));
		}else
			row=i->second;
		
		return IStat(row,newRow);
	}
	inline ~ExonGroupPostMatrix()
	{
		for(I i=begin();i!=end();i++)
		{
			delete i->second;
		}
	}	
};


#define LIVE_FOREVER INT_MAX

class SpliceTraversalGraph
{
public:
	vector<AESPathThread*> exonThreads;
	vector<AESPathThreadExonGroup*> exonGroupThreads;
	list<AESPathThread*> *activeExonThreads;
	ExonPostMatrix exonPostMatrix;
	ExonGroupPostMatrix exonGroupPostMatrix;
	string locusName;
	bool checkLocusName;
	/* some virtual function here for overloading allowing for specific actions
	 */
	
	virtual void startGraph()
	{
		
	}
	virtual void discoverExon(const ExonPostItem& curItem)
	{
		
	}
	virtual void discoverExonGroup(const ExonGroupPostItem& curItem)
	{
		
	}
	virtual void arriveVisitedExon(const ExonPostItem&curItem, ExonPostRow* prevItems)
	{
		
	}
	virtual void arriveVisitedExonGroup(const ExonGroupPostItem&curItem, ExonGroupPostRow* prevItems)
	{
		
	}
	virtual void endGraph()
	{
		
	}
	virtual void endPath(AESPathThread* exonThread)
	{
		
	}
	
	int lifeSpan;
	int curStep;
	//initialize the splice traversal graph starting at neg with specified lifespan

	
	inline void visitExon(GffEntry::ExonPtr exon, AESPathThread* exonThread,int indexOnThread,int time)
	{
		//record in the ExonPost
		ExonPostMatrix::IStat iStat=exonPostMatrix.getRow(exon);

		ExonPostItem curItem(exonThread,indexOnThread,time);
		
		if(iStat.second) //new Item
		{
			this->discoverExon(curItem);
		}
		else
		{
			this->arriveVisitedExon(curItem,iStat.first);
		}
		
		ExonPostRow* row=iStat.first;
		row->postItem(time,curItem);		
		//trigger action the virtual functions
		
	}
	
	inline void visitExonGroup(NExonGroup::NExonGroupPtr exongroup,AESPathThreadExonGroup* exonGroupThread,int indexOnThread,int time)
	{
		//Debugerr<<"veg a"<<endl;
		ExonGroupPostMatrix::IStat iStat=exonGroupPostMatrix.getRow(exongroup);
		//Debugerr<<"veg b"<<endl;
		ExonGroupPostItem curItem(exonGroupThread,indexOnThread,time);
		//Debugerr<<"veg c"<<endl;
		if(iStat.second) //new item
		{
			this->discoverExonGroup(curItem);
		}
		else
		{
			this->arriveVisitedExonGroup(curItem,iStat.first);
		}
		//Debugerr<<"veg d"<<endl;
		ExonGroupPostRow* row=iStat.first;
		//Debugerr<<"veg e"<<endl;
		row->postItem(time,curItem);
		//Debugerr<<"veg f"<<endl;
	}
	
	
	typedef vector<AESPathThread*> ThreadV;
	typedef map<NExonGroup::NExonGroupPtr, ThreadV*>   OGSMap1;
	typedef map<AESPathThreadExonGroup*, OGSMap1*> OGSMap2;
	
	
	class OutGroupStat: public OGSMap2
	{
	public:
		
		inline void registerThread(NExonGroup::NExonGroupPtr toGroup,AESPathThread* exonThread)
		{
			//Debugerr<<"rt a"<<endl;
			OGSMap1* innerMap;
			ThreadV* threadv;
			//Debugerr<<"rt b"<<endl;
			AESPathThreadExonGroup* egt=exonThread->parentExonGroupThread;
			//Debugerr<<"rt c"<<endl;
			iterator i=find(egt);
			//Debugerr<<"rt d"<<endl;
			if(i==end())
			{
				innerMap=new OGSMap1;
				insert(value_type(egt,innerMap));
			}else
				innerMap=i->second;
			//Debugerr<<"rt e:"<<toGroup<<endl;
			//Debugerr<<"rt e1:"<<toGroup->getBound()<<endl;
			OGSMap1::iterator j=innerMap->find(toGroup);
			if(j==innerMap->end())
			{
			    threadv=new ThreadV;
				innerMap->insert(OGSMap1::value_type(toGroup,threadv));
			}else
				threadv=j->second;
			//Debugerr<<"rt f"<<endl;
			threadv->push_back(exonThread);
			//Debugerr<<"rt z"<<endl;
		}
		
		inline ~OutGroupStat()
		{
			for(iterator i=begin();i!=end();i++)
			{
				OGSMap1* map1=i->second;
				for(OGSMap1::iterator j=map1->begin();j!=map1->end();j++)
				{
					ThreadV* threadv=j->second;
					delete threadv;
				}
				
				delete map1;
			}
		}
	};
	
	inline bool nextStep()
	{
		if(curStep>=lifeSpan || !this->activeExonThreads || this->activeExonThreads->size()==0)
			return false;
		
		curStep++;
		
		//cerr<<"ns a:"<<curStep<<endl;
		
		//construct new list for next round;
		list<AESPathThread*> *newThreads=new list<AESPathThread*>;
		OutGroupStat outGroupStat;
	
		//retrieve the next task from the active exon thread
		for(list<AESPathThread*>::iterator li=this->activeExonThreads->begin();li!=this->activeExonThreads->end();li++)
		{
			//Debugerr<<"ns b"<<endl;		
			AESPathThread* thr=*li;
			GffEntry::ExonPtr fromExon=thr->getLastExon();
			//Debugerr<<"ns c"<<endl;	
			//now the branch out!
			
			//cerr<<"trav: from"<<fromExon->egsid<<endl;
			int numOutJnx=(fromExon->outJnxs)?fromExon->outJnxs->size():0;
			//cerr<<"njnx: "<<numOutJnx<<endl;
			
			if(numOutJnx==0)
			{
				//cerr<<"ns e"<<endl;
				endPath(thr);
				//end of transcript in this path
			}
			else if(numOutJnx==1)
			{
				//Debugerr<<"ns f>"<<endl;
				map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator iJ=fromExon->outJnxs->begin();
				//no need to branch exonThread
				//Debugerr<<"ns g"<<endl;
				GffEntry::ExonPtr toExon=iJ->first;
				GffEntry::Jnx* outEdge=iJ->second;
				if(!toExon->exonGroup)
				{
					
					//cerr<<"warning: exon group not valid: could be same fromExon used by different Locus"<<endl;

					endPath(thr);
					
					continue;
				}				
				if(this->checkLocusName && toExon->exonGroup->locusName!=this->locusName) //stupid same fromExon used by different locus
				{
					cerr<<"locus name not equal"<<endl;
					endPath(thr);
					continue;
				}
				
				//cerr<<"\tto Exon:"<<toExon->egsid<<endl;
				//Debugerr<<"ns h"<<endl;
				int index=thr->append(InEdgeAndVertex(outEdge,toExon));
				//add the thr again to the new active list
				//Debugerr<<"ns i"<<endl;
				newThreads->push_back(thr);
				//now register this ExonGroup to the next ExonGroup transition
				//Debugerr<<"ns j:"<<toExon->getBound()<<endl;
				
				outGroupStat.registerThread(toExon->exonGroup,thr);
				// visiting the exon post and also triggering the exon visit event
				//Debugerr<<"ns k"<<endl;
				this->visitExon(toExon,thr,index,curStep);
				//Debugerr<<"<ns l"<<endl;
				
			}
			else //numOutJnx>1 -> goes to different exons, need to branch off
			{
				//cerr<<"ns m>"<<endl;

				bool nothingActually=true;
				for(map<GffEntry::ExonPtr,GffEntry::Jnx*>::iterator ji=fromExon->outJnxs->begin();ji!=fromExon->outJnxs->end();ji++)
				{
					
					GffEntry::ExonPtr toExon=ji->first;
					GffEntry::Jnx* outEdge=ji->second;
					
					if(!toExon->exonGroup)
					{
						//cerr<<"warning: exon group not valid: could be same fromExon used by different Locus"<<endl;
						continue;
					}
					if(this->checkLocusName && toExon->exonGroup->locusName!=this->locusName) //stupid same fromExon used by different locus
					{
						//endPath(thr);
						continue;
					}
					
					//cerr<<"\tto Exon:"<<toExon->egsid<<endl;
					
					nothingActually=false;
					//branch exonThread
					//Debugerr<<"ns n"<<endl;
					AESPathThread* forkThread=thr->fork(InEdgeAndVertex(outEdge,toExon));
					//add the new branched thread to the new active threads list
					//Debugerr<<"ns o"<<endl;
					newThreads->push_back(forkThread);
					//now register this ExonGruop to the next ExonGroup transition
					//Debugerr<<"ns p:";
					//Debugerr<<toExon->egsid<<" "<<toExon->getBound()<<endl;
					outGroupStat.registerThread(toExon->exonGroup,forkThread);
					//to add: vistiing the exon post and also triggering the exon visit event
					//Debugerr<<"ns q"<<endl;
					this->visitExon(toExon,forkThread,0,curStep);
					
				}
				
				
				if(nothingActually)
				{
					//cerr<<"nothing actually"<<endl;
					endPath(thr);
					continue;
				}
				
				//Debugerr<<"<ns r"<<endl;
			}
			
			
			
		}
		
		
		//Debugerr<<"s"<<endl;
		//now see how many exon group threads we need
		//go through the outGroupStat, and create the appropriate exon group threads (by forking? or appending?)
		
		//vistting the exongroup post and also triggering the exon group visit event
		for(OutGroupStat::iterator i=outGroupStat.begin();i!=outGroupStat.end();i++)
		{
			AESPathThreadExonGroup* egt=i->first;
			OGSMap1* innerMap=i->second;
			int numOGroup=innerMap->size();
			//Debugerr<<"t:"<<numOGroup<<endl;
			if(numOGroup==1) //numOGroup = 1 -> no need to fork
			{
				//Debugerr<<"u"<<endl;
				NExonGroup::NExonGroupPtr toExonGroup=innerMap->begin()->first;
				//that thread is already in the exon group thread, no need to add to the exongroup thread
				//Debugerr<<"v"<<endl;
				int index=egt->append(toExonGroup);
				ThreadV* threadv=(innerMap->begin())->second;
				egt->registerExonThreads(*threadv);
				//register exon group post and trigger visit event
				//Debugerr<<"w"<<endl;
				this->visitExonGroup(toExonGroup,egt,index,curStep);
				//Debugerr<<"x"<<endl;
			}
			else //numOGroup >1 need to fork!
			{
				for(OGSMap1::iterator j=innerMap->begin();j!=innerMap->end();j++)
				{
					//Debugerr<<"y"<<endl;
					NExonGroup::NExonGroupPtr toExonGroup=j->first;
					ThreadV* threadv=j->second;
					//Debugerr<<"z"<<endl;
					AESPathThreadExonGroup* negt=egt->fork(toExonGroup,*threadv);
					//Debugerr<<"za"<<endl;
					this->exonGroupThreads.push_back(negt);
					//register exon group post and trigger visit event
					//Debugerr<<"zb"<<endl;
					this->visitExonGroup(toExonGroup,negt,0,curStep);
					//Debugerr<<"zc"<<endl;
				}
			}
		}
		
		//destroy the outGroupStat map (auto by class destructor)
		
		//Debugerr<<"zd"<<endl;
		delete this->activeExonThreads;
		this->activeExonThreads=newThreads;
		
		//Debugerr<<"ze"<<endl;
		if(curStep>=lifeSpan || !this->activeExonThreads || this->activeExonThreads->size()==0)
		{
			this->endGraph();
		}
		
		return !(curStep==lifeSpan);
	}
	
	
	inline void enterGraphLoop()
	{
		while(this->nextStep()){}
	}

	inline SpliceTraversalGraph(NESuperGroup* _sgroup,string _locusName,bool _checkLocusName,int _lifeSpan=LIVE_FOREVER):lifeSpan(_lifeSpan),curStep(0),locusName(_locusName),checkLocusName(_checkLocusName)
	{
			//cerr<<"sg a"<<endl;
			//AESPathThreadExonGroup *egt=new AESPathThreadExonGroup(neg);
			
		
			AESPathThreadExonGroup *egt;
			map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*> mexonGroupThreads;
			
			activeExonThreads=new list<AESPathThread*>;
			for(NESuperGroup::iterator i=_sgroup->begin();i!=_sgroup->end();i++)
			{
				NExonGroup* exonGroup=*i;
				
				cerr<<_locusName<<" adding in exongroup "<<exonGroup->sid<<" from sgroup "<<_sgroup->getID()<<endl;
				
				for(NExonGroup::ExonI exonI=exonGroup->levelExons.begin();exonI!=exonGroup->levelExons.end();exonI++)
				{
					
					

					GffEntry::ExonPtr exon=*exonI;

					cerr<<"\t"<<_locusName<<" add in exon "<<exon->egsid<<endl;

					//cerr<<"adding exon from startgroup :"<<exon->egsid<<endl;
					AESPathThread* et=new AESPathThread(egt,InEdgeAndVertex(NULL,exon),NULL,NULL);
				
					NExonGroup::NExonGroupPtr egp=exon->exonGroup;
				
					map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::iterator egti=mexonGroupThreads.find(egp);
				
					if(egti==mexonGroupThreads.end())
					{
						egt=new AESPathThreadExonGroup(egp);
						mexonGroupThreads.insert(map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::value_type(egp,egt));
						exonGroupThreads.push_back(egt);
					}
					else
						egt=egti->second;
					
					egt->addExonThread(et);
					
					exonThreads.push_back(et);
					activeExonThreads->push_back(et);
					this->visitExon(exon,et,0,curStep);
				}
			}
			
			for(map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::iterator i=mexonGroupThreads.begin();i!=mexonGroupThreads.end();i++)
			{
				this->visitExonGroup(i->first,i->second,0,curStep);
			}
			
			//cerr<<"sg b"<<endl;
			this->startGraph();
			//cerr<<"sg z"<<endl;
			
	}		
	inline SpliceTraversalGraph(NExonGroup::IntExonPtrMMapDI range,string _locusName,bool _checkLocusName,int _lifeSpan=LIVE_FOREVER):lifeSpan(_lifeSpan),curStep(0),locusName(_locusName),checkLocusName(_checkLocusName)
	{
			//Debugerr<<"sg a"<<endl;
			//AESPathThreadExonGroup *egt=new AESPathThreadExonGroup(neg);
			
		
			AESPathThreadExonGroup *egt;
			map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*> mexonGroupThreads;
			
			activeExonThreads=new list<AESPathThread*>;
			
			for(NExonGroup::IntExonPtrMMapI exonI=range.first;exonI!=range.second;exonI++)
			{
				
				GffEntry::ExonPtr exon=exonI->second;
				AESPathThread* et=new AESPathThread(egt,InEdgeAndVertex(NULL,exon),NULL,NULL);
				
				NExonGroup::NExonGroupPtr egp=exon->exonGroup;
				
				map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::iterator egti=mexonGroupThreads.find(egp);
				
				if(egti==mexonGroupThreads.end())
				{
					egt=new AESPathThreadExonGroup(egp);
					mexonGroupThreads.insert(map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::value_type(egp,egt));
					exonGroupThreads.push_back(egt);
				}
				else
					egt=egti->second;
				
				egt->addExonThread(et);
				
				exonThreads.push_back(et);
				activeExonThreads->push_back(et);
				this->visitExon(exon,et,0,curStep);
			}
			
			for(map<NExonGroup::NExonGroupPtr,AESPathThreadExonGroup*>::iterator i=mexonGroupThreads.begin();i!=mexonGroupThreads.end();i++)
			{
				this->visitExonGroup(i->first,i->second,0,curStep);
			}
			
			
			this->startGraph();
			//Debugerr<<"sg z"<<endl;
			
	}	
	
	inline SpliceTraversalGraph(NExonGroup::NExonGroupPtr neg,string _locusName,bool _checkLocusName,int _lifeSpan=LIVE_FOREVER):lifeSpan(_lifeSpan),curStep(0),locusName(_locusName),checkLocusName(_checkLocusName)
	{
		//Debugerr<<"sg a"<<endl;
		AESPathThreadExonGroup *egt=new AESPathThreadExonGroup(neg);
		
		activeExonThreads=new list<AESPathThread*>;
		
		for(NExonGroup::ExonI exonI=neg->levelExons.begin();exonI!=neg->levelExons.end();exonI++)
		{
			GffEntry::ExonPtr exon=*exonI;
			AESPathThread* et=new AESPathThread(egt,InEdgeAndVertex(NULL,exon),NULL,NULL);
			egt->addExonThread(et);
			exonThreads.push_back(et);
			activeExonThreads->push_back(et);
			this->visitExon(exon,et,0,curStep);
		}
		
		exonGroupThreads.push_back(egt);
		this->visitExonGroup(neg,egt,0,curStep);
	
		
		
		this->startGraph();
		//Debugerr<<"sg z"<<endl;
		
	}
	
	/*inline*/ virtual ~SpliceTraversalGraph()
	{
		for(vector<AESPathThread*>::iterator i=exonThreads.begin();i!=exonThreads.end();i++)
		{
			AESPathThread* thr=*i;
			delete thr;
		}
		
		for(vector<AESPathThreadExonGroup*>::iterator i=exonGroupThreads.begin();i!=exonGroupThreads.end();i++)
		{
			AESPathThreadExonGroup* thr=*i;
			delete thr;
		}
		
		if(this->activeExonThreads)
		{
			delete this->activeExonThreads;
		}
	}
};


#endif /*SPLICEGRAPH_H_*/
