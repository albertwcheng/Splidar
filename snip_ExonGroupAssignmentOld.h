/*
 * snip_ExonGroupAssignmentOld.h
 *
 *  Created on: Jan 26, 2010
 *      Author: awcheng
 */

#ifndef SNIP_EXONGROUPASSIGNMENTOLD_H_
#define SNIP_EXONGROUPASSIGNMENTOLD_H_

NExonGroup::ChildI collapseGroups(NExonGroup* parent,NExonGroup::ChildI startOverlap,NExonGroup::ChildI endOverlap1, GffEntry::ExonPtr exon);

bool assignExonGroupPerExon(NExonGroup* parentOfSubtree, NExonGroup* subtree, NExonGroup::ChildI& childI,GffEntry::ExonPtr exon,int level);

void assignExonGroupFirstTranscript(NExonGroup* root, GffEntry* transcript,int jobID);

void assignExonGroupPerTranscript(NExonGroup* root, GffEntry* transcript,int JobID);




bool validateNExonSubtree(NExonGroup * subtree)
{
	if(!::isValid(subtree->bound))
	{
		cerr<<"NEGError: Invalid bound of ExonGroup"<<endl;
		return false;

	}




	KeyPair<int,int> exonOverlapBound(INT_MIN,INT_MAX); //going to take max,min
	KeyPair<int,int> subtreeBound(INT_MAX,INT_MIN);     //going to take min,max
	//now validate levelExons that they overlaps each other and they report the true bound of this subtree


	KeyPair<int,int> levelExonUnionBound=subtree->bound;

	if(subtree->levelExons.size()>0)
	{

		ExonBoundGroup *leftexonboundgroup=NULL;
		ExonBoundGroup *rightexonboundgroup=NULL;

		for(NExonGroup::ExonI exonI=subtree->levelExons.begin();exonI!=subtree->levelExons.end();exonI++)
		{
			GffEntry::Exon* curExon=*exonI;
			overlapBoundMutFirst(exonOverlapBound,curExon->getBound());  //all 01
			expandBoundMutFirst(subtreeBound,curExon->getBound());

			//validate whether all level Exons have the same exonbound group

			ExonBoundGroup* curLeftExonBoundGroup=curExon->leftExonBoundGroup;
			ExonBoundGroup* curRightExonBoundGroup=curExon->rightExonBoundGroup;

			if(!curLeftExonBoundGroup)
			{
				cerr<<"NEGError: left exon bound group not set for exon:"<<subtree->locusName<<"/"<<subtree->sid<<endl;
				return false;
			}

			if(!curRightExonBoundGroup)
			{
				cerr<<"NEGError: right exon bound group not set for exon:"<<subtree->locusName<<"/"<<subtree->sid<<endl;
			}

			if(!leftexonboundgroup)
			{
				leftexonboundgroup=curLeftExonBoundGroup;
				rightexonboundgroup=curRightExonBoundGroup;
			}
			else
			{
				if(leftexonboundgroup!=curLeftExonBoundGroup)
				{
					cerr<<"NEGError: left exon bound group not consistent for exon:"<<subtree->locusName<<"/"<<subtree->sid<<endl;
					return false;
				}

				if(rightexonboundgroup!=curRightExonBoundGroup)
				{
					cerr<<"NEGError: right exon bound group not consistent for exon:"<<subtree->locusName<<"/"<<subtree->sid<<endl;
					return false;
				}

			}
		}

		levelExonUnionBound=subtreeBound;
	}

	//check if overlapBound is still ok
	if(exonOverlapBound.k1<exonOverlapBound.k2)
	{
		//correct
	}
	else
	{
		cerr<<"NEGError: some exons in levelExons don't overlap "<<exonOverlapBound<<endl;
		return false;
	}



	//now validate children trees
	if(subtree->hasChildren())
	{

		NExonGroup* neg;
		NExonGroup::ChildI childI=subtree->children->begin();
		neg=*childI;
		KeyPair<int,int> prevBoundChild=neg->bound;

		if(!areOverlapping(levelExonUnionBound,prevBoundChild))
		{
			cerr<<"NEGError: Child do not overlap with leveExons "<<levelExonUnionBound<<" vs "<<prevBoundChild<<endl;
			return false;
		}
		//overlapBoundMutFirst(exonOverlapBound,prevBoundChild);  //all 01
		expandBoundMutFirst(subtreeBound,prevBoundChild);
		childI++;


		if(!validateNExonSubtree(neg))
			return false;

		for(;childI!=subtree->children->end();childI++)
		{
			neg=*childI;
			KeyPair<int,int> curBound=neg->bound;

			if(!areOverlapping(levelExonUnionBound,curBound))
			{
				cerr<<"NEGError: Child do not overlap with leveExons "<<levelExonUnionBound<<" vs "<<curBound<<endl;
				return false;
			}
			//overlapBoundMutFirst(exonOverlapBound,curBound);  //all 01
			expandBoundMutFirst(subtreeBound,curBound);

			if(prevBoundChild.k1<=curBound.k1 && prevBoundChild.k2<=curBound.k1)
			{
				//correct
				if(!validateNExonSubtree(neg))
					return false;

				prevBoundChild=curBound;
			}
			else
			{
				//wrong
				cerr<<"NEGError: child bound overlap: "<<prevBoundChild<<" vs "<<curBound<<endl;
				return false;
			}



		}



	}


	//now check if the expanded bound is the reported bound

	if(subtreeBound!=subtree->bound)
	{
		cerr<<"NEGError: expanded bound is not the reported bound "<<subtreeBound<<" vs "<<subtree->bound<<endl;
		return false;
	}


	//pass every test
	return true;


}
bool validateNExonGroupAssignment(GffEntry::Locus* locus,int jobID)
{
	if(!locus->root)
	{
		cerr<<"NEGError: No Exon Groups assigned"<<endl;
		return false;
	}

	//each exon in locus has to be visited with the flag jobID!

	for(int i=0;i<locus->exonCount;i++)
	{
		GffEntry::ExonPtr exon=locus->exons[i];

		if(exon->usageFlag!=jobID)
		{
			cerr<<"NEGError: Exon not visited by job"<<endl;
			return false;
		}
		if(exon->egsid=="")
		{
			cerr<<"NEGError: Exon not finalized with an egsid"<<endl;
			return false;
		}
		if(!exon->exonGroup)
		{
			cerr<<"NEGError: Exon group not given to exon"<<endl;
		}
	}

	//now validate the tree

	return validateNExonSubtree(locus->root);



}

void assignExonGroupPerLocus(GffEntry::Locus* locus,int jobID, bool validateAssignment)
{
	//foreach transcript in locus
//	cerr<<"a"<<endl;
	if(locus->root) //already has exon tree
		return;
//	cerr<<"b"<<endl;
	locus->root= new NExonGroup;

	bool firstTranscript=true;


	for(vector<GffEntry*>::iterator tI=locus->transcripts.begin();tI!=locus->transcripts.end();tI++)
	{
//		cerr<<"c"<<endl;



		GffEntry* transcript=(*tI);

		//ignore single exon transcript
		//if(transcript->exonCount<2) //why?
		//	continue;

		for(int ei=0;ei<transcript->exonCount;ei++)
		{
			GffEntry::ExonPtr exon=transcript->exons[ei];
			if(exon->visited)
					continue; // counted already ignore;

			exon->visited=true;

			for(vector<GffEntry::GBlockPtr>::iterator bi=exon->blocks->begin();bi!=exon->blocks->end();bi++)
			{
				GffEntry::GBlockPtr block=(*bi);
				block->usageFrequency++;
			}
		}

//		cerr<<"d"<<endl;
	//	cerr<<"process transcript: "<<transcript->name<<"\t"<<transcript->chrom<<":"<<transcript->exons[0]->getStart1()<<"-"<<transcript->exons[transcript->exonCount-1]->getEnd1()<<endl;

		if(firstTranscript)
		{
//			cerr<<"e1"<<endl;
			assignExonGroupFirstTranscript(locus->root,transcript,jobID);
			firstTranscript=false;
		}
		else
		{
//			cerr<<"e2"<<endl;
			assignExonGroupPerTranscript(locus->root,transcript,jobID);
		}

	//	printExonTree(cerr,locus);
	//	cerr<<"***************";

	}

	//finish for this locus, now finalize to label the exon groups and exons with ids
	locus->root->finalizeExonTree(locus,*locus->names.begin());

	//now finalize the new exon bound group
	locus->assignExonBoundGroups();




	if(validateAssignment)
	{

		if(::validateNExonGroupAssignment(locus,jobID))
		{
			cerr<<"Validation returns success"<<endl;
		}
		else
		{
			cerr<<"Validation failed"<<endl;
		}
	}
}

void assignExonGroupPerTranscript(NExonGroup* root, GffEntry* transcript,int jobID)
{
//	cerr<<"f2"<<endl;
	if(!root->hasChildren())
	{
		cerr<<"Has not Children!"<<endl;
		return;
	}
	NExonGroup::ChildI topLevelI=root->children->begin();

//	cerr<<"g2"<<endl;
	for(int ei=0;ei<transcript->exonCount;ei++)
	{
		GffEntry::ExonPtr exon=transcript->exons[ei];
		if(exon->usageFlag==jobID)
			continue;




		exon->usageFlag=jobID;

	//	cerr<<"process exon "<<exon->chr<<":"<<exon->getStart1()<<"-"<<exon->getEnd1()<<endl;

		assignExonGroupPerExon(NULL,root,topLevelI,exon,0);
	}
//	cerr<<"h2"<<endl;
}

void assignExonGroupFirstTranscript(NExonGroup* root, GffEntry* transcript,int jobID)
{
//	cerr<<"f1"<<endl;
	//just add everything as single groups into root
	for(int ei=0;ei<transcript->exonCount;ei++)
	{

		GffEntry::ExonPtr exon=transcript->exons[ei];
		root->insertChild(new NExonGroup(exon));
		if(exon->usageFlag==jobID)
		{
			cerr<<"STRANGE!!! exon has already been used in assigning group: "<<transcript->name<<"/"<<transcript->name2<<"//"<<exon->chr<<":"<<exon->getStart1()<<"-"<<exon->getEnd1()<<endl;
			continue;
		}
		exon->usageFlag=jobID;
	}
//	cerr<<"g1"<<endl;
}

//return whether the childI of the upper thread need to be updated.
bool assignExonGroupPerExon(NExonGroup* parentOfSubtree, NExonGroup* subtree, NExonGroup::ChildI& childI,GffEntry::ExonPtr exon,int level)
{
	//cerr<<"level "<<level<<endl;
//	cerr<<"h"<<endl;

	if(subtree->levelExons.size()>0)
	{
		if(!areOverlapping(exon->getBound(),subtree->getBound()))
		{

			if(parentOfSubtree)
			{
				cerr<<"Strange: should not happen with non-root"<<endl;
			}
			else
			{
			//	cerr<<"aaa"<<endl;
				NExonGroup *oldRoot=new NExonGroup(*subtree);
				NExonGroup *newGroup=new NExonGroup(exon);
				//just apply for a new children because the whole children has been passed to oldRoot
				subtree->children=new set<NExonGroup::NExonGroupPtr>;
				subtree->levelExons.clear();
				subtree->bound.k1=INT_MAX;
				subtree->bound.k2=INT_MIN;
				subtree->insertChild(oldRoot);
				childI=subtree->insertChild(newGroup);
			//	cerr<<"bbb"<<endl;
			/*?*/ return true;
			}
		}
	}

	subtree->expandBoundBy(exon);

	if(subtree->hasChildren()) //is a nodex
	{
		/* go through each children
		 *  case 1: overlaps only 1 children: continue to search from that children
		 *  case 2: overlaps no children: add a new child
		 *  case 3: overlaps > 1 children but not all: collapse children together
		 *  case 4: overlaps all children: add to level Exons
		 */

		unsigned int nOverlaps=0; //CHANGED unsigned-signed
//		cerr<<"i"<<endl;
		while(childI!=subtree->children->end() && isSmaller((*childI)->getBound(),exon->getBound()))
		{
			childI++;
		}
//		cerr<<"j"<<endl;
		if(childI==subtree->children->end()) //overlaps no children
		{
			NExonGroup::ChildI newChildI=subtree->insertChild(new NExonGroup(exon));
			childI=newChildI;
			return false;
		}

		//now ChildI should points to something that overlap or pass.
//		cerr<<"k"<<endl;
		if(areOverlapping((*childI)->getBound(),exon->getBound()))
		{
			nOverlaps++;
			NExonGroup::ChildI startOverlap=childI;
			NExonGroup::ChildI itChild=childI;
//			cerr<<"l1"<<endl;

			itChild++; //new

			while(itChild!=subtree->children->end() && areOverlapping((*itChild)->getBound(),exon->getBound()))
			{
				itChild++;
				nOverlaps++; //new bug 10/11/08
			}
//			cerr<<"l2"<<endl;

			if(nOverlaps==1)
			{
				//search Deeper (recursively)
					NExonGroup::ChildI itChildR;
					NExonGroup* childTree=(*childI);
					if(childTree->hasChildren())
						itChildR=childTree->children->begin();

					bool transmitChildI=assignExonGroupPerExon(subtree, childTree, itChildR, exon, level+1);
					if(transmitChildI)
						childI=itChildR;
			}
			else if(nOverlaps==subtree->children->size()) //overlap all children
			{
				//cerr<<"insert higher level"<<endl;
				subtree->insertExon(exon);

			}
			else if(nOverlaps>1) // overlap >1 but not all children
 			{
				//cerr<<"collapse Groups"<<endl;
				childI=collapseGroups(subtree,startOverlap,itChild,exon);
			}else //if(nOverlaps==1)
			{
				cerr<<"Should not happen at all"<<endl;
			}
//			cerr<<"l3"<<endl;
		}
		else

		{ //in between groups
//			cerr<<"l1a"<<endl;
			childI=subtree->insertChild(new NExonGroup(exon));
		}
//		cerr<<"m"<<endl;
	}
	else // is a leaf, only has levelExons but not children NExonGroups
	{
	//	cerr<<"n"<<endl;
		NExonGroup *NOE=NULL;
		bool allOverlap=true;
		vector<GffEntry::ExonPtr> OE; //overlap set;
	//	cerr<<"UsingSubTree:"<<subtree<<endl;
		for(NExonGroup::ExonI ei=subtree->levelExons.begin();ei!=subtree->levelExons.end();ei++)
		{
	//		cerr<<"n1"<<endl;

			GffEntry::Exon* thisExon=(*ei);

	//		cerr<<"n1a"<<endl;


			if(!thisExon)
				cerr<<"ei=NULL"<<endl;

	//		cerr<<"n2a:"<<thisExon<<endl;

	//		cerr<<"n3p"<<endl;
	//		cerr<<"exonID="<<thisExon->exonID<<endl;
			KeyPair<int,int> b2=thisExon->getBound();
	//		cerr<<"n3a"<<endl;

			KeyPair<int,int> b1=exon->getBound();
	//		cerr<<"n4a"<<endl;

			if(areOverlapping(b1,b2))
			{
	//			cerr<<"n2"<<endl;
				OE.push_back(thisExon);

	//			cerr<<"n3"<<endl;
			}else
			{
				allOverlap=false;
				if(NOE==NULL)
				{

	//				cerr<<"n4"<<endl;
					NOE=new NExonGroup;
	//				cerr<<"n5"<<endl;
				}

			//	cerr<<"n6p"<<endl;
				NOE->insertExon(thisExon);
			//	cerr<<"n6"<<endl;
			}
		}
	//	cerr<<"o"<<endl;
		if(allOverlap)
		{
			subtree->insertExon(exon);

		}
		else
		{
			NExonGroup* DOE=new NExonGroup;
			NExonGroup* OEG=new NExonGroup(exon);

			for(vector<GffEntry::ExonPtr>::iterator i=OE.begin();i!=OE.end();i++)
			{
				if(areOverlapping((*i)->getBound(),NOE->getBound()))
				{
					DOE->insertExon(*i);
				}else
				{
					OEG->insertExon(*i);
				}
			}

			DOE->insertChild(OEG);
			DOE->insertChild(NOE);
			parentOfSubtree->eraseChild(subtree);
			//cerr<<"deleteSubTree "<<subtree->bound.k1<<"-"<<subtree->bound.k2<<endl;
			delete subtree;
			childI=parentOfSubtree->insertChild(DOE);
			return true;

		}
	//	cerr<<"p"<<endl;
	}

	return false;
}

NExonGroup::ChildI collapseGroups(NExonGroup* parent,NExonGroup::ChildI startOverlap,NExonGroup::ChildI endOverlap1, GffEntry::ExonPtr exon)
{
	//cerr<<"collapse group"<<endl;
//	cerr<<"q"<<endl;
	//printExonTree(cerr,parent);
	//cerr<<"*******"<<endl;
	NExonGroup* newNEG=new NExonGroup(exon);
	newNEG->insertChildren(startOverlap,endOverlap1);
//	cerr<<"r"<<endl;
	parent->eraseChildren(startOverlap,endOverlap1);

	NExonGroup::ChildI newI=parent->insertChild(newNEG);
	//cerr<<"collapsed"<<endl;
	//printExonTree(cerr,parent);
	return newI;
}


void printExonTree(ostream& os,NExonGroup* root,string ident,string prefix,string chr,string geneName)
{
	os<<prefix<<geneName<<" "<<root->sid<<",superIds="<<((root->leftSuperGroup)?root->leftSuperGroup->getID():"null")<<","<<((root->rightSuperGroup)?root->rightSuperGroup->getID():"null")<<"@"<<chr<<":"<<(root->bound.k1+1)<<"-"<<root->bound.k2<<":[";
	for(NExonGroup::ExonI i=root->levelExons.begin();i!=root->levelExons.end();i++)
	{
		GffEntry::ExonPtr exon=(*i);
		os<<exon->exonGroup->relRootString<<";"<<exon->leftExonBoundGroup->name<<","<<exon->rightExonBoundGroup->name<<"="<<exon->egsid<<":"<<exon->getStart1()<<"-"<<exon->getEnd1()<<"|";

	}

	os<<"]"<<endl;

	if(root->hasChildren())
	{
		for(NExonGroup::ChildI i=root->children->begin();i!=root->children->end();i++)
		{
			printExonTree(os,(*i),ident,prefix+ident,chr,geneName);
		}
	}


}

void printExonTree(ostream& os,GffEntry::Locus* locus,string ident)
{
	printExonTree(os,locus->root,ident,"",locus->chr,(*locus->names.begin()));
}


#endif /* SNIP_EXONGROUPASSIGNMENTOLD_H_ */
