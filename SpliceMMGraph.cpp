/*
 * SpliceMMGraph.cpp
 *
 *  Created on: Dec 1, 2009
 *      Author: awcheng
 */

#include "SpliceMMGraph.h"

 const int DenObj::Junction=DENOBJ_JUNCTION;
 const int DenObj::Gblock=DENOBJ_GBLOCK;
 const int DenObj::Terminal=DENOBJ_TERMINAL;
 const int DenObj::Collapsed=DENOBJ_COLLAPSED;

DenObjPtr DenObj::addRightEnsured(DenGraph*_parent,DenObjPtr _right)
{
	DenObjPtr working=_parent->ensureDenObj(_right);
	this->addRightDirectly(working);
	return working;
}

DenObjPtr DenObj::addLeftEnsured(DenGraph*_parent,DenObjPtr _left)
{
	DenObjPtr working=_parent->ensureDenObj(_left);
	this->addLeftDirectly(working);
	return working;
}


