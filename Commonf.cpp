/*
 * Commonf.cpp
 *
 *  Created on: Jan 18, 2010
 *      Author: awcheng
 */

#include "snipIncludes.h"

void die_clean_up()
{
//	cerr<<"Die::Clean Up"<<endl;
//	cerr<<"GffEntry::resetEverything()"<<endl;
//	GffEntry::resetEverything();
//	cerr<<"Done Clean Up"<<endl;
}

void die_exit_clean_up()
{
	GffEntry::resetEverything();
	PyConfigFile::deinit_Python();
}
