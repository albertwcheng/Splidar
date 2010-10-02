/*
 * PyConfigFile.h
 *
 *  Created on: Jan 7, 2010
 *      Author: awcheng
 */

#ifndef PYCONFIGFILE_H_
#define PYCONFIGFILE_H_


#include <Python.h>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "StringUtil.h"
using namespace std;

#define PYCONFIGFILE_TMP_VAR_NAME  "____TMP____XAXSDFSD"

class PyConfigFile
{
private:
	PyObject* pyDict;
	inline void Py_RunStringSys(string command)
	{
		PyRun_String(command.c_str(),Py_file_input,pyDict,pyDict);
	}

	//static int refCount;




	inline void init()
	{
		//Py_Initialize();

		PyObject* main_module =
		   PyImport_AddModule("__main__");

		pyDict =
		   PyModule_GetDict(main_module);


		//now import all popular modules
		//PyRun_String("__sys_E_T=exc_info()",Py_file_input,pyDict,pyDict);
		//string getExceptionCode="def ___get_ex(clearEx=True):\n\tglobal __sys_E_T\n\t__sys_E_T=exc_info()\n\tif clearEx:\n\t\texc_clear()";

		//cerr<<getExceptionCode<<endl;
		//PyRun_String(getExceptionCode,Py_file_input,pyDict,pyDict);

		string initializerCode=
		string("from sys import *\n")+
		"from os import *\n"+
		"from glob import *\n"+
		"import os\n"+
		"import types\n"+
		"def ___X_X_X_toStr(x):\n\treturn str(x)\n";

		Py_RunStringSys(initializerCode);


		string printObjCode="def ___Print_Obj(x,padding=''):\n\ttry:\n\t\tretstr=''\n\t\tif type(x).__name__=='dict':\n\t\t\tretstr+='dict('+str(len(x))+')'+'\\n'\n\t\t\tfor k,v in x.items():\n\t\t\t\tretstr+=padding+'\\t'+str(k)+'='+___Print_Obj(v,padding+'\\t')+'\\n'\n\t\telif type(x).__name__=='list':\n\t\t\tretstr+='list('+str(len(x))+')'+'\\n'\n\t\t\tfor v in x:\n\t\t\t\tretstr+=padding+'\\t'+___Print_Obj(v,padding+'\\t')+'\\n'\n\t\telse:\n\t\t\tretstr+=str(x)\n\t\treturn retstr\n\texcept:\n\t\tprint >> stderr,'error occured',exc_info()";
		//string tryRecursive="def _re(x):\n\tif x<10:\n\t\t_re(x+1)\n\tprint >> stderr,x";
		//cerr<<tryRecursive<<endl;

		//runString(tryRecursive);
		//runString("_re(0)");
		//cerr<<printObjCode<<endl;
		Py_RunStringSys(printObjCode);
		Py_RunStringSys("def sources(L):\n\tfor s in L:\n\t\tprint >>stderr,'loading',s\n\t\texecfile(s,globals())");
		Py_RunStringSys("____GGG=''");
		Py_RunStringSys("___fancy___X_X_X=''");

		Py_RunStringSys("Yes=True\nNo=False\nYES=Yes\nNO=No\nyes=Yes\nno=No\nT=True\nF=False");
		pyFunc_Str=getObject("___X_X_X_toStr");

		string getExceptionCode="___EX_STR=''\ndef putExcepStr():\n\tglobal ___EX_STR\n\t___EX_STR=str(exc_info())\n\texc_clear()\n";
		//cerr<<getExceptionCode<<endl;
		Py_RunStringSys(getExceptionCode);


		//cerr<<"prev"<<endl;
		//pyFunc_Print_Obj=getObject("___Print_Obj");
		//cerr<<"done init"<<endl;
		//this->addDeref(main_module);
		//this->addDeref(pyDict);




		snapShotSysVars();
	}
	set<PyObject* > derefee;
	PyObject* pyFunc_Str;
	PyObject* pyFunc_Print_Obj;
	set<string> sysVars;
	vector<string> general_argvs;

	inline string getExceptionString()
	{
		Py_RunStringSys("putExcepStr();");
		string excepstring=getStringDirect("___EX_STR");
		if(excepstring!="(None, None, None)")
			return excepstring;
		else
			return "";

	}


	inline vector<string> getGlobals()
	{
		runString("____GGG=globals().keys()\n____GGG.sort()");

		PyObject* o1=PyDict_GetItemString(pyDict,"____GGG");
		int n=PyList_Size(o1);
		vector<string> tmp;
		for(int i=0;i<n;i++)
		{
			tmp.push_back(PyString_AS_STRING(PyList_GET_ITEM(o1,i)));
		}
		//runString("del ____GGG");
		return tmp;
	}
	inline void snapShotSysVars()
	{
		/*sysVars.clear();
		PyObject* items=NEW(PyDict_Items(pyDict));

		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			PyObject* po=NEW(PyTuple_GET_ITEM(item,0));
			string key=PyString_AS_STRING(po);
			sysVars.insert(key);
		}*/
		vector<string> tmp=this->getGlobals();
		for(vector<string>::iterator i=tmp.begin();i!=tmp.end();i++)
			sysVars.insert(*i);
	}
public:

	static void initPython()
	{
		Py_Initialize();
	}

	static void deinit_Python()
	{
		Py_Finalize();
	}

	struct file_not_found {
		string filename;
		file_not_found( const string& filename_ = string() )
			: filename(filename_) {} };
	struct key_not_found {  // thrown only by T read(key) variant of read()
		string key;
		key_not_found( const string& key_ = string() )
			: key(key_) {} };
	struct python_error{
		string errlog;
		python_error(const string & errlog_=string())
			:errlog(errlog_){}
	};

	inline PyConfigFile()
	{


		init();

	}
	inline PyConfigFile(int argc,const char**argv)
	{
		init();
		this->loadFromArgs(argc,argv);
	}

	inline PyConfigFile(string filename)
	{
		init();
		loadPythonFile(filename);
	}


	inline void runFile(FILE* exp_file,string filename)
	{
		//cerr<<"****run file"<<filename<<endl;
		PyRun_File(exp_file,filename.c_str(),Py_file_input,pyDict,pyDict);
		string excString=this->getExceptionString();
		//cerr<<"***RUN FILE "<<filename<<" got exc string "<<excString<<endl;
		if(excString!="")
			throw python_error(excString);
	}

	inline void loadPythonFile(string filename)
	{
		FILE*       exp_file;
		exp_file = fopen(filename.c_str(), "r");
		if (!exp_file)
			throw file_not_found(filename.c_str());

		runFile(exp_file,filename.c_str());
		fclose(exp_file);
	}




	inline string escape(string value)
	{
		return StringUtil::escape(value);
	}

	inline void addObj(string name,PyObject* value)
	{
		PyDict_SetItemString(pyDict,name.c_str(),value);
	}
	inline PyObject* NEW(PyObject* obj)
	{
		this->addDeref(obj);
		return obj;
	}
	inline void setString(string name,string value)
	{


		addObj(name,NEW(PyString_FromString(value.c_str())));
	}
	inline void setFloat(string name,double value)
	{
		addObj(name,NEW(PyFloat_FromDouble(value)));
	}
	inline void setInt(string name,long value)
	{
		addObj(name,NEW(PyInt_FromLong(value)));
	}

	inline void addArgItem(string name,string type,string value)
	{
		type=StringUtil::toLower(type);
		if(name=="@externs")
		{
			loadPythonFile(value);
		}else
		{
			if(type=="b" || type=="boolean" || type=="bool")
			{

				value=StringUtil::toLower(value);
				if(value=="1" || value=="true" || value=="yes" || value=="y" || value=="t")
				{
					setInt(name,1);
				}else
				{
					setInt(name,0);
				}
			}
			else if (type=="float" || type=="double" || type=="d" || type=="f")
			{
				setFloat(name,StringUtil::atof(value));

			}
			else if(type=="int" || type=="i" )
			{
				setInt(name,StringUtil::atoi(value));
			}else

			{
								string escaped=escape(value);
								setString(name,value);

			}
		}
	}
	//return args
	inline vector<string> loadFromArgs(int argc,const char** argv)
	{
		general_argvs.clear();
		bool acceptValue=false;
		string key="";
		string type="";
		vector<string> splits;
		//vector<string> general_argvs;
		for(int i=0;i<argc;i++)
		{
			if(!acceptValue)
			{
				if(argv[i][0]=='-')
				{
					acceptValue=true;
					const char*curarg=argv[i];
					key=string(curarg+1);
					type="";
					StringUtil::split(key,":",splits);
					key=splits[0];
					if(splits.size()==2)
					{
						type=splits[1];
					}

				}
				else
				{
					general_argvs.push_back(argv[i]);
				}
			}
			else
			{
				this->addArgItem(key,type,argv[i]);
				acceptValue=false;
			}

		}

		return general_argvs;

	}

	inline void runString(string str)
	{
		PyRun_String(str.c_str(),Py_file_input,pyDict,pyDict);
		string excString=this->getExceptionString();
		if(excString!="")
			throw python_error(excString);
	}
	inline PyObject* getObject(string name,bool addToDeref=true)
	{
		PyObject* obj=PyDict_GetItemString(pyDict,name.c_str());
		if(!obj)
			throw key_not_found(name);

		if(addToDeref)
			this->addDeref(obj);

		return obj;
	}
	inline string getStringDirect(string name)
	{
		return PyString_AS_STRING(getObject(name.c_str()));
	}
	inline void addDeref(PyObject* obj)
	{
		this->derefee.insert(obj);
	}
	inline string getString(string name)
	{
		return getString(getObject(name));
	}
	inline string getString(PyObject* rawObj)
	{



		PyObject* pArgs=NEW(PyTuple_New(1));

		PyTuple_SetItem(pArgs,0,rawObj);

		PyObject* pvalue=NEW(PyObject_CallObject(pyFunc_Str,pArgs));


		return PyString_AS_STRING(pvalue);


	}
	inline ostream& listAllItemsAsSuperString(ostream& os)
	{
		return this->listAllVariables(os);
	}
	inline vector<string> readArray(string name)
	{
		return this->getListOfStrings(name);
	}
	inline string readSuperString(string name)
	{
		return getString(name);
	}
	inline bool readBool(string name)
	{
		return this->getLongDirect(name);
	}

	inline bool readBoolWithDefaultValue(string name,bool defaultvalue)
	{
		try{
			return this->readBool(name);
		}catch(key_not_found & knf)
		{
			return defaultvalue;
		}
	}

	inline string readSuperStringWithDefaultValue(string name,const string& defaultValue)
	{
		try{
			return this->readSuperString(name);
		}catch(key_not_found& knf)
		{
			return defaultValue;
		}
	}

	inline int readInt(string name)
	{
		return this->getLongDirect(name);
	}
	inline string getObjectInFancyString(string name)
	{


		/*
		PyObject* pArgs=NEW(PyTuple_New(1));

		PyTuple_SetItem(pArgs,0,rawObj);

		PyObject* pvalue=NEW(PyObject_CallObject(pyFunc_Print_Obj,pArgs));


		return PyString_AS_STRING(pvalue);
		*/
		runString("___fancy___X_X_X=___Print_Obj("+name+")");
		//runString("print >> stderr,'print direct',___fancy___X_X_X");
		return getStringDirect("___fancy___X_X_X");

	}

	inline long getLongDirect(string name)
	{
		return PyInt_AS_LONG(getObject(name.c_str()));
	}
	inline double getFloatDirect(string name)
	{
		return PyFloat_AS_DOUBLE(getObject(name.c_str()));
	}

	inline double getFloatDirectWithDefaultValue(string name,double defaultvalue)
	{
		try{
			return getFloatDirect(name);
		}catch(key_not_found &knf)
		{
			return defaultvalue;
		}
	}
	inline string getListItemString(PyObject* listObj,int i)
	{
		return PyString_AS_STRING(NEW(PyList_GET_ITEM(listObj,i)));
	}
	inline long getListItemInt(PyObject* listObj,int i)
	{
		return PyInt_AS_LONG(NEW(PyList_GET_ITEM(listObj,i)));
	}
	inline double getListItemDouble(PyObject* listObj,int i)
	{
		return PyFloat_AS_DOUBLE(NEW(PyList_GET_ITEM(listObj,i)));
	}
	inline vector<string> getListOfStrings(PyObject* listObj)
	{

		vector<string> result;
		if(!listObj)
			return result;
		int size=PyList_Size(listObj);
		for(int i=0;i<size;i++)
		{
			result.push_back(getListItemString(listObj,i));
		}
		return result;
	}
	inline void getStringMapOfStringList(string name,map<string, vector<string> >& stringMapVectorString,bool clearMap=true)
	{
		map<string,PyObject*> tmpMap=getStringDictOfPyObjects(name);

		if(clearMap)
		{
			stringMapVectorString.clear();
		}

		for(map<string,PyObject*>::iterator i=tmpMap.begin();i!=tmpMap.end();i++)
		{
			string chr=i->first;
			PyObject* list=i->second;
			stringMapVectorString.insert(map<string,vector<string> >::value_type(chr,getListOfStrings(list)));

		}
	}
	inline vector<string> getListOfStrings(string name)
	{

		PyObject* listObj=getObject(name.c_str());
		return getListOfStrings(listObj);
	}
	inline vector<int> getListOfLongs(string name)
	{
		vector<int> result;
		PyObject* listObj=getObject(name.c_str());
		int size=PyList_Size(listObj);
		for(int i=0;i<size;i++)
		{
			result.push_back(getListItemInt(listObj,i));
		}
		return result;
	}

	inline vector<double> getListOfDoubles(string name)
	{
		vector<double> result;
		PyObject* listObj=getObject(name.c_str());
		int size=PyList_Size(listObj);
		for(int i=0;i<size;i++)
		{
			result.push_back(getListItemDouble(listObj,i));
		}
		return result;
	}

	inline map<string,string> getStringDictOfStrings(PyObject* dictObj)
	{
		map<string,string> result;

		if(!dictObj)
			return result;
		PyObject* items=NEW(PyDict_Items(dictObj));
		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			string key=getString(NEW(PyTuple_GET_ITEM(item,0)));
			string value=getString(NEW(PyTuple_GET_ITEM(item,1)));
			result.insert(map<string,string>::value_type(key,value));
		}

		return result;

	}
	inline int readIntWithDefaultValue(string name,int defaultValue)
	{
		try
		{
			return readInt(name);
		}catch(PyConfigFile::key_not_found & knf)
		{
			return defaultValue;
		}
	}

	inline map<string,string> getStringDictOfStrings(string name)
	{
		return getStringDictOfStrings(getObject(name));

	}

	inline map<string,long> getStringDictOfLongs(string name)
	{
		map<string,long> result;
		PyObject* dictObj=getObject(name.c_str());
		PyObject* items=NEW(PyDict_Items(dictObj));
		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			string key=getString(NEW(PyTuple_GET_ITEM(item,0)));
			long value=PyInt_AS_LONG(NEW(PyTuple_GET_ITEM(item,1)));
			result.insert(map<string,long>::value_type(key,value));
		}

		return result;

	}

	inline map<string,PyObject*> getStringDictOfPyObjects(string name)
	{
		map<string,PyObject*> result;
		PyObject* dictObj=getObject(name.c_str());
		PyObject* items=NEW(PyDict_Items(dictObj));
		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			string key=getString(NEW(PyTuple_GET_ITEM(item,0)));
			PyObject* value=(NEW(PyTuple_GET_ITEM(item,1)));
			result.insert(map<string,PyObject*>::value_type(key,value));
		}

		return result;

	}
	inline map<string,double> getStringDictOfDoubles(string name)
	{
		map<string,double> result;
		PyObject* dictObj=getObject(name.c_str());
		PyObject* items=NEW(PyDict_Items(dictObj));
		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			string key=getString(NEW(PyTuple_GET_ITEM(item,0)));
			double value=PyFloat_AS_DOUBLE(NEW(PyTuple_GET_ITEM(item,1)));
			result.insert(map<string,double>::value_type(key,value));
		}

		return result;

	}


	inline void finalize()
	{
		///cerr<<"fstart"<<endl;
		if(pyDict)
		{

			//cerr<<"A"<<endl;
			//for(set<PyObject*>::iterator i=this->derefee.begin();i!=this->derefee.end();i++)
			//	Py_XDECREF(*i);
			//cerr<<"B"<<endl;
			//Py_Finalize();
			//cerr<<"C"<<endl;
			pyDict=NULL;

		}
		//cerr<<"fend"<<endl;

	}

	inline void listD()
	{


		cerr<<getObjectInFancyString("D")<<endl;
	}

	inline ostream& listAllVariables(ostream& os)
	{
		/*PyObject* items=NEW(PyDict_Items(pyDict));

		int nItems=PyList_Size(items);
		for(int i=0;i<nItems;i++)
		{
			PyObject* item=NEW(PyList_GET_ITEM(items,i)); //item is now a (key,value) tuple
			PyObject* po=NEW(PyTuple_GET_ITEM(item,0));
			string key=PyString_AS_STRING(po);

			if(sysVars.find(key)!=sysVars.end())
				continue;

			//if(key.length()<2)
				//this->derefee.erase(po);
			//cerr<<"getting key="<<key<<endl;
			string value=this->getObjectInFancyString(key);

			os<<key<<"="<<value<<endl;
		}*/
		vector<string> curGlobals=this->getGlobals();
		for(vector<string>::iterator i=curGlobals.begin();i!=curGlobals.end();i++)
		{
			if(sysVars.find(*i)!=sysVars.end())
				continue;
			string key=*i;
			//os<<key<<endl;
			os<<key<<"="<<this->getObjectInFancyString(key)<<endl;
		}



		return os;
	}
	inline ~PyConfigFile()
	{
		finalize();
	}

};

#endif /* PYCONFIGFILE_H_ */
