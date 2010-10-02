#ifndef NUCLEICBIT_H_
#define NUCLEICBIT_H_

#include "Nucleic.h"
#include <inttypes.h>

template < typename T >
inline T highbit(T& t)
{
   return t = (((T)(-1)) >> 1) + 1;
}

template < typename T >
std::ostream& bin(const T& value, std::ostream &o)
{
   for ( T bit = highbit(bit); bit; bit >>= 1 )
   {
	  o << ( ( value & bit ) ? '1' : '0' );
   }
   return o;
}



typedef  unsigned long Nuc16Key;
typedef  uint64_t NucKey3b;
#define NULL_KEY 0x00;
#define KEY_CAPACITY 21
#define SHIFT(x,y) (x<<(y*3))
#define ENCROUND(y) { \
		if(ps>=pe) \
		{      \
			if(y<21)\
				key|=SHIFT(uint64_t(7),y); \
			return key; \
		} \
		switch(*ps) \
		{ \
			case A_SYM:case A_SYML: \
				break; \
			case C_SYM:case C_SYML: \
				key|=SHIFT(uint64_t(1),y); \
				break; \
			case T_SYM:case T_SYML:case U_SYM:case U_SYML: \
				key|=SHIFT(uint64_t(2),y); \
				break; \
			case G_SYM:case G_SYML: \
				key|=SHIFT(uint64_t(3),y); \
				break; \
			case R_SYM:case R_SYML: \
				key|=SHIFT(uint64_t(5),y); \
				break; \
			case Y_SYM:case Y_SYML: \
				key|=SHIFT(uint64_t(6),y); \
				break;\
			default:\
				key|=SHIFT(uint64_t(4),y); \
				break;\
		} \
		ps++;\
} 

inline NucKey3b Nuc2Key3b(const char* nuc,int length)
{
	NucKey3b key=NULL_KEY;
	register const char* ps=nuc;
	register const char* pe=nuc+length;
	//register unsigned char* kw=((unsigned char*)((void*)&key));
	
	if(length>KEY_CAPACITY)
	{
		cerr<<"error: try to encode string of > "<<KEY_CAPACITY<<" in length"<<endl;
		return key;
	}
	
	
	ENCROUND(0);
	ENCROUND(1);	
	ENCROUND(2);	
	ENCROUND(3);	
	ENCROUND(4);	
	ENCROUND(5);
	ENCROUND(6);	
	ENCROUND(7);	
	ENCROUND(8);	
	ENCROUND(9);	
	ENCROUND(10);
	ENCROUND(11);	
	ENCROUND(12);	
	ENCROUND(13);	
	ENCROUND(14);	
	ENCROUND(15);
	ENCROUND(16);	
	ENCROUND(17);	
	ENCROUND(18);	
	ENCROUND(19);		
	ENCROUND(20);
	ENCROUND(21);	
	
	return key;
}

inline Nuc16Key Nuc2Key(const char *nuc) //
{

	Nuc16Key key=0x00;
	register const char* ps=nuc;
	register const char* pe=nuc+16;
	register unsigned char* kw=((unsigned char*)((void*)&key));

	
	while(ps<pe)
	{
		switch(*ps)
		{
		case A_SYM:case A_SYML:
		
		break;
		case C_SYM:case C_SYML:
		*kw|=64;
		break;
		case T_SYM:case T_SYML:case U_SYM:case U_SYML:
		*kw|=128;
		break;
		case G_SYM:case G_SYML:
		*kw|=192;
		break;
		
		}
		
		ps++;
		
		switch(*ps)
		{
		case A_SYM:case A_SYML:
		
		break;
		case C_SYM:case C_SYML:
		*kw|=16;
		break;
		case T_SYM:case T_SYML:case U_SYM:case U_SYML:
		*kw|=32;
		break;
		case G_SYM:case G_SYML:
		*kw|=48;
		break;
		
		}
		
		ps++;
		
		switch(*ps)
		{
		case A_SYM:case A_SYML:
		
		break;
		case C_SYM:case C_SYML:
		*kw|=4;
		break;
		case T_SYM:case T_SYML:case U_SYM:case U_SYML:
		*kw|=8;
		break;
		case G_SYM:case G_SYML:
		*kw|=12;
		break;
		
		}
		
		
		ps++;
		
		switch(*ps)
		{
		case A_SYM:case A_SYML:
		
		break;
		case C_SYM:case C_SYML:
		*kw|=1;
		break;
		case T_SYM:case T_SYML:case U_SYM:case U_SYML:
		*kw|=2;
		break;
		case G_SYM:case G_SYML:
		*kw|=3;
		break;
		
		}
		
		
		ps++;
		kw++;
		
		
	}
	
	
	
	
	return key;
	
}

#endif /*NUCLEICBIT_H_*/
