
#ifndef MONKEY_H_INCLUDED_
#define MONKEY_H_INCLUDED_

void ilPSPinit();
void ilPSPshutdown();
void* Infrastructure_Addr(char* name);
void InitBinding();

#ifdef DEFINE_MONKEY_INTERNALS
#ifndef DECLARE_MONKEY
#define BINDING_
#else
#define BINDING_ extern
#endif 

BINDING_ MonoAssembly *assembly;
BINDING_ MonoDomain *domain;
BINDING_ MonoImage *img;

#endif

#include "binding.h"

#endif


