
#include <assert.h>

//#include <glib-2.0/glib.h>
#include <mono/jit/jit.h>
#include <mono/metadata/assembly.h>
#include <mono/metadata/debug-helpers.h>

#define DECLARE_MONKEY
#define DEFINE_MONKEY_INTERNALS
#include "monkey.h"

//MonoMethod *Managed_Infrastructure_Addr;
int initdone = 0;

void ilPSPinit() {
    MonoClass *klass;
    MonoMethodDesc* desc_X;

    assert(initdone == 0);
    initdone = 1;

    domain = mono_jit_init_version("domain_name","v2.0.50727"); 
    assert(domain);

    assembly = mono_domain_assembly_open (domain, "ilPSP.ExternalBinding.dll");
    assert(assembly);

    img = mono_assembly_get_image(assembly);
    assert(img);

    klass = mono_class_from_name (img, "ilPSP.ExternalBinding", "Infrastructure");
    assert(klass);

    desc_X = mono_method_desc_new ("ilPSP.ExternalBinding.Infrastructure:Addr(string)", 1);
    assert(desc_X);

    //Managed_Infrastructure_Addr = mono_method_desc_search_in_class (desc_X, klass);
    //assert(Managed_Infrastructure_Addr);

    InitBinding();

    Common_ilPSPInitialize();
}
/* all kinds of name mangling */
void _ilPSPinit() { ilPSPinit(); }
void ilPSPinit_() { ilPSPinit(); }
void _ilpspinit() { ilPSPinit(); }
void ilpspinit_() { ilPSPinit(); }
void _ILPSPINIT() { ilPSPinit(); }
void ILPSPINIT_() { ilPSPinit(); }
void  ILPSPINIT() { ilPSPinit(); }

/*
void* Infrastructure_Addr(char* name) {
    void* args[1];
    MonoObject *retval;
    void* funcPointer;

    args [0] = mono_string_new (domain, name);
    retval = mono_runtime_invoke (Managed_Infrastructure_Addr, NULL, args, NULL);

    funcPointer = *((void**)mono_object_unbox (retval));
    return funcPointer;
}
*/

void ilPSPshutdown() {
    assert(initdone);
    Common_ilPSPFinalize();
    mono_jit_cleanup (domain);
}

void _ilPSPshutdown() { ilPSPshutdown(); }
void ilPSPshutdown_() { ilPSPshutdown(); }
void _ilpspshutdown() { ilPSPshutdown(); }
void ilpspshutdown_() { ilPSPshutdown(); }
void _ILPSPSHUTDOWN() { ilPSPshutdown(); }
void ILPSPSHUTDOWN_() { ilPSPshutdown(); }

