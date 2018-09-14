#include <mpi.h>
#include "DllExportPreProc.h"

/*
 * Within this file, we define the MPI functions that should be exported by the Platform_Native - library;
 * We add a BoSSS-prefix to avoid naming confusion.
 * (On Windows, it would be easier to us an .def - export file, but that is not supported by fu**ing Linux,
 * so I have to write this extremely silly piece of code.)
 */

int DLL_EXPORT BoSSS_Get_MPI_COMM_WORLD()   { return MPI_Comm_c2f(MPI_COMM_WORLD); }
int DLL_EXPORT BoSSS_Get_MPI_COMM_SELF()    { return MPI_Comm_c2f(MPI_COMM_SELF); }
    
int DLL_EXPORT BoSSS_Get_MPI_Datatype_CHAR()           { return MPI_Type_c2f(MPI_CHAR); }
//int BoSSS_Get_MPI_Datatype_SIGNED_CHAR() { return MPI_SIGNED_CHAR; } // only MS-MPI
int DLL_EXPORT BoSSS_Get_MPI_Datatype_UNSIGNED_CHAR()  { return MPI_Type_c2f(MPI_UNSIGNED_CHAR); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_BYTE()           { return MPI_Type_c2f(MPI_BYTE); }
//int BoSSS_Get_MPI_Datatype_WCHAR() { return MPI_WCHAR; } // only MS-MPI
int DLL_EXPORT BoSSS_Get_MPI_Datatype_SHORT()          { return MPI_Type_c2f(MPI_SHORT); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_UNSIGNED_SHORT() { return MPI_Type_c2f(MPI_UNSIGNED_SHORT); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_INT()            { return MPI_Type_c2f(MPI_INT); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_UNSIGNED()       { return MPI_Type_c2f(MPI_UNSIGNED); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_LONG()           { return MPI_Type_c2f(MPI_LONG); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_UNSIGNED_LONG()  { return MPI_Type_c2f(MPI_UNSIGNED_LONG); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_FLOAT()          { return MPI_Type_c2f(MPI_FLOAT); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_DOUBLE()         { return MPI_Type_c2f(MPI_DOUBLE); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_LONG_DOUBLE()    { return MPI_Type_c2f(MPI_LONG_DOUBLE); }
int DLL_EXPORT BoSSS_Get_MPI_Datatype_LONG_LONG_INT()  { return MPI_Type_c2f(MPI_LONG_LONG_INT); }
//int BoSSS_Get_MPI_Datatype_UNSIGNED_LONG_LONG() { return MPI_UNSIGNED_LONG_LONG; } // only MS-MPI

int DLL_EXPORT BoSSS_Get_MPI_MAX()    { return MPI_Op_c2f(MPI_MAX); }
int DLL_EXPORT BoSSS_Get_MPI_MIN()    { return MPI_Op_c2f(MPI_MIN); }
int DLL_EXPORT BoSSS_Get_MPI_SUM()    { return MPI_Op_c2f(MPI_SUM); }
int DLL_EXPORT BoSSS_Get_MPI_PROD()   { return MPI_Op_c2f(MPI_PROD); }
int DLL_EXPORT BoSSS_Get_MPI_LAND()   { return MPI_Op_c2f(MPI_LAND); }
int DLL_EXPORT BoSSS_Get_MPI_BAND()   { return MPI_Op_c2f(MPI_BAND); }
int DLL_EXPORT BoSSS_Get_MPI_LOR()    { return MPI_Op_c2f(MPI_LOR); }
int DLL_EXPORT BoSSS_Get_MPI_BOR()    { return MPI_Op_c2f(MPI_BOR); }
int DLL_EXPORT BoSSS_Get_MPI_LXOR()   { return MPI_Op_c2f(MPI_LXOR); }
int DLL_EXPORT BoSSS_Get_MPI_BXOR()   { return MPI_Op_c2f(MPI_BXOR); }
int DLL_EXPORT BoSSS_Get_MPI_MINLOC() { return MPI_Op_c2f(MPI_MINLOC); }
int DLL_EXPORT BoSSS_Get_MPI_MAXLOC() { return MPI_Op_c2f(MPI_MAXLOC); }
//int BoSSS_Get_MPI_REPLACE() { return MPI_REPLACE; }



int DLL_EXPORT BoSSS_Get_MPI_ANY_SOURCE() {  return MPI_ANY_SOURCE; }       
int DLL_EXPORT BoSSS_Get_MPI_PROC_NULL() { return MPI_PROC_NULL; }  
int DLL_EXPORT BoSSS_Get_MPI_ROOT() { return MPI_ROOT; } 
int DLL_EXPORT BoSSS_Get_MPI_ANY_TAG() { return MPI_ANY_TAG; } 
int DLL_EXPORT BoSSS_Get_MPI_MAX_PROCESSOR_NAME() { return MPI_MAX_PROCESSOR_NAME; }
int DLL_EXPORT BoSSS_Get_MPI_MAX_ERROR_STRING() { return MPI_MAX_ERROR_STRING; }  
int DLL_EXPORT BoSSS_Get_MPI_MAX_OBJECT_NAME() {  return MPI_MAX_OBJECT_NAME; } 
int DLL_EXPORT BoSSS_Get_MPI_UNDEFINED() {  return MPI_UNDEFINED; } 
int DLL_EXPORT BoSSS_Get_MPI_CART() {  return MPI_CART; } 
int DLL_EXPORT BoSSS_Get_MPI_GRAPH() {  return MPI_GRAPH; } 
int DLL_EXPORT BoSSS_Get_MPI_KEYVAL_INVALID() {  return MPI_KEYVAL_INVALID; }  
int DLL_EXPORT BoSSS_Get_MPI_REQUEST_NULL() { return MPI_Request_c2f(MPI_REQUEST_NULL); }

int DLL_EXPORT BoSSS_Get_MPI_Status_Size() { return sizeof(MPI_Status); }

int DLL_EXPORT  BoSSS_MPI_Init(int* argc, char*** argv) {
    return MPI_Init(argc,argv);
}

int DLL_EXPORT  BoSSS_MPI_Finalize() {
    return MPI_Finalize();
}

int DLL_EXPORT BoSSS_MPI_Status_C2f(MPI_Status* cSt, int fSt[]) {
	MPI_Status_c2f(cSt, fSt);
}

int DLL_EXPORT BoSSS_MPI_Status_f2c(int fSt[], MPI_Status* cSt) {
	MPI_Status_f2c(fSt, cSt);
}

int DLL_EXPORT BoSSS_MPI_Request_C2f(MPI_Request cRq, int* fRq) {
	*fRq = MPI_Request_f2c(cRq);
	return 0;
}

int DLL_EXPORT BoSSS_MPI_Request_f2c(int* fRq, MPI_Request* cRq) {
	*cRq = MPI_Request_c2f(fRq);
	return 0;
}
