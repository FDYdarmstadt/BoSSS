#define PARAM1    void* _1
#define CALL1     _1
#define PARAM2    void* _1, void* _2
#define CALL2     _1, _2
#define PARAM3    void* _1, void* _2, void* _3
#define CALL3     _1, _2, _3
#define PARAM4    void* _1, void* _2, void* _3, void* _4
#define CALL4     _1, _2, _3, _4
#define PARAM5    void* _1, void* _2, void* _3, void* _4, void* _5
#define CALL5     _1, _2, _3, _4, _5
#define PARAM6    void* _1, void* _2, void* _3, void* _4, void* _5, void* _6
#define CALL6     _1, _2, _3, _4, _5, _6
#define PARAM7    void* _1, void* _2, void* _3, void* _4, void* _5, void* _6, void* _7
#define CALL7     _1, _2, _3, _4, _5, _6, _7
#define PARAM8    void* _1, void* _2, void* _3, void* _4, void* _5, void* _6, void* _7, void* _8
#define CALL8     _1, _2, _3, _4, _5, _6, _7, _8


#define MAKE_MPIF(funcname,param,call)         \
void MPI_API funcname(param);                     \
DLL_EXPORT void BoSSS_##funcname(param) {      \
    funcname(call);                  \
}

//#ifdef WIN32
//#define DLL_EXPORT __declspec(dllexport)
//#define MPI_API
//#else
#define DLL_EXPORT 
#define MPI_API
//#endif