// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME aDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "bitset256.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *vectorlEbitsetlE256gRsPgR_Dictionary();
   static void vectorlEbitsetlE256gRsPgR_TClassManip(TClass*);
   static void *new_vectorlEbitsetlE256gRsPgR(void *p = 0);
   static void *newArray_vectorlEbitsetlE256gRsPgR(Long_t size, void *p);
   static void delete_vectorlEbitsetlE256gRsPgR(void *p);
   static void deleteArray_vectorlEbitsetlE256gRsPgR(void *p);
   static void destruct_vectorlEbitsetlE256gRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bitset<256> >*)
   {
      vector<bitset<256> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bitset<256> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bitset<256> >", -2, "vector", 214,
                  typeid(vector<bitset<256> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEbitsetlE256gRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<bitset<256> >) );
      instance.SetNew(&new_vectorlEbitsetlE256gRsPgR);
      instance.SetNewArray(&newArray_vectorlEbitsetlE256gRsPgR);
      instance.SetDelete(&delete_vectorlEbitsetlE256gRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEbitsetlE256gRsPgR);
      instance.SetDestructor(&destruct_vectorlEbitsetlE256gRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bitset<256> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<bitset<256> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEbitsetlE256gRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bitset<256> >*)0x0)->GetClass();
      vectorlEbitsetlE256gRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEbitsetlE256gRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEbitsetlE256gRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bitset<256> > : new vector<bitset<256> >;
   }
   static void *newArray_vectorlEbitsetlE256gRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bitset<256> >[nElements] : new vector<bitset<256> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEbitsetlE256gRsPgR(void *p) {
      delete ((vector<bitset<256> >*)p);
   }
   static void deleteArray_vectorlEbitsetlE256gRsPgR(void *p) {
      delete [] ((vector<bitset<256> >*)p);
   }
   static void destruct_vectorlEbitsetlE256gRsPgR(void *p) {
      typedef vector<bitset<256> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<bitset<256> >

namespace ROOT {
   static TClass *bitsetlE256gR_Dictionary();
   static void bitsetlE256gR_TClassManip(TClass*);
   static void *new_bitsetlE256gR(void *p = 0);
   static void *newArray_bitsetlE256gR(Long_t size, void *p);
   static void delete_bitsetlE256gR(void *p);
   static void deleteArray_bitsetlE256gR(void *p);
   static void destruct_bitsetlE256gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const bitset<256>*)
   {
      bitset<256> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(bitset<256>));
      static ::ROOT::TGenericClassInfo 
         instance("bitset<256>", 2, "bitset", 749,
                  typeid(bitset<256>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &bitsetlE256gR_Dictionary, isa_proxy, 4,
                  sizeof(bitset<256>) );
      instance.SetNew(&new_bitsetlE256gR);
      instance.SetNewArray(&newArray_bitsetlE256gR);
      instance.SetDelete(&delete_bitsetlE256gR);
      instance.SetDeleteArray(&deleteArray_bitsetlE256gR);
      instance.SetDestructor(&destruct_bitsetlE256gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback<Internal::TStdBitsetHelper< bitset<256> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const bitset<256>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *bitsetlE256gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const bitset<256>*)0x0)->GetClass();
      bitsetlE256gR_TClassManip(theClass);
   return theClass;
   }

   static void bitsetlE256gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_bitsetlE256gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<256> : new bitset<256>;
   }
   static void *newArray_bitsetlE256gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) bitset<256>[nElements] : new bitset<256>[nElements];
   }
   // Wrapper around operator delete
   static void delete_bitsetlE256gR(void *p) {
      delete ((bitset<256>*)p);
   }
   static void deleteArray_bitsetlE256gR(void *p) {
      delete [] ((bitset<256>*)p);
   }
   static void destruct_bitsetlE256gR(void *p) {
      typedef bitset<256> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class bitset<256>

namespace {
  void TriggerDictionaryInitialization_aDict_Impl() {
    static const char* headers[] = {
"bitset256.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include",
"/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/BaconProd/Ntupler/plugins/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "aDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <size_t _Nb> class __attribute__((annotate("$clingAutoload$bitset256.h")))  bitset;
}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "aDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "bitset256.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("aDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_aDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_aDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_aDict() {
  TriggerDictionaryInitialization_aDict_Impl();
}
