#ifndef __LINKDOS_TYPES_H__
#define __LINKDOS_TYPES_H__

#include "LinkCConfig.h"

typedef int* matpop;
typedef float*** matal;
typedef int** matald;
typedef float** matdl;
typedef int** matind;
typedef int* matnumal;
typedef float** matdes;
typedef float** matstat;
typedef float** matst;

typedef matald* ptald;
typedef matdl* ptdl;
typedef matdes* ptdes;
typedef matstat* ptstat;
typedef matst* ptst;

matpop createMatpop(LinkCConfiguration config);
matal createMatal(LinkCConfiguration config);
matald createMatald(LinkCConfiguration config);
matdl createMatdl(LinkCConfiguration config);
matind createMatind(LinkCConfiguration config);
matnumal createMatnumal(LinkCConfiguration config);
matdes createMatdes(LinkCConfiguration config);
matstat createMatstat(LinkCConfiguration config);
matst createMatst(LinkCConfiguration config);

void destroyMatpop(LinkCConfiguration config, matpop value);
void destroyMatal(LinkCConfiguration config, matal value);
void destroyMatald(LinkCConfiguration config, matald value);
void destroyMatdl(LinkCConfiguration config, matdl value);
void destroyMatind(LinkCConfiguration config, matind value);
void destroyMatnumal(LinkCConfiguration config, matnumal value);
void destroyMatdes(LinkCConfiguration config, matdes value);
void destroyMatstat(LinkCConfiguration config, matstat value);
void destroyMatst(LinkCConfiguration config, matst value);

#endif
