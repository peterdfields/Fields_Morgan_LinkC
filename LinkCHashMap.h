#ifndef __LINKDOS_HASH_MAP_H__
#define __LINKDOS_HASH_MAP_H__

#include <stdio.h>

typedef struct LinkCHashMap {
	void (*put)(struct LinkCHashMap*, const char *key, const char *value);
	const char* (*get)(struct LinkCHashMap*, const char *key);
	size_t (*length)(struct LinkCHashMap*);

	void (*print)(struct LinkCHashMap*, FILE*);

	void (*destroy)(struct LinkCHashMap*);
} LinkCHashMap;

LinkCHashMap* createLinkCHashMap(unsigned int initialCapacity);

#endif
