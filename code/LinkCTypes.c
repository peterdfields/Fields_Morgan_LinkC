#include <stdio.h>
#include <stdlib.h>

#include "LinkCConfig.h"
#include "LinkCTypes.h"

#ifdef LINKDOS_MEMORY_DEBUG
static void* debugMalloc(size_t size);

#define MALLOC debugMalloc
#else
#define MALLOC malloc
#endif

matpop createMatpop(LinkCConfiguration config) {
	return (matpop) MALLOC(sizeof(int) * config.npopm);
}

matal createMatal(LinkCConfiguration config) {
	float*** ret;
	int i, j;

	ret = (float***) MALLOC(sizeof(float**) * config.npopm);

	for (i = 0; i < config.npopm; i++) {
		ret[i] = (float**) MALLOC(sizeof(float*) * config.nlocm);

		for (j = 0; j < config.nlocm; j++)
			ret[i][j] = (float*) MALLOC(sizeof(float) * config.alm);
	}

	return ret;
}

matald createMatald(LinkCConfiguration config) {
	int** ret;
	int i;

	ret = (int**) MALLOC(sizeof(int*) * config.alm);

	for (i = 0; i < config.alm; i++)
		ret[i] = (int*) MALLOC(sizeof(int) * config.nlocm);

	return ret;
}

matdl createMatdl(LinkCConfiguration config) {
	float** ret;
	int i;

	ret = (float**) MALLOC(sizeof(float*) * config.alm);

	for (i = 0; i < config.alm; i++)
		ret[i] = (float*) MALLOC(sizeof(float) * config.alm);

	return ret;
}

matind createMatind(LinkCConfiguration config) {
	int** ret;
	int i;

	ret = (int**) MALLOC(sizeof(int*) * config.npopm);

	for (i = 0; i < config.npopm; i++)
		ret[i] = (int*) MALLOC(sizeof(int) * config.nlocm);

	return ret;
}

matnumal createMatnumal(LinkCConfiguration config) {
	return (int*) MALLOC(sizeof(int) * config.nlocm);
}

matdes createMatdes(LinkCConfiguration config) {
	float** ret;
	int i;

	ret = (float**) MALLOC(sizeof(float*) * config.alm);

	for (i = 0; i < config.alm; i++)
		ret[i] = (float*) MALLOC(sizeof(float) * config.alm);

	return ret;
}

matstat createMatstat(LinkCConfiguration config) {
	float** ret;
	int i;

	ret = (float**) MALLOC(sizeof(float*) * config.npopm);

	for (i = 0; i < config.npopm; i++)
		ret[i] = (float*) MALLOC(sizeof(float) * config.nlocm);

	return ret;
}

matst createMatst(LinkCConfiguration config) {
	float** ret;
	int i;

	ret = (float**) MALLOC(sizeof(float*) * config.nlocm);

	for (i = 0; i < config.nlocm; i++)
		ret[i] = (float*) MALLOC(sizeof(float) * config.alm);

	return ret;
}

void destroyMatpop(LinkCConfiguration config, matpop value) {
	free(value);
}

void destroyMatal(LinkCConfiguration config, matal value) {
	int i, j;

	for (i = 0; i < config.npopm; i++) {
		for (j = 0; j < config.nlocm; j++)
			free(value[i][j]);

		free(value[i]);
	}

	free(value);
}

void destroyMatald(LinkCConfiguration config, matald value) {
	int i;

	for (i = 0; i < config.alm; i++)
		free(value[i]);

	free(value);
}

void destroyMatdl(LinkCConfiguration config, matdl value) {
	int i;

	for (i = 0; i < config.alm; i++)
		free(value[i]);

	free(value);
}

void destroyMatind(LinkCConfiguration config, matind value) {
	int i;

	for (i = 0; i < config.npopm; i++)
		free(value[i]);

	free(value);
}

void destroyMatnumal(LinkCConfiguration config, matnumal value) {
	free(value);
}

void destroyMatdes(LinkCConfiguration config, matdes value) {
	int i;

	for (i = 0; i < config.alm; i++)
		free(value[i]);

	free(value);
}

void destroyMatstat(LinkCConfiguration config, matstat value) {
	int i;

	for (i = 0; i < config.npopm; i++)
		free(value[i]);

	free(value);
}

void destroyMatst(LinkCConfiguration config, matst value) {
	int i;

	for (i = 0; i < config.nlocm; i++)
		free(value[i]);

	free(value);
}

#ifdef LINKDOS_MEMORY_DEBUG
void* debugMalloc(size_t size)
{
	void *ret = malloc(size);
	if (!ret)
	{
		fprintf(stderr, "Unable to allocate memory.\n");
		exit(-1);
	}

	return ret;
}
#endif
