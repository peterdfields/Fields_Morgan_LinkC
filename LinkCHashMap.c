#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LinkCUtilities.h"
#include "LinkCHashMap.h"

typedef struct LinkCHashMapNode {
	char *_key;
	char *_data;
	struct LinkCHashMapNode *_next;
} LinkCHashMapNode;

static LinkCHashMapNode* createNode(const char *key, const char *value,
		LinkCHashMapNode *next);
static void destroyNode(LinkCHashMapNode *node);

LinkCHashMapNode* createNode(const char *key, const char *value,
		LinkCHashMapNode *next) {
	LinkCHashMapNode *ret = (LinkCHashMapNode*) malloc(
			sizeof(LinkCHashMapNode));

	ret->_key = copyString(key);
	ret->_data = copyString(value);
	ret->_next = next;

	return ret;
}

void destroyNode(LinkCHashMapNode *node) {
	if (node->_next)
		destroyNode(node->_next);

	free(node->_key);
	free(node->_data);
}

typedef struct LinkCHashMapImpl {
	LinkCHashMap _interface;

	size_t _capacity;
	size_t _length;
	LinkCHashMapNode **_table;
} LinkCHashMapImpl;

static size_t hashString(const char *string);
static void growTable(LinkCHashMapImpl*);
static void putImpl(struct LinkCHashMap*, const char *key, const char *value);
static const char* getImpl(struct LinkCHashMap*, const char *key);
static size_t lengthImpl(struct LinkCHashMap*);
static void printImpl(struct LinkCHashMap*, FILE*);
static void destroyImpl(struct LinkCHashMap*);

LinkCHashMap* createLinkCHashMap(unsigned int initialCapacity) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) malloc(
			sizeof(LinkCHashMapImpl));

	impl->_interface.put = putImpl;
	impl->_interface.get = getImpl;
	impl->_interface.length = lengthImpl;
	impl->_interface.print = printImpl;
	impl->_interface.destroy = destroyImpl;

	impl->_capacity = initialCapacity;
	impl->_length = 0;
	impl->_table = (LinkCHashMapNode**) malloc(
			sizeof(LinkCHashMapNode*) * initialCapacity);
	memset(impl->_table, 0, sizeof(LinkCHashMapNode*) * initialCapacity);

	return (LinkCHashMap*) impl;
}

size_t hashString(const char *string) {
	size_t lcv;
	size_t ret = 0x0;

	for (lcv = 0; string[lcv]; lcv++) {
		ret <<= 1;
		ret ^= (size_t) (string[lcv]);
	}

	return ret;
}

void growTable(LinkCHashMapImpl *impl) {
	unsigned int lcv;
	LinkCHashMapNode *node;
	LinkCHashMapNode **tmpTable;
	LinkCHashMapImpl *tmpImpl = (LinkCHashMapImpl*) createLinkCHashMap(
			impl->_capacity * 2);

	for (lcv = 0; lcv < impl->_capacity; lcv++) {
		for (node = impl->_table[lcv]; node; node = node->_next) {
			putImpl((LinkCHashMap*) tmpImpl, node->_key, node->_data);
		}
	}

	tmpTable = tmpImpl->_table;
	tmpImpl->_table = impl->_table;
	impl->_table = tmpTable;
	tmpImpl->_capacity = impl->_capacity;
	impl->_capacity *= 2;

	destroyImpl((LinkCHashMap*) tmpImpl);
}

void putImpl(struct LinkCHashMap *map, const char *key, const char *value) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) map;
	LinkCHashMapNode *node;
	size_t index;

	if (impl->_length >= (impl->_capacity * 2 / 3))
		growTable(impl);

	index = hashString(key) % impl->_capacity;
	for (node = impl->_table[index]; node; node = node->_next) {
		if (strcmp(node->_key, key) == 0) {
			free(node->_data);
			node->_data = copyString(value);
			return;
		}
	}

	node = createNode(key, value, impl->_table[index]);
	impl->_table[index] = node;
	impl->_length++;
}

const char* getImpl(struct LinkCHashMap *map, const char *key) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) map;
	LinkCHashMapNode *node;
	size_t index;

	index = hashString(key) % impl->_capacity;
	for (node = impl->_table[index]; node; node = node->_next) {
		if (strcmp(node->_key, key) == 0)
			return node->_data;
	}

	return NULL ;
}

size_t lengthImpl(struct LinkCHashMap *map) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) map;
	return impl->_length;
}

void printImpl(struct LinkCHashMap *map, FILE *out) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) map;
	int lcv;
	LinkCHashMapNode *node;

	fprintf(out, "LinkCHashMap(length = %lu, capacity = %lu):\n", impl->_length,
			impl->_capacity);
	for (lcv = 0; lcv < impl->_capacity; lcv++) {
		for (node = impl->_table[lcv]; node; node = node->_next)
			fprintf(out, "\t%s -> %s\n", node->_key, node->_data);
	}
}

void destroyImpl(struct LinkCHashMap *map) {
	LinkCHashMapImpl *impl = (LinkCHashMapImpl*) map;
	int lcv;

	for (lcv = 0; lcv < impl->_capacity; lcv++) {
		if (impl->_table[lcv])
			destroyNode(impl->_table[lcv]);
	}

	free(impl->_table);
	free(impl);
}
