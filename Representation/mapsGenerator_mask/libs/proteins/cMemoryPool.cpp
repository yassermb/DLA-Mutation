/*************************************************************************\
 
 Copyright 2010 Sergei Grudinin, INRIA.
 All Rights Reserved.
 
 The author may be contacted via:
 
 Mail:					Sergei Grudinin
 INRIA Rhone-Alpes Research Unit
 Zirst - 655 avenue de l'Europe - Montbonnot
 38334 Saint Ismier Cedex - France
 
 Phone:				+33 4 76 61 53 24
 
 EMail:				sergei.grudinin@inria.fr
 
 \**************************************************************************/


#include "cMemoryPool.hpp"
#include <iostream>

cMemoryPool::cMemoryPool(size_t nbytes) {
	
	currentChunk = firstChunk = new cMemoryChunk(nbytes);	
	
}

cMemoryChunk::cMemoryChunk(size_t nbytes): length(nbytes), used(0), nextChunk(NULL) {

	data = malloc(nbytes);	
}

cMemoryPool::~cMemoryPool() {

	cMemoryChunk *chunk = firstChunk;
	cMemoryChunk *next;
	while (chunk) {
		next = chunk->nextChunk;
		delete(chunk);
		chunk=next;
	}
}

cMemoryChunk::~cMemoryChunk() {
	
	free(data);
}

void* cMemoryPool::alloc(size_t nbytes) {
	
	void *p;
	currentChunk = currentChunk->alloc(nbytes, p);
	return p;
}

cMemoryChunk* cMemoryChunk::alloc(size_t nbytes, void *&p) {

	if (nbytes>length) {
#if DEBUG
		std::cerr<< "memory requested higher than chunk size of our pool!!!"<<std::endl;
#endif
		if (!used) {
			length = nbytes;
			free(data);
			data = malloc(nbytes);	
			p = (char*) data;
			used += nbytes;
			return this;
		}
		cMemoryChunk *newChunk = new cMemoryChunk(nbytes);
		this->nextChunk = newChunk;
		p = (char*)newChunk->data;
		newChunk->used += nbytes;
		return newChunk;
	}
	
	if (nbytes > length-used) {
		cMemoryChunk *newChunk = new cMemoryChunk(length);
		this->nextChunk = newChunk;
		p = (char*)newChunk->data;
		newChunk->used += nbytes;
		return newChunk;
	}

	p = (char*)data + used;
	used += nbytes;

	return this;
}

void cMemoryPool::dealloc(void* p) {
}

/*
void* cMemoryPool::operator new(size_t nbytes)
{
	if (nbytes == 0)
		nbytes = 1;                    // so all alloc's get a distinct address
	void* ans = malloc(nbytes + 4);  // overallocate by 4 bytes
	*(Pool**)ans = NULL;             // use NULL in the global new
	return (char*)ans + 4;           // don't let users see the Pool*
}

void* cMemoryPool::operator new(size_t nbytes, Pool& pool)
{
	if (nbytes == 0)
		nbytes = 1;                    // so all alloc's get a distinct address
	void* ans = pool.alloc(nbytes + 4); // overallocate by 4 bytes
	*(Pool**)ans = &pool;            // put the Pool* here
	return (char*)ans + 4;           // don't let users see the Pool*
}

void cMemoryPool::operator delete(void* p)
{
	if (p != NULL) {
		p = (char*)p - 4;              // back off to the Pool*
		Pool* pool = *(Pool**)p;
		if (pool == null)
			free(p);                     // note: 4 bytes left of the original p
		else
			pool->dealloc(p);            // note: 4 bytes left of the original p
	}
} 
*/