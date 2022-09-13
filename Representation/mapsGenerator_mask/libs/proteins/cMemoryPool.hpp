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


#pragma once

#include <list>
//#include <stdio.h>
#include <cstdlib>

class cMemoryChunk;

class cMemoryPool {
public:
	cMemoryPool(size_t nbytes);
	~cMemoryPool();
	void* alloc(size_t nbytes);
	void dealloc(void* p);
	
private:
	cMemoryChunk* firstChunk;
	cMemoryChunk* currentChunk;
	
};

class  cMemoryChunk {

	friend class cMemoryPool;
	public:

		cMemoryChunk(size_t nbytes);
		~cMemoryChunk();
		
		cMemoryChunk* alloc(size_t nbytes, void *&);
		//cMemoryChunk *getNextChunk() const {return nextChunk;}
		//void* operator new(size_t nbytes);
		//void* operator new(size_t nbytes, Pool& pool);
		//void operator delete(void* p);

	private:

		size_t length;
		size_t used;
		void *data;
		cMemoryChunk *nextChunk;
		

};

inline void* operator new(size_t nbytes, cMemoryPool* pool)
{
	return pool->alloc(nbytes);
} 

inline void* operator new[](size_t nbytes, cMemoryPool* pool) {
	return pool->alloc(nbytes); 
}

inline void* operator new(size_t nbytes, cMemoryPool &pool)
{
	return pool.alloc(nbytes);
} 

inline void* operator new[](size_t nbytes, cMemoryPool& pool) {
	return pool.alloc(nbytes); 
}

inline void operator delete(void* p, cMemoryPool& pool) // to handle an exception
{
	pool.dealloc(p);
} 		

inline void operator delete[](void* p, cMemoryPool& pool) // to handle an exception
{
	pool.dealloc(p);
} 		

inline void operator delete(void* p, cMemoryPool *pool) // to handle an exception
{
	pool->dealloc(p);
} 		

inline void operator delete[](void* p, cMemoryPool *pool) // to handle an exception
{
	pool->dealloc(p);
} 		

