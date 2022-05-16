#ifndef _THREAD_H
#define _THREAD_H

#include <iostream>
#include <string>
#include <vector>

#include <pthread.h>
#include <unistd.h>
#include "Block.h"
#include "varCand.h"
#include "blatAlnTra.h"

using namespace std;

class Thread
{
	private:
		pthread_t tid;
		size_t user_tid;   // 0-based: 0, 1, 2, 3, ...
		static void* runDetect0(void* pVoid);  // the pointer to executing function
		void* runDetect1();  // inner executing method

	public:
//		size_t user_tid;   // 0-based: 0, 1, 2, 3, ...

		vector<Block*> *blockVector;
		Thread();
		virtual ~Thread();
		virtual void runDetect() = 0;  // thread running entity
		bool startDetect();  // start the thread for detection
		pthread_t getThreadID();
		void setUserThreadID(size_t user_tid); // set user tid
		void setBlockVec(vector<Block*> *vec);
		size_t getUserThreadID();
		bool join();
};

class MultiThread: public Thread {
	public:
		size_t num_threads = 0;
		void runDetect();
		void setNumThreads(size_t n);
};


#endif /* _THREAD_H */
