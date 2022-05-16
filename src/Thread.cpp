#include "Thread.h"

Thread::Thread()
{
	tid = 0;
	user_tid = 0;
	blockVector = NULL;
}

Thread::~Thread(){
}

void* Thread::runDetect0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->runDetect1();
	return p;
}

void* Thread::runDetect1()
{
	runDetect();
	pthread_exit(NULL);
}

bool Thread::startDetect()
{
	return pthread_create(&tid, NULL, runDetect0, this) == 0;
}
pthread_t Thread::getThreadID()
{
	return tid;
}

void Thread::setUserThreadID(size_t user_tid){
	this->user_tid = user_tid;
}

void Thread::setBlockVec(vector<Block*> *blockVector){
	this->blockVector = blockVector;
}

size_t Thread::getUserThreadID(){
	return user_tid;
}

bool Thread::join()
{
	return pthread_join(tid, NULL) == 0;
}

void MultiThread::setNumThreads(size_t n){
	num_threads = n;
}

void MultiThread::runDetect(){
	size_t user_tid_tmp = getUserThreadID();

	Block *bloc;
	for (size_t i=0; i < blockVector->size(); i++){
		if(i%num_threads==user_tid_tmp){

			bloc = blockVector->at(i);
			if(bloc->process_flag) {
				bloc->blockIlluminaDetect();
			}
		}
	}
}
