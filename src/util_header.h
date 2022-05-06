#ifndef SRC_UTIL_HEADER_H_
#define SRC_UTIL_HEADER_H_

#include <iostream>
#include <string>
#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace std;

const char* getTargetName(const bam_hdr_t *h, int tid);
int64_t getTargetLen(const bam_hdr_t *h, int tid);

#endif /* SRC_UTIL_HEADER_H_ */
