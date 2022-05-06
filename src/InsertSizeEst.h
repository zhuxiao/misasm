#ifndef SRC_INSERTSIZEEST_H_
#define SRC_INSERTSIZEEST_H_

#include <iostream>
#include <string>
#include <cmath>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "Paras.h"
#include "util_header.h"
#include "alnDataLoader.h"

using namespace std;

#define MAX_PAIR_NUM		200000
#define MAX_ISIZE			2000
#define MIN_CHR_LEN			1000

class InsertSizeEst {

	public:
		string inBamFile;
		double isize, sdev, sum, x2;
		int64_t count;
	public:
		InsertSizeEst(string &inBamFile);
		virtual ~InsertSizeEst();
		void estInsertSize(double &isize_1, double &sdev_1, double &max, double &min);

	private:
		void estInsertSizeFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string &reg);
};

#endif /* SRC_INSERTSIZEEST_H_ */
