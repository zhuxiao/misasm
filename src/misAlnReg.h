#ifndef SRC_MISALNREG_H_
#define SRC_MISALNREG_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "Base.h"

using namespace std;

#define REF_END_SKIP_SIZE					1000

#define SUB_MIS_ALN_REG_SIZE				10
#define SUB_MIS_ALN_REG_RATIO_THRES			(0.05f)
#define HIGH_SUB_MIS_ALN_REG_RATIO_THRES	(0.3f)
#define GAPPED_MIS_ALN_REG_NUM_THRES		1
#define MIN_MIS_ALN_REG_NUM_THRES			2
#define SUB_REG_MAX_MIS_ALN_NUM_THRES		2

#define HIGH_CLIP_RATIO_THRES				(0.1f)
#define CLIP_SUPPORT_READS_NUM_THRES		3  // should be input parameter

class misAlnReg {
	public:
		int64_t startPos, endPos, chrlen;
		Base *misAlnRegBaseArr;
		uint16_t disagrNum;
		uint16_t misAlnSubregNum, subRegNum, highClipBaseNum, zeroCovBaseNum;		// 10-bp sub-region
		float disagrRegRatio;
		bool misAlnFlag;

	public:
		misAlnReg(size_t startPos, size_t endPos, size_t chrlen, Base *misAlnRegBaseArr);
		void computeDisagrSubreg();
		void computeHighClipBaseNum();
		void computeZeroCovBaseNum();
};

#endif /* SRC_MISALNREG_H_ */
