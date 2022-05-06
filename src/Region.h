#ifndef SRC_PARSER_REGION_H_
#define SRC_PARSER_REGION_H_

#include <iostream>
#include <fstream>
#include <htslib/faidx.h>
#include <time.h>

#include "structures.h"
#include "Paras.h"
#include "Base.h"
#include "misAlnReg.h"

using namespace std;

#define HEAD_REGION		0
#define INNER_REGION	1
#define TAIL_REGION		2

#define REG_DIF_NONE		0
#define REG_DIF_SNV			1
#define REG_DIF_INDEL		2
#define REG_DIF_BOTH		3
#define REG_DIF_UNCERTAIN	4

#define SUB_REG_SIZE				10
#define SUB_CLIP_REG_SIZE			100
//#define SUB_MISJOIN_REG_SIZE			100
//#define SUB_EXTEND_REG_SIZE			0//500//100//500
//#define GAP_EXTEND_REG_SIZE		500	//20220502

#define LARGE_INDEL_RATIO_THRES				(0.1f)
#define HIGH_INDEL_CLIP_RATIO_THRES			(0.6f)
#define SECOND_INDEL_CLIP_RATIO_THRES		(0.3f)
#define HIGH_INDEL_CLIP_BASE_RATIO_THRES	(0.1f)

#define MIN_COVRATIO		0.6
#define MAX_COVRATIO		2

#define MIN_CLIP_NUM		5
#define MIN_GAP_SIZE		5

#define ILLUMINA_HIGH_CLIP_THRES			0.8	//0.7
#define ILLUMINA_HIGH_CLIP_BASE_THRES		0.8	//0.5

#define MAX_PAIR_THRES		0.3	//
#define MAX_PAIR_THRES_POS		0.1	//20220426

// 20220421
// #define MISASSEMBLY_REGION_START_THRES		501

//#define NOISE_REG			500	//20220421

class Region {
	public:
		Paras *paras;
		string chrname;
		int64_t startRPos, endRPos, startMidPartPos, endMidPartPos, chrlen, minRPos, maxRPos;
		Base *regBaseArr;  // the local region base array
		size_t regFlag, subRegSize;
		bool wholeRefGapFlag; // true -- if all the bases in the region of reference are 'N'; false -- otherwise

		double meanBlockCov, localRegCov;
		size_t highIndelSubRegNum;

		// candidate flag
		bool indelCandFlag;
		size_t difType;

		// output file
		string out_dir_assemble;

		double fragmentsize, ratio_exisizemax, ratio_lessisizemin, ratio_abdirection, ratio_abref;

		int64_t abpairstartpos, abpairendpos;	//20220424

		string misType;	//20220504

	private:
		vector<int64_t> disagrePosVector;  // element: the reference position with disagreement
		vector<int64_t> zeroCovPosVector;  // element: the reference position with zero coverage
		vector<reg_t*> abCovRegVector;    // element: the reference region with abnormal coverage, format: pos1-pos2

		// SNV and indel
		vector<int64_t> snvVector;
		vector<reg_t*> indelVector;
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

	public:
		Region(string& chrname, int64_t startRpos, int64_t endRPos, int64_t chrlen, int64_t minRPos, int64_t maxRPos, Base *regBaseArr, size_t regFlag, Paras *paras);
		virtual ~Region();
		void setOutputDir(string& out_dir_assemble_prefix);

	public:
		void setMeanBlockCov(double meanBlockCov);
		int computeRegAbSigs();
		int computeDisagreements();
		void printRegAbSigs();
		void detectSNV();
		void detectIndelReg();
		void detectIlluminaIndelReg();
		void detectHighClipReg();
		void detectIlluminaHighClipReg(vector<bam1_t*> &slideAlnVector);
		void detectIlluminaMisjoinReg(vector<bam1_t*> &slideAlnVector);	//20220423
		void detectIlluminaGapMisjoinReg(vector<bam1_t*> &slideAlnVector);	//20220423
		void detectIlluminaBaseClipReg();	//20220422
		size_t getDisagreeNum();
		size_t getSNVNum();
		size_t getIndelNum();
		void determineDifType();
		vector<int64_t> getSnvVector();
		vector<reg_t*> getIndelVector();
		vector<reg_t*> getClipRegVector();
		vector<int64_t> getZeroCovPosVector();



	private:
		bool IsWholeRefGap();
		void addDisagrePos(int64_t pos);
		void addZeroCovPos(int64_t pos);
		void destroyDisagrePosVector();
		void destroyZeroCovPosVector();
		void destroyAbCovRegVector();
		void destroySnvVector();
		void destroyIndelVector();
		void destroyClipRegVector();
		int computeAbCovReg();
		reg_t* allocateReg(string &chrname, int64_t startPosReg, int64_t endPosReg);
		reg_t* allocateMisReg(string &chrname, int64_t startPosReg, int64_t endPosReg, string misType);	//20220504
		double computeMeanCovReg(int64_t startPosReg, int64_t endPosReg);
		int computeHighIndelEventRegNum();
		int32_t computeReadIndelEventNumReg(int64_t startPosReg, int64_t endPosReg);
		bool computeSNVFlag(int64_t pos, int64_t startCheckPos, int64_t endCheckPos);
		bool haveMuchShortIndelsAround(int64_t startCheckPos, int64_t endCheckPos);
		reg_t* getIndelReg(int64_t startCheckPos);
		reg_t* getIlluminaIndelReg(int64_t startCheckPos);
		bool haveNoAbSigs(Base *base, int64_t pos);
		int32_t getMismatchBasesAround(int64_t pos1, int64_t pos2);
		int32_t getDisZeroCovNum(int64_t startPos, int64_t endPos);
		int32_t getLargeIndelBaseNum(int64_t startPos, int64_t endPos);
		int32_t getLargeIndelNum(int64_t startPos, int64_t endPos);
		int32_t getHighConIndelNum(int64_t startPos, int64_t endPos, float threshold, float polymer_ignore_ratio_thres);
		int32_t getIlluminaHighConIndelNum(int64_t startPos, int64_t endPos, float threshold);

		// duplication and inversion
		reg_t* getClipReg(int64_t startCheckPos);
		bool haveNoClipSig(int64_t startPos, int64_t endPos, double clip_ratio_thres);

		reg_t* getIlluminaClipReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);
		bool isIlluminaClipReg(int64_t startPos, int64_t endPos, double clip_ratio_thres);
		bool isIlluminaBaseClipReg(int64_t startPos, int64_t endPos, double clip_ratio_thres);
		bool isIlluminaAbisize(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos);
		bool isIlluminaAbPairing(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos);

		//	20220421
		//	get Misjoin position
		reg_t* getIlluminaBaseClipReg(int64_t startCheckPos);
		reg_t* getIlluminaMisjoinReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);
		reg_t* getIlluminaGapMisjoinReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);
		bool isIlluminaAbPair(vector<bam1_t*> &slideAlnVector, float threshold);
		int estIlluminaAbPairStartPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);
		int estIlluminaAbPairEndPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);
		bool isIlluminaAbPairPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos);

		bool isIlluminaGapClipReg(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos);
		bool isIlluminaGapReg(int64_t startPos, int64_t endPos);

		void estIlluminaPairing(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos);
		void estIlluminaRegInsertsize(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos);
};

#endif /* SRC_PARSER_REGION_H_ */
