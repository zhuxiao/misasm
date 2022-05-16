#ifndef SRC_BLOCK_H_
#define SRC_BLOCK_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <iterator>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "Base.h"
#include "structures.h"
#include "Paras.h"
#include "alnDataLoader.h"
#include "events.h"
#include "Region.h"
#include "varCand.h"
#include "clipReg.h"

//#define MIN_ADJACENT_REG_DIST		50
#define MIN_CLIP_NUM			5
#define NOISE_REG				500

using namespace std;

class Block{
	public:
		Paras *paras;
		string chrname;
		int64_t chrlen, startPos, endPos, winSize;      // 1-based position
		string workdir, outCovFile;
		vector<simpleReg_t*> sub_limit_reg_vec;
		bool process_flag;

		Base *baseArr;
		vector<bam1_t*> alnDataVector;
		faidx_t *fai;

		// SNV and indel
		vector<size_t> snvVector;
		vector<reg_t*> indelVector;
		vector<misAlnReg> misAlnRegVector;
		vector<reg_t*> zeroCovRegVector;	// for long deletions

		// clip regions
		vector<reg_t*> clipRegVector;
		vector<mateClipReg_t*> mateClipRegVector;

		// assembly information
		vector<varCand*> assembled_indel_vec;
		vector<varCand*> *assembled_chr_clipReg_vec;

		// detect abnormal signatures
		bool headIgnFlag, tailIgnFlag;

		double meanCov;  // excluding the gap regions filled by Ns

		// output file
		string out_dir_detect, out_dir_assemble, out_dir_call;
		string snvDetectPrefix_local, indelDetectPrefix_local, clipRegDetectPrefix_local, snvFilenameDetect, indelFilenameDetect, clipRegFilenameDetect;

		string indelAssemblePrefix_local, indelFilenameAssemble;
		ofstream *var_cand_indel_file, *misAln_reg_file, *var_cand_clipReg_file;

		int64_t TotalBaseNum;

		vector<bam1_t*> slideAlnVector;
		size_t tmp_startpos=0;

	public:
		Block(string chrname, size_t startPos, size_t endPos, faidx_t *fai, Paras *paras);
		virtual ~Block();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix);
		void setLimitRegs(vector<simpleReg_t*> &sub_limit_reg_vec);
		void setProcessFlag(bool process_flag);
		void setRegIngFlag(bool headIgnFlag, bool tailIgnFlag);
		void blockFillDataEst(size_t op_est);
		void blockDetect();
		void blockGenerateLocalAssembleWorkOpt();
		void setVarCandFiles(ofstream *var_cand_indel_file, ofstream *var_cand_clipReg_file);
		void resetVarCandFiles();
		void setMisAlnRegFile(ofstream *misAln_reg_file);
		void resetMisAlnRegFile();
		void setIlluminaMisAlnRegFile(ofstream *misAln_reg_file);
		void resetIlluminaMisAlnRegFile();
		void setAssembledChrClipRegVec(vector<varCand*> *assembled_chr_clipReg_vec);
		void blockAnalyse();
		void blockIlluminaDetect();
		void checkpos(int64_t startpos,  int64_t endpos);
		void checkposClip(int64_t startpos,  int64_t endpos);
		void estTotalBaseNum(int64_t startpos,  int64_t endpos);
		void estcoveragevalue(int64_t startpos,  int64_t endpos);
		void estdisagreement(int64_t startpos,  int64_t endpos);
		void estzerocov(int64_t startpos,  int64_t endpos);
		void estEvent(int64_t startpos, int64_t endpos);

		bool isHighIndelReg(int64_t startpos, int64_t endpos, float threshold, float ignore_polymer_ratio_threshold);
		bool isIlluminaHighIndelReg(int64_t startpos, int64_t endpos, float threshold);
		bool isHighInsReg(int64_t startpos, int64_t endpos, float threshold, float ignore_polymer_ratio_threshold);
		bool isHighDelReg(int64_t startpos, int64_t endpos, float threshold, float ignore_polymer_ratio_threshold);

	private:
		void destroyBaseArray();
		void destoryAlnData();
		void destroySnvVector();
		void destroyIndelVector();
		void destroyClipRegVector();
		void destroyMisAlnRegVector();
		Base *initBaseArray();
		int loadAlnData();
		int computeBlockBaseInfo();
		int computeBlockBaseInfoIllumina();
		void computeBlockMeanCov();
		void fillDataEst(size_t op_est);
		void AddReadLenEstInfo();
		int outputCovFile();
		void maskMisAlnRegs();
		int computeMisAlnDisagrReg();
		void extractMisAlnRegions();
		void saveMisAlnRegToFile();
		void computeDisagrNumSingleRegion(size_t startRpos, size_t endRPos, size_t regFlag);
		bool isMisAlnReg(Region &reg);
		int computeAbSigs();
		int computeAbSigsIllumina(vector<bam1_t*> &alnDataVector);
		int processSingleRegion(int64_t startRpos, int64_t endRPos, int64_t regFlag);
		int processSingleRegionIllumina(vector<bam1_t*> &alnDataVector, int64_t startRpos, int64_t endRPos, int64_t regFlag);
		void copySVEvents(Region &reg);
		void computeZeroCovReg(Region &reg);
		void updateZeroCovRegUsingIndelReg(vector<reg_t*> &zeroCovRegVec, vector<reg_t*> &indelVec);
		void updateSVRegUsingLongZeroCov();
		void updateClipRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void updateIndelRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec);
		void removeFalseIndel();
		void removeOverlapIndel();
		void removeOverlapClip();
		void removeFalseSNV();
		void sortRegVec(vector<reg_t*> &regVector);
		void blockGenerateLocalAssembleWorkOpt_Indel();
		void blockGenerateLocalAssembleWorkOpt_ClipReg();
		void saveSV2File();
		void saveMisassembly2File();
		vector<simpleReg_t*> computeLimitRegsForAssembleWork(vector<reg_t*> &varVec, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void generateAssembleWork(vector<reg_t*> &varVec, bool limit_reg_process_flag, vector<simpleReg_t*> &sub_limit_reg_vec_work, bool clip_reg_flag);
		bool getPrevAssembledDoneFlag(string &contigfilename, bool clipReg_flag);
		bool getPrevAssembledDoneFlag2(string &contigfilename, vector<varCand*> *assembled_varcand_vec);
		bool getPrevAssembledDoneFlag2(string &contigfilename, vector<string> &assembled_filename_vec);

		// duplication and inversion
		void mergeOverlappedClipReg();

		// reads aln data of slide window
		int detectIlluminaSlideAlnData(int64_t startRpos, int64_t endRPos);
		//get the accurate clipping regions from clipRegVector
//		int getAccurateClipReg(vector<reg_t*> &regVector);
};

#endif /* SRC_BLOCK_H_ */
