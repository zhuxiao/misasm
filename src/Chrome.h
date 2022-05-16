#ifndef SRC_CHROME_H_
#define SRC_CHROME_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <htslib/faidx.h>

#include "structures.h"
#include "Paras.h"
#include "Block.h"
#include "util.h"
#include "varCand.h"
#include "clipReg.h"

using namespace std;

#define BLOCK_SIZE_EST		20000
#define DIST_CHR_END		20000
#define MIN_CHR_SIZE_EST	(DIST_CHR_END*2+BLOCK_SIZE_EST)

#define REFSEQ_PATTERN		"refseq"
#define CLIPREG_PATTERN		"clipReg_refseq"

#define MIN_CHR_LEN			1000

class Chrome{
	public:
		Paras *paras;
		vector<Chrome*> *chr_vec;
		string chrname;
		int64_t chrlen;
		//bool process_flag;
		bool print_flag;

		// output directory
		string out_dir_detect, out_dir_assemble, out_dir_call;
		string out_filename_detect_snv, out_filename_detect_indel, out_filename_detect_clipReg;
		string out_filename_call_snv, out_filename_call_indel, out_filename_call_clipReg;

		vector<mateClipReg_t*> mateClipRegVector;
		vector<varCand*> var_cand_vec;
		vector<varCand*> var_cand_clipReg_vec;

		vector<Block*> blockVector;

		vector<varCand*> assembled_chr_clipReg_vec;		// previously assembled file names
		vector<varCand*> blat_aligned_chr_varCand_vec;	// previously blat aligned information
		vector<varCand*> blat_aligned_chr_clipReg_varCand_vec;	// previously blat aligned information


	//private:
		int32_t blockNum, process_block_num;
		string blocks_out_file;
		faidx_t *fai;

		// clip regions
		vector<reg_t*> clipRegVector;  // only for duplication and inversion

		// call
		vector<reg_t*> mis_aln_vec;

		// assembly info and misAln region
		string var_cand_indel_filename, misAln_reg_filename, var_cand_clipReg_filename;
		ofstream var_cand_indel_file, misAln_reg_file, var_cand_clipReg_file;

		// blat align information
		string blat_var_cand_indel_filename, blat_var_cand_clipReg_filename;
		ofstream blat_var_cand_indel_file, blat_var_cand_clipReg_file;


	public:
		Chrome(string& chrname, int chrlen, faidx_t *fai, Paras *paras, vector<Chrome*> *chr_vec);
		virtual ~Chrome();
		void setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix);
		void setMisasmOutputDir(string& out_dir_detect_prefix);
		string getVarcandIndelFilename();
		string getVarcandClipregFilename();
		int generateChrBlocks();
		void saveChrBlocksToFile();
		void chrFillDataEst(size_t op_est);
		int chrDetect();
		void removeFPIndelSnvInClipReg(vector<mateClipReg_t*> &mate_clipReg_vec);;
		void chrMergeDetectResultToFile();
		void chrMergeDetectResultToFileIllumina();
		void chrMergeMisasmResultToFile();
		void loadPrevAssembledInfo();
		void chrLoadDataAssemble();
		int chrGenerateLocalAssembleWorkOpt();
		void chrResetAssembleData();
		void chrLoadDataCall();
		int chrCall();
		void mergeSameRegTRA();
		void removeVarCandNodeIndel(varCand *var_cand);
		void removeVarCandNodeClipReg(varCand *var_cand);
		void chrFillVarseq();
		void saveCallSV2File();
		void chrLoadMateClipRegData();

		int chrIlluminaMisIdentify();


	private:
		void init();
		Block *allocateBlock(string& chrname, size_t chrlen, int begPos, int endPos, faidx_t *fai, bool headIgnFlag, bool tailIgnFlag, bool block_process_flag);
		void destroyBlockVector();
		void destroyVarCandVector(vector<varCand*> &var_cand_vec);
		void destroyMisAlnVector();
		void destroyClipRegVector();
		void destroyMateClipRegVector();
		void chrLoadIndelDataAssemble();
		void chrLoadIndelData(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void chrLoadClipRegDataAssemble();
		int chrDetect_st();
		int chrDetect_mt();
		void removeRedundantIndelDetect();
		void removeRedundantIndelItemDetect(reg_t *reg, size_t bloc_idx, size_t indel_vec_idx);
		void chrSetVarCandFiles();
		void chrResetVarCandFiles();
		void chrSetMisAlnRegFile();
		void chrResetMisAlnRegFile();
		void chrSetIlluminaMisAlnRegFile();
		void chrResetIlluminaMisAlnRegFile();
		Block *computeBlocByPos(int64_t begPos, vector<Block*> &block_vec);
		int32_t computeBlocID(int64_t begPos, vector<Block*> &block_vec);
		int chrGenerateLocalAssembleWorkOpt_st();
		int chrGenerateLocalAssembleWorkOpt_mt();
		void outputAssemDataToFile(string &filename);
		void removeVarCandNode(varCand *var_cand, vector<varCand*> &var_cand_vec);
		void loadPrevAssembledInfo2(bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		string getAssembleFileHeaderLine();

		void chrCall_st();
		void chrCall_mt();
		void chrCallVariants(vector<varCand*> &var_cand_vec);
		void loadVarCandData();
		void loadVarCandDataFromFile(vector<varCand*> &var_cand_vec, string &var_cand_filename, bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void loadClipRegCandData();
		void sortVarCandData(vector<varCand*> &var_cand_vec);
		bool isVarCandDataSorted(vector<varCand*> &var_cand_vec);
		void loadMisAlnRegData();
		void sortMisAlnRegData();
		void loadPrevBlatAlnItems(bool clipReg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		void setBlatVarcandFiles();
		string getBlatVarcandFileHeaderLine();
		void resetBlatVarcandFiles();
		void chrFillVarseqSingleVec(vector<varCand*> &var_cand_vec);
		void removeRedundantVar();
		void removeRedundantIndel(vector<varCand*> &var_cand_vec);
		//void removeRedundantClipReg(vector<varCand*> &var_cand_clipReg_vec);
		bool isRedundantVarItemBinSearch(reg_t *reg, vector<varCand*> &var_cand_vec);
		void removeRedundantVarItemsInNewCalledVarvec(reg_t *reg_idx, int32_t idx, vector<varCand*> &var_cand_vec);
		void removeFPNewVarVec();
		void removeFPNewVarVecIndel(vector<varCand*> &var_cand_vec);
		bool isIndelInClipReg(reg_t *reg, vector<mateClipReg_t*> &mate_clipReg_vec);
		bool isSnvInClipReg(size_t pos, vector<mateClipReg_t*> &mate_clipReg_vec);

		// DUP, INV, TRA
		void chrComputeMateClipReg();
		void removeRedundantItemsClipReg(vector<reg_t*> &clipReg_vec, vector<bool> &clip_processed_flag_vec);
		void processClipRegs(size_t idx, vector<bool> &clip_processed_flag_vec, mateClipReg_t &mate_clip_reg, reg_t *reg);
		mateClipReg_t *getMateRegEndSameClipReg(mateClipReg_t *clip_reg, vector<mateClipReg_t*> &mate_clipReg_vec);
		void removeFPClipRegsDupInv();
		vector<size_t> getOverlapClipReg(reg_t *given_reg);
		void chrLoadMateClipRegDataOp(bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec);
		string getDirnameCall(string &chrname_given);
		mateClipReg_t* getMateClipReg(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, string &chrname_mate_clip_reg);
		mateClipReg_t* getMateClipRegOp(reg_t *reg1, reg_t *reg2, int32_t *clip_reg_idx_tra, vector<mateClipReg_t*> &mate_clipReg_vec);

		// save SV to file
		void saveCallIndelClipReg2File(string &outfilename_indel, string &outfilename_clipReg);

		// identify
		int chrIlluminaMisIdentify_st();
		int chrIlluminaMisIdentify_mt();
};

#endif /* SRC_CHROME_H_ */
