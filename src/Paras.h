#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include "structures.h"

using namespace std;

// program variables
//#define PROG_NAME					"ASVCLR"
//#define PROG_DESC					"Accurate Structural Variant Caller for Long Reads"
//#define PROG_VERSION				"0.7.0"
#define PROG_NAME					"misasm"
#define PROG_DESC					"misassembly detector"
#define PROG_VERSION				"0.1.0"

#define SIZE_EST_OP					0
#define NUM_EST_OP					1
#define SNV_EST_OP					2

#define LARGE_INDEL_SIZE_FACTOR		3

#define MD_MISMATCH   				10
#define BASE_INDEL_CON_UNUSED		(-1)

// default parameter values
//#define BLOCKSIZE  					1000000
//#define SLIDESIZE  					500
//#define ASSEM_SLIDE_SIZE			5000
#define BLOCKSIZE  					1000000
#define SLIDESIZE  					1000
#define SUBBLOCKSIZE  				100	//20220502
#define GAP_EXTEND_REG_SIZE  		500	//20220502
#define ASSEM_SLIDE_SIZE			5000

#define MIN_INDEL_EVENT_SIZE		2
#define MIN_INDEL_EVENT_NUM			5
#define MIN_SV_SIZE_USR				2

#define MIN_CLIP_READS_NUM_THRES	7	// to be parameterized

#define MAX_CLIP_REG_SIZE			10000  // to be parameterized

#define SIZE_PERCENTILE_EST			0.95
#define NUM_PERCENTILE_EST			0.99995
#define AUX_ARR_SIZE				1001

#define MAX_REG_SUM_SIZE_EST		500000

#define ASSEMBLY_SUCCESS			"ASS_SUCCESS"
#define ASSEMBLY_FAILURE			"ASS_FAILURE"
#define ALIGN_SUCCESS				"ALN_SUCCESS"
#define ALIGN_FAILURE				"ALN_FAILURE"
#define CALL_SUCCESS				"CALL_SUCCESS"
#define CALL_FAILURE				"CALL_FAILURE"

#define DONE_STR					"DONE"

#define SAMPLED_STR					"SAMPLED"
#define UNSAMPLED_STR				"UNSAMPLED"

#define LIMIT_REG_ALL_STR			"ALL"

#define EXPECTED_COV_ASSEMBLE		30.0f
#define MASK_VAL_DEFAULT			0

#define NUM_PARTS_PROGRESS			50
#define NUM_THREADS_PER_ASSEM_WORK	0  // 0: unspecify the limited number of threads for each Canu work

#define OUT_DIR						"output"

#define MAX_ULTRA_HIGH_COV_THRES	300		// maximal coverage threshold for ultra-high coverage

#define MIN_ADJACENT_REG_DIST		50

#define MIN_HIGH_CONSENSUS_INS_RATIO		0.3f
#define MIN_HIGH_CONSENSUS_DEL_RATIO		0.4f

//20220502
#define MIN_ABNORMAL_PAIRING_RATIO			0.3
#define HIGH_CLIP_BASE_RATIO				0.8
#define MIN_INDEL_BASE_RATIO				0.3

// program parameters
class Paras
{
	public:
		// user/system defined parameters
		string command, refFile, inBamFile, outFilePrefix; //, canu_version;
		string outDir;
		string out_dir_detect = "misassemblies";	//20220429
//		string out_dir_detect = "1_candidates";    // "1_candidates"
		string out_dir_assemble = "2_assemble";  // "2_assemble"
		string out_dir_call = "3_call";      // "3_call"
		string out_dir_tra = out_dir_call + "/" + "tra";
		string out_dir_result = "4_results";	// "4_results"
		size_t blockSize, slideSize, assemSlideSize, min_sv_size_usr, num_threads, large_indel_size_thres;

		size_t subblockSize, gapextendSize;	//20220502
		float abpairRatio, highclipRatio, indelRatio;//20220502

		bool maskMisAlnRegFlag, load_from_file_flag;
		size_t misAlnRegLenSum = 0;
		size_t minClipReadsNumSupportSV;

		// limit SV regions, item format: CHR | CHR:START-END
		vector<simpleReg_t*> limit_reg_vec;
		bool limit_reg_process_flag = false;	// true for limit regions; default is false for disable limit regions (i.e. process all regions)
		string limit_reg_filename = "limit_regions.bed";

		size_t maxClipRegSize;

		size_t mean_read_len, total_read_num_est;
		size_t reg_sum_size_est, max_reg_sum_size_est;

		// reads sampling parameters
		double expected_cov_assemble;
		bool delete_reads_flag;

		// clipping reads sampling parameters
		double max_ultra_high_cov;

		// estimated parameters
		size_t min_ins_size_filt, min_del_size_filt, min_clip_size_filt;//equal 0
		size_t min_ins_num_filt, min_del_num_filt, min_clip_num_filt;//equal 0

		// auxiliary data for estimation
		size_t insSizeEstArr[AUX_ARR_SIZE], delSizeEstArr[AUX_ARR_SIZE], clipSizeEstArr[AUX_ARR_SIZE];
		size_t insNumEstArr[AUX_ARR_SIZE], delNumEstArr[AUX_ARR_SIZE], clipNumEstArr[AUX_ARR_SIZE];

		// assembly regions for thread pool
		vector<assembleWork_opt*> assem_work_vec;
		size_t assemble_reg_preDone_num, assemble_reg_work_total, assemble_reg_workDone_num, num_parts_progress;
		pthread_mutex_t mtx_assemble_reg_workDone_num;
		size_t num_threads_per_assem_work;

		// previously assembled regions
		vector<string> assembled_clipReg_filename_vec;

		double insert_isize, insert_sdev, insert_max, insert_min;

	public:
		Paras();
		Paras(int argc, char **argv);
		virtual ~Paras();
		void initEst();
		void estimate(size_t op_est);
		void outputParas();
		void outputEstParas(string info);
		int checkBamFile();

	private:
		void init();
		//string getCanuVersion();
		int parseParas(int argc, char **argv);
		int misasmParas(int argc, char **argv);	//20220502
		int parseMisasmParas(int argc, char **argv);//20220502
		int parseDetectParas(int argc, char **argv);
		int parseAssembleParas(int argc, char **argv);
		int parseCallParas(int argc, char **argv);
		int parseAllParas(int argc, char **argv);
		void showUsage();
		void showDetectUsage();
		void showAssembleUsage();
		void showCallUsage();
		void showAllUsage();
		size_t estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val);
};

#endif /* SRC_PARAS_H_ */
