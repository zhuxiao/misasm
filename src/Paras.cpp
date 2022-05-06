#include <string.h>
#include <vector>
#include <htslib/hts.h>
#include <cmath>
#include "Paras.h"
#include "util.h"

// Constructor with parameters
Paras::Paras(){
	init();
}

// Constructor with parameters
/*Paras::Paras(int argc, char **argv){
	init();
	if(parseParas(argc, argv)!=0) exit(1);

	// check bam file
	if(checkBamFile()!=0) exit(1);
}
*/

// 20220502
Paras::Paras(int argc, char **argv){
	init();
//	if(parseParas(argc, argv)!=0) exit(1);
	if(misasmParas(argc, argv)!=0) exit(1);

	// check bam file
	if(checkBamFile()!=0) exit(1);
}

//Destructor
Paras::~Paras(){
	if(!limit_reg_vec.empty()) destroyLimitRegVector(limit_reg_vec);
}

// initialization
void Paras::init(){
	command = "";
	inBamFile = "";
	outFilePrefix = "";
	outDir = OUT_DIR;
	num_threads = 0;

	min_ins_size_filt = 0;
	min_del_size_filt = 0;
	min_clip_size_filt = 0;
	min_ins_num_filt = 0;
	min_del_num_filt = 0;
	min_clip_num_filt = 0;

	mean_read_len = total_read_num_est = 0;

	reg_sum_size_est = 0;
	max_reg_sum_size_est = MAX_REG_SUM_SIZE_EST;

	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	max_ultra_high_cov = MAX_ULTRA_HIGH_COV_THRES;

	assemble_reg_preDone_num = assemble_reg_work_total = assemble_reg_workDone_num = 0;
	num_parts_progress = NUM_PARTS_PROGRESS;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;

	//canu_version = getCanuVersion();
}

// check Bam file, and generate the BAM index if it is unavailable
int Paras::checkBamFile(){
	samFile *in = 0;
	string idx_filename;
	int ret = 1;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile.c_str() << " for reading" << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == NULL) { // index is unavailable, then generate it
		if(in) sam_close(in);

		cout << __func__ << ": BAM index is unavailable, now generate it, please wait ...\n" << endl;

		// construct the index file name
		if(inBamFile.substr(inBamFile.size()-4).compare(".bam")==0)
			idx_filename = inBamFile.substr(0, inBamFile.size()-4) + ".bai";
		else
			idx_filename = inBamFile + ".bai";

		ret = sam_index_build3(inBamFile.c_str(), idx_filename.c_str(), 0, num_threads);
		switch (ret) {
			case 0:
				break;
			case -1:
				cerr << __func__ << ", error occurred while opening " << inBamFile << endl;
				exit(-1);
			case -2:
				cerr << __func__ << ", failed to open " << inBamFile << endl;
				exit(-2);
			case -3:
				cerr << __func__ << inBamFile << " is in a format that cannot be usefully indexed" << endl;
				exit(-3);
			case -4:
				cerr << __func__ << ", failed to create or write index" << endl;
				exit(-4);
			default:
				cerr << __func__ << ", failed to create index for " << inBamFile << endl;
				exit(1);
		}
	}else{
		ret = 0;
		hts_idx_destroy(idx); // destroy the BAM index
		if(in) sam_close(in);
	}

	return ret;
}

// parse the parameters
int Paras::parseParas(int argc, char **argv){
	if (argc < 2) { showUsage(); return 1; }

    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { showUsage(); exit(0); }
        argv++;
        argc = 2;
    }

    if (strcmp(argv[1], "detect")==0){
    	if(argc==2){ showDetectUsage(); exit(0); }
    	command = "detect";
    	return parseDetectParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "assemble")==0){
    	if(argc==2){ showAssembleUsage(); exit(0); }
    	command = "assemble";
    	return parseAssembleParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "call")==0){
    	if(argc==2){ showCallUsage(); exit(0); }
    	command = "call";
    	return parseCallParas(argc-1, argv+1);
    }else if(strcmp(argv[1], "all")==0){
    	if(argc==2){ showAllUsage(); exit(0); }
    	command = "all";
    	return parseAllParas(argc-1, argv+1);
    }else{
    	cerr << "Error: invalid command " << argv[1] << endl << endl;
    	showUsage(); return 1;
    }
}

//20220502
int Paras::misasmParas(int argc, char **argv){
	if (argc < 2) { showUsage(); return 1; }

    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { showUsage(); exit(0); }
        argv++;
        argc = 1;
    }

    return parseMisasmParas(argc-1, argv);
}

// parse the parameters for misasm	20220503
int Paras::parseMisasmParas(int argc, char **argv){
	int opt, threadNum_tmp = 1, mask_val;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;

	subblockSize = SUBBLOCKSIZE;	//20220502
	gapextendSize = GAP_EXTEND_REG_SIZE; //20220502
	abpairRatio = MIN_ABNORMAL_PAIRING_RATIO; //20220502
	highclipRatio = HIGH_CLIP_BASE_RATIO; //20220502
	indelRatio = MIN_INDEL_BASE_RATIO; //20220502

	min_sv_size_usr = MIN_SV_SIZE_USR;
	minClipReadsNumSupportSV = MIN_CLIP_READS_NUM_THRES;
	maxClipRegSize = MAX_CLIP_REG_SIZE;
	mask_val = MASK_VAL_DEFAULT;
	outDir = OUT_DIR;
	simpleReg_t *simple_reg;
	string simple_reg_str;

	while( (opt = getopt(argc, argv, ":b:s:u:g:a:l:i:c:o:p:t:M:h")) != -1 ){	//20220502
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			case 'u': subblockSize = stoi(optarg); break;	//20220502
			case 'g': gapextendSize = stoi(optarg); break;	//20220502
			case 'a': abpairRatio = stof(optarg); break;	//20220502
			case 'l': highclipRatio = stof(optarg); break;	//20220502
			case 'i': indelRatio = stof(optarg); break;		//20220502
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'n': minClipReadsNumSupportSV = stoi(optarg); break;
			case 'c': maxClipRegSize = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': mask_val = stoi(optarg); break;
			case 'h': showDetectUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	load_from_file_flag = true;
	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	outDir = deleteTailPathChar(outDir);

	if(mask_val==1) maskMisAlnRegFlag = true;
	else if(mask_val==0) maskMisAlnRegFlag = false;
	else{
		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
		showUsage();	//20220503
		return 1;
	}

	delete_reads_flag = true;

	opt = argc - optind; // the reference and BAM file on the command line
	if(opt>=1) {	//20220503
		refFile = argv[optind];
		inBamFile = argv[optind+1];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for detect command
int Paras::parseDetectParas(int argc, char **argv){
	int opt, threadNum_tmp = 1, mask_val;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;

	subblockSize = SUBBLOCKSIZE;	//20220502
	gapextendSize = GAP_EXTEND_REG_SIZE; //20220502
	abpairRatio = MIN_ABNORMAL_PAIRING_RATIO; //20220502
	highclipRatio = HIGH_CLIP_BASE_RATIO; //20220502
	indelRatio = MIN_INDEL_BASE_RATIO; //20220502

	min_sv_size_usr = MIN_SV_SIZE_USR;
	minClipReadsNumSupportSV = MIN_CLIP_READS_NUM_THRES;
	maxClipRegSize = MAX_CLIP_REG_SIZE;
	mask_val = MASK_VAL_DEFAULT;
	outDir = OUT_DIR;
	simpleReg_t *simple_reg;
	string simple_reg_str;

//	while( (opt = getopt(argc, argv, ":b:s:c:o:p:t:M:h")) != -1 ){
	while( (opt = getopt(argc, argv, ":b:s:u:g:a:l:i:c:o:p:t:M:h")) != -1 ){	//20220502
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			case 'u': subblockSize = stoi(optarg); break;	//20220502
			case 'g': gapextendSize = stoi(optarg); break;	//20220502
			case 'a': abpairRatio = stof(optarg); break;	//20220502
			case 'l': highclipRatio = stof(optarg); break;	//20220502
			case 'i': indelRatio = stof(optarg); break;		//20220502
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'n': minClipReadsNumSupportSV = stoi(optarg); break;
			case 'c': maxClipRegSize = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': mask_val = stoi(optarg); break;
			case 'h': showDetectUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	load_from_file_flag = true;
	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	outDir = deleteTailPathChar(outDir);

	if(mask_val==1) maskMisAlnRegFlag = true;
	else if(mask_val==0) maskMisAlnRegFlag = false;
	else{
		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
		showDetectUsage();
		return 1;
	}

	delete_reads_flag = true;

	opt = argc - optind; // the reference and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for assemble command
int Paras::parseAssembleParas(int argc, char **argv){
	int opt, threadNum_tmp = 1, mask_val, delete_reads_val;
	blockSize = BLOCKSIZE;
	slideSize = ASSEM_SLIDE_SIZE;
	assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	minClipReadsNumSupportSV = MIN_CLIP_READS_NUM_THRES;
	maxClipRegSize = MAX_CLIP_REG_SIZE;
	mask_val = MASK_VAL_DEFAULT;
	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	delete_reads_val = 1;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;
	outDir = OUT_DIR;

	while( (opt = getopt(argc, argv, ":b:S:m:n:c:x:o:p:t:T:M:R:h")) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'n': minClipReadsNumSupportSV = stoi(optarg); break;
			case 'c': maxClipRegSize = stoi(optarg); break;
			case 'x': expected_cov_assemble = stod(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'T': num_threads_per_assem_work = stoi(optarg); break;
			case 'M': mask_val = stoi(optarg); break;
			case 'R': delete_reads_val = stoi(optarg); break;
			case 'h': showAssembleUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	load_from_file_flag = true;
	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	if(num_threads*num_threads_per_assem_work>(size_t)sysconf(_SC_NPROCESSORS_ONLN)){ // warning
		cout << "Warning: the user-specified total number of concurrent assemble work is " << num_threads << ", and the user-specified number of threads for each assemble work is " << num_threads_per_assem_work << ", which exceeds the total number of available processors on the machine (" << sysconf(_SC_NPROCESSORS_ONLN) << ")." << endl;
	}

	outDir = deleteTailPathChar(outDir);

	if(mask_val==1) maskMisAlnRegFlag = true;
	else if(mask_val==0) maskMisAlnRegFlag = false;
	else{
		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
		showAssembleUsage();
		return 1;
	}
	if(delete_reads_val==1) delete_reads_flag = true;
	else if(delete_reads_val==0) delete_reads_flag = false;
	else{
		cout << "Error: Please specify the correct reads delete flag for local assembly using -R option." << endl << endl;
		showAllUsage();
		return 1;
	}

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showAssembleUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for call command
int Paras::parseCallParas(int argc, char **argv){
	int opt, threadNum_tmp = 1, mask_val, delete_reads_val;
	blockSize = BLOCKSIZE;
	slideSize = ASSEM_SLIDE_SIZE;
	assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	minClipReadsNumSupportSV = MIN_CLIP_READS_NUM_THRES;
	maxClipRegSize = MAX_CLIP_REG_SIZE;
	mask_val = MASK_VAL_DEFAULT;
	delete_reads_val = 1;
	outDir = OUT_DIR;

	while( (opt = getopt(argc, argv, ":b:S:m:n:c:o:p:t:M:R:h")) != -1 ){
		switch(opt){
			case 'b': blockSize = stoi(optarg); break;
			case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'n': minClipReadsNumSupportSV = stoi(optarg); break;
			case 'c': maxClipRegSize = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'M': mask_val = stoi(optarg); break;
			case 'R': delete_reads_val = stoi(optarg); break;
			case 'h': showCallUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	load_from_file_flag = true;
	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	outDir = deleteTailPathChar(outDir);

	if(mask_val==1) maskMisAlnRegFlag = true;
	else if(mask_val==0) maskMisAlnRegFlag = false;
	else{
		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
		showCallUsage();
		return 1;
	}
	if(delete_reads_val==1) delete_reads_flag = true;
	else if(delete_reads_val==0) delete_reads_flag = false;
	else{
		cout << "Error: Please specify the correct reads delete flag for local assembly using -R option." << endl << endl;
		showAllUsage();
		return 1;
	}

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showCallUsage();
		return 1;
	}

	return 0;
}

// parse the parameters for 'all' command
int Paras::parseAllParas(int argc, char **argv){
	int opt, threadNum_tmp = 1, mask_val, delete_reads_val;
	blockSize = BLOCKSIZE;
	slideSize = SLIDESIZE;
	assemSlideSize = ASSEM_SLIDE_SIZE;
	min_sv_size_usr = MIN_SV_SIZE_USR;
	minClipReadsNumSupportSV = MIN_CLIP_READS_NUM_THRES;
	maxClipRegSize = MAX_CLIP_REG_SIZE;
	mask_val = MASK_VAL_DEFAULT;
	expected_cov_assemble = EXPECTED_COV_ASSEMBLE;
	delete_reads_val = 1;
	num_threads_per_assem_work = NUM_THREADS_PER_ASSEM_WORK;
	outDir = OUT_DIR;
	simpleReg_t *simple_reg;
	string simple_reg_str;

	while( (opt = getopt(argc, argv, ":b:s:S:m:n:c:x:o:p:t:T:M:R:h")) != -1 ){
		switch(opt){
			case 'b': assemSlideSize = stoi(optarg); break;
			case 's': slideSize = stoi(optarg); break;
			case 'S': assemSlideSize = stoi(optarg); break;
			case 'm': min_sv_size_usr = stoi(optarg); break;
			case 'n': minClipReadsNumSupportSV = stoi(optarg); break;
			case 'c': maxClipRegSize = stoi(optarg); break;
			case 'x': expected_cov_assemble = stod(optarg); break;
			case 'o': outDir = optarg; break;
			case 'p': outFilePrefix = optarg; break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'T': num_threads_per_assem_work = stoi(optarg); break;
			case 'M': mask_val = stoi(optarg); break;
			case 'R': delete_reads_val = stoi(optarg); break;
			case 'h': showAllUsage(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	load_from_file_flag = false;
	num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	if(num_threads*num_threads_per_assem_work>(size_t)sysconf(_SC_NPROCESSORS_ONLN)){ // warning
		cout << "Warning: the user-specified total number of concurrent assemble work is " << num_threads << ", and the user-specified number of threads for each assemble work is " << num_threads_per_assem_work << ", which exceeds the total number of available processors on the machine (" << sysconf(_SC_NPROCESSORS_ONLN) << ")." << endl;
	}

	outDir = deleteTailPathChar(outDir);

	if(mask_val==1) maskMisAlnRegFlag = true;
	else if(mask_val==0) maskMisAlnRegFlag = false;
	else{
		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
		showAllUsage();
		return 1;
	}
	if(delete_reads_val==1) delete_reads_flag = true;
	else if(delete_reads_val==0) delete_reads_flag = false;
	else{
		cout << "Error: Please specify the correct reads delete flag for local assembly using -R option." << endl << endl;
		showAllUsage();
		return 1;
	}

	opt = argc - optind; // the reference file and BAM file on the command line
	if(opt>=2) {
		refFile = argv[optind];
		inBamFile = argv[optind+1];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showAllUsage();
		return 1;
	}

	if(refFile.size()==0){
		cout << "Error: Please specify the reference" << endl << endl;
		showAllUsage();
		return 1;
	}

	return 0;
}

// show the usage
/*void Paras::showUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib V" << hts_version() << ")" << endl << endl;
	cout << "Usage: misasm <command> [options] <ASM_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "     ASM_FILE     Assembled sequence file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted BAM file (required)" << endl;
	cout << "     Region       Limit assembled sequence region to process: CHR|CHR:START-END." << endl;
	cout << "                  If unspecified, all assembled sequence regions will be " << endl;
	cout << "                  processed (optional)" << endl << endl;

	cout << "Commands:" << endl;
	cout << "     detect       detect misjoin and indel signatures in aligned reads" << endl;
//	cout << "     assemble     assemble candidate regions" << endl;
//	cout << "     call         call indels by alignments of local genome assemblies" << endl;
//	cout << "     all          run the above commands in turn" << endl;
}
*/

// 20220502 show the usage
void Paras::showUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
//	cout << "Usage: misasm [options] <ASM_FILE> <BAM_FILE> [Region ...]" << endl << endl;
	cout << "Usage: misasm [options] <ASM_FILE> <BAM_FILE>" << endl << endl;//20220504

	cout << "Description:" << endl;
	cout << "     ASM_FILE     Assembled sequence file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted BAM file (required)" << endl;
//	cout << "     Region       Limit assembled sequence region to process: CHR|CHR:START-END." << endl;
//	cout << "                  If unspecified, all assembled sequence regions will be " << endl;
//	cout << "                  processed (optional)" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -s INT       detect slide size [1000]" << endl;
	cout << "     -u INT       detect subblock size [100]" << endl;	//20220502
	cout << "     -g INT       gap extend region size [500]" << endl;	//20220502
	cout << "     -a FLOAT     minimum proportion of abnormal pairing reads in misassemblies [0.3]" << endl;	//20220502
	cout << "     -l FLOAT     minimum proportion of high-clipping reads in misassemblies [0.8]" << endl;	//20220502
	cout << "     -i FLOAT     minimum proportion of indel reads in misassemblies [0.3]" << endl;	//20220502

	cout << "     -n INT       minimal clipping reads supporting a misassembly [7]" << endl;
	cout << "     -c INT       maximal clipping region size to detect [10000]" << endl;
//	cout << "     -r FILE      limit reference regions to process [null]: CHR|CHR:START-END" << endl;
	cout << "     -o DIR       output directory [output]" << endl;
	cout << "     -p STR       prefix of output result files [null]" << endl;
//	cout << "     -t INT       number of concurrent work [1]" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for detect command
void Paras::showDetectUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: misasm detect [options] <ASM_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "     ASM_FILE     Assembled sequence file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted BAM file (required)" << endl;
	cout << "     Region       Limit assembled sequence region to process: CHR|CHR:START-END." << endl;
	cout << "                  If unspecified, all assembled sequence regions will be " << endl;
	cout << "                  processed (optional)" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -s INT       detect slide size [1000]" << endl;
	cout << "     -u INT       detect subblock size [100]" << endl;	//20220502
	cout << "     -g INT       gap extend region size [500]" << endl;	//20220502
	cout << "     -a FLOAT     minimum proportion of abnormal pairing reads in misassemblies [0.3]" << endl;	//20220502
	cout << "     -l FLOAT     minimum proportion of high-clipping reads in misassemblies [0.8]" << endl;	//20220502
	cout << "     -i FLOAT     minimum proportion of indel reads in misassemblies [0.3]" << endl;	//20220502

//	cout << "     -m INT       minimal SV size to detect [2]" << endl;
	cout << "     -n INT       minimal clipping reads supporting a misassembly [7]" << endl;
	cout << "     -c INT       maximal clipping region size to detect [10000]" << endl;
	cout << "     -r FILE      limit reference regions to process [null]: CHR|CHR:START-END" << endl;
	cout << "     -o DIR       output directory [output]" << endl;
	cout << "     -p STR       prefix of output result files [null]" << endl;
//	cout << "     -t INT       number of concurrent work [1]" << endl;
//	cout << "     -M INT       Mask mis-aligned regions [0]: 1 for yes, 0 for no" << endl; //20220429
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for assemble command
void Paras::showAssembleUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr assemble [options] <REF_FILE> <BAM_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "     REF_FILE     Reference file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted BAM file (required)" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -S INT       assemble slide size [10000]" << endl;
	cout << "     -c INT       maximal clipping region size [10000]" << endl;
	cout << "     -x FLOAT     expected sampling coverage for local assemble [" << EXPECTED_COV_ASSEMBLE << "], " << endl;
	cout << "                  0 for no coverage sampling" << endl;
	cout << "     -o DIR       output directory [output]" << endl;
	cout << "     -p STR       prefix of output result files [null]" << endl;
	cout << "     -t INT       number of concurrent work [1]" << endl;
	cout << "     -T INT       limited number of threads for each assemble work [0]:" << endl;
	cout << "                  0 for unlimited, and positive INT for the limited" << endl;
	cout << "                  number of threads for each assemble work" << endl;
	cout << "     -M INT       Mask mis-aligned regions [0]: 1 for yes, 0 for no" << endl;
	cout << "     -R INT       Delete temporary reads during local assembly [1]:" << endl;
	cout << "                  1 for yes, 0 for no" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for call command
void Paras::showCallUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr call [options] <REF_FILE> <BAM_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "     REF_FILE     Reference file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted BAM file (required)" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -S INT       assemble slide size used in 'assemble' command [10000]" << endl;
	cout << "     -c INT       maximal clipping region size [10000]" << endl;
	cout << "     -o DIR       output directory [output]" << endl;
	cout << "     -p STR       prefix of output result files [null]" << endl;
	cout << "     -t INT       number of concurrent work [1]" << endl;
	cout << "     -M INT       Mask mis-aligned regions [0]: 1 for yes, 0 for no" << endl;
	cout << "     -R INT       Delete temporary reads during local assembly [1]:" << endl;
	cout << "                  1 for yes, 0 for no" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show the usage for all command
void Paras::showAllUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;
	cout << "Usage: asvclr all [options] <REF_FILE> <BAM_FILE> [Region ...]" << endl << endl;

	cout << "Description:" << endl;
	cout << "     REF_FILE     Reference file (required)" << endl;
	cout << "     BAM_FILE     Coordinate-sorted file (required)" << endl;
	cout << "     Region       Limit reference region to process: CHR|CHR:START-END." << endl;
	cout << "                  If unspecified, all reference regions will be " << endl;
	cout << "                  processed (optional)" << endl << endl;

	cout << "Options: " << endl;
	cout << "     -b INT       block size [1000000]" << endl;
	cout << "     -s INT       detect slide size [500]" << endl;
	cout << "     -S INT       assemble slide size [10000]" << endl;
	cout << "     -m INT       minimal SV size to detect [2]" << endl;
	cout << "     -n INT       minimal clipping reads supporting a SV [7]" << endl;
	cout << "     -c INT       maximal clipping region size [10000]" << endl;
	cout << "     -x FLOAT     expected sampling coverage for local assemble [" << EXPECTED_COV_ASSEMBLE << "], " << endl;
	cout << "                  0 for no coverage sampling" << endl;
	cout << "     -o DIR       output directory [output]" << endl;
	cout << "     -p STR       prefix of output result files [null]" << endl;
	cout << "     -t INT       number of concurrent work [1]" << endl;
	cout << "     -T INT       limited number of threads for each assemble work [0]:" << endl;
	cout << "                  0 for unlimited, and positive INT for the limited" << endl;
	cout << "                  number of threads for each assemble work" << endl;
	cout << "     -M INT       Mask mis-aligned regions [0]: 1 for yes, 0 for no" << endl;
	cout << "     -R INT       Delete temporary reads during local assembly [1]:" << endl;
	cout << "                  1 for yes, 0 for no" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// output parameters
void Paras::outputParas(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;

	if(refFile.size()) cout << "Assemble file: " << refFile << endl;
	if(inBamFile.size()) cout << "Align file: " << inBamFile << endl;
	if(outDir.size()) cout << "Output directory: " << outDir << endl;
	if(outFilePrefix.size()) cout << "Output result file prefix: " << outFilePrefix << endl;

	cout << "Clipping number supporting misassembly: " << minClipReadsNumSupportSV << endl;
	cout << "Block size: " << blockSize << " bp" << endl;
	cout << "Slide size: " << slideSize << " bp" << endl;
	cout << "Subblock size: " << subblockSize << " bp" << endl;	//20220502
	cout << "Gap extend region size: " << gapextendSize << " bp" << endl;	//20220502
	cout << "Minimum proportion of abnormal pairing reads in misassemblies: " << abpairRatio << endl;	//20220502
	cout << "Minimum proportion of high-clipping reads in misassemblies: " << highclipRatio << endl;	//20220502
	cout << "Minimum proportion of indel reads in misassemblies: " << indelRatio << endl;	//20220502
	cout << "Maximal clipping region size: " << maxClipRegSize << " bp" << endl << endl;
//	cout << "Expected sampling coverage: " << expected_cov_assemble << endl;	//20220504
//	cout << "Number of concurrent works: " << num_threads << endl << endl;


//	20220429
//	cout << "Limited number of threads for each assemble work: " << num_threads_per_assem_work << endl;
//	if(maskMisAlnRegFlag) cout << "Mask mis-aligned regions: yes" << endl;
//	else cout << "Mask mis-aligned regions: no" << endl;
//	if(delete_reads_flag) cout << "Delete local temporary reads: yes" << endl << endl;
//	else cout << "Delete local temporary reads: no" << endl << endl;

	string desc = "Limit regions to process:";
	printLimitRegs(limit_reg_vec, desc);
}

// output the estimation parameters
void Paras::outputEstParas(string info){
	cout << info << endl;
	cout << "Mean read length: " << mean_read_len << " bp" << endl;
	cout << "min_ins_size_filt: " << min_ins_size_filt << " bp" << endl;
	cout << "min_del_size_filt: " << min_del_size_filt << " bp" << endl;
//	cout << "min_clip_size_filt: " << min_clip_size_filt << " bp" << endl;
	cout << "min_ins_num_filt: " << min_ins_num_filt << endl;
	cout << "min_del_num_filt: " << min_del_num_filt << endl;
	cout << "min_clip_num_filt: " << min_clip_num_filt << endl;
//	cout << "large_indel_size_thres: " << large_indel_size_thres << " bp" << endl << endl;
	cout << "large_indel_size_thres: " << large_indel_size_thres << " bp" << endl;
}

// initialize the estimation auxiliary data
void Paras::initEst(){
	size_t i;
	for(i=0; i<AUX_ARR_SIZE; i++) insSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) delSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) clipSizeEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) insNumEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) delNumEstArr[i] = 0;
	for(i=0; i<AUX_ARR_SIZE; i++) clipNumEstArr[i] = 0;
}

// estimate parameters
void Paras::estimate(size_t op_est){
	if(op_est==SIZE_EST_OP){
		// min_ins_size_filt
//		cout << "min_ins_size_filt:" << endl;
		min_ins_size_filt = estimateSinglePara(insSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		// min_del_size_filt
//		cout << "min_del_size_filt:" << endl;
		min_del_size_filt = estimateSinglePara(delSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		// min_clip_size_filt
//		min_clip_size_filt = estimateSinglePara(clipSizeEstArr, AUX_ARR_SIZE, SIZE_PERCENTILE_EST, MIN_INDEL_EVENT_SIZE);

		if(min_ins_size_filt>min_del_size_filt) large_indel_size_thres = min_ins_size_filt * LARGE_INDEL_SIZE_FACTOR;
		else large_indel_size_thres = min_del_size_filt * LARGE_INDEL_SIZE_FACTOR;

		// compute mean read length
		if(total_read_num_est>0) mean_read_len = round((double)mean_read_len / total_read_num_est);
		else mean_read_len = total_read_num_est = 0;

	}else if(op_est==NUM_EST_OP){
		// min_ins_num_filt
//		cout << "min_ins_num_filt:" << endl;
		min_ins_num_filt = estimateSinglePara(insNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);

		// min_del_num_filt
//		cout << "min_del_num_filt:" << endl;
		min_del_num_filt = estimateSinglePara(delNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);

		// min_clip_num_filt
//		cout << "min_clip_num_filt:" << endl;
		min_clip_num_filt = estimateSinglePara(clipNumEstArr, AUX_ARR_SIZE, NUM_PERCENTILE_EST, MIN_INDEL_EVENT_NUM);
	}else if(op_est==SNV_EST_OP){

	}else{
		cerr << __func__ << ", line=" << __LINE__ << ", invalid estimation op_flag: " << op_est << endl;
		exit(1);
	}
}

// estimate single parameter
size_t Paras::estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val){
	size_t i, total, val = 1;
	double total1;

	for(i=0, total=0; i<AUX_ARR_SIZE; i++) total += arr[i];
	if(total>0)
//		cout << "total=" << total << endl;
		for(i=0, total1=0; i<AUX_ARR_SIZE; i++){
			total1 += arr[i];
//			cout << "\t" << i << ": " << arr[i] << ", " << total1/total << endl;
			if(total1/total>=threshold){
				val = i + 1;
				break;
			}
		}
	if(val<min_val) val = min_val;

	return val;
}
