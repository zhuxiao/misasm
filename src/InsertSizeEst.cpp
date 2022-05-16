#include "InsertSizeEst.h"

InsertSizeEst::InsertSizeEst(string &inBamFile) {
	this->inBamFile = inBamFile;
	count = 0;
	isize = sdev = sum = x2 = 0;
}

InsertSizeEst::~InsertSizeEst() {
}

// get the mean insert size and standard deviation from the BAM
void InsertSizeEst::estInsertSize(double &isize_1, double &sdev_1, double &max, double &min){
	samFile *in = NULL;
	bam_hdr_t *header;
	string reg_str_bam, chrname;
	double s;
	int64_t chrlen;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
		exit(1);
	}

	if ((header = sam_hdr_read(in)) == 0) {
		cerr << __func__ << ": fail to read the header from " << inBamFile << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == 0) { // index is unavailable
		cerr << __func__ << ": random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}

	for(int32_t i=0; i<header->n_targets; i++){
		chrname = getTargetName(header, i); // get the reference name
		chrlen = getTargetLen(header, i); // get the reference length
		reg_str_bam = chrname + ":1-" + to_string(chrlen);

		if(chrlen>MIN_CHR_LEN){

			hts_itr_t *iter = sam_itr_querys(idx, header, reg_str_bam.c_str()); // parse a region in the format like `chr2:100-200'
			if (iter == NULL) { // region invalid or reference name not found
				int beg, end;
				if (hts_parse_reg(reg_str_bam.c_str(), &beg, &end))
					cerr <<  __func__ << ": region " << reg_str_bam << " specifies an unknown reference name." << endl;
				else
					cerr <<  __func__ << ": region " << reg_str_bam << " could not be parsed." << endl;
				exit(1);
			}

			estInsertSizeFromIter(in, header, iter, reg_str_bam);
			hts_itr_destroy(iter);
			if(count>=MAX_PAIR_NUM) break;
		}
	}
	hts_idx_destroy(idx); // destroy the BAM index
	bam_hdr_destroy(header); // destroy the BAM header
	if(in) sam_close(in);

	isize = (double) sum / count; // the mean insert size
	s = x2 + count*isize*isize - 2*isize*sum;
	sdev = (double) sqrt(s / (count-1)); // the standard deviation
	isize_1 = isize;
	sdev_1 = sdev;
	max = isize + 3*sdev;
	min= isize - 3*sdev;

//	cout << "isize:" << isize << ", sdev:" << sdev << endl;
}

void InsertSizeEst::estInsertSizeFromIter(samFile *in, bam_hdr_t *header, hts_itr_t *iter, string &reg){
	int result;
	bam1_t *b;

	b = bam_init1();
	while ((result = sam_itr_next(in, iter, b)) >= 0) {

		if(b->core.tid==b->core.mtid){ // the query and the query's mate are on the reference
			if(b->core.pos < b->core.mpos and !bam_is_rev(b) and bam_is_mrev(b) and ((b)->core.flag&BAM_FSUPPLEMENTARY)==0){ // the query pos less the query's mate pos, the query is not on the reverse strand,the query's mate is on the reverse strand,flag<2048
				if(b->core.isize > 0 and b->core.isize < MAX_ISIZE ){
					sum += b->core.isize; // get the sum of insert size from BAM
					x2 += b->core.isize*b->core.isize;
					count++; // get the paired reads count
					if(count>=MAX_PAIR_NUM) break;
				}
			}
		}

		bam_destroy1(b);
		b = bam_init1();
	}

	bam_destroy1(b);
	if (result < -1) {
		cerr <<  __func__ << ": retrieval of region " << reg << " failed due to truncated file or corrupt BAM index file." << endl;
		exit(1);
	}
}
