#include <algorithm>
#include "util.h"
#include "clipReg.h"
#include "Region.h"

//Constructor
Region::Region(string& chrname, int64_t startRpos, int64_t endRPos, int64_t chrlen, int64_t minRPos, int64_t maxRPos, Base *regBaseArr, size_t regFlag, Paras *paras) {
	this->paras = paras;
	this->chrname = chrname;
	this->startRPos = startRpos;
	this->endRPos = endRPos;
	this->chrlen = chrlen;
	this->minRPos = minRPos;
	this->maxRPos = maxRPos;
	this->regBaseArr = regBaseArr;
	this->regFlag = regFlag;

	switch(regFlag){
		case HEAD_REGION:
			startMidPartPos = startRpos;
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) endMidPartPos = endRPos;
			break;
		case INNER_REGION:
			startMidPartPos = startRpos + paras->slideSize;
			if(startMidPartPos>endRPos) {
				cerr << __func__ << "(), line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) {
				cerr << __func__ << "(), line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			break;
		case TAIL_REGION:
			startMidPartPos = startRpos + paras->slideSize;
//			if(startMidPartPos>endRPos) {
			if(startMidPartPos>endRPos+SLIDESIZE) {	//20220429
//			if(startMidPartPos>endRPos+paras->slideSize) {	//20220504
				cerr << __func__ << "(), line=" << __LINE__ << ": error!" << endl; exit(1);
			}
			endMidPartPos = startMidPartPos + paras->slideSize - 1;
			if(endMidPartPos>endRPos) endMidPartPos = endRPos;
			break;
	}

	subRegSize = SUB_REG_SIZE;
	meanBlockCov = 0;

	// determine if all the base in the region are 'N' bases in reference
	wholeRefGapFlag = IsWholeRefGap();

	// compute the local region coverage
	localRegCov = computeMeanCovReg(startMidPartPos, endMidPartPos);
	highIndelSubRegNum = 0;

	indelCandFlag = false;
	difType = REG_DIF_NONE;

	fragmentsize = ratio_exisizemax = ratio_lessisizemin = ratio_abdirection = ratio_abref = 0;

	abpairstartpos = abpairendpos = -1;	//20220424

	misType = "";	//20220504
}

//Destructor
Region::~Region(){
//	fai_destroy(fai);
	if(!disagrePosVector.empty()) destroyDisagrePosVector();
	if(!zeroCovPosVector.empty()) destroyZeroCovPosVector();
	if(!abCovRegVector.empty()) destroyAbCovRegVector();
	if(!snvVector.empty()) destroySnvVector();
	if(!indelVector.empty()) destroyIndelVector();
	if(!clipRegVector.empty()) destroyClipRegVector();
}

// set the output directory
void Region::setOutputDir(string& out_dir_assemble_prefix){
	out_dir_assemble = out_dir_assemble_prefix + "/" + chrname;
}

// determine if all the base in the region are 'N' bases in reference
bool Region::IsWholeRefGap(){
	bool flag = true;
	for(int64_t i=startMidPartPos; i<=endMidPartPos; i++)
		if(regBaseArr[i-startRPos].coverage.idx_RefBase!=4){ // excluding 'N'
			flag = false;
			break;
		}
	return flag;
}

// set the block mean coverage
void Region::setMeanBlockCov(double meanBlockCov){
	this->meanBlockCov = meanBlockCov;
}

// compute abnormal signatures in region
int Region::computeRegAbSigs(){
	// compute disagreements
	computeDisagreements();

	// compute abnormal coverage, including high/low coverage
	computeAbCovReg();

	// compute the number sub-regions of high indel events in the region
	computeHighIndelEventRegNum();

	return 0;
}

// compute disagreements, only check the middle part sub-region
int Region::computeDisagreements(){
	int64_t pos, regIdx;
	for(pos=startMidPartPos; pos<endMidPartPos; pos++){
		regIdx = pos - startRPos;
		if(regBaseArr[regIdx].coverage.idx_RefBase!=4){ // A, C, G, T, but N
			if(regBaseArr[regIdx].isDisagreeBase()) addDisagrePos(pos);
			if(regBaseArr[regIdx].isZeroCovBase()) addZeroCovPos(pos);
		}
	}
	disagrePosVector.shrink_to_fit();
	zeroCovPosVector.shrink_to_fit();
	return 0;
}

// add base position having disagreement to vector
void Region::addDisagrePos(int64_t  pos){
	disagrePosVector.push_back(pos);
}

// add base position having zero coverage to vector
void Region::addZeroCovPos(int64_t pos){
	zeroCovPosVector.push_back(pos);
}

// destroy the disagreements vector
void Region::destroyDisagrePosVector(){
	vector<int64_t>().swap(disagrePosVector);
}

// destroy the zero coverage vector
void Region::destroyZeroCovPosVector(){
	vector<int64_t>().swap(zeroCovPosVector);
}

// destroy the abnormal coverage vector
void Region::destroyAbCovRegVector(){
	vector<reg_t*>::iterator it;
	for(it=abCovRegVector.begin(); it!=abCovRegVector.end(); it++)
		delete (*it);
	vector<reg_t*>().swap(abCovRegVector);
}

// destroy the SNV vector
void Region::destroySnvVector(){
	vector<int64_t>().swap(snvVector);
}

// destroy the indel vector
void Region::destroyIndelVector(){
	vector<reg_t*>().swap(indelVector);
}

// destroy the clip region vector
void Region::destroyClipRegVector(){
	vector<reg_t*>().swap(clipRegVector);
}

// output abnormal signatures in region
void Region::printRegAbSigs(){
	vector<int64_t>::iterator it;

	// output SNVs
	cout << "region: " << startMidPartPos << "-" << endMidPartPos << endl;
	cout << "    Number of SNVs: " << getSNVNum() << endl;
	for(it=snvVector.begin(); it!=snvVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output Indels
	cout << "    Number of Indels: " << getIndelNum() << endl;
	vector<reg_t*>::iterator reg;
	for(reg=indelVector.begin(); reg!=indelVector.end(); reg++)
		cout << "        " << (*reg)->startRefPos << "-" << (*reg)->endRefPos << endl;

	// output position of disagreements
	cout << "    Number of disagreements: " << getDisagreeNum() << endl;
	for(it=disagrePosVector.begin(); it!=disagrePosVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output position of zero coverage
	cout << "    Number of zero coverage: " << zeroCovPosVector.size() << endl;
	for(it=zeroCovPosVector.begin(); it!=zeroCovPosVector.end(); it++)
		cout << "        " << (*it) << endl;

	// output sub-regions of abnormal coverage
	cout << "    Number of sub-regions with abnormal coverage: " << abCovRegVector.size() << endl;
	for(reg=abCovRegVector.begin(); reg!=abCovRegVector.end(); reg++)
		cout << "        " << (*reg)->startRefPos << "-" << (*reg)->endRefPos << endl;

	// output the number of high indel sub-regions
	cout << "    Number of high indel sub-regions: " << highIndelSubRegNum << endl;
}

// compute the abnormal coverage region
int Region::computeAbCovReg(){
	int64_t pos, startPosSubReg, endPosSubReg;
	double tmp_cov, covRatio2block, covRatio2local;
	reg_t *reg;

	for(pos=startMidPartPos; pos<=endMidPartPos; pos+=subRegSize){
		startPosSubReg = pos;
		endPosSubReg = startPosSubReg + subRegSize - 1;
		if(endPosSubReg>endMidPartPos) endPosSubReg = endMidPartPos;
		tmp_cov = computeMeanCovReg(startPosSubReg, endPosSubReg);
		covRatio2block = tmp_cov / meanBlockCov;
		covRatio2local = tmp_cov / localRegCov;
		if((covRatio2block<MIN_COVRATIO and covRatio2local<MIN_COVRATIO) or (covRatio2block>MAX_COVRATIO and covRatio2local>MAX_COVRATIO)) {  // abnormal coverage condition
			reg = allocateReg(chrname, startPosSubReg, endPosSubReg);
			abCovRegVector.push_back(reg);
		}
	}
	abCovRegVector.shrink_to_fit();

	return 0;
}

// allocate reg
reg_t* Region::allocateReg(string &chrname, int64_t startPosReg, int64_t endPosReg){
	reg_t *reg = new reg_t();
	if(!reg){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	reg->chrname = chrname;
	reg->startRefPos = startPosReg;
	reg->endRefPos = endPosReg;
	reg->startLocalRefPos = reg->endLocalRefPos = 0;
	reg->startQueryPos = reg->endQueryPos = 0;
	reg->var_type = VAR_UNC;
	reg->sv_len = 0;
	reg->query_id = -1;
	reg->blat_aln_id = -1;
	reg->call_success_status = false;
	reg->short_sv_flag = false;
	reg->zero_cov_flag = false;
	reg->aln_seg_end_flag = false;
	return reg;
}

// allocate reg for misasm
reg_t* Region::allocateMisReg(string &chrname, int64_t startPosReg, int64_t endPosReg, string misType){
	reg_t *reg = new reg_t();
	if(!reg){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot allocate memory" << endl;
		exit(1);
	}
	reg->chrname = chrname;
	reg->startRefPos = startPosReg;
	reg->endRefPos = endPosReg;
	reg->startLocalRefPos = reg->endLocalRefPos = 0;
	reg->startQueryPos = reg->endQueryPos = 0;
	reg->var_type = VAR_UNC;
	reg->sv_len = 0;
	reg->query_id = -1;
	reg->blat_aln_id = -1;
	reg->call_success_status = false;
	reg->short_sv_flag = false;
	reg->zero_cov_flag = false;
	reg->aln_seg_end_flag = false;

	reg->misType = misType;
	return reg;
}

// compute the mean coverage of the sub-region, excluding the gap region
double Region::computeMeanCovReg(int64_t startPosReg, int64_t endPosReg){
	int64_t i, totalReadBeseNum = 0, totalRefBaseNum = 0;
	double mean_cov;
	for(i=startPosReg; i<=endPosReg; i++)
		if(regBaseArr[i-startRPos].coverage.idx_RefBase!=4){ // excluding 'N'
			totalReadBeseNum += regBaseArr[i-startRPos].coverage.num_bases[5];
			totalRefBaseNum ++;
		}
	if(totalRefBaseNum) mean_cov = (double)totalReadBeseNum/totalRefBaseNum;
	else mean_cov = 0;
	return mean_cov;
}

// compute the number of high indel events sub-regions
int Region::computeHighIndelEventRegNum(){
	int64_t pos, startPosSubReg, endPosSubReg;
	highIndelSubRegNum = 0;
	for(pos=startMidPartPos; pos<=endMidPartPos; pos+=subRegSize){
		startPosSubReg = pos;
		endPosSubReg = startPosSubReg + subRegSize - 1;
		if(endPosSubReg>endMidPartPos) endPosSubReg = endMidPartPos;
		if(computeReadIndelEventNumReg(startPosSubReg, endPosSubReg)>=5) highIndelSubRegNum ++;
	}
	return 0;
}

// compute the total indel number in reads, including insertions, deletions and clippings
// return: the total number of the above indel events
int32_t Region::computeReadIndelEventNumReg(int64_t startPosReg, int64_t endPosReg){
	int64_t i, total = 0;
	Base *base;
	for(i=startPosReg; i<=endPosReg; i++){
		base = regBaseArr + i - startRPos;
		total += base->insVector.size() + base->delVector.size() + base->clipVector.size();
	}
	return total;
}

// detect SNVs for the region
void Region::detectSNV(){
	int64_t startCheckPos, endCheckPos;
	bool SNV_flag;
	vector<int64_t>::iterator disagr;
	for(disagr=disagrePosVector.begin(); disagr!=disagrePosVector.end();){
		startCheckPos = (*disagr) - subRegSize;
		endCheckPos = (*disagr) + subRegSize;
		if(startCheckPos<startRPos) startCheckPos = startRPos;
		if(endCheckPos>endRPos) endCheckPos = endRPos;

//		if(*disagr==18486270)
//			cout << *disagr << endl;

		SNV_flag = computeSNVFlag(*disagr, startCheckPos, endCheckPos);
		if(SNV_flag){
			if(haveMuchShortIndelsAround(startCheckPos, endCheckPos))  // further check around
				indelVector.push_back(allocateReg(chrname, *disagr, *disagr));
			else{
				snvVector.push_back(*disagr); // add the SNV
				disagr = disagrePosVector.erase(disagr);
			}
		}else disagr ++;
	}
	snvVector.shrink_to_fit();
}

// compute the SNV flag for a base position
bool Region::computeSNVFlag(int64_t pos, int64_t startCheckPos, int64_t endCheckPos){
	bool SNV_flag = true;
	double cov_tmp;

	if(isInReg(pos, indelVector)) SNV_flag = false;

	if(SNV_flag){ // check the maxBase and idx_ref
		if(regBaseArr[pos-startRPos].coverage.idx_max==regBaseArr[pos-startRPos].coverage.idx_RefBase or (double)regBaseArr[pos-startRPos].coverage.num_max/regBaseArr[pos-startRPos].coverage.num_bases[5]<MIN_RATIO_SNV){
			SNV_flag = false;
		}else if(regBaseArr[pos-startRPos].coverage.idx_RefBase==5){ // mixed base
			switch(regBaseArr[pos-startRPos].coverage.refBase){
				case 'M':
				case 'm':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==1) SNV_flag = false;
					break;
				case 'R':
				case 'r':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==2) SNV_flag = false;
					break;
				case 'S':
				case 's':
					if(regBaseArr[pos-startRPos].coverage.idx_max==1 or regBaseArr[pos-startRPos].coverage.idx_max==2) SNV_flag = false;
					break;
				case 'V':
				case 'v':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==1 or regBaseArr[pos-startRPos].coverage.idx_max==2) SNV_flag = false;
					break;
				case 'W':
				case 'w':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				case 'Y':
				case 'y':
					if(regBaseArr[pos-startRPos].coverage.idx_max==1 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				case 'H':
				case 'h':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==1 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				case 'K':
				case 'k':
					if(regBaseArr[pos-startRPos].coverage.idx_max==2 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				case 'D':
				case 'd':
					if(regBaseArr[pos-startRPos].coverage.idx_max==0 or regBaseArr[pos-startRPos].coverage.idx_max==2 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				case 'B':
				case 'b':
					if(regBaseArr[pos-startRPos].coverage.idx_max==1 or regBaseArr[pos-startRPos].coverage.idx_max==2 or regBaseArr[pos-startRPos].coverage.idx_max==3) SNV_flag = false;
					break;
				default: cerr << __func__ << ": unknown base: " << regBaseArr[pos-startRPos].coverage.refBase << endl; exit(1);
			}
		}
	}

	// check the coverage
	if(SNV_flag){
		cov_tmp = computeMeanCovReg(startCheckPos, endCheckPos);
		if((cov_tmp/meanBlockCov<0.6 and cov_tmp/localRegCov<0.6) or (cov_tmp/meanBlockCov>2 and cov_tmp/localRegCov>2))
			SNV_flag = false;
	}
	return SNV_flag;
}

// determine whether there are much short indel events around
bool Region::haveMuchShortIndelsAround(int64_t startCheckPos, int64_t endCheckPos){
	bool flag = false;
	Base *base;

	for(int64_t i=startCheckPos; i<=endCheckPos; i++){
		base = regBaseArr + i - startRPos;
		//if(base->insVector.size()>=paras->min_ins_num_filt or base->delVector.size()>=paras->min_del_num_filt or base->clipVector.size()>=paras->min_clip_num_filt
		if(base->insVector.size()>=paras->min_ins_num_filt or base->del_num_from_del_vec>=(int32_t)paras->min_del_num_filt or base->clipVector.size()>=paras->min_clip_num_filt
			or base->num_shortIns>=(int32_t)paras->min_ins_num_filt or base->num_shortdel>=(int32_t)paras->min_del_num_filt or base->num_shortClip>=(int32_t)paras->min_clip_num_filt
			/*or base->getLargerInsNum(paras->min_ins_size_filt)>0 or base->getLargerDelNum(paras->min_del_size_filt)>0 or base->getLargerClipNum(paras->min_clip_size_filt)>0*/){
			flag = true;
			break;
		}
	}

	return flag;
}

// get the number of disagreements excluding the SNV items
size_t Region::getDisagreeNum(){
	return disagrePosVector.size();
}

// get the number of SNVs
size_t Region::getSNVNum(){
	return snvVector.size();
}

// get the number of Indels
size_t Region::getIndelNum(){
	return indelVector.size();
}

// detect the indel regions
void Region::detectIndelReg(){
	reg_t *reg = NULL;
	int64_t i = startMidPartPos - subRegSize;
	if(i<minRPos) i = minRPos;
	while(i<endMidPartPos){
//		if(i>7163000)  //109500, 5851000, 11812500, 12601500, 14319500, 14868000, 18343500, <18710000>
//			cout << i << endl;

		reg = getIndelReg(i);
		if(reg){
			//cout << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", localRef: " << reg->startLocalRefPos << "-" << reg->endLocalRefPos << ", Query :" << reg->startQueryPos << "-" << reg->endQueryPos << endl;
			indelVector.push_back(reg);
			i = reg->endRefPos + subRegSize;
		} else break;
	}
	indelVector.shrink_to_fit();
}

// get the indel region given the start checking chromosome position
reg_t* Region::getIndelReg(int64_t startCheckPos){
	reg_t *reg = NULL;
	int32_t reg_size1, reg_size2, num1, num3, num4, extendSize, high_con_indel_base_num;
	vector<double> num_vec;
	double high_indel_clip_ratio;
	int64_t i, checkPos, startPos1, endPos1, startPos2;
	int64_t startPos_indel = -1, endPos_indel = -1;
	bool indel_beg_flag, indel_end_flag;

	checkPos = startCheckPos;
	while(checkPos<=endMidPartPos){
		startPos_indel = -1;
		endPos_indel = -1;
		indel_beg_flag = false;
		indel_end_flag = false;

		// get normal region1
		startPos1 = -1;
		reg_size1 = 0;
		for(i=checkPos; i<=endMidPartPos; i++){
			if(regBaseArr[i-startRPos].coverage.idx_RefBase==4){ // skip the Ns region
				startPos1 = -1;
				reg_size1 = 0;
				continue;
			}

			if(haveNoAbSigs(regBaseArr+i-startRPos, i)){
				if(startPos1==-1) startPos1 = i;
				reg_size1 ++;
			}else{
				if(reg_size1>=(int64_t)subRegSize) {
					indel_beg_flag = true;
					break;
				}else{
					startPos1 = -1;
					reg_size1 = 0;
				}
			}
		}
		if(indel_beg_flag){
			endPos1 = startPos1 + reg_size1 - 1;

			// get normal region2
			startPos2 = -1;
			reg_size2 = 0;
			extendSize = 0;
			for(i=endPos1+1; i<=endRPos+extendSize; i++){
				if(i>chrlen) break;

				if(regBaseArr[i-startRPos].coverage.idx_RefBase==4){ // skip the Ns region
					startPos2 = -1;
					reg_size2 = 0;
					continue;
				}

				if(haveNoAbSigs(regBaseArr+i-startRPos, i)){
					if(startPos2==-1)
						startPos2 = i;
					++reg_size2;
					if(reg_size2>=2*(int64_t)subRegSize){
						num3 = getLargeIndelNum(startPos2, i);
						if(num3<3){
							indel_end_flag = true;
							break;
						}else{
							startPos2 = -1;
							reg_size2 = 0;
						}
					}
				}else{
					startPos2 = -1;
					reg_size2 = 0;
				}

				if(i==endRPos and indel_end_flag==false){
					extendSize = paras->slideSize * 2;
					if(extendSize+endRPos>maxRPos)
						extendSize = maxRPos - endRPos;
				}
			}

			if(indel_end_flag){
				startPos_indel = endPos1 + 1;
				endPos_indel = startPos2 - 1;
			}else if(extendSize>0){
				indel_end_flag = true;
				startPos_indel = endPos1 + 1;
				endPos_indel = endRPos + extendSize - 1;
			}
		}
		// allocate the indel region
		if(indel_beg_flag and indel_end_flag){
			high_con_indel_base_num = getHighConIndelNum(startPos_indel, endPos_indel, MIN_HIGH_INDEL_BASE_RATIO, IGNORE_POLYMER_RATIO_THRES);
			if(endPos_indel-startPos_indel+1>=(int32_t)paras->min_sv_size_usr or high_con_indel_base_num>0){
				num1 = getDisZeroCovNum(startPos_indel, endPos_indel);
				//num2 = getMismatchBasesAround(startPos_indel, endPos_indel);
				num3 = getLargeIndelBaseNum(startPos_indel, endPos_indel);
				num_vec = getTotalHighIndelClipRatioBaseNum(regBaseArr+startPos_indel-startRPos, endPos_indel-startPos_indel+1);
				num4 = num_vec.at(0);
				high_indel_clip_ratio = num_vec.at(1);
				//if(num1>0 or num2>=DISAGREE_NUM_THRES_REG or num3>0) {
				if(num1>0 or num3>0 or num4>0 or high_indel_clip_ratio>=HIGH_INDEL_CLIP_BASE_RATIO_THRES) {
					reg = allocateReg(chrname, startPos_indel, endPos_indel);
					break;
				}else checkPos = endPos_indel + 1;
			}else checkPos = endPos_indel + 1;
		}else break;
	}

	return reg;
}

// detect the Illumina indel regions
void Region::detectIlluminaIndelReg(){
	reg_t *reg = NULL;
	int64_t i = startMidPartPos;

	while(i<endMidPartPos){
		reg = getIlluminaIndelReg(i);
		if(reg){
			indelVector.push_back(reg);
			i = reg->endRefPos + 1;
		} else break;
	}
	indelVector.shrink_to_fit();
}

// get the Illumina indel region given the start checking chromosome position
reg_t* Region::getIlluminaIndelReg(int64_t startCheckPos){
	reg_t *reg = NULL;
	int64_t i, checkPos, startPos_indel, endPos_indel;
	string misType_indel;

	checkPos = startCheckPos;
	startPos_indel = endPos_indel = -1;
	misType_indel = "";

	while(checkPos<=endMidPartPos){
		if(regBaseArr[checkPos-startRPos].isIlluminaHighInsBase(paras->indelRatio)){
			startPos_indel = checkPos;
			for(i=startPos_indel+1; i<=endMidPartPos; i++){
				if(regBaseArr[i-startRPos].isIlluminaHighInsBase(paras->indelRatio)){
					continue;
				}else{
					endPos_indel = i-1;
					misType_indel = "Insertion";
					reg = allocateMisReg(chrname, startPos_indel, endPos_indel, misType_indel);
					break;
				}
			}

			break;

		}

		else if(regBaseArr[checkPos-startRPos].isIlluminaHighDelBase(paras->indelRatio)){
			startPos_indel = checkPos;
			for(i=startPos_indel+1; i<=endMidPartPos; i++){
				if(regBaseArr[i-startRPos].isIlluminaHighDelBase(paras->indelRatio)){
					continue;
				}else{
					endPos_indel = i-1;
					misType_indel = "Deletion";
					reg = allocateMisReg(chrname, startPos_indel, endPos_indel, misType_indel);
					break;
				}
			}

			break;

		}else{
			checkPos++;
			startPos_indel = endPos_indel = -1;
		}

	}
	return reg;
}

// determine whether the base have no abnormal signatures
bool Region::haveNoAbSigs(Base *base, int64_t pos){
	if(base->isDisagreeBase())
		if(find(snvVector.begin(), snvVector.end(), pos)==snvVector.end()) return false;
	if(base->isZeroCovBase() or base->insVector.size()>=paras->min_ins_num_filt or base->del_num_from_del_vec>=(int32_t)paras->min_del_num_filt or base->clipVector.size()>=paras->min_clip_num_filt)
		return false;
	else if(base->getLargeIndelNum(paras->large_indel_size_thres)>=3 or base->getLargeIndelNum(paras->large_indel_size_thres*2)>=2 or (double)(base->getTotalIndelNum()+base->getTotalClipNum())/base->getTotalCovNum()>=HIGH_INDEL_CLIP_RATIO_THRES)
		return false;
	return true;
}

// check [-2, 2] region around
int32_t Region::getMismatchBasesAround(int64_t pos1, int64_t pos2){
	int64_t i, num, startPos, endPos;
	Base *base;

	startPos = pos1;
	if(startPos<startMidPartPos) startPos = startMidPartPos;
	endPos = pos2;
	if(endPos>endMidPartPos) endPos = endMidPartPos;

	for(num=0, i=startPos; i<=endPos; i++){
		base = regBaseArr + i - startRPos;
		if(base->coverage.idx_max!=base->coverage.idx_RefBase or (double)base->coverage.num_max/base->coverage.num_bases[5]<=DISAGREE_THRES_REG)
			num ++;
	}
	return num;
}

int32_t Region::getDisZeroCovNum(int64_t startPos, int64_t endPos){
	int64_t i, total = 0;
	for(i=startPos; i<=endPos; i++)
		if(regBaseArr[i-startRPos].isDisagreeBase() or regBaseArr[i-startRPos].isZeroCovBase())
			total ++;
	return total;
}

// get the number of bases with long indels
int32_t Region::getLargeIndelBaseNum(int64_t startPos, int64_t endPos){
	int64_t i, large_indel_num, total;
	double ratio;

	total = 0;
	for(i=startPos; i<=endPos; i++){
		large_indel_num = regBaseArr[i-startRPos].getLargeIndelNum(paras->large_indel_size_thres);
		ratio = (double)large_indel_num / regBaseArr[i-startRPos].coverage.num_bases[5];
		if(ratio>=LARGE_INDEL_RATIO_THRES)
			total ++;
		else{
			large_indel_num = regBaseArr[i-startRPos].getLargeIndelNum(paras->large_indel_size_thres*2);
			ratio = (double)large_indel_num / regBaseArr[i-startRPos].coverage.num_bases[5];
			if(ratio>=0.5*LARGE_INDEL_RATIO_THRES)
				total ++;
		}
	}

	return total;
}

// get the number of bases with long indels
int32_t Region::getLargeIndelNum(int64_t startPos, int64_t endPos){
	int64_t i, large_indel_num;
	large_indel_num = 0;
	for(i=startPos; i<=endPos; i++){
		large_indel_num += regBaseArr[i-startRPos].getLargeIndelNum(paras->large_indel_size_thres);
	}
	return large_indel_num;
}

// get the number of bases with high consensus indels
int32_t Region::getHighConIndelNum(int64_t startPos, int64_t endPos, float threshold, float polymer_ignore_ratio_thres){
	int64_t i, high_con_indel_base_num = 0;
	bool flag;
	for(i=startPos; i<=endPos; i++){
		flag = regBaseArr[i-startRPos].isHighConIndelBase(threshold, polymer_ignore_ratio_thres);
		if(flag) high_con_indel_base_num ++;
	}
	return high_con_indel_base_num;
}

// get the number of bases with Illumina high consensus indels
int32_t Region::getIlluminaHighConIndelNum(int64_t startPos, int64_t endPos, float threshold){
	int64_t i, high_con_indel_base_num = 0;
	bool flag;
	for(i=startPos; i<=endPos; i++){
		flag = regBaseArr[i-startRPos].isIlluminaHighConIndelBase(threshold);
		if(flag) high_con_indel_base_num ++;
	}
	return high_con_indel_base_num;
}

// determine the dif type for indel candidate
void Region::determineDifType(){
	if(getDisagreeNum()>0 or zeroCovPosVector.size()>0 or abCovRegVector.size()>0 or highIndelSubRegNum>0)
		indelCandFlag = true;
}

// get SNV vector
vector<int64_t> Region::getSnvVector(){
	return snvVector;
}

// get indel vector
vector<reg_t*> Region::getIndelVector(){
	return indelVector;
}

// get indel vector
vector<reg_t*> Region::getClipRegVector(){
	return clipRegVector;
}

// get indel vector
vector<int64_t> Region::getZeroCovPosVector(){
	return zeroCovPosVector;
}

//
void Region::detectHighClipReg(){
	int64_t i = startMidPartPos - SUB_CLIP_REG_SIZE;
	if(i<1) i = 1;
	while(i<endMidPartPos){
//		if(i>2932500)
//			cout << i << endl;

		reg_t *reg = getClipReg(i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + SUB_CLIP_REG_SIZE;
		} else break;
	}
	clipRegVector.shrink_to_fit();
}

// get the clip region
reg_t* Region::getClipReg(int64_t startCheckPos){
	reg_t *reg = NULL;
	int64_t checkPos, subclipreg_size, startPos_tmp, endPos_tmp;
	int64_t startPos_clip, endPos_clip;
	bool normal_reg_flag, clip_reg_flag;

	//cout << chrname << ":" << to_string(startRPos) << "-" << to_string(endRPos) << ", localRegCov=" << localRegCov << endl;

	subclipreg_size = SUB_CLIP_REG_SIZE;
	clip_reg_flag = false;

	startPos_clip = endPos_clip = -1;
	checkPos = startCheckPos;
	while(checkPos<=endMidPartPos){
		startPos_tmp = checkPos;
		endPos_tmp = checkPos + subclipreg_size - 1;
		if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;
		normal_reg_flag = haveNoClipSig(startPos_tmp, endPos_tmp, HIGH_CLIP_RATIO_THRES);
		if(normal_reg_flag==false){ // clip region
			if(clip_reg_flag==false){ // first clip region
				clip_reg_flag = true;
				startPos_clip = startPos_tmp;
			}
			endPos_clip = endPos_tmp;
		}else { // normal region
			if(clip_reg_flag==true){
				break;
			}
		}

		if(startPos_clip!=-1 and endPos_clip!=-1){ // check disagreements
			// compute the number of disagreements
			int32_t disagreeNum = computeDisagreeNum(regBaseArr+startPos_clip-startRPos, endPos_clip-startPos_clip+1);
			normal_reg_flag = haveNoClipSig(startPos_clip, endPos_clip, HIGH_CLIP_RATIO_THRES*3);

			if(disagreeNum>=1 or normal_reg_flag==false) {
				reg = allocateReg(chrname, startPos_clip, endPos_clip);
				break;
			}else{
				clip_reg_flag = false;
				startPos_clip = endPos_clip = -1;
			}
		}

		checkPos = endPos_tmp + 1;
	}

	return reg;
}

bool Region::haveNoClipSig(int64_t startPos, int64_t endPos, double clip_ratio_thres){
	bool flag = true;
	int64_t i;
	size_t j, clip_num;
	vector<clipEvent_t*> clip_vec;
	double ratio;

	clip_num = 0;
	for(i=startPos; i<=endPos; i++){
		clip_vec = regBaseArr[i-startRPos].clipVector;
		for(j=0; j<clip_vec.size(); j++){
			if(stoi(clip_vec.at(j)->seq)>=MIN_CLIP_END_SIZE)
				clip_num ++;
		}
	}

	if(localRegCov>0){
		ratio = clip_num / localRegCov;
		if(ratio>=clip_ratio_thres) flag = false;
	}

	if(flag){
		for(i=startPos; i<=endPos; i++)
			if((double)regBaseArr[i-startRPos].clipVector.size()/regBaseArr[i-startRPos].coverage.num_bases[5]>=clip_ratio_thres){
				flag = false;
				break;
			}
	}

	return flag;
}

void Region::detectIlluminaHighClipReg(vector<bam1_t*> &slideAlnVector){
	int64_t i = startMidPartPos;
	if(i<1) i = 1;
	while(i<endMidPartPos){

		reg_t *reg = getIlluminaClipReg(slideAlnVector, i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + 1;
		} else break;
	}
	clipRegVector.shrink_to_fit();
}

// get the Illumina clip region
reg_t* Region::getIlluminaClipReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){

	reg_t *reg = NULL;
	int64_t checkPos, subclipreg_size, startPos_tmp, endPos_tmp;
	int64_t startPos_clip, endPos_clip;
	bool gap_flag, regclip_flag, baseclip_flag, abisize_flag, gapclipreg_flag;

	subclipreg_size = SUB_CLIP_REG_SIZE;
	gap_flag = regclip_flag = baseclip_flag = abisize_flag = gapclipreg_flag = false;
	startPos_clip = endPos_clip = -1;
	checkPos = startCheckPos;
	while(checkPos<=endMidPartPos){
		startPos_tmp = checkPos;
		endPos_tmp = checkPos + subclipreg_size - 1;
		if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;

		gap_flag = isIlluminaGapReg(startPos_tmp, endPos_tmp);	//0s

		if(gap_flag == false){  //no gap
			regclip_flag = isIlluminaClipReg(startPos_tmp, endPos_tmp, ILLUMINA_HIGH_CLIP_THRES);	//0s
			baseclip_flag = isIlluminaBaseClipReg(startPos_tmp, endPos_tmp, ILLUMINA_HIGH_CLIP_BASE_THRES);	//0s
//			abisize_flag = isIlluminaAbisize(slideAlnVector, startPos_tmp, endPos_tmp);
			abisize_flag = isIlluminaAbPairing(slideAlnVector, startPos_tmp, endPos_tmp);	//time=0.1s

			if(regclip_flag==true or baseclip_flag==true or abisize_flag==true){
//			if((regclip_flag==true or baseclip_flag==true or abisize_flag==true) and startPos_tmp>=MISASSEMBLY_REGION_START_THRES and endPos_tmp<=chrlen-MISASSEMBLY_REGION_START_THRES+1){//20210625
				startPos_clip = startPos_tmp;
				endPos_clip = endPos_tmp;
				reg = allocateReg(chrname, startPos_clip, endPos_clip);
				break;
			}else startPos_clip = endPos_clip = -1;

		}
		else{  //include gap
			gapclipreg_flag = isIlluminaGapClipReg(slideAlnVector, startPos_tmp, endPos_tmp);	//time=0.18s
			if(gapclipreg_flag==true){
//			if(gapclipreg_flag==true and startPos_tmp>=MISASSEMBLY_REGION_START_THRES and endPos_tmp<=chrlen-MISASSEMBLY_REGION_START_THRES+1){//20210625
				startPos_clip = startPos_tmp;
				endPos_clip = endPos_tmp;
				reg = allocateReg(chrname, startPos_clip, endPos_clip);
				break;
			}else startPos_clip = endPos_clip = -1;
		}

		checkPos = endPos_tmp + 1;
	}
	return reg;
}

// get Misjoin
void Region::detectIlluminaMisjoinReg(vector<bam1_t*> &slideAlnVector){
	reg_t *reg = NULL;
	int64_t i = startMidPartPos;
	if(i<1) i = 1;
	while(i<endMidPartPos){

		reg = getIlluminaMisjoinReg(slideAlnVector, i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + 1;
		} else break;

		break;

	}
	clipRegVector.shrink_to_fit();
}

// get the Illumina misjoin region
reg_t* Region::getIlluminaMisjoinReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){
	reg_t *reg = NULL;
	int64_t checkPos,startPos_tmp, endPos_tmp, startPos_misjoin, endPos_misjoin, misjoinPos_tmp;
	bool gap_flag, abpair_flag;
	string misType_misjoin;

	gap_flag = abpair_flag = false;
	startPos_misjoin = endPos_misjoin = misjoinPos_tmp = -1;
	checkPos = startCheckPos;
	misType_misjoin = "";

	while(checkPos<=endMidPartPos){
		startPos_tmp = checkPos;
		endPos_tmp = checkPos + paras->subblockSize - 1;
		if(endPos_tmp>endMidPartPos) endPos_tmp = endMidPartPos;

		gap_flag = isIlluminaGapReg(startPos_tmp, endPos_tmp);

		if(gap_flag == false){  //no gap
			abpair_flag = isIlluminaAbPairing(slideAlnVector, startPos_tmp, endPos_tmp);

			if(abpair_flag==true){
				estIlluminaAbPairStartPos(slideAlnVector, startPos_tmp);
				estIlluminaAbPairEndPos(slideAlnVector, endPos_tmp);
				startPos_misjoin = abpairstartpos;
				endPos_misjoin = abpairendpos;

				if(startPos_misjoin>endPos_misjoin){
					misjoinPos_tmp = endPos_misjoin;
					endPos_misjoin = startPos_misjoin;
					startPos_misjoin = misjoinPos_tmp;
				}
				misType_misjoin = "Abnormal_pairing";
				reg = allocateMisReg(chrname, startPos_misjoin, endPos_misjoin, misType_misjoin);
				break;
			}else startPos_misjoin = endPos_misjoin = -1;

		}

		checkPos = endPos_tmp + 1;
	}
	return reg;
}

int Region::estIlluminaAbPairStartPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){
	int64_t checkPos;
	bool flag_abpairpos= false;

	abpairstartpos = -1;
	checkPos = startCheckPos;
	flag_abpairpos = isIlluminaAbPairPos(slideAlnVector, checkPos);
	if(flag_abpairpos==true){

		for(int64_t i=checkPos-1; i>=startMidPartPos; i--){
			if(!isIlluminaAbPairPos(slideAlnVector, i)){
				abpairstartpos = i+1;
				break;
			}
		}

		if(abpairstartpos==-1) abpairstartpos = startMidPartPos;

	}else{
		for(int64_t j=checkPos+1; j<=endMidPartPos; j++){
			if(isIlluminaAbPairPos(slideAlnVector, j)){
				abpairstartpos = j;
				break;
			}
		}

		if(abpairstartpos==-1) abpairstartpos = endMidPartPos;
	}
	return 0;
}

int Region::estIlluminaAbPairEndPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){
	int64_t checkPos;
	bool flag_abpairpos= false;

	abpairendpos = -1;
	checkPos = startCheckPos;
	flag_abpairpos = isIlluminaAbPairPos(slideAlnVector, checkPos);
	if(flag_abpairpos==true){

		for(int64_t i=checkPos+1; i<=endMidPartPos; i++){
			if(!isIlluminaAbPairPos(slideAlnVector, i)){
				abpairendpos = i-1;
				break;
			}
		}

		if(abpairendpos==-1) abpairendpos = endMidPartPos;
	}else{
		for(int64_t j=checkPos-1; j>=startMidPartPos; j--){
			if(isIlluminaAbPairPos(slideAlnVector, j)){
				abpairendpos = j;
				break;
			}
		}

		if(abpairendpos==-1) abpairendpos = startMidPartPos;
	}
	return 0;
}

bool Region::isIlluminaAbPairPos(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){
	bool flag = false;
	int64_t checkPos;
	int64_t count_abref, count_abdirection, count_paired_reads, count_exisizemax, count_lessisizemin;
	double ratio_baseexisizemax, ratio_baselessisizemin, ratio_baseabdirection, ratio_baseabref;

	checkPos = startCheckPos;
	count_abref = count_abdirection = count_paired_reads = count_exisizemax = count_lessisizemin = 0;
	ratio_baseexisizemax = ratio_baselessisizemin = ratio_baseabdirection = ratio_baseabref = 0;

	for(size_t i=0; i<slideAlnVector.size(); i++){
		if(bam_endpos(slideAlnVector[i])>=checkPos){
			if(slideAlnVector[i]->core.pos<=checkPos){
				if(slideAlnVector[i]->core.tid == slideAlnVector[i]->core.mtid){ // the query and the query's mate are on the reference
					if(slideAlnVector[i]->core.pos < slideAlnVector[i]->core.mpos and (slideAlnVector[i]->core.flag&BAM_FSUPPLEMENTARY) == 0 and slideAlnVector[i]->core.isize > 0){

						if(!bam_is_rev(slideAlnVector[i]) and bam_is_mrev(slideAlnVector[i])){
							count_paired_reads++; // get the paired reads count

							if(slideAlnVector[i]->core.isize >= paras->insert_max) count_exisizemax++; // the paired reads count of the insert size more than the insert size maximum
							if(slideAlnVector[i]->core.isize <= paras->insert_min) count_lessisizemin++; // the paired reads count of the insert size lower than the insert size minimum
						}else
							count_abdirection++; // the paired reads count of the error pairing direction
					}
				}else
					count_abref ++;
			}
			else if(slideAlnVector[i]->core.pos>checkPos) break;
		}

	}

	if(regBaseArr[checkPos-startRPos].coverage.num_bases[5]>0){
		ratio_baseexisizemax = (double) count_exisizemax / regBaseArr[checkPos-startRPos].coverage.num_bases[5];
		ratio_baselessisizemin = (double) count_lessisizemin / regBaseArr[checkPos-startRPos].coverage.num_bases[5];
		ratio_baseabdirection = (double) count_abdirection / regBaseArr[checkPos-startRPos].coverage.num_bases[5];
		ratio_baseabref = (double) count_abref / regBaseArr[checkPos-startRPos].coverage.num_bases[5];
		if(ratio_baseexisizemax>MAX_PAIR_THRES_POS or ratio_baselessisizemin>=MAX_PAIR_THRES_POS or ratio_baseabdirection>=MAX_PAIR_THRES_POS or ratio_baseabref>=MAX_PAIR_THRES_POS*2){
			flag = true;
		}
	}

	return flag;
}

// get clipping
void Region::detectIlluminaBaseClipReg(){
	reg_t *reg = NULL;
	int64_t i = startMidPartPos;
	if(i<1) i = 1;
	while(i<endMidPartPos){

		reg = getIlluminaBaseClipReg(i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + 1;
		} else break;
	}
	clipRegVector.shrink_to_fit();
}


// clip region without gap
reg_t* Region::getIlluminaBaseClipReg(int64_t startCheckPos){
	reg_t *reg = NULL;
	int64_t i, checkPos, startPos_clip, endPos_clip;
	string misType_clip;

	checkPos = startCheckPos;
	startPos_clip = endPos_clip = -1;
	misType_clip = "";

	if(!isIlluminaGapReg(checkPos, endMidPartPos)){
		while(checkPos<=endMidPartPos){
			if(regBaseArr[checkPos-startRPos].isIlluminaHighClipBase(paras->highclipRatio)){
				startPos_clip = checkPos;
				for(i=startPos_clip+1; i<=endMidPartPos; i++){
					if(regBaseArr[i-startRPos].isIlluminaHighClipBase(paras->highclipRatio)){
						continue;
					}else{
						endPos_clip = i-1;
						misType_clip = "High_clipping";
						reg = allocateMisReg(chrname, startPos_clip, endPos_clip, misType_clip);
						break;
					}
				}
				break;

			}else{
				checkPos++;
				startPos_clip = endPos_clip = -1;
			}
		}
	}

	return reg;
}

void Region::detectIlluminaGapMisjoinReg(vector<bam1_t*> &slideAlnVector){
	reg_t *reg = NULL;
	int64_t i = startMidPartPos;
	if(i<1) i = 1;
	while(i<endMidPartPos){

		reg = getIlluminaGapMisjoinReg(slideAlnVector, i);
		if(reg) {
			clipRegVector.push_back(reg);
			i = reg->endRefPos + 1;
		} else break;
	}
	clipRegVector.shrink_to_fit();
}

reg_t* Region::getIlluminaGapMisjoinReg(vector<bam1_t*> &slideAlnVector, int64_t startCheckPos){
	reg_t *reg = NULL;
	int64_t checkPos, i, startPos_gap, endPos_gap;
	size_t N_num = 0;
	string misType_gap;

	checkPos = startCheckPos;
	startPos_gap = endPos_gap = -1;
	misType_gap = "";

	while(checkPos<=endMidPartPos){
		if(regBaseArr[checkPos-startRPos].coverage.idx_RefBase==4){
			startPos_gap = checkPos;
			N_num++;
			for(i=checkPos+1; i<=endRPos; i++){
				if(regBaseArr[i-startRPos].coverage.idx_RefBase==4) N_num++;
				else break;
			}
			if(N_num>=MIN_GAP_SIZE){
				endPos_gap = i-1;
				estIlluminaRegInsertsize(slideAlnVector, startPos_gap-paras->gapextendSize, endPos_gap+paras->gapextendSize);
				if(fragmentsize>=paras->insert_max or fragmentsize<=paras->insert_min){
				misType_gap = "Abnormal_gap";
				reg = allocateMisReg(chrname, startPos_gap, endPos_gap, misType_gap);	//20220504
				}
				break;
			}
		}else{
			checkPos++;
			N_num = 0;
			startPos_gap = endPos_gap = -1;
		}
	}
	return reg;
}


bool Region::isIlluminaClipReg(int64_t startPos, int64_t endPos, double clip_ratio_thres){
	bool flag = false;
	int64_t i;
	size_t clip_num;
	double ratio1;

	clip_num = 0;
	for(i=startPos; i<=endPos; i++){
		clip_num += regBaseArr[i-startRPos].clipVector.size();
	}

	if(localRegCov>0){
		ratio1 = clip_num / localRegCov;
		if(ratio1>=clip_ratio_thres){
			flag = true;
			cout << chrname << ":" << startPos << "-" << endPos << "  clip_num:" << clip_num << ", ratio1:" << ratio1 << endl;
		}
	}

	return flag;
}

bool Region::isIlluminaBaseClipReg(int64_t startPos, int64_t endPos, double clip_ratio_thres){
	bool flag = false;
	int64_t i;
	double ratio2;

	for(i=startPos; i<=endPos; i++){
		if(regBaseArr[i-startRPos].coverage.num_bases[5]>0){
			ratio2 = (double) regBaseArr[i-startRPos].clipVector.size()/regBaseArr[i-startRPos].coverage.num_bases[5];
			if(ratio2>=clip_ratio_thres){
				flag = true;
				cout << chrname << ":" << i << ", clip_num:" << regBaseArr[i-startRPos].clipVector.size() << ", ratio2:" << ratio2 << endl;
			}
		}
	}
	return flag;
}

bool Region::isIlluminaAbisize(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos){
	bool flag = false;

	estIlluminaRegInsertsize(slideAlnVector, startPos, endPos);

	if(fragmentsize>=paras->insert_max or fragmentsize<=paras->insert_min){
		flag = true;
		cout << chrname << ":" << startPos << "-" << endPos << " , abnormal insert size:" << fragmentsize << endl;
	}
	return flag;
}

bool Region::isIlluminaAbPairing(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos){
	bool flag = false;

	estIlluminaPairing(slideAlnVector, startPos, endPos);
//	if(ratio_exisizemax>=MAX_PAIR_THRES or ratio_lessisizemin>=MAX_PAIR_THRES or ratio_abdirection>=MAX_PAIR_THRES){
//	if(ratio_exisizemax>=MAX_PAIR_THRES or ratio_lessisizemin>=MAX_PAIR_THRES or ratio_abdirection>=MAX_PAIR_THRES or ratio_abref>=MAX_PAIR_THRES*2){
	if(ratio_exisizemax>=paras->abpairRatio or ratio_lessisizemin>=paras->abpairRatio or ratio_abdirection>=paras->abpairRatio or ratio_abref>=paras->abpairRatio*2){
		flag = true;
//		cout << chrname << ":" << startPos << "-" << endPos << " , abnormal pairing, ratio_exisizemax:" << ratio_exisizemax << ", ratio_lessisizemin:" << ratio_lessisizemin << ", ratio_abdirection:" << ratio_abdirection << ", ratio_abref:" << ratio_abref << endl;
	}
	return flag;
}

void Region::estIlluminaPairing(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos){
	int64_t count_reads, count_abref, count_abdirection, count_paired_reads, count_exisizemax, count_lessisizemin;
	count_reads = count_abref = count_abdirection = count_paired_reads = count_exisizemax = count_lessisizemin = 0;

	for(size_t i=0; i<slideAlnVector.size(); i++){
		if(bam_endpos(slideAlnVector[i])>=startPos){
			if(bam_endpos(slideAlnVector[i])<=endPos or (slideAlnVector[i]->core.pos<=endPos and bam_endpos(slideAlnVector[i])>endPos)){
				if(slideAlnVector[i]->core.tid == slideAlnVector[i]->core.mtid){ // the query and the query's mate are on the reference
					if(slideAlnVector[i]->core.pos < slideAlnVector[i]->core.mpos and (slideAlnVector[i]->core.flag&BAM_FSUPPLEMENTARY) == 0 and slideAlnVector[i]->core.isize > 0){

						if(!bam_is_rev(slideAlnVector[i]) and bam_is_mrev(slideAlnVector[i])){
							count_paired_reads++; // get the paired reads count

							if(slideAlnVector[i]->core.isize >= paras->insert_max) count_exisizemax++; // the paired reads count of the insert size more than the insert size maximum
							if(slideAlnVector[i]->core.isize <= paras->insert_min) count_lessisizemin++; // the paired reads count of the insert size lower than the insert size minimum
						}else
							count_abdirection++; // the paired reads count of the error pairing direction
					}
				}else
					count_abref ++;
			}
			else if(slideAlnVector[i]->core.pos>endPos) break;
		}

	}
	count_reads = count_abref + count_abdirection + count_paired_reads;

	ratio_exisizemax = (double) count_exisizemax / count_reads;
	ratio_lessisizemin = (double) count_lessisizemin / count_reads;
	ratio_abdirection = (double) count_abdirection / count_reads;
	ratio_abref = (double) count_abref / count_reads;
}

bool Region::isIlluminaGapClipReg(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos){
	bool flag = false;
	int64_t startPos_tmp, endPos_tmp;
//	20220502
//	startPos_tmp = startPos-GAP_EXTEND_REG_SIZE;
//	endPos_tmp = endPos+GAP_EXTEND_REG_SIZE;
	startPos_tmp = startPos-paras->gapextendSize;
	endPos_tmp = endPos+paras->gapextendSize;

	estIlluminaRegInsertsize(slideAlnVector, startPos_tmp, endPos_tmp);

	if(fragmentsize>=paras->insert_max or fragmentsize<=paras->insert_min) {
		flag = true;
		cout << chrname << ":" << startPos << "-" << endPos << " , include gap, abnormal insert size: " << fragmentsize << endl;
	}

	return flag;
}

void Region::estIlluminaRegInsertsize(vector<bam1_t*> &slideAlnVector, int64_t startPos, int64_t endPos){
	size_t i;
	int64_t count = 0;
	double sum = 0;

	for(i=0; i<slideAlnVector.size(); i++){
		if(bam_endpos(slideAlnVector[i])>=startPos){
			if(bam_endpos(slideAlnVector[i])<=endPos or (slideAlnVector[i]->core.pos<=endPos and bam_endpos(slideAlnVector[i])>endPos)){
				if(slideAlnVector[i]->core.tid==slideAlnVector[i]->core.mtid and slideAlnVector[i]->core.pos < slideAlnVector[i]->core.mpos and slideAlnVector[i]->core.isize > 0){ // the query and the query's mate are on the reference
							sum += slideAlnVector[i]->core.isize; // get the sum of the insert size
							count++; // get the paired reads count of the query and the query's mate are on the reference
				}
			}
			else if(slideAlnVector[i]->core.pos>endPos) break;
		}
	}

	fragmentsize = (double) sum/count;
}

bool Region::isIlluminaGapReg(int64_t startPos, int64_t endPos){
	bool flag = false;
	int64_t i, j;
	size_t N_num = 0;

	i=startPos;
	while(i<=endPos){
		if(regBaseArr[i-startRPos].coverage.idx_RefBase==4){
			N_num++;
			for(j=i+1; j<=endRPos; j++){
				if(regBaseArr[j-startRPos].coverage.idx_RefBase==4) N_num++;
				else break;
			}
			if(N_num>=MIN_GAP_SIZE){
				flag = true;
//				cout << chrname << ":" << i << "-" << j-1 << ", gapsize:" << N_num << endl;
				break;
			}
		}else{
			i++;
			N_num = 0;
		}
	}

	return flag;
}


