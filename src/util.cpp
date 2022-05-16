#include <iostream>
#include <cstring>
#include <unistd.h>
#include <pthread.h>
#include <htslib/thread_pool.h>

#include "util.h"
#include "covLoader.h"
#include "LocalAssembly.h"

// string split function
vector<string> split(const string& s, const string& delim)
{
    vector<string> elems;
    size_t pos = 0;
    size_t len = s.length();
    size_t delim_len = delim.length();
    if (delim_len == 0) return elems;
    while (pos < len)
    {
        int find_pos = s.find(delim, pos);
        if (find_pos < 0)
        {
            elems.push_back(s.substr(pos, len - pos));
            break;
        }
        elems.push_back(s.substr(pos, find_pos - pos));
        pos = find_pos + delim_len;
    }
    return elems;
}

// determine whether the str exist in string vector
bool isExistStr(string &str, vector<string> &vec){
	bool exist_flag = false;
	for(size_t i=0; i<vec.size(); i++)
		if(str.compare(vec.at(i))==0){
			exist_flag = true;
			break;
		}
	return exist_flag;
}

// copy single file
int copySingleFile(string &infilename, ofstream &outfile){
	ifstream infile;
	string line;

	if(isFileExist(infilename)==false) return 0;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}
	// read each line and save to the output file
	while(getline(infile, line))
		if(line.size()>0 and line.at(0)!='#'){
			outfile << line << endl;
		}
	infile.close();

	return 0;
}

char getBase(string &seq, size_t pos, size_t orient){
	char ch;

	if(pos<1 or pos>seq.size()){
		cerr << __func__ << ", line=" << __LINE__ << ": invalid pos: " << pos << endl;
		exit(1);
	}

	ch = seq[pos-1];
	if(orient==ALN_MINUS_ORIENT){
		switch(ch){
			case 'A': ch = 'T'; break;
			case 'C': ch = 'G'; break;
			case 'G': ch = 'C'; break;
			case 'T': ch = 'A'; break;
			case 'a': ch = 't'; break;
			case 'c': ch = 'g'; break;
			case 'g': ch = 'c'; break;
			case 't': ch = 'a'; break;
			case 'N':
			case 'n': break;
			default: ch = 'N'; break;
		}
	}

	return ch;
}

// reverse the sequence
void reverseSeq(string &seq){
	char tmp;
	int32_t len = seq.size();
	for(int32_t i=0; i<len/2; i++){
		tmp = seq[i];
		seq[i] = seq[len-1-i];
		seq[len-1-i] = tmp;
	}
}

// reverse complement the sequence
void reverseComplement(string &seq){
	reverseSeq(seq);

	int32_t len = seq.size();
	for(int32_t i=0; i<len; i++){
		switch(seq[i]){
			case 'a': seq[i] = 't'; break;
			case 'c': seq[i] = 'g'; break;
			case 'g': seq[i] = 'c'; break;
			case 't': seq[i] = 'a'; break;
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
			case 'N': seq[i] = 'N'; break;
			case 'n': seq[i] = 'n'; break;

			// mixed bases
			case 'M': seq[i] = 'K'; break;
			case 'm': seq[i] = 'k'; break;
			case 'R': seq[i] = 'Y'; break;
			case 'r': seq[i] = 'y'; break;
			case 'S': seq[i] = 'S'; break;
			case 's': seq[i] = 's'; break;
			case 'V': seq[i] = 'B'; break;
			case 'v': seq[i] = 'b'; break;
			case 'W': seq[i] = 'W'; break;
			case 'w': seq[i] = 'w'; break;
			case 'Y': seq[i] = 'R'; break;
			case 'y': seq[i] = 'r'; break;
			case 'H': seq[i] = 'D'; break;
			case 'h': seq[i] = 'd'; break;
			case 'K': seq[i] = 'M'; break;
			case 'k': seq[i] = 'm'; break;
			case 'D': seq[i] = 'H'; break;
			case 'd': seq[i] = 'h'; break;
			case 'B': seq[i] = 'V'; break;
			case 'b': seq[i] = 'v'; break;
			default: cerr << __func__ << ": unknown base: " << seq[i] << endl; exit(1);
		}
	}
}

// upper the sequence
void upperSeq(string &seq){
	for(size_t i=0; i<seq.size(); i++)
		if(seq[i]>='a' and seq[i]<='z') seq[i] -= 32;
}

// get the number of contigs
size_t getCtgCount(string &contigfilename){
	size_t ctg_num;
	string line;
	ifstream infile;
	bool flag = false;

	flag = isFileExist(contigfilename);
	ctg_num = 0;
	if(flag){
		infile.open(contigfilename);
		if(!infile.is_open()){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << contigfilename << endl;
			exit(1);
		}

		while(getline(infile, line))
			if(line.size())
				if(line[0]=='>')	ctg_num ++;
		infile.close();
	}

	return ctg_num;
}

// find the vector item
reg_t* findVarvecItem(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec){
	reg_t *reg, *target_reg = NULL;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos>=reg->startRefPos and startRefPos<=reg->endRefPos) or (endRefPos>=reg->startRefPos and endRefPos<=reg->endRefPos)
			or (reg->startRefPos>=startRefPos and reg->startRefPos<=endRefPos) or (reg->endRefPos>=startRefPos and reg->endRefPos<=endRefPos)){
			// overlap
			target_reg = reg;
			break;
		}
	}
	return target_reg;
}

// find the vector item 20210424
reg_t* findOverlapvecItem(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec){
	reg_t *reg, *target_reg = NULL;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos<=reg->endRefPos and reg->endRefPos<=endRefPos) or (reg->endRefPos>endRefPos and reg->startRefPos<=endRefPos)){
			// overlap
			target_reg = reg;
			break;
		}
	}
	return target_reg;
}

// find all the vector item
vector<reg_t*> findVarvecItemAll(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec){
	reg_t *reg;
	vector<reg_t*> reg_vec_ret;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos>=reg->startRefPos and startRefPos<=reg->endRefPos) or (endRefPos>=reg->startRefPos and endRefPos<=reg->endRefPos)
			or (reg->startRefPos>=startRefPos and reg->startRefPos<=endRefPos) or (reg->endRefPos>=startRefPos and reg->endRefPos<=endRefPos)){
			// overlap
			reg_vec_ret.push_back(reg);
		}
	}
	return reg_vec_ret;
}

// find the vector item according to extended margin size
reg_t* findVarvecItemExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize){
	reg_t *reg, *target_reg = NULL;
	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos+leftExtSize>=reg->startRefPos and startRefPos<=reg->endRefPos+rightExtSize) or (endRefPos+leftExtSize>=reg->startRefPos and endRefPos<=reg->endRefPos+rightExtSize)
			or (reg->startRefPos+leftExtSize>=startRefPos and reg->startRefPos<=endRefPos+rightExtSize) or (reg->endRefPos+leftExtSize>=startRefPos and reg->endRefPos<=endRefPos+rightExtSize)){
			// overlap
			target_reg = reg;
			break;
		}
	}
	return target_reg;
}

// find all the vector item according to extended margin size
vector<reg_t*> findVarvecItemAllExtSize(int32_t startRefPos, int32_t endRefPos, vector<reg_t*> &varVec, int32_t leftExtSize, int32_t rightExtSize){
	reg_t *reg;
	vector<reg_t*> reg_vec_ret;

	for(size_t i=0; i<varVec.size(); i++){
		reg = varVec[i];
		if((startRefPos+leftExtSize>=reg->startRefPos and startRefPos<=reg->endRefPos+rightExtSize) or (endRefPos+leftExtSize>=reg->startRefPos and endRefPos<=reg->endRefPos+rightExtSize)
			or (reg->startRefPos+leftExtSize>=startRefPos and reg->startRefPos<=endRefPos+rightExtSize) or (reg->endRefPos+leftExtSize>=startRefPos and reg->endRefPos<=endRefPos+rightExtSize)){
			// overlap
			reg_vec_ret.push_back(reg);
		}
	}
	return reg_vec_ret;
}

// get the vector index
int32_t getVectorIdx(reg_t *reg, vector<reg_t*> &varVec){
	int32_t idx = -1;
	for(size_t i=0; i<varVec.size(); i++)
		if(varVec[i]==reg){
			idx = i;
			break;
		}
	return idx;
}

// get overlapped region
int32_t getOverlappedRegIdxByCallFlag(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_tmp;
	bool flag;
	int32_t idx_ret = -1;
	for(size_t i=0; i<varVec.size(); i++){
		reg_tmp = varVec.at(i);
		if(reg_tmp->call_success_status){
			flag = isOverlappedReg(reg, reg_tmp);
			if(flag){
				idx_ret = i;
				break;
			}
		}
	}
	return idx_ret;
}

reg_t* getOverlappedRegByCallFlag(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_ret = NULL;
	int32_t idx = getOverlappedRegIdxByCallFlag(reg, varVec);
	if(idx!=-1) reg_ret = varVec.at(idx);
	return reg_ret;
}

reg_t* getOverlappedReg(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_ret = NULL;
	int32_t idx = getOverlappedRegIdx(reg, varVec);
	if(idx!=-1) reg_ret = varVec.at(idx);
	return reg_ret;
}


// get overlapped region
int32_t getOverlappedRegIdx(reg_t *reg, vector<reg_t*> &varVec){
	reg_t *reg_tmp;
	bool flag;
	int32_t idx_ret = -1;
	for(size_t i=0; i<varVec.size(); i++){
		reg_tmp = varVec.at(i);
		flag = isOverlappedReg(reg, reg_tmp);
		if(flag){
			idx_ret = i;
			break;
		}
	}
	return idx_ret;
}

// determine whether two regions are overlapped
bool isOverlappedReg(reg_t* reg1, reg_t* reg2){
	bool flag = false;
	if(reg1->chrname.compare(reg2->chrname)==0){
		if(isOverlappedPos(reg1->startRefPos, reg1->endRefPos, reg2->startRefPos, reg2->endRefPos))
			flag = true;
	}
	return flag;
}

// determine whether the given positions of two regions are overlapped
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	bool flag = false;
	if((startPos1>=startPos2 and startPos1<=endPos2)
		or (endPos2>=startPos1 and endPos2<=endPos1)
		or (startPos2>=startPos1 and startPos2<=endPos1)
		or (endPos1>=startPos2 and endPos1<=endPos2))
			flag = true;
	return flag;
}

int32_t getOverlapSize(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	int32_t overlap_size;

	if(startPos1<=startPos2) overlap_size = endPos1 - startPos2 + 1;
	else overlap_size = endPos2 - startPos1 + 1;

	return overlap_size;
}

bool isAdjacent(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2, int32_t dist_thres){
	bool flag;
	int32_t overlap_size;

	if(dist_thres<=0){
		cerr << __func__ << ", line=" << __LINE__ << ", invalid distance threshold: " << dist_thres << endl;
		exit(1);
	}

	overlap_size = getOverlapSize(startPos1, endPos1, startPos2, endPos2);
	//if(overlap_size>=-dist_thres and overlap_size<=dist_thres) flag = true;
	if(overlap_size>=-dist_thres) flag = true;
	else flag = false;

	return flag;
}

// determine whether the two mate chip regions are overlapped
bool isOverlappedMateClipReg(mateClipReg_t *mate_clip_reg1, mateClipReg_t *mate_clip_reg2){
	bool overlap_flag, overlap_flag1, overlap_flag2;
	reg_t *left_reg1, *right_reg1, *left_reg2, *right_reg2;
	size_t mate_end_flag, start_pos1, end_pos1, start_pos2, end_pos2;
	string chrname1, chrname2;

	overlap_flag = overlap_flag1 = overlap_flag2 = false;
	if(mate_clip_reg1->leftClipRegNum==1 and mate_clip_reg1->rightClipRegNum==1 and mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==1){
		left_reg1 = mate_clip_reg1->leftClipReg ? mate_clip_reg1->leftClipReg : mate_clip_reg1->leftClipReg2;
		right_reg1 = mate_clip_reg1->rightClipReg ? mate_clip_reg1->rightClipReg : mate_clip_reg1->rightClipReg2;
		left_reg2 = mate_clip_reg2->leftClipReg ? mate_clip_reg2->leftClipReg : mate_clip_reg2->leftClipReg2;
		right_reg2 = mate_clip_reg2->rightClipReg ? mate_clip_reg2->rightClipReg : mate_clip_reg2->rightClipReg2;

		if((left_reg1->chrname.compare(left_reg2->chrname)==0 and right_reg1->chrname.compare(right_reg2->chrname)==0 and isOverlappedReg(left_reg1, left_reg2) and isOverlappedReg(right_reg1, right_reg2))
			or (left_reg1->chrname.compare(right_reg2->chrname)==0 and right_reg1->chrname.compare(left_reg2->chrname)==0 and isOverlappedReg(left_reg1, right_reg2) and isOverlappedReg(right_reg1, left_reg2)))
			overlap_flag = true;

	}else if(mate_clip_reg1->leftClipRegNum==2 and mate_clip_reg1->rightClipRegNum==2 and mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==2){
		chrname1 = chrname2 = "";
		start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

		if(mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg1->leftClipReg2->chrname)==0 and mate_clip_reg1->rightClipReg->chrname.compare(mate_clip_reg1->rightClipReg2->chrname)==0
			and mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0 and mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
			// check left region
			chrname1 = mate_clip_reg1->leftClipReg->chrname;
			if(mate_clip_reg1->leftMeanClipPos<mate_clip_reg1->leftMeanClipPos2){
				start_pos1 = mate_clip_reg1->leftMeanClipPos;
				end_pos1 = mate_clip_reg1->leftMeanClipPos2;
			}else{
				start_pos1 = mate_clip_reg1->leftMeanClipPos2;
				end_pos1 = mate_clip_reg1->leftMeanClipPos;
			}

			chrname2 = mate_clip_reg2->leftClipReg->chrname;
			if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
				start_pos2 = mate_clip_reg2->leftMeanClipPos;
				end_pos2 = mate_clip_reg2->leftMeanClipPos2;
			}else{
				start_pos2 = mate_clip_reg2->leftMeanClipPos2;
				end_pos2 = mate_clip_reg2->leftMeanClipPos;
			}

			if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
			else overlap_flag1 = false;

			// check right region
			chrname1 = mate_clip_reg1->rightClipReg->chrname;
			if(mate_clip_reg1->rightMeanClipPos<mate_clip_reg1->rightMeanClipPos2){
				start_pos1 = mate_clip_reg1->rightMeanClipPos;
				end_pos1 = mate_clip_reg1->rightMeanClipPos2;
			}else{
				start_pos1 = mate_clip_reg1->rightMeanClipPos2;
				end_pos1 = mate_clip_reg1->rightMeanClipPos;
			}

			chrname2 = mate_clip_reg2->rightClipReg->chrname;
			if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
				start_pos2 = mate_clip_reg2->rightMeanClipPos;
				end_pos2 = mate_clip_reg2->rightMeanClipPos2;
			}else{
				start_pos2 = mate_clip_reg2->rightMeanClipPos2;
				end_pos2 = mate_clip_reg2->rightMeanClipPos;
			}

			if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
			else overlap_flag2 = false;

			if(overlap_flag1 and overlap_flag2) overlap_flag = true;
		}
	}else{
		if(mate_clip_reg1->leftClipRegNum==1 and mate_clip_reg1->rightClipRegNum==2){
			// check left region
			chrname1 = chrname2 = "";
			start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

			if(mate_clip_reg1->leftClipReg){
				chrname1 = mate_clip_reg1->leftClipReg->chrname;
				start_pos1 = mate_clip_reg1->leftClipReg->startRefPos;
				end_pos1 = mate_clip_reg1->leftClipReg->endRefPos;
			}else{
				chrname1 = mate_clip_reg1->leftClipReg2->chrname;
				start_pos1 = mate_clip_reg1->leftClipReg2->startRefPos;
				end_pos1 = mate_clip_reg1->leftClipReg2->endRefPos;
			}

			// get the mate end
			mate_end_flag = 0;
			if(mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==2){
				mate_end_flag = 1;
			}else if(mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==1){
				mate_end_flag = 2;
			}

			if(mate_end_flag==1){
				if(mate_clip_reg2->leftClipReg){
					chrname2 = mate_clip_reg2->leftClipReg->chrname;
					start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
					end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
				}else{
					chrname2 = mate_clip_reg2->leftClipReg2->chrname;
					start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
					end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
				}
			}else if(mate_end_flag==2){
				if(mate_clip_reg2->rightClipReg){
					chrname2 = mate_clip_reg2->rightClipReg->chrname;
					start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
					end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
				}else{
					chrname2 = mate_clip_reg2->rightClipReg2->chrname;
					start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
					end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
				}
			}

			if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
			else overlap_flag1 = false;

			// check right regions
			chrname1 = chrname2 = "";
			start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

			if(mate_end_flag==1) mate_end_flag = 2;
			else if(mate_end_flag==2) mate_end_flag = 1;

			if(mate_clip_reg1->rightClipReg->chrname.compare(mate_clip_reg1->rightClipReg2->chrname)==0){
				chrname1 = mate_clip_reg1->rightClipReg->chrname;
				if(mate_clip_reg1->rightMeanClipPos<mate_clip_reg1->rightMeanClipPos2){
					start_pos1 = mate_clip_reg1->rightMeanClipPos;
					end_pos1 = mate_clip_reg1->rightMeanClipPos2;
				}else{
					start_pos1 = mate_clip_reg1->rightMeanClipPos2;
					end_pos1 = mate_clip_reg1->rightMeanClipPos;
				}

				if(mate_end_flag==1){
					if(mate_clip_reg2->leftClipRegNum==1){
						if(mate_clip_reg2->leftClipReg){
							chrname2 = mate_clip_reg2->leftClipReg->chrname;
							start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
							end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
						}else{
							chrname2 = mate_clip_reg2->leftClipReg2->chrname;
							start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
							end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
						}
					}else if(mate_clip_reg2->leftClipRegNum==2){
						if(mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0){
							chrname2 = mate_clip_reg2->leftClipReg->chrname;
							if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
								start_pos2 = mate_clip_reg2->leftMeanClipPos;
								end_pos2 = mate_clip_reg2->leftMeanClipPos2;
							}else{
								start_pos2 = mate_clip_reg2->leftMeanClipPos2;
								end_pos2 = mate_clip_reg2->leftMeanClipPos;
							}
						}
					}
				}else if(mate_end_flag==2){
					if(mate_clip_reg2->rightClipRegNum==1){
						if(mate_clip_reg2->rightClipReg){
							chrname2 = mate_clip_reg2->rightClipReg->chrname;
							start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
							end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
						}else{
							chrname2 = mate_clip_reg2->rightClipReg2->chrname;
							start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
							end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
						}
					}else if(mate_clip_reg2->rightClipRegNum==2){
						if(mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
							chrname2 = mate_clip_reg2->rightClipReg->chrname;
							if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
								start_pos2 = mate_clip_reg2->rightMeanClipPos;
								end_pos2 = mate_clip_reg2->rightMeanClipPos2;
							}else{
								start_pos2 = mate_clip_reg2->rightMeanClipPos2;
								end_pos2 = mate_clip_reg2->rightMeanClipPos;
							}
						}
					}
				}
				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag2 = false;
			}

			if(overlap_flag1 and overlap_flag2) overlap_flag = true;
		}else if(mate_clip_reg1->leftClipRegNum==2 and mate_clip_reg1->rightClipRegNum==1){
			// check right region
			chrname1 = chrname2 = "";
			start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

			if(mate_clip_reg1->rightClipReg){
				chrname1 = mate_clip_reg1->rightClipReg->chrname;
				start_pos1 = mate_clip_reg1->rightClipReg->startRefPos;
				end_pos1 = mate_clip_reg1->rightClipReg->endRefPos;
			}else{
				chrname1 = mate_clip_reg1->rightClipReg2->chrname;
				start_pos1 = mate_clip_reg1->rightClipReg2->startRefPos;
				end_pos1 = mate_clip_reg1->rightClipReg2->endRefPos;
			}

			// get the mate end
			mate_end_flag = 0;
			if(mate_clip_reg2->leftClipRegNum==1 and mate_clip_reg2->rightClipRegNum==2){
				mate_end_flag = 1;
			}else if(mate_clip_reg2->leftClipRegNum==2 and mate_clip_reg2->rightClipRegNum==1){
				mate_end_flag = 2;
			}

			if(mate_end_flag==1){
				if(mate_clip_reg2->leftClipReg){
					chrname2 = mate_clip_reg2->leftClipReg->chrname;
					start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
					end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
				}else{
					chrname2 = mate_clip_reg2->leftClipReg2->chrname;
					start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
					end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
				}
			}else if(mate_end_flag==2){
				if(mate_clip_reg2->rightClipReg){
					chrname2 = mate_clip_reg2->rightClipReg->chrname;
					start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
					end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
				}else{
					chrname2 = mate_clip_reg2->rightClipReg2->chrname;
					start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
					end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
				}
			}

			if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag1 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
			else overlap_flag1 = false;

			// check left regions
			chrname1 = chrname2 = "";
			start_pos1 = end_pos1 = start_pos2 = end_pos2 = 0;

			if(mate_end_flag==1) mate_end_flag = 2;
			else if(mate_end_flag==2) mate_end_flag = 1;

			if(mate_clip_reg1->leftClipReg->chrname.compare(mate_clip_reg1->leftClipReg2->chrname)==0){
				chrname1 = mate_clip_reg1->leftClipReg->chrname;
				if(mate_clip_reg1->leftMeanClipPos<mate_clip_reg1->leftMeanClipPos2){
					start_pos1 = mate_clip_reg1->leftMeanClipPos;
					end_pos1 = mate_clip_reg1->leftMeanClipPos2;
				}else{
					start_pos1 = mate_clip_reg1->leftMeanClipPos2;
					end_pos1 = mate_clip_reg1->leftMeanClipPos;
				}

				if(mate_end_flag==1){
					if(mate_clip_reg2->leftClipRegNum==1){
						if(mate_clip_reg2->leftClipReg){
							chrname2 = mate_clip_reg2->leftClipReg->chrname;
							start_pos2 = mate_clip_reg2->leftClipReg->startRefPos;
							end_pos2 = mate_clip_reg2->leftClipReg->endRefPos;
						}else{
							chrname2 = mate_clip_reg2->leftClipReg2->chrname;
							start_pos2 = mate_clip_reg2->leftClipReg2->startRefPos;
							end_pos2 = mate_clip_reg2->leftClipReg2->endRefPos;
						}
					}else if(mate_clip_reg2->leftClipRegNum==2){
						if(mate_clip_reg2->leftClipReg->chrname.compare(mate_clip_reg2->leftClipReg2->chrname)==0){
							chrname2 = mate_clip_reg2->leftClipReg->chrname;
							if(mate_clip_reg2->leftMeanClipPos<mate_clip_reg2->leftMeanClipPos2){
								start_pos2 = mate_clip_reg2->leftMeanClipPos;
								end_pos2 = mate_clip_reg2->leftMeanClipPos2;
							}else{
								start_pos2 = mate_clip_reg2->leftMeanClipPos2;
								end_pos2 = mate_clip_reg2->leftMeanClipPos;
							}
						}
					}
				}else if(mate_end_flag==2){
					if(mate_clip_reg2->rightClipRegNum==1){
						if(mate_clip_reg2->rightClipReg){
							chrname2 = mate_clip_reg2->rightClipReg->chrname;
							start_pos2 = mate_clip_reg2->rightClipReg->startRefPos;
							end_pos2 = mate_clip_reg2->rightClipReg->endRefPos;
						}else{
							chrname2 = mate_clip_reg2->rightClipReg2->chrname;
							start_pos2 = mate_clip_reg2->rightClipReg2->startRefPos;
							end_pos2 = mate_clip_reg2->rightClipReg2->endRefPos;
						}
					}else if(mate_clip_reg2->rightClipRegNum==2){
						if(mate_clip_reg2->rightClipReg->chrname.compare(mate_clip_reg2->rightClipReg2->chrname)==0){
							chrname2 = mate_clip_reg2->rightClipReg->chrname;
							if(mate_clip_reg2->rightMeanClipPos<mate_clip_reg2->rightMeanClipPos2){
								start_pos2 = mate_clip_reg2->rightMeanClipPos;
								end_pos2 = mate_clip_reg2->rightMeanClipPos2;
							}else{
								start_pos2 = mate_clip_reg2->rightMeanClipPos2;
								end_pos2 = mate_clip_reg2->rightMeanClipPos;
							}
						}
					}
				}
				if(chrname1.size()>0 and chrname1.compare(chrname2)==0) overlap_flag2 = isOverlappedPos(start_pos1, end_pos1, start_pos2, end_pos2); // same chrome
				else overlap_flag2 = false;
			}

			if(overlap_flag1 and overlap_flag2) overlap_flag = true;
		}
	}

	return overlap_flag;
}

// get the overlapped mate clip region
mateClipReg_t* getOverlappedMateClipReg(mateClipReg_t *mate_clip_reg_given, vector<mateClipReg_t*> &mateClipRegVec){
	mateClipReg_t *mate_clip_reg_overlapped = NULL, *clip_reg;

	for(size_t i=0; i<mateClipRegVec.size(); i++){
		clip_reg = mateClipRegVec.at(i);
		if(clip_reg!=mate_clip_reg_given and clip_reg->valid_flag and clip_reg->reg_mated_flag and clip_reg->sv_type==mate_clip_reg_given->sv_type){
			if(isOverlappedMateClipReg(mate_clip_reg_given, clip_reg)==true){ // same chrome
				mate_clip_reg_overlapped = clip_reg;
				break;
			}
		}
	}

	return mate_clip_reg_overlapped;
}

// load sam/bam header
bam_hdr_t* loadSamHeader(string &inBamFile)
{
	bam_hdr_t *header = NULL;
	samFile *in = 0;

    if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
        cerr << __func__ << ": failed to open " << inBamFile << " for reading" << endl;
        exit(1);
    }

    if ((header = sam_hdr_read(in)) == 0) {
        cerr << "fail to read the header from " << inBamFile << endl;
        exit(1);
    }
    sam_close(in);

	return header;
}

// determine whether a position in an region or not.
bool isInReg(int32_t pos, vector<reg_t*> &vec){
	bool flag = false;
	reg_t* reg;
	for(size_t i=0; i<vec.size(); i++){
		reg = vec.at(i);
		if(pos>=reg->startRefPos and pos<=reg->endRefPos) { flag = true; break; }
	}
	return flag;
}

// compute the number of disagreements
int32_t computeDisagreeNum(Base *baseArray, int32_t arr_size){
	int32_t i, disagreeNum = 0;
	for(i=0; i<arr_size; i++)
		if(baseArray[i].isZeroCovBase() or baseArray[i].isDisagreeBase())
			disagreeNum ++;
	return disagreeNum;
}

// merge overlapped indels
void mergeOverlappedReg(vector<reg_t*> &regVector){
	size_t i, j;
	reg_t *reg1, *reg2;

	for(i=0; i<regVector.size(); i++){
		reg1 = regVector.at(i);
		if(reg1->call_success_status){
			for(j=i+1; j<regVector.size(); ){
				reg2 = regVector.at(j);
				if(reg2->call_success_status and isOverlappedReg(reg1, reg2) and reg1->query_id!=-1 and reg1->query_id==reg2->query_id){ // same reference region and query
					updateReg(reg1, reg2);
					delete reg2;
					regVector.erase(regVector.begin()+j);
				}else j++;
			}
		}
	}
	regVector.shrink_to_fit();
}

// merge two overlapped indels
void updateReg(reg_t* reg1, reg_t* reg2){
	size_t new_startRefPos, new_endRefPos, new_startLocalRefPos, new_endLocalRefPos, new_startQueryPos, new_endQueryPos;
	if(reg1->startRefPos<=reg2->startRefPos) {
		new_startRefPos = reg1->startRefPos;
		new_startLocalRefPos = reg1->startLocalRefPos;
		new_startQueryPos = reg1->startQueryPos;
	}else{
		new_startRefPos = reg2->startRefPos;
		new_startLocalRefPos = reg2->startLocalRefPos;
		new_startQueryPos = reg2->startQueryPos;
	}
	if(reg1->endRefPos>=reg2->endRefPos){
		new_endRefPos = reg1->endRefPos;
		new_endLocalRefPos = reg1->endLocalRefPos;
		new_endQueryPos = reg1->endQueryPos;
	}else{
		new_endRefPos = reg2->endRefPos;
		new_endLocalRefPos = reg2->endLocalRefPos;
		new_endQueryPos = reg2->endQueryPos;
	}
	reg1->startRefPos = new_startRefPos;
	reg1->endRefPos = new_endRefPos;
	reg1->startLocalRefPos = new_startLocalRefPos;
	reg1->endLocalRefPos = new_endLocalRefPos;
	reg1->startQueryPos = new_startQueryPos;
	reg1->endQueryPos = new_endQueryPos;
}

// merge adjacent regions
void mergeAdjacentReg(vector<reg_t*> &regVec, size_t dist_thres){
	size_t i;
	reg_t *reg_tmp, *reg_tmp2;
	bool flag;
	for(i=1; i<regVec.size(); ){
		reg_tmp = regVec.at(i-1);
		reg_tmp2 = regVec.at(i);
		flag = isAdjacent(reg_tmp->startRefPos, reg_tmp->endRefPos, reg_tmp2->startRefPos, reg_tmp2->endRefPos, dist_thres);
		if(flag){ // merge
			if(reg_tmp->startRefPos>reg_tmp2->startRefPos) reg_tmp->startRefPos = reg_tmp2->startRefPos;
			if(reg_tmp->endRefPos<reg_tmp2->endRefPos) reg_tmp->endRefPos = reg_tmp2->endRefPos;
			delete reg_tmp2;
			regVec.erase(regVec.begin()+i);
		}else i ++;
	}
}

void printRegVec(vector<reg_t*> &regVec, string header){
	reg_t *reg;
	cout << header << " region size: " << regVec.size() << endl;
	for(size_t i=0; i<regVec.size(); i++){
		reg = regVec.at(i);
		cout << "[" << i << "]: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", localRef: " << reg->startLocalRefPos << "-" << reg->endLocalRefPos << ", Query :" << reg->startQueryPos << "-" << reg->endQueryPos << endl;
	}
}

vector<string> getLeftRightPartChrname(mateClipReg_t *mate_clip_reg){
	vector<string> chrname_vec;
	string chrname_left = "", chrname_right = "";

	// left part
	if(mate_clip_reg->leftClipRegNum==2){
		chrname_left = mate_clip_reg->leftClipReg->chrname;
	}else if(mate_clip_reg->leftClipRegNum==1){
		if(mate_clip_reg->leftClipReg)
			chrname_left = mate_clip_reg->leftClipReg->chrname;
		else if(mate_clip_reg->leftClipReg2)
			chrname_left = mate_clip_reg->leftClipReg2->chrname;
	}

	// right part
	if(mate_clip_reg->rightClipRegNum==2){
		chrname_right = mate_clip_reg->rightClipReg->chrname;
	}else if(mate_clip_reg->rightClipRegNum==1){
		if(mate_clip_reg->rightClipReg)
			chrname_right = mate_clip_reg->rightClipReg->chrname;
		else if(mate_clip_reg->rightClipReg2)
			chrname_right = mate_clip_reg->rightClipReg2->chrname;
	}

	chrname_vec.push_back(chrname_left);
	chrname_vec.push_back(chrname_right);

	return chrname_vec;
}

// preprocess pipe command characters
string preprocessPipeChar(string &cmd_str){
	char ch;
	string ret_str = "";

	for(size_t i=0; i<cmd_str.size(); i++){
		ch = cmd_str.at(i);
		if(ch=='|') // pipe character
			ret_str += "_";
		else
			ret_str += ch;
	}
	return ret_str;
}

bool isFileExist(string &filename){
	bool flag = false;
	struct stat fileStat;
	if (stat(filename.c_str(), &fileStat) == 0)
		if(fileStat.st_size>0)
			flag = true;
	return flag;
}

void removeRedundantItems(vector<reg_t*> &reg_vec){
	size_t i, j;
	reg_t *reg, *reg_tmp;

	for(i=0; i<reg_vec.size(); i++){
		reg = reg_vec.at(i);
		for(j=i+1; j<reg_vec.size(); ){
			reg_tmp = reg_vec.at(j);
			if(reg_tmp->chrname.compare(reg->chrname)==0 and reg_tmp->startRefPos==reg->startRefPos and reg_tmp->endRefPos==reg->endRefPos){ // redundant item
				delete reg_tmp;
				reg_vec.erase(reg_vec.begin()+j);
			}else j++;
		}
	}
}

// get line count of file
int32_t getLineCount(string &filename){
	int32_t num;
	ifstream infile;
	string line;

	infile.open(filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	num = 0;
	while(getline(infile, line)) if(line.size()>0 and line.at(0)!='#') num ++;
	infile.close();

	return num;
}

bool isBaseMatch(char ctgBase, char refBase){
	bool match_flag = false;

	// Upper case
	if(ctgBase>='a' and ctgBase<='z') ctgBase -= 32;
	if(refBase>='a' and refBase<='z') refBase -= 32;

	if(ctgBase==refBase){
		match_flag = true;
	}else if(ctgBase=='-' or refBase=='-'){
		match_flag = false;
	}else{
		switch(refBase){
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
			case 'N':
			case 'n':
				match_flag = false;
				break;
			case 'M':
			case 'm':
				if(ctgBase=='A' or ctgBase=='C') match_flag = true;
				break;
			case 'R':
			case 'r':
				if(ctgBase=='A' or ctgBase=='G') match_flag = true;
				break;
			case 'S':
			case 's':
				if(ctgBase=='C' or ctgBase=='G') match_flag = true;
				break;
			case 'V':
			case 'v':
				if(ctgBase=='A' or ctgBase=='C' or ctgBase=='G') match_flag = true;
				break;
			case 'W':
			case 'w':
				if(ctgBase=='A' or ctgBase=='T') match_flag = true;
				break;
			case 'Y':
			case 'y':
				if(ctgBase=='C' or ctgBase=='T') match_flag = true;
				break;
			case 'H':
			case 'h':
				if(ctgBase=='A' or ctgBase=='C' or ctgBase=='T') match_flag = true;
				break;
			case 'K':
			case 'k':
				if(ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			case 'D':
			case 'd':
				if(ctgBase=='A' or ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			case 'B':
			case 'b':
				if(ctgBase=='C' or ctgBase=='G' or ctgBase=='T') match_flag = true;
				break;
			default: cerr << __func__ << ": unknown base: " << refBase << endl; exit(1);
		}
	}

	return match_flag;
}

bool isRegValid(reg_t *reg, int32_t min_size){
	bool flag = false;
	if(reg->startRefPos<=reg->endRefPos and reg->startLocalRefPos<=reg->endLocalRefPos and reg->startQueryPos<=reg->endQueryPos){
		if(reg->endRefPos-reg->startRefPos+1>=min_size or reg->endQueryPos-reg->startQueryPos+1>=min_size)
			flag = true;
	}
//	if(flag==false){
//		cout << "========= Invalid region: " << reg->chrname << ":" << reg->startRefPos << "-" << reg->endRefPos << ", local_Loc: " << reg->startLocalRefPos << "-" << reg->endLocalRefPos << ", query_Loc: " << reg->startQueryPos << "-" << reg->endQueryPos << endl;
//	}
	return flag;
}


void exchangeRegLoc(reg_t *reg){
	int64_t tmp;
	if(reg->startRefPos>reg->endRefPos){
		//cout << "startRefPos<-->endRefPos and startLocalRefPos<-->endLocalRefPos has been exchanged." << endl;
		tmp = reg->startRefPos;
		reg->startRefPos = reg->endRefPos;
		reg->endRefPos = tmp;
		tmp = reg->startLocalRefPos;
		reg->startLocalRefPos = reg->endLocalRefPos;
		reg->endLocalRefPos = tmp;
	}

	if(reg->startQueryPos>reg->endQueryPos){
		tmp = reg->startQueryPos;
		reg->startQueryPos = reg->endQueryPos;
		reg->endQueryPos = tmp;
	}
}

// BLAT alignment, and the output is in sim4 format
void blatAln(string &alnfilename, string &contigfilename, string &refseqfilename){
	string blat_cmd, out_opt;
	int i, ret, sleep_sec;

	out_opt = "-out=sim4 " + alnfilename;
	blat_cmd = "blat " + refseqfilename + " " + contigfilename + " " + out_opt + " > /dev/null 2>&1";

	//cout << "blat_cmd: " + blat_cmd << endl;

	for(i=0; i<3; i++){
		ret = system(blat_cmd.c_str());
		if(ret!=0){
			if(i<2){
				sleep_sec = (i + 1) * (i + 1) * 10;
				sleep(sleep_sec);
				cout << __func__ << ": retry aligning " << contigfilename << endl;
			}else{
				cerr << "Please run the correct blat command or check whether blat was correctly installed when aligning " << contigfilename <<"." << endl;
				exit(1);
			}
		}else break;
	}
}

// check whether the blat alignment is matched to the contig
bool isBlatAlnResultMatch(string &contigfilename, string &alnfilename){
	bool flag;
	vector<string> query_name_vec_ctg, query_name_vec_blat, line_vec, line_vec1, len_vec;
	vector<size_t> query_len_vec_blat;
	ifstream infile;
	string line, query_name_blat, query_name_ctg;
	size_t i, query_len_blat, query_len_ctg;
	int32_t query_loc;

	infile.open(alnfilename);
	if(!infile.is_open()){
		cerr << __func__ << "(), line=" << __LINE__ << ": cannot open file:" << alnfilename << endl;
		exit(1);
	}

	// load query names
	FastaSeqLoader fa_loader(contigfilename);
	query_name_vec_ctg = fa_loader.getFastaSeqNames();

	// get blat alignment information
	while(getline(infile, line)){
		if(line.size()){
			if(line.substr(0, 3).compare("seq")==0){  // query and subject titles
				line_vec = split(line, "=");
				line_vec1 = split(line_vec[1].substr(1), ",");
				if(line_vec[0][3]=='1') {
					// get the query name and length
					query_name_vec_blat.push_back(line_vec1.at(0));
					len_vec = split(line_vec1[line_vec1.size()-1].substr(1), " ");
					query_len_vec_blat.push_back(stoi(len_vec[0]));
				}
			}
		}
	}
	infile.close();

	// check the match items
	flag = true;
	for(i=0; i<query_name_vec_blat.size(); i++){
		query_name_blat = query_name_vec_blat.at(i);
		query_len_blat = query_len_vec_blat.at(i);

		// search the query from contig vector
		query_loc = getQueryNameLoc(query_name_blat, query_name_vec_ctg);
		if(query_loc!=-1){ // query name matched
			query_len_ctg = fa_loader.getFastaSeqLen(query_loc);
			if(query_len_blat!=query_len_ctg) {
				flag = false;
				break;
			}
		}else{ // query name unmatched
			flag = false;
			break;
		}
	}

	if(flag==false){
		cerr << __func__ << ", line=" << __LINE__ << ", it is not matched between the blat alignment and its corresponding contig, contigfile=" << contigfilename << ", alnfile=" << alnfilename << ", error!" << endl;
	}

	return flag;
}

// get query name location
int32_t getQueryNameLoc(string &query_name, vector<string> &query_name_vec){
	int32_t i, loc = -1, query_name_len;
	string query_name_tmp;
	query_name_len = query_name.size();
	if(query_name_vec.size()){
		for(i=query_name_vec.size()-1; i>=0; i--){
			query_name_tmp = query_name_vec.at(i).substr(0,query_name_len);
			if(query_name_tmp.compare(query_name)==0) { loc = i; break; }
		}
	}
	return loc;
}

// clean assemble temporary folders
void cleanPrevAssembledTmpDir(const string &assem_dir_str, const string &dir_prefix)
{
	DIR *dp;
	struct dirent *entry;
	struct stat statbuf;
	string cmd_str, dir_str;
	char *path;

	// get current work directory
	path = getcwd(NULL, 0);

	if((dp=opendir(assem_dir_str.c_str()))==NULL){
		cerr << "cannot open directory: " << assem_dir_str << endl;
		exit(1);
	}
	chdir(assem_dir_str.c_str());
	while((entry = readdir(dp)) != NULL){
		lstat(entry->d_name, &statbuf);
		if(S_ISDIR(statbuf.st_mode)){
			if(strcmp(".", entry->d_name) == 0 || strcmp("..", entry->d_name) == 0) continue;

			if(strlen(entry->d_name)>4){
				dir_str = entry->d_name;
				if(dir_str.substr(0, 4).compare(dir_prefix)==0){
					cmd_str = "rm -rf ";
					cmd_str += entry->d_name;
					system(cmd_str.c_str());
				}
			}
		}
	}
	closedir(dp);
	chdir(path);
	free(path);
}

// get call file header line for INDEL which starts with '#'
string getCallFileHeaderBed(){
	string header_line;
	header_line = "#Chr\tStart\tEnd\tSVType\tSVLen\tDupNum\tRef\tAlt";
	return header_line;
}

// get call file header line for INDEL which starts with '#'
string getCallFileHeaderBedpe(){
	string header_line;
	header_line = "#Chr1\tStart1\tEnd1\tChr2\tStart2\tEnd2\tSVType\tSVLen1\tSVLen2\tRef1\tAlt1\tRef2\tAlt2";
	return header_line;
}

assembleWork_opt* allocateAssemWorkOpt(string &chrname, string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, vector<reg_t*> &varVec, bool clip_reg_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){
	assembleWork_opt *assem_work_opt;
	assem_work_opt = new assembleWork_opt();
	assem_work_opt->chrname = chrname;
	assem_work_opt->readsfilename = readsfilename;
	assem_work_opt->contigfilename = contigfilename;
	assem_work_opt->refseqfilename = refseqfilename;
	assem_work_opt->tmpdir = tmpdir;
	assem_work_opt->clip_reg_flag = clip_reg_flag;

	// sub-regions
	assem_work_opt->var_array = (reg_t**)malloc(varVec.size()*sizeof(reg_t*));
	assem_work_opt->arr_size = varVec.size();
	for(size_t i=0; i<varVec.size(); i++) assem_work_opt->var_array[i] = varVec.at(i);

	// limit process regions
	assem_work_opt->limit_reg_process_flag = limit_reg_process_flag;
	assem_work_opt->limit_reg_array = NULL;
	assem_work_opt->limit_reg_array_size = 0;
	if(limit_reg_process_flag and limit_reg_vec.size()){
		assem_work_opt->limit_reg_array = (simpleReg_t**)malloc(limit_reg_vec.size()*sizeof(simpleReg_t*));
		assem_work_opt->limit_reg_array_size = limit_reg_vec.size();
		for(size_t i=0; i<limit_reg_vec.size(); i++) assem_work_opt->limit_reg_array[i] = limit_reg_vec.at(i);
	}

	return assem_work_opt;
}

void releaseAssemWorkOpt(assembleWork_opt *assem_work_opt){
	free(assem_work_opt->var_array); assem_work_opt->var_array = NULL;
	if(assem_work_opt->limit_reg_array_size>0) { free(assem_work_opt->limit_reg_array); assem_work_opt->limit_reg_array = NULL; assem_work_opt->limit_reg_array_size = 0; }
	delete assem_work_opt;
}

void destroyAssembleWorkOptVec(vector<assembleWork_opt*> &assem_work_vec){
	assembleWork_opt *assem_work_opt;
	for(size_t i=0; i<assem_work_vec.size(); i++){
		assem_work_opt = assem_work_vec.at(i);
		releaseAssemWorkOpt(assem_work_opt);
	}
	vector<assembleWork_opt*>().swap(assem_work_vec);
}

//time canu1.8 -p assembly -d out_1.8 genomeSize=30000 -pacbio-raw clipReg_reads_hs37d5_21480275-21480297.fq
void *doit_canu(void *arg) {
    char *cmd_job = (char *)arg;

    //usleep(random() % 100000); // to coerce job completion out of order


    cout << cmd_job << endl;
    system(cmd_job);

    free(arg);
    return NULL;
}

int test_canu(int n, vector<string> &cmd_vec){

    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 1);
    //string *ip;
    string str;

    // Dispatch jobs
    for (size_t i = 0; i < cmd_vec.size(); i++) {
        //int *ip = (int*)malloc(sizeof(*ip));
    	//int *ip = (int*)malloc(sizeof(int));
        //*ip = i;
    	str = cmd_vec.at(i);
    	char *ip = (char*)malloc((str.size()+1) * sizeof(*ip));
        strcpy(ip, str.c_str());
        hts_tpool_dispatch(p, q, doit_canu, ip);
    }

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    return 0;
}

// process single assemble work
void* processSingleAssembleWork(void *arg){
	assembleWork *assem_work = (assembleWork *)arg;
	assembleWork_opt *assem_work_opt = assem_work->assem_work_opt;
	vector<reg_t*> varVec;
	vector<simpleReg_t*> sub_limit_reg_vec;
	size_t i, num_done, num_work, num_work_per_ten_percent;
	double percentage;
	Time time;

//	pthread_mutex_lock(assem_work->p_mtx_assemble_reg_workDone_num);
//	cout << "assemble region [" << assem_work->work_id << "]: " << assem_work_opt->readsfilename << endl;
//	pthread_mutex_unlock(assem_work->p_mtx_assemble_reg_workDone_num);

	for(i=0; i<assem_work_opt->arr_size; i++) varVec.push_back(assem_work_opt->var_array[i]);
	for(i=0; i<assem_work_opt->limit_reg_array_size; i++) sub_limit_reg_vec.push_back(assem_work_opt->limit_reg_array[i]);

	performLocalAssembly(assem_work_opt->readsfilename, assem_work_opt->contigfilename, assem_work_opt->refseqfilename, assem_work_opt->tmpdir, assem_work->num_threads_per_assem_work, varVec, assem_work_opt->chrname, assem_work->inBamFile, assem_work->fai, *(assem_work->var_cand_file), assem_work->expected_cov_assemble, assem_work->delete_reads_flag, assem_work_opt->limit_reg_process_flag, sub_limit_reg_vec);

	// release memory
	sub_limit_reg_vec.clear();
	if(assem_work_opt->arr_size>0) { free(assem_work_opt->var_array); assem_work_opt->var_array = NULL; assem_work_opt->arr_size = 0; }
	if(assem_work_opt->limit_reg_array_size>0) { free(assem_work_opt->limit_reg_array); assem_work_opt->limit_reg_array = NULL; assem_work_opt->limit_reg_array_size = 0; }

	// output progress information
	pthread_mutex_lock(assem_work->p_mtx_assemble_reg_workDone_num);
	(*assem_work->p_assemble_reg_workDone_num) ++;
	num_done = *assem_work->p_assemble_reg_workDone_num;

	num_work = assem_work->num_work;
	num_work_per_ten_percent = assem_work->num_work_per_ten_percent;
	if(num_done==1 and num_done!=num_work){
		num_done = 0;
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}else if(num_done%num_work_per_ten_percent==0 or num_done==num_work){
		percentage = (double)num_done / num_work * 100;
		cout << "[" << time.getTime() << "]: processed regions: " << num_done << "/" << num_work << " (" << percentage << "%)" << endl;
	}
	pthread_mutex_unlock(assem_work->p_mtx_assemble_reg_workDone_num);

	delete (assembleWork *)arg;

	return NULL;
}

void performLocalAssembly(string &readsfilename, string &contigfilename, string &refseqfilename, string &tmpdir, size_t num_threads_per_assem_work, vector<reg_t*> &varVec, string &chrname, string &inBamFile, faidx_t *fai, ofstream &assembly_info_file, double expected_cov_assemble, bool delete_reads_flag, bool limit_reg_process_flag, vector<simpleReg_t*> &limit_reg_vec){

	LocalAssembly local_assembly(readsfilename, contigfilename, refseqfilename, tmpdir, num_threads_per_assem_work, varVec, chrname, inBamFile, fai, 0, expected_cov_assemble, delete_reads_flag);

	local_assembly.setLimitRegs(limit_reg_process_flag, limit_reg_vec);
	if(local_assembly.assem_success_preDone_flag==false){
		// extract the corresponding refseq from reference
		local_assembly.extractRefseq();

		// extract the reads data from BAM file
		local_assembly.extractReadsDataFromBAM();

		// local assembly using Canu
		local_assembly.localAssembleCanu();
	}

	// record assembly information
	local_assembly.recordAssemblyInfo(assembly_info_file);

	// empty the varVec
	varVec.clear();
}

// get old output directory name
string getOldOutDirname(string &filename, string &sub_work_dir){
	string old_dir = "";
	size_t pos = filename.find(sub_work_dir);
	if(pos!=filename.npos){
		if(pos==0) old_dir = "";
		else{
			old_dir = filename.substr(0, pos);
			old_dir = deleteTailPathChar(old_dir);
		}
	}else{
		cerr << __func__ << ", line=" << __LINE__ << ", filename=" << filename << ", cannot find old output directory for " << sub_work_dir << ", error!" << endl;
		exit(1);
	}
	return old_dir;
}

// get updated item file name
string getUpdatedItemFilename(string &filename, string &out_dir, string &old_out_dir){
	string new_filename;

	if(old_out_dir.size()==0)
		new_filename = out_dir + "/" + filename;
	else if(filename.at(old_out_dir.size())=='/')
		new_filename = out_dir + filename.substr(old_out_dir.size());
	else
		new_filename = out_dir + "/" + filename.substr(old_out_dir.size());

	return new_filename;
}

// get chromosome name by given file name
string getChrnameByFilename(string &filename){
	string chrname;
	int32_t start_idx, end_idx;

	if(filename.size()>0){
		start_idx = end_idx = -1;
		end_idx = filename.find_last_of("/") - 1;
		for(int32_t i=end_idx; i>=0; i--){
			if(filename.at(i)=='/'){
				start_idx = i + 1;
				break;
			}
		}
		if(start_idx!=-1 and end_idx!=-1)
			chrname = filename.substr(start_idx, end_idx-start_idx+1);
	}

	return chrname;
}

// delete the tail '/' path character
string deleteTailPathChar(string &dirname){
	if(dirname.size()>0 and dirname.at(dirname.size()-1)=='/'){
		return dirname.substr(0, dirname.size()-1);
	}else{
		return dirname;
	}
}

// get overlapped simple regions
vector<simpleReg_t*> getOverlappedSimpleRegsExt(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec, int32_t ext_size){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;
	int64_t new_startRegPos, new_endRegPos;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) or (simple_reg->startPos==-1 and simple_reg->endPos==-1)) flag = true;
				else{
					new_startRegPos = simple_reg->startPos - ext_size;
					if(new_startRegPos<1) new_startRegPos = 1;
					new_endRegPos = simple_reg->endPos + ext_size;
					if(isOverlappedPos(begPos, endPos, new_startRegPos, new_endRegPos)) flag = true;
				}

				if(flag) sub_simple_reg_vec.push_back(simple_reg);
			}
		}
	}

	return sub_simple_reg_vec;
}

// get overlapped simple regions
vector<simpleReg_t*> getOverlappedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) or (simple_reg->startPos==-1 and simple_reg->endPos==-1)) flag = true;
				else if(isOverlappedPos(begPos, endPos, simple_reg->startPos, simple_reg->endPos)) flag = true;

				if(flag) sub_simple_reg_vec.push_back(simple_reg);
			}
		}
	}

	return sub_simple_reg_vec;
}

// get fully contained simple regions
vector<simpleReg_t*> getFullyContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if((begPos==-1 and endPos==-1) and (simple_reg->startPos==-1 and simple_reg->endPos==-1)){ // chr, chr --> true
					flag = true;
				}else if(simple_reg->startPos==-1 and simple_reg->endPos==-1){ // chr:start-end, chr --> true
					flag = true;
				}else{ // CHR:START-END fully contained in 'simple_reg'
					if(isFullyContainedReg(chrname, begPos, endPos, simple_reg->chrname, simple_reg->startPos, simple_reg->endPos))
						flag = true;
				}
			}
			if(flag) sub_simple_reg_vec.push_back(simple_reg);
		}
	}
	return sub_simple_reg_vec;
}

// determine whether region1 is fully contained in region2
bool isFullyContainedReg(string &chrname1, int64_t begPos1,  int64_t endPos1, string &chrname2, int64_t startPos2, int64_t endPos2){
	bool flag = false;
	if(chrname1.compare(chrname2)==0 and begPos1>=startPos2 and endPos1<=endPos2) flag = true;
	return flag;
}

// get simple regions whose region contain the given start or end position
vector<simpleReg_t*> getPosContainedSimpleRegs(string &chrname, int64_t begPos, int64_t endPos, vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *simple_reg;
	vector<simpleReg_t*> sub_simple_reg_vec;
	bool flag;

	if(chrname.size()){
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			simple_reg = limit_reg_vec.at(i);
			flag = false;
			if(simple_reg->chrname.compare(chrname)==0){
				if(simple_reg->startPos==-1 and simple_reg->endPos==-1){ // pos, chr --> true
					flag = true;
				}else if((begPos>=simple_reg->startPos and begPos<=simple_reg->endPos) or (endPos>=simple_reg->startPos and endPos<=simple_reg->endPos)){
					flag = true;
				}
			}
			if(flag) sub_simple_reg_vec.push_back(simple_reg);
		}
	}

	return sub_simple_reg_vec;
}

// determine whether the start position is contained in region2
bool isPosContained(string &chrname1, int64_t begPos1, string &chrname2, int64_t startPos2, int64_t endPos2){
	bool flag = false;
	if(chrname1.compare(chrname2)==0 and begPos1>=startPos2 and begPos1<=endPos2) flag = true;
	return flag;
}

// allocate simple region node
simpleReg_t* allocateSimpleReg(string &simple_reg_str){
	simpleReg_t *simple_reg = NULL;
	vector<string> str_vec, pos_vec;
	string chrname_tmp;
	int64_t pos1, pos2;

	str_vec = split(simple_reg_str, ":");
	if(str_vec.size()){
		chrname_tmp = str_vec.at(0);
		if(str_vec.size()==1){
			pos1 = pos2 = -1;
		}else if(str_vec.size()==2){
			pos_vec = split(str_vec.at(1), "-");
			if(pos_vec.size()==2){
				pos1 = stoi(pos_vec.at(0));
				pos2 = stoi(pos_vec.at(1));
			}else goto fail;
		}else goto fail;
	}else goto fail;

	simple_reg = new simpleReg_t();
	simple_reg->chrname = chrname_tmp;
	simple_reg->startPos = pos1;
	simple_reg->endPos = pos2;

	return simple_reg;

fail:
	cout << "Skipped invalid region: " << simple_reg_str << endl;
	return NULL;
}

// free the memory of limited regions
void destroyLimitRegVector(vector<simpleReg_t*> &limit_reg_vec){
	vector<simpleReg_t*>::iterator s_reg;
	for(s_reg=limit_reg_vec.begin(); s_reg!=limit_reg_vec.end(); s_reg++)
		delete (*s_reg);   // free each node
	vector<simpleReg_t*>().swap(limit_reg_vec);
}

// print limit regions
void printLimitRegs(vector<simpleReg_t*> &limit_reg_vec, string &description){
	simpleReg_t *limit_reg;
	string limit_reg_str;

	if(limit_reg_vec.size()){
		cout << description << endl;
		for(size_t i=0; i<limit_reg_vec.size(); i++){
			limit_reg = limit_reg_vec.at(i);
			limit_reg_str = limit_reg->chrname;
			if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
			cout << "region [" << i << "]: " << limit_reg_str << endl;
		}
		cout << endl;
	}
}

string getLimitRegStr(vector<simpleReg_t*> &limit_reg_vec){
	simpleReg_t *limit_reg;
	string limit_reg_str;

	if(limit_reg_vec.size()){
		limit_reg = limit_reg_vec.at(0);
		limit_reg_str = limit_reg->chrname;
		if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
		for(size_t i=1; i<limit_reg_vec.size(); i++){
			limit_reg = limit_reg_vec.at(i);
			limit_reg_str = limit_reg->chrname;
			if(limit_reg->startPos!=-1 and limit_reg->endPos!=-1) limit_reg_str += ":" + to_string(limit_reg->startPos) + "-" + to_string(limit_reg->endPos);
			limit_reg_str += ":" + limit_reg_str;
		}
	}

	return limit_reg_str;
}

// get region by given file name
void getRegByFilename(simpleReg_t *reg, string &filename, string &pattern_str){
	string simple_filename, chr_pos_str, chr_str, pos_str;
	vector<string> pos_vec;
	size_t pos, pattern_size;

	reg->chrname = "";
	reg->startPos = reg->endPos = -1;
	pattern_size = pattern_str.size();

	pos = filename.find_last_of("/");			// path/pattern_chr_start-end
	simple_filename = filename.substr(pos+1);	// pattern_chr_start-end
	pos = simple_filename.find(pattern_str);
	chr_pos_str = simple_filename.substr(pos+pattern_size+1); // chr_start-end
	pos = chr_pos_str.find_last_of("_");
	chr_str = chr_pos_str.substr(0, pos);
	pos_str = chr_pos_str.substr(pos+1);
	pos_vec = split(pos_str, "-");

	reg->chrname = chr_str;
	reg->startPos = stoi(pos_vec.at(0));
	reg->endPos = stoi(pos_vec.at(1));
}

// extract simple regions by string formated regions
vector<simpleReg_t*> extractSimpleRegsByStr(string &regs_str){
	vector<simpleReg_t*> limit_reg_vec;
	simpleReg_t *simple_reg;
	vector<string> reg_vec, chr_pos_vec, pos_vec;
	string reg_str, chrname_reg;
	int64_t start_pos, end_pos;

	reg_vec = split(regs_str, ";");
	for(size_t i=0; i<reg_vec.size(); i++){
		reg_str = reg_vec.at(i);	// CHR or CHR:START-END

		chr_pos_vec = split(reg_str, ":");
		chrname_reg = chr_pos_vec.at(0);
		if(chr_pos_vec.size()==1)
			start_pos = end_pos = -1;
		else if(chr_pos_vec.size()==2){
			pos_vec = split(chr_pos_vec.at(1), "-");
			start_pos = stoi(pos_vec.at(0));
			end_pos = stoi(pos_vec.at(1));
		}else{
			cerr << __func__ << ", line=" << __LINE__ << ", invalid simple region: " << reg_str << ", error!" << endl;
			exit(1);
		}

		simple_reg = new simpleReg_t();
		simple_reg->chrname = chrname_reg;
		simple_reg->startPos = start_pos;
		simple_reg->endPos = end_pos;

		limit_reg_vec.push_back(simple_reg);
	}

	return limit_reg_vec;
}

void adjustVarLocByMismatchRegs(reg_t *reg, vector<mismatchReg_t*> &misReg_vec, int32_t start_aln_idx_var, int32_t end_aln_idx_var){
	int32_t i, max_reg_idx, maxValue, secMax_reg_idx, secMaxValue, left_reg_idx, right_reg_idx;
	mismatchReg_t *mis_reg, *mis_reg2;

	// get maximum and second maximum mismatch region
	max_reg_idx = secMax_reg_idx = -1;
	maxValue = secMaxValue = 0;
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		if(mis_reg->reg_size>maxValue){
			max_reg_idx = i;
			maxValue = mis_reg->reg_size;
		}
		if(mis_reg->reg_size>=secMaxValue){
			secMax_reg_idx = i;
			secMaxValue = mis_reg->reg_size;
		}
	}

	// check its neighbors
	if(max_reg_idx!=-1){
		left_reg_idx = right_reg_idx = -1;

		// left side extend
		mis_reg = misReg_vec.at(max_reg_idx);
		for(i=max_reg_idx-1; i>=0; i--){
			mis_reg2 = misReg_vec.at(i);
			if((mis_reg->start_aln_idx-mis_reg2->end_aln_idx<=MIN_ADJACENT_REG_DIST) or (mis_reg->start_aln_idx>=end_aln_idx_var and mis_reg2->end_aln_idx<=start_aln_idx_var)){ // distance less than threshold
				left_reg_idx = i;
				mis_reg = mis_reg2;
				mis_reg2 = NULL;
			}else break;
		}
		if(left_reg_idx==-1) left_reg_idx = max_reg_idx;

		// right side extend
		mis_reg = misReg_vec.at(max_reg_idx);
		for(i=max_reg_idx+1; i<(int32_t)misReg_vec.size(); i++){
			mis_reg2 = misReg_vec.at(i);
			if((mis_reg2->start_aln_idx-mis_reg->end_aln_idx<=MIN_ADJACENT_REG_DIST) or (mis_reg->end_aln_idx<=start_aln_idx_var and mis_reg2->start_aln_idx>=end_aln_idx_var)){ // distance less than threshold
				right_reg_idx = i;
				mis_reg = mis_reg2;
				mis_reg2 = NULL;
			}else break;
		}
		if(right_reg_idx==-1) right_reg_idx = max_reg_idx;

		mis_reg = misReg_vec.at(left_reg_idx);
		mis_reg2 = misReg_vec.at(right_reg_idx);
		reg->startRefPos = mis_reg->startRefPos;
		reg->endRefPos = mis_reg2->endRefPos;
		reg->startLocalRefPos = mis_reg->startLocalRefPos;
		reg->endLocalRefPos = mis_reg2->endLocalRefPos;
		reg->startQueryPos = mis_reg->startQueryPos;
		reg->endQueryPos = mis_reg2->endQueryPos;
	}
}

void releaseMismatchRegVec(vector<mismatchReg_t*> &misReg_vec){
	// release memory
	for(size_t i=0; i<misReg_vec.size(); i++) delete misReg_vec.at(i);
	vector<mismatchReg_t*>().swap(misReg_vec);
}

mismatchReg_t *getMismatchReg(int32_t aln_idx, vector<mismatchReg_t*> &misReg_vec){
	mismatchReg_t *mis_item_ret = NULL, *mis_item;
	for(size_t i=0; i<misReg_vec.size(); i++){
		mis_item = misReg_vec.at(i);
		if(mis_item->valid_flag){
			if(aln_idx>=mis_item->start_aln_idx and aln_idx<=mis_item->end_aln_idx){
				mis_item_ret = mis_item;
				break;
			}
		}
	}
	return mis_item_ret;
}

bool isPolymerSeq(string &seq){
	bool flag = true;
	for(size_t m=1; m<seq.size(); m++){
		if(seq.at(m-1)!=seq.at(m)){
			flag = false;
			break;
		}
	}
	return flag;
}

// create directory
void createDir(string &dirname){
	string dirname_tmp;
	vector<string> str_vec;
	int32_t ret;

	str_vec = split(dirname, "/");
	dirname_tmp = str_vec.at(0);
	for(size_t i=0; i<str_vec.size(); i++){
		ret = mkdir(dirname_tmp.c_str(), S_IRWXU | S_IROTH);  // create the output directory
		if(ret==-1 and errno!=EEXIST){
			cerr << __func__ << ", line=" << __LINE__ << ": cannot create directory: " << dirname_tmp << endl;
			exit(1);
		}
		if(i+1<str_vec.size()) dirname_tmp = dirname_tmp + "/" + str_vec.at(i+1);
	}
}

// get the number of high ratio indel bases
vector<double> getTotalHighIndelClipRatioBaseNum(Base *regBaseArr, int64_t arr_size){
	int32_t i, indel_num, clip_num, total_cov;
	double ratio, total, total2;
	Base *base;
	vector<double> base_num_vec;

	total = total2 = 0;
	for(i=0; i<arr_size; i++){
		base  = regBaseArr + i;
		indel_num = base->getTotalIndelNum();
		clip_num = base->getTotalClipNum();
		total_cov = base->getTotalCovNum();
		ratio = (double)(indel_num + clip_num) / total_cov;
		if(ratio>=HIGH_INDEL_CLIP_RATIO_THRES) total ++;
		if(ratio>=SECOND_INDEL_CLIP_RATIO_THRES) total2 ++;
	}

	ratio = (double)total2 / arr_size * 100;

	base_num_vec.push_back(total);
	base_num_vec.push_back(ratio);

	return base_num_vec;
}

vector<mismatchReg_t*> getMismatchRegVec(localAln_t *local_aln){
	vector<mismatchReg_t*> misReg_vec;	// all the mismatch regions including gap regions
	int64_t ref_pos, local_ref_pos, query_pos, tmp, startMisRefPos, endMisRefPos, startMisLocalRefPos, endMisLocalRefPos, startMisQueryPos, endMisQueryPos;
	int32_t i, j, start_aln_idx, start_check_idx, end_check_idx, mis_reg_size, begin_mismatch_aln_idx, end_mismatch_aln_idx;
	string ctgseq_aln, midseq_aln, refseq_aln;
	mismatchReg_t *mis_reg, *mis_reg2;
	reg_t *reg;
	bool gap_flag, extend_flag, flag;
	int16_t mis_loc_flag;	// -1 for unused, 0 for gap in query, 1 for gap in reference, 2 for mismatch

	if(local_aln->reg==NULL or local_aln->start_aln_idx_var==-1 or local_aln->end_aln_idx_var==-1) return misReg_vec;

	ctgseq_aln = local_aln->alignResultVec[0];
	midseq_aln = local_aln->alignResultVec[1];
	refseq_aln = local_aln->alignResultVec[2];

	// compute refPos and queryPos at start_aln_idx
	start_check_idx = local_aln->start_aln_idx_var - EXT_SIZE_CHK_VAR_LOC;
	end_check_idx = local_aln->end_aln_idx_var + EXT_SIZE_CHK_VAR_LOC;
	if(start_check_idx<0) start_check_idx = 0;
	if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
	extend_flag = false;

	reg = local_aln->reg;
	ref_pos = reg->startRefPos;
	local_ref_pos = reg->startLocalRefPos;
	query_pos = reg->startQueryPos;

	// initialize the start positions
	start_aln_idx = local_aln->start_aln_idx_var;
	if(ctgseq_aln[start_aln_idx]=='-'){ // gap in query
		query_pos --;
		while(ctgseq_aln[start_aln_idx]=='-' and start_aln_idx>0){
			ref_pos --;
			local_ref_pos --;
			start_aln_idx --;
		}
	}else if(refseq_aln[start_aln_idx]=='-'){ // gap in reference
		while(refseq_aln[start_aln_idx]=='-' and start_aln_idx>0){
			query_pos --;
			start_aln_idx --;
		}
		ref_pos --;
		local_ref_pos --;
	}

	for(i=start_aln_idx; i>=start_check_idx; i--){
		if(midseq_aln[i]==' '){ // mismatch including gap
			if(ctgseq_aln[i]=='-'){ // gap in query
				ref_pos --;
				local_ref_pos --;
			}else if(refseq_aln[i]=='-'){ // gap in reference
				query_pos --;
			}else{ // mismatch
				query_pos --;
				ref_pos --;
				local_ref_pos --;
			}

			if(i==end_check_idx and extend_flag==false){
				end_check_idx -= EXT_SIZE_CHK_VAR_LOC;
				if(end_check_idx<0) end_check_idx = 0;
				extend_flag = true;
				//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
			}
		}else{ // match
			if(extend_flag) break;

			query_pos --;
			ref_pos --;
			local_ref_pos --;
		}
	}
	if(i>=start_check_idx) {
		start_check_idx = i;
	}else{ // adjust to correct value
		start_check_idx = i + 1;
		query_pos ++;
		ref_pos ++;
		local_ref_pos ++;
	}

	startMisRefPos = endMisRefPos = startMisLocalRefPos = endMisLocalRefPos = startMisQueryPos = endMisQueryPos = -1;
	begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
	mis_reg_size = 0;
	gap_flag = extend_flag = false;
	mis_loc_flag = -1;
	for(i=start_check_idx; i<=end_check_idx; i++){
		if(midseq_aln[i]==' '){ // mismatch including gap
			if(begin_mismatch_aln_idx==-1){
				begin_mismatch_aln_idx = i;
				startMisRefPos = ref_pos;
				startMisLocalRefPos = local_ref_pos;
				startMisQueryPos = query_pos;

				if(ctgseq_aln[i]=='-'){ // gap in query
					mis_loc_flag = 0;
				}else if(refseq_aln[i]=='-'){ // gap in reference
					mis_loc_flag = 1;
				}else{ // mismatch
					mis_loc_flag = 2;
				}
				mis_reg_size = 1;
			}else{
				mis_reg_size ++;
			}

			if(ctgseq_aln[i]=='-'){ // gap in query
				ref_pos ++;
				local_ref_pos ++;
				gap_flag = true;
			}else if(refseq_aln[i]=='-'){ // gap in reference
				query_pos ++;
				gap_flag = true;
			}else{ // mismatch
				query_pos ++;
				ref_pos ++;
				local_ref_pos ++;
			}

			if(i==end_check_idx){
				if(extend_flag==false){
					end_check_idx += EXT_SIZE_CHK_VAR_LOC;
					if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
					extend_flag = true;
					//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
				}else if(begin_mismatch_aln_idx!=-1){ // mismatch region not closed
					end_check_idx += EXT_SIZE_CHK_VAR_LOC;
					if(end_check_idx>(int32_t)midseq_aln.size()-1) end_check_idx = midseq_aln.size() - 1;
					//cout << "line=" << __LINE__ << ", end_check_idx=" << end_check_idx << endl;
				}
			}
		}else{ // match
			if(begin_mismatch_aln_idx!=-1){
				end_mismatch_aln_idx = i - 1;
				switch(mis_loc_flag){
					case 0:  // gap in query
						endMisRefPos = ref_pos - 1;
						endMisLocalRefPos = local_ref_pos - 1;
						endMisQueryPos = query_pos;
						break;
					case 1:  // gap in reference
						endMisRefPos = ref_pos;
						endMisLocalRefPos = local_ref_pos;
						endMisQueryPos = query_pos - 1;
						break;
					case 2:  // mismatch
						endMisRefPos = ref_pos - 1;
						endMisLocalRefPos = local_ref_pos - 1;
						endMisQueryPos = query_pos - 1;
						break;
					default: // error
						cerr << __func__ << ", line=" << __LINE__ << ", invalid mis_loc_flag: " << mis_loc_flag << ", error!" << endl;
						exit(1);
				}

				mis_reg = new mismatchReg_t();
				mis_reg->start_aln_idx = begin_mismatch_aln_idx;
				mis_reg->end_aln_idx = end_mismatch_aln_idx;
				mis_reg->reg_size = mis_reg_size;
				mis_reg->startRefPos = startMisRefPos;
				mis_reg->endRefPos = endMisRefPos;
				mis_reg->startLocalRefPos = startMisLocalRefPos;
				mis_reg->endLocalRefPos = endMisLocalRefPos;
				mis_reg->startQueryPos = startMisQueryPos;
				mis_reg->endQueryPos = endMisQueryPos;
				mis_reg->valid_flag = true;
				mis_reg->gap_flag = gap_flag;
				misReg_vec.push_back(mis_reg);

				begin_mismatch_aln_idx = end_mismatch_aln_idx = -1;
				startMisRefPos = endMisRefPos = startMisLocalRefPos = endMisLocalRefPos = startMisQueryPos = endMisQueryPos = -1;
				mis_reg_size = 0;
				mis_loc_flag = -1;

				if(extend_flag) break;
			}
			if(extend_flag) break;

			query_pos ++;
			ref_pos ++;
			local_ref_pos ++;
			gap_flag = false;
		}
	}

	// sort
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		for(j=i+1; j<(int32_t)misReg_vec.size(); j++){
			mis_reg2 = misReg_vec.at(j);
			if(mis_reg->start_aln_idx>mis_reg2->start_aln_idx){ // exchange
				tmp = mis_reg->start_aln_idx; mis_reg->start_aln_idx = mis_reg2->start_aln_idx; mis_reg2->start_aln_idx = tmp;
				tmp = mis_reg->end_aln_idx; mis_reg->end_aln_idx = mis_reg2->end_aln_idx; mis_reg2->end_aln_idx = tmp;
				tmp = mis_reg->reg_size; mis_reg->reg_size = mis_reg2->reg_size; mis_reg2->reg_size = tmp;
				tmp = mis_reg->startRefPos; mis_reg->startRefPos = mis_reg2->startRefPos; mis_reg2->startRefPos = tmp;
				tmp = mis_reg->endRefPos; mis_reg->endRefPos = mis_reg2->endRefPos; mis_reg2->endRefPos = tmp;
				tmp = mis_reg->startLocalRefPos; mis_reg->startLocalRefPos = mis_reg2->startLocalRefPos; mis_reg2->startLocalRefPos = tmp;
				tmp = mis_reg->endLocalRefPos; mis_reg->endLocalRefPos = mis_reg2->endLocalRefPos; mis_reg2->endLocalRefPos = tmp;
				tmp = mis_reg->startQueryPos; mis_reg->startQueryPos = mis_reg2->startQueryPos; mis_reg2->startQueryPos = tmp;
				tmp = mis_reg->endQueryPos; mis_reg->endQueryPos = mis_reg2->endQueryPos; mis_reg2->endQueryPos = tmp;
				flag = mis_reg->valid_flag; mis_reg->valid_flag = mis_reg2->valid_flag; mis_reg2->valid_flag = flag;
				flag = mis_reg->gap_flag; mis_reg->gap_flag = mis_reg2->gap_flag; mis_reg2->gap_flag = flag;
			}
		}
	}

	return misReg_vec;
}

void removeShortPolymerMismatchRegItems(localAln_t *local_aln, vector<mismatchReg_t*> &misReg_vec, string &inBamFile, faidx_t *fai){
	size_t m, start_ref_pos, end_ref_pos, chrlen_tmp;
	int32_t i, j, start_check_idx, end_check_idx;
	string chrname_tmp, ctgseq_aln, refseq_aln, query_substr, subject_substr, substr, substr_gap;
	mismatchReg_t *mis_reg, *mis_reg0, *mis_reg2, *tmp_mis_reg;
	char ch_left, ch_right;
	bool flag, large_neighbor_dist_flag, high_con_ratio_flag, disagree_flag;
	vector<bam1_t*> alnDataVector;
	Base *start_base, *end_base;

	if(misReg_vec.empty()) return;

	if(local_aln->reg) chrname_tmp = local_aln->reg->chrname;
	else if(local_aln->cand_reg) chrname_tmp = local_aln->cand_reg->chrname;
	else return;
	chrlen_tmp = faidx_seq_len(fai, chrname_tmp.c_str()); // get the reference length

	start_ref_pos = misReg_vec.at(0)->startRefPos - EXT_SIZE_CHK_VAR_LOC;
	end_ref_pos = misReg_vec.at(misReg_vec.size()-1)->endRefPos + EXT_SIZE_CHK_VAR_LOC;
	if(start_ref_pos<1) start_ref_pos = 1;
	if(end_ref_pos>chrlen_tmp) end_ref_pos = chrlen_tmp;

	// load align data
	alnDataLoader data_loader(chrname_tmp, start_ref_pos, end_ref_pos, inBamFile);
	data_loader.loadAlnData(alnDataVector);

	// load coverage
	covLoader cov_loader(chrname_tmp, start_ref_pos, end_ref_pos, fai);
	Base *baseArray = cov_loader.initBaseArray();
	cov_loader.generateBaseCoverage(baseArray, alnDataVector);

	ctgseq_aln = local_aln->alignResultVec[0];
	refseq_aln = local_aln->alignResultVec[2];

	// filter out short (e.g. <= 1 bp) mismatched polymer items
	for(i=0; i<(int32_t)misReg_vec.size(); i++){
		mis_reg = misReg_vec.at(i);
		start_check_idx = mis_reg->start_aln_idx;
		end_check_idx = mis_reg->end_aln_idx;

		// determine whether the base contains many insertions or deletions
		start_base = baseArray + mis_reg->startRefPos - start_ref_pos;
		end_base = baseArray + mis_reg->endRefPos - start_ref_pos;
		high_con_ratio_flag = disagree_flag = false;

		if(((start_base->max_con_type==BAM_CINS and start_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_INS_RATIO) or (start_base->max_con_type==BAM_CDEL and start_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_DEL_RATIO)) or ((end_base->max_con_type==BAM_CINS and end_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_INS_RATIO) or (end_base->max_con_type==BAM_CDEL and end_base->maxConIndelEventRatio>=MIN_HIGH_CONSENSUS_DEL_RATIO)))
			high_con_ratio_flag = true;
		if(start_base->isDisagreeBase() or end_base->isDisagreeBase())
			disagree_flag = true;

		if(high_con_ratio_flag or disagree_flag) continue; // skip below operations

		if(mis_reg->gap_flag){
			query_substr = ctgseq_aln.substr(start_check_idx, end_check_idx-start_check_idx+1);
			subject_substr = refseq_aln.substr(start_check_idx, end_check_idx-start_check_idx+1);

			// check polymer
			if(query_substr.find("-")!=string::npos){ // gap in query
				substr_gap = query_substr;
				substr = subject_substr;
				ch_left = refseq_aln.at(start_check_idx-1);
				ch_right = refseq_aln.at(end_check_idx+1);
			}else{ // gap in subject
				substr = query_substr;
				substr_gap = subject_substr;
				ch_left = ctgseq_aln.at(start_check_idx-1);
				ch_right = ctgseq_aln.at(end_check_idx+1);
			}

			flag = isPolymerSeq(substr);
			if(flag){ // check
				if(high_con_ratio_flag==false and disagree_flag==false and (substr.at(0)==ch_left or substr.at(substr.size()-1)==ch_right) and substr.size()<MIN_VALID_POLYMER_SIZE)
					mis_reg->valid_flag = false;
			}
		}

		// further to check the distance of its neighbors
		large_neighbor_dist_flag = false;
		if(mis_reg->valid_flag and mis_reg->reg_size<=MAX_SHORT_MIS_REG_SIZE){
			if(i==0){ // first item
				mis_reg2 = NULL;
				for(m=i+1; m<misReg_vec.size(); m++){ // search right side
					tmp_mis_reg = misReg_vec.at(m);
					if(tmp_mis_reg->valid_flag){
						mis_reg2 = tmp_mis_reg;
						break;
					}
				}
				if(mis_reg2==NULL or (mis_reg2->start_aln_idx-mis_reg->end_aln_idx>MIN_ADJACENT_REG_DIST))
					large_neighbor_dist_flag = true;
			}else if(i==(int32_t)misReg_vec.size()-1){ // last item
				mis_reg0 = NULL;
				for(j=i-1; j>=0; j--){ // search left side
					tmp_mis_reg = misReg_vec.at(j);
					if(tmp_mis_reg->valid_flag){
						mis_reg0 = tmp_mis_reg;
						break;
					}
				}
				if(mis_reg0==NULL or (mis_reg->start_aln_idx-mis_reg0->end_aln_idx>MIN_ADJACENT_REG_DIST))
					large_neighbor_dist_flag = true;
			}else{ // middle items
				mis_reg0 = NULL;
				for(j=i-1; j>=0; j--){ // search left side
					tmp_mis_reg = misReg_vec.at(j);
					if(tmp_mis_reg->valid_flag){
						mis_reg0 = tmp_mis_reg;
						break;
					}
				}
				mis_reg2 = NULL;
				for(m=i+1; m<misReg_vec.size(); m++){ // search right side
					tmp_mis_reg = misReg_vec.at(m);
					if(tmp_mis_reg->valid_flag){
						mis_reg2 = tmp_mis_reg;
						break;
					}
				}
				if((mis_reg0==NULL or (mis_reg->start_aln_idx-mis_reg0->end_aln_idx>MIN_ADJACENT_REG_DIST)) and (mis_reg2==NULL or (mis_reg2->start_aln_idx-mis_reg->end_aln_idx>MIN_ADJACENT_REG_DIST)))
					large_neighbor_dist_flag = true;
			}

			if(large_neighbor_dist_flag==true and high_con_ratio_flag==false and disagree_flag==false)
				mis_reg->valid_flag = false;
		}
	}

	// remove invalid items
	for(m=0; m<misReg_vec.size(); ){
		mis_reg = misReg_vec.at(m);
		if(mis_reg->valid_flag==false){
			delete mis_reg;
			misReg_vec.erase(misReg_vec.begin()+m);
		}else m++;
	}

	// release memory
	data_loader.freeAlnData(alnDataVector);
	cov_loader.freeBaseArray(baseArray);
}


Time::Time() {
	timestamp = start = end = start_all = end_all = 0;
	time_tm = NULL;
	time(&start_all);
}

Time::~Time() {
}

string Time::getTime(){
	string time_str;
	time(&timestamp);	// get the Unix time stamp
	time_tm = localtime(&timestamp); // convert the Unix time to Year-Mon-Day...
	time_str = ctime(&timestamp);
	time_str.erase(time_str.find("\n"), 1);
	return time_str;
}

void Time::printTime(){
	cout << getTime() << endl;
}

void Time::setStartTime(){
	time(&start);
}

void Time::printSubCmdElapsedTime(){
	double cost_sec, cost_hour, cost_min, cost_day;

	if(start!=0){
		time(&end);
		cost_sec = difftime(end, start);
		cost_min = cost_sec / 60;
		cost_hour = cost_min / 60;
		cost_day = cost_hour / 24;

		start = end = 0; // reset time

		if(cost_day>1)
			cout << "Sub-command elapsed time: " << cost_day << " days = " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_hour>1)
			cout << "Sub-command elapsed time: " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_min>1)
			cout << "Sub-command elapsed time: " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else
			cout << "Sub-command elapsed time: " << cost_sec << " seconds" << endl << endl;
	}else{
		cerr << "Please set the start time before calculating the running time" << endl;
		exit(1);
	}
}

void Time::printOverallElapsedTime(){
	double cost_sec, cost_hour, cost_min, cost_day;

	if(start_all!=0){
		time(&end_all);
		cost_sec = difftime(end_all, start_all);
		cost_min = cost_sec / 60;
		cost_hour = cost_min / 60;
		cost_day = cost_hour / 24;

		if(cost_day>1)
			cout << "Overall elapsed time: " << cost_day << " days = " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_hour>1)
			cout << "Overall elapsed time: " << cost_hour << " hours = " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else if(cost_min>1)
			cout << "Overall elapsed time: " << cost_min << " minutes = " << cost_sec << " seconds" << endl << endl;
		else
			cout << "Overall elapsed time: " << cost_sec << " seconds" << endl << endl;
	}else{
		cerr << "Please set the overall start time before calculating the overall running time" << endl;
		exit(1);
	}
}

// Merge continuous regions
void mergeContinuousReg(vector<reg_t*> &regVector){
	size_t i, j;
	reg_t *reg1, *reg2;

	for(i=0; i<regVector.size(); i++){
		reg1 = regVector.at(i);
		for(j=i+1; j<regVector.size(); ){
			reg2 = regVector.at(j);
			if(isContinuousReg(reg1, reg2)){
				updateContinuousReg(reg1, reg2);
				delete reg2;
				regVector.erase(regVector.begin()+j);
			}else j++;
		}
	}
	regVector.shrink_to_fit();
}

// determine whether two regions are Continuous
bool isContinuousReg(reg_t* reg1, reg_t* reg2){
	bool flag = false;
	if(reg1->chrname.compare(reg2->chrname)==0){
		if(isContinuousPos(reg1->startRefPos, reg1->endRefPos, reg2->startRefPos, reg2->endRefPos))
			flag = true;
	}
	return flag;
}

// determine whether the given positions of two regions are Continuous
bool isContinuousPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	bool flag = false;
	if(startPos2-endPos1>=0 and startPos2-endPos1<=MISASSEMBLY_CONTINUOUS_REGION_THRES)
		flag = true;
	return flag;
}

void updateContinuousReg(reg_t* reg1, reg_t* reg2){
	size_t new_startRefPos, new_endRefPos;
	if(reg1->startRefPos<=reg2->startRefPos) {
		new_startRefPos = reg1->startRefPos;
	}else{
		new_startRefPos = reg2->startRefPos;
	}
	if(reg1->endRefPos>=reg2->endRefPos){
		new_endRefPos = reg1->endRefPos;
	}else{
		new_endRefPos = reg2->endRefPos;
	}
	reg1->startRefPos = new_startRefPos;
	reg1->endRefPos = new_endRefPos;
}

// Remove overlapped regions
void removeOverlappedReg(vector<reg_t*> &regVector){
	size_t i, j;
	reg_t *reg1, *reg2;

	for(i=0; i<regVector.size(); i++){
		reg1 = regVector.at(i);
		for(j=i+1; j<regVector.size(); ){
			reg2 = regVector.at(j);
			if(isOverlappedReg(reg1, reg2)){
				updateOverlappedReg(reg1, reg2);
				delete reg2;
				regVector.erase(regVector.begin()+j);
			}else j++;
		}
	}
	regVector.shrink_to_fit();
}

void updateOverlappedReg(reg_t* reg1, reg_t* reg2){
	size_t new_startRefPos, new_endRefPos;
	if(reg1->startRefPos<=reg2->startRefPos) {
		new_startRefPos = reg1->startRefPos;
	}else{
		new_startRefPos = reg2->startRefPos;
	}
	if(reg1->endRefPos>=reg2->endRefPos){
		new_endRefPos = reg1->endRefPos;
	}else{
		new_endRefPos = reg2->endRefPos;
	}
	reg1->startRefPos = new_startRefPos;
	reg1->endRefPos = new_endRefPos;
}
