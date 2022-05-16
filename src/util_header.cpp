#include "util_header.h"

const char* getTargetName(const bam_hdr_t *h, int tid){
#ifdef HTS_VERSION
#if HTS_VERSION >= 101000
	return sam_hdr_tid2name(h, tid);
#endif
#else
	if (!h || tid < 0)
		return 0;

	if (tid < h->n_targets)
		return h->target_name[tid];
#endif
	return NULL;
}

int64_t getTargetLen(const bam_hdr_t *h, int tid){
#ifdef HTS_VERSION
#if HTS_VERSION >= 101000
	return sam_hdr_tid2len(h, tid);
#endif
#else
	if (!h || tid < 0)
		return 0;

	if (tid < h->n_targets) {
		if (h->target_len[tid] < UINT32_MAX || !h->sdict) {
			return h->target_len[tid];
		}
	}
#endif
	return 0;
}
