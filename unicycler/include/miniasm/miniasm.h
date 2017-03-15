#ifndef MINIASM_H
#define MINIASM_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <string>

#pragma GCC diagnostic ignored "-Wsign-compare"

#include "miniasm/sdict.h"
#include "miniasm/asg.h"

extern int ma_verbose;

typedef struct {
	int min_span;
	int min_match;
	int min_dp;
	float min_iden;

	int max_hang;
	int min_ovlp;
	float int_frac;

	int gap_fuzz;
	int n_rounds;
	int bub_dist;
	int max_ext;
	float min_ovlp_drop_ratio, max_ovlp_drop_ratio, final_ovlp_drop_ratio;
} ma_opt_t;

typedef struct {
	uint64_t qns;
	uint32_t qe, tn, ts, te;
	uint32_t ml:31, rev:1;
	uint32_t bl:31, del:1;
} ma_hit_t;

typedef struct { size_t n, m; ma_hit_t *a; } ma_hit_v;

typedef struct {
	uint32_t s:31, del:1, e;
} ma_sub_t;

typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;

typedef struct { size_t n, m; ma_utg_t *a; } ma_utg_v;

typedef struct {
	ma_utg_v u;
	asg_t *g;
} ma_ug_t;

void ma_opt_init(ma_opt_t *opt);
sdict_t *prefilter_contained_reads(const char *fn, int min_span, int min_match, int max_hang, float int_frac);
ma_hit_t *read_hits_file(const char *fn, int min_span, int min_match, sdict_t *d, size_t *n, int bi_dir, const sdict_t *excl);
ma_sub_t *filter_reads_using_depth(int min_dp, float min_iden, int end_clip, size_t n, const ma_hit_t *a, const sdict_t *read_dict);
size_t filter_hits_using_span(const ma_sub_t *reg, int min_span, size_t n, ma_hit_t *a);
size_t filter_hits_using_overhang(const ma_sub_t *sub, int max_hang, int min_ovlp, size_t n, ma_hit_t *a, float *cov);
void merge_subreads(size_t n_sub, ma_sub_t *a, const ma_sub_t *b);
size_t remove_chimeric_reads(int max_hang, int min_dp, size_t n, const ma_hit_t *a, const sdict_t *d, ma_sub_t *sub, std::string chimeric_read_list);
size_t remove_contained_reads(int max_hang, float int_frac, int min_ovlp, sdict_t *d, ma_sub_t *sub, size_t n, ma_hit_t *a, std::string contained_read_list);
bool is_read_illumina_contig(const sdict_t *read_dict, int id);
void ma_hit_mark_unused(sdict_t *read_dict, int n, const ma_hit_t *a);

asg_t *make_string_graph(int max_hang, float int_frac, int min_ovlp, sdict_t const *d, ma_sub_t const *sub, unsigned long n_hits, ma_hit_t const *hit);
void save_string_graph(const asg_t *g, const sdict_t *d, const ma_sub_t *sub, std::string graph_filename, const char *reads_filename);
ma_ug_t *make_unitig_graph(asg_t *g);
int generate_unitig_seqs(ma_ug_t *g, const sdict_t *d, const ma_sub_t *sub, const char *fn);
void save_unitig_graph(const ma_ug_t *ug, const sdict_t *d, const ma_sub_t *sub, std::string graph_filename);
void destroy_unitig_graph(ma_ug_t *ug);

#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

static inline int ma_hit2arc(const ma_hit_t *h, int ql, int tl, int max_hang, float int_frac, int min_ovlp, asg_arc_t *p)
{
	int32_t tl5, tl3, ext5, ext3, qs = (int32_t)h->qns;
	uint32_t u, v, l; // u: query end; v: target end; l: length from u to v

    // tl5: 5'-end overhang (on the query strand); tl3: similar
	if (h->rev)
        tl5 = tl - h->te, tl3 = h->ts;
	else
        tl5 = h->ts, tl3 = tl - h->te;

	ext5 = qs < tl5? qs : tl5;
	ext3 = ql - h->qe < tl3? ql - h->qe : tl3;

	if (ext5 > max_hang || ext3 > max_hang || h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac)
		return MA_HT_INT;
	if (qs <= tl5 && ql - h->qe <= tl3) // query contained
        return MA_HT_QCONT;
	else if (qs >= tl5 && ql - h->qe >= tl3) // target contained
        return MA_HT_TCONT;
	else if (qs > tl5)
        u = 0, v = !!h->rev, l = qs - tl5;
	else
        u = 1, v = !h->rev, l = (ql - h->qe) - tl3;

	if (h->qe - qs + ext5 + ext3 < min_ovlp || h->te - h->ts + ext5 + ext3 < min_ovlp) // short overlap
        return MA_HT_SHORT_OVLP;

	u |= h->qns>>32<<1, v |= h->tn<<1;
	p->ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0, p->ml = h->ml, p->mr = (float)h->ml / h->bl;
	return l;
}

#endif
