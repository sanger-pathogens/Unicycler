#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <limits>

#pragma GCC diagnostic ignored "-Wvla-extension"

#include "miniasm/sdict.h"
#include "miniasm/paf.h"
#include "miniasm/kvec.h"
#include "miniasm/sys.h"
#include "miniasm/miniasm.h"
#include "miniasm/ksort.h"

#define ma_hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, ma_hit_t, ma_hit_key, 8)

KSORT_INIT_GENERIC(uint32_t)

using namespace std;

typedef kvec_t(uint32_t) uint32_v;

void ma_hit_sort(size_t n, ma_hit_t *a)
{
    radix_sort_hit(a, a + n);
}

void ma_hit_mark_unused(sdict_t *read_dict, int n, const ma_hit_t *a)
{
    size_t i;
    for (i = 0; i < read_dict->n_seq; ++i)
        read_dict->seq[i].aux = 0;
    for (i = 0; i < n; ++i)
        read_dict->seq[a[i].qns>>32].aux = read_dict->seq[a[i].tn].aux = 1;
    for (i = 0; i < read_dict->n_seq; ++i) {
        sd_seq_t *s = &read_dict->seq[i];
        if (!s->aux) s->del = 1;
        else s->aux = 0;
    }
}

sdict_t *prefilter_contained_reads(const char *fn, int min_span, int min_match, int max_hang, float int_frac)
{
    paf_file_t *fp;
    paf_rec_t r;
    sdict_t *d;

    fp = paf_open(fn);
    d = init_seq_dict();
    while (paf_read(fp, &r) >= 0) {
        int l5, l3;
        if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue;
        l5 = r.rev? r.tl - r.te : r.ts;
        l3 = r.rev? r.ts : r.tl - r.te;
        if (r.ql>>1 > r.tl) {
            if (l5 > max_hang>>2 || l3 > max_hang>>2 || r.te - r.ts < r.tl * int_frac) continue; // internal match
            if ((int)r.qs - l5 > max_hang<<1 && (int)(r.ql - r.qe) - l3 > max_hang<<1)
                sd_put(d, r.tn, r.tl);
        } else if (r.ql < r.tl>>1) {
            if (r.qs > max_hang>>2 || r.ql - r.qe > max_hang>>2 || r.qe - r.qs < r.ql * int_frac) continue; // internal
            if (l5 - (int)r.qs > max_hang<<1 && l3 - (int)(r.ql - r.qe) > max_hang<<1)
                sd_put(d, r.qn, r.ql);
        }
    }
    paf_close(fp);
//    fprintf(stderr, "[M::%s::%s] dropped %d contained reads\n", __func__, sys_timestamp(), d->n_seq);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] dropped " << d->n_seq << " contained reads\n";
    return d;
}

ma_hit_t *read_hits_file(const char *fn, int min_span, int min_match, sdict_t *read_dict, size_t *n, int bi_dir, const sdict_t *excl)
{
    paf_file_t *fp;
    paf_rec_t r;
    ma_hit_v h = {0,0,0};
    size_t i, tot = 0, tot_len = 0;

    fp = paf_open(fn);
    while (paf_read(fp, &r) >= 0) {
        ma_hit_t *p;
        ++tot;
        if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue;
        if (excl && (sd_get(excl, r.qn) >= 0 || sd_get(excl, r.tn) >= 0)) continue;
        kv_pushp(ma_hit_t, h, &p);
        p->qns = (uint64_t)sd_put(read_dict, r.qn, r.ql)<<32 | r.qs;
        p->qe = r.qe;
        p->tn = sd_put(read_dict, r.tn, r.tl);
        p->ts = r.ts, p->te = r.te, p->rev = r.rev, p->ml = r.ml, p->bl = r.bl;
        if (bi_dir && p->qns>>32 != p->tn) {
            kv_pushp(ma_hit_t, h, &p);
            p->qns = (uint64_t)sd_put(read_dict, r.tn, r.tl)<<32 | r.ts;
            p->qe = r.te;
            p->tn = sd_put(read_dict, r.qn, r.ql);
            p->ts = r.qs, p->te = r.qe, p->rev = r.rev, p->ml = r.ml, p->bl = r.bl;
        }
    }
    paf_close(fp);
    for (i = 0; i < read_dict->n_seq; ++i)
        tot_len += read_dict->seq[i].len;
//    fprintf(stderr, "[M::%s::%s] read %ld hits; stored %ld hits and %d sequences (%ld bp)\n", __func__, sys_timestamp(), tot, h.n, read_dict->n_seq, tot_len);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] read " << tot << " hits; stored " << h.n << " hits and " << read_dict->n_seq << " sequences (" << tot_len << " bp)\n";
    ma_hit_sort(h.n, h.a);
    *n = h.n;
    return h.a;
}


string get_read_name(const sdict_t *read_dict, int id) {
    return read_dict->seq[id].name;
}

bool is_read_illumina_contig(const sdict_t *read_dict, int id) {
    return get_read_name(read_dict, id).find("CONTIG_") == 0;
}

ma_sub_t *filter_reads_using_depth(int min_dp, float min_iden, int end_clip, size_t n, const ma_hit_t *a, const sdict_t *read_dict)
{
    size_t num_reads = read_dict->n_seq;

    size_t j, n_remained = 0;
    kvec_t(uint32_t) b = {0,0,0};
    ma_sub_t *subreads = 0;

    subreads = (ma_sub_t*)calloc(num_reads, sizeof(ma_sub_t));
    for (size_t i = 1, last = 0; i <= n; ++i) {
        if (i == n || a[i].qns>>32 != a[i-1].qns>>32) { // we come to a new query sequence
            size_t start = 0;
            int dp, qid = int(a[i-1].qns>>32);

            ma_sub_t max, max2;
            kv_resize(uint32_t, b, i - last);
            b.n = 0;
            for (j = last; j < i; ++j) { // collect all starts and ends
                uint32_t qs, qe;
                if (a[j].tn == qid || a[j].ml < a[j].bl * min_iden) continue; // skip self match
                qs = (uint32_t)a[j].qns + end_clip, qe = a[j].qe - end_clip;
                if (qe > qs) {
                    kv_push(uint32_t, b, qs<<1);
                    kv_push(uint32_t, b, qe<<1|1);
                }
            }

            // If the read is an Illumina contig, then we may not have alignments to the middle of
            // the read. So we don't do the normal miniasm stuff but instead only clip off unaligned
            // parts from its ends.
            int min_depth;
            if (is_read_illumina_contig(read_dict, qid)) {

                // If the Illumina contig has no alignments, then we just include the whole thing.
                if (b.n == 0) {
                    subreads[qid].s = 0;
                    subreads[qid].e = read_dict->seq[qid].len;
                }
                // If the read has alignments, we use those to clip off unaligned parts.
                else {
                    uint32_t min_start = numeric_limits<uint32_t>::max();
                    uint32_t max_end = 0;
                    for (j = 0; j < b.n; ++j) {
                        if (b.a[j] & 1) {  // is an end position
                            uint32_t read_end = b.a[j] >> 1;
                            max_end = std::max(read_end, max_end);
                        } else {  // is a start position
                            uint32_t read_start = b.a[j] >> 1;
                            min_start = std::min(read_start, min_start);
                        }
                    }
                    subreads[qid].s = min_start;
                    subreads[qid].e = max_end;
                }
                subreads[qid].del = 0;
                ++n_remained;
            }

            // If the read isn't a contig (i.e. it's a normal read) then we do the standard miniasm
            // behaviour: clip the read to the best depth region.
            else {
                min_depth = min_dp;

                ks_introsort_uint32_t(b.n, b.a);
                max.s = max.e = max.del = max2.s = max2.e = max2.del = 0;
                for (j = 0, dp = 0; j < b.n; ++j) {
                    int old_dp = dp;
                    if (b.a[j]&1) --dp;
                    else ++dp;
                    if (old_dp < min_depth && dp >= min_depth) {
                        start = b.a[j]>>1;
                    } else if (old_dp >= min_depth && dp < min_depth) {
                        int len = int((b.a[j]>>1) - start);
                        if (len > max.e - max.s) max2 = max, max.s = u_int32_t(start), max.e = b.a[j]>>1;
                        else if (len > max2.e - max2.s) max2.s = u_int32_t(start), max2.e = b.a[j]>>1;
                    }
                }
                if (max.e - max.s > 0) {
                    assert(qid < num_reads);
                    subreads[qid].s = max.s - end_clip;
                    subreads[qid].e = max.e + end_clip;
                    subreads[qid].del = 0;
                    ++n_remained;
                } else subreads[qid].del = 1;
            }


            last = i;
        }
    }
    free(b.a);
//    fprintf(stderr, "[M::%s::%s] %ld query sequences remain after sub\n", __func__, sys_timestamp(), n_remained);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] " << n_remained << " query sequences remain after sub\n";
    return subreads;
}

size_t filter_hits_using_span(const ma_sub_t *subreads, int min_span, size_t n, ma_hit_t *a)
{
    size_t i, m;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *p = &a[i];
        const ma_sub_t *rq = &subreads[p->qns>>32], *rt = &subreads[p->tn];
        int qs, qe, ts, te;
        if (rq->del || rt->del) continue;
        if (p->rev) {
            qs = p->te < rt->e? (uint32_t)p->qns : (uint32_t)p->qns + (p->te - rt->e);
            qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts);
            ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e);
            te = (uint32_t)p->qns > rq->s? p->te : p->te - (rq->s - (uint32_t)p->qns);
        } else {
            qs = p->ts > rt->s? (uint32_t)p->qns : (uint32_t)p->qns + (rt->s - p->ts);
            qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
            ts = (uint32_t)p->qns > rq->s? p->ts : p->ts + (rq->s - (uint32_t)p->qns);
            te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
        }
        qs = (qs > rq->s? qs : rq->s) - rq->s;
        qe = (qe < rq->e? qe : rq->e) - rq->s;
        ts = (ts > rt->s? ts : rt->s) - rt->s;
        te = (te < rt->e? te : rt->e) - rt->s;
        if (qe - qs >= min_span && te - ts >= min_span) {
            double r = (double)((qe - qs) + (te - ts)) / ((p->qe - (uint32_t)p->qns) + (p->te - p->ts));
            p->bl = (int)(p->bl * r + .499);
            p->ml = (int)(p->ml * r + .499);
            p->qns = p->qns>>32<<32 | qs, p->qe = qe, p->ts = ts, p->te = te;
            a[m++] = *p;
        }
    }
//    fprintf(stderr, "[M::%s::%s] %ld hits remain after cut\n", __func__, sys_timestamp(), m);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] " << m << " hits remain after cut\n";
    return m;
}

size_t filter_hits_using_overhang(const ma_sub_t *subreads, int max_hang, int min_ovlp, size_t n, ma_hit_t *a, float *cov)
{
    size_t i, m;
    asg_arc_t t;
    uint64_t tot_dp = 0, tot_len = 0;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *h = &a[i];
        const ma_sub_t *sq = &subreads[h->qns>>32], *st = &subreads[h->tn];
        int r;
        if (sq->del || st->del) continue;
        r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, .5, min_ovlp, &t);
        if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
            a[m++] = *h, tot_dp += r >= 0? r : r == MA_HT_QCONT? sq->e - sq->s : st->e - st->s;
    }
    for (i = 1; i <= m; ++i)
        if (i == m || a[i].qns>>32 != a[i-1].qns>>32)
            tot_len += subreads[a[i-1].qns>>32].e - subreads[a[i-1].qns>>32].s;
    *cov = (double)tot_dp / tot_len;
//    fprintf(stderr, "[M::%s::%s] %ld hits remain after filtering; crude coverage after filtering: %.2f\n", __func__, sys_timestamp(), m, *cov);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] " << m << " hits remain after filtering; crude coverage after filtering: " << *cov << "\n";
    return m;
}

void merge_subreads(size_t n_sub, ma_sub_t *a, const ma_sub_t *b)
{
    size_t i;
    for (i = 0; i < n_sub; ++i)
        a[i].e = a[i].s + b[i].e, a[i].s += b[i].s;
}

static inline int is_chimeric(int max_hang, int min_dp, size_t st, size_t en, const ma_hit_t *a, const ma_sub_t *sub, uint32_v c[2])
{
    size_t i;
    int k, chi[2];
    c[0].n = c[1].n = 0;
    for (i = st; i < en; ++i) {
        const ma_hit_t *h = &a[i];
        int ql = sub[h->qns>>32].e - sub[h->qns>>32].s, tl = sub[h->tn].e - sub[h->tn].s;
        int tl5, tl3, ql5 = (uint32_t)h->qns, ql3 = ql - h->qe;
        tl5 = h->rev? tl - h->te : h->ts;
        tl3 = h->rev? h->ts : tl - h->te;
        if (ql5 < max_hang && ql5 < tl5) {
            if (ql3 > max_hang && tl3 > max_hang) {
                kv_push(uint32_t, c[0], (ql - h->qe) << 1 | 1);
            } else if (ql3 > tl3 && tl3 < max_hang) {
                kv_push(uint32_t, c[0], (ql - h->qe) << 1 | 0);
            }
        } else if (ql3 < max_hang && ql3 < tl3) {
            if (ql5 > max_hang && tl5 > max_hang) {
                kv_push(uint32_t, c[1], (uint32_t)h->qns << 1 | 1);
            } else if (ql5 > tl5 && tl5 < max_hang) {
                kv_push(uint32_t, c[1], (uint32_t)h->qns << 1 | 0);
            }
        }
    }
    if (c[0].n < min_dp || c[1].n < min_dp) return 0;
    chi[0] = chi[1] = -1;
    for (k = 0; k < 2; ++k) {
        int cnt[2], max = 0;
        ks_introsort_uint32_t(c[k].n, c[k].a);
        cnt[0] = cnt[1] = 0;
        for (i = 0; i < c[k].n; ++i) {
            ++cnt[c[k].a[i]&1];
            max = max > cnt[1] - cnt[0]? max : cnt[1] - cnt[0];
        }
        if (max >= min_dp) chi[k] = max;
    }
    return (chi[0] > 0 || chi[1] > 0);
}

size_t remove_chimeric_reads(int max_hang, int min_dp, size_t n, const ma_hit_t *a, const sdict_t *read_dict, ma_sub_t *subreads)
{
    size_t i, start = 0, n_chi = 0;
    uint32_v c[2] = {{0,0,0}, {0,0,0}};
    for (i = 1; i <= n; ++i) {
        if (i == n || a[i].qns>>32 != a[start].qns>>32) {
            if (is_chimeric(max_hang, min_dp, start, i, a, subreads, c)) {

                // Illumina contigs are exempt from being labelled chimeric.
                int id = int(a[start].qns>>32);
                if (!is_read_illumina_contig(read_dict, id)) {
                    subreads[id].del = 1;
                    ++n_chi;
                }
            }
            start = i;
        }
    }
    free(c[0].a); free(c[1].a);
//    fprintf(stderr, "[M::%s::%s] identified %ld chimeric reads\n", __func__, sys_timestamp(), n_chi);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] identified " << n_chi << " chimeric reads\n";
    return n_chi;
}

size_t remove_contained_reads(int max_hang, float int_frac, int min_ovlp, sdict_t *read_dict, ma_sub_t *subreads, size_t n, ma_hit_t *a)
{
    int32_t *map, r;
    size_t i, m, old_n_seq = read_dict->n_seq;
    asg_arc_t t;
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *h = &a[i];

        int query_i = int(h->qns>>32);
        int target_i = h->tn;
        ma_sub_t *query_subread = &subreads[query_i];
        ma_sub_t *target_subread = &subreads[target_i];

        r = ma_hit2arc(h, query_subread->e - query_subread->s, target_subread->e - target_subread->s, max_hang, int_frac, min_ovlp, &t);

        if (r == MA_HT_QCONT) {  // If the query is contained in the target
            query_subread->del = 1;
        }
        else if (r == MA_HT_TCONT) {  // If the target is contained in the query
            target_subread->del = 1;
        }
    }

    for (i = 0; i < read_dict->n_seq; ++i) {

        // Illumina contigs can't be deleted!
        if (is_read_illumina_contig(read_dict, i))
            subreads[i].del = 0;

        if (subreads[i].del)
            read_dict->seq[i].del = 1;
    }

    ma_hit_mark_unused(read_dict, n, a);
    map = sd_squeeze(read_dict);
    for (i = 0; i < old_n_seq; ++i)
        if (map[i] >= 0) subreads[map[i]] = subreads[i];
    for (i = m = 0; i < n; ++i) {
        ma_hit_t *h = &a[i];
        int32_t qn = map[h->qns>>32], tn = map[h->tn];
        if (qn >= 0 && tn >= 0) {
            a[i].qns = (uint64_t)qn<<32 | (uint32_t)a[i].qns;
            a[i].tn = tn;
            a[m++] = a[i];
        }
    }
    free(map);
//    fprintf(stderr, "[M::%s::%s] %d sequences and %ld hits remain after containment removal\n", __func__, sys_timestamp(), read_dict->n_seq, m);
    std::cerr << "[M::" << __func__ << "::" << sys_timestamp() << "] " << read_dict->n_seq << " sequences and " << m << " hits remain after containment removal\n";
    return m;
}
