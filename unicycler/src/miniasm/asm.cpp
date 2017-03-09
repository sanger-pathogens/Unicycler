#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include <iostream>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <miniasm/sdict.h>

#include "miniasm/miniasm.h"
#include "miniasm/kvec.h"
#include "miniasm/sdict.h"
#include "miniasm/kseq.h"

KSEQ_INIT(gzFile, gzread)


using namespace std;


asg_t *make_string_graph(int max_hang, float int_frac, int min_ovlp, sdict_t const *read_dict, ma_sub_t const *sub, unsigned long n_hits, ma_hit_t const *hit)
{
    size_t i;
    asg_t *g;
    g = asg_init();
    for (i = 0; i < read_dict->n_seq; ++i) {
        if (sub)
            asg_seq_set(g, i, sub[i].e - sub[i].s, (sub[i].del || read_dict->seq[i].del));
        else
            asg_seq_set(g, i, read_dict->seq[i].len, read_dict->seq[i].del);
    }
    for (i = 0; i < n_hits; ++i) {
        int r;
        asg_arc_t t, *p;
        const ma_hit_t *h = &hit[i];
        uint32_t qn = h->qns>>32;
        int ql = sub? sub[qn].e - sub[qn].s : read_dict->seq[qn].len;
        int tl = sub? sub[h->tn].e - sub[h->tn].s : read_dict->seq[h->tn].len;
        r = ma_hit2arc(h, ql, tl, max_hang, int_frac, min_ovlp, &t);
        if (r >= 0) {
            if (qn == h->tn) { // self match
                if ((uint32_t)h->qns == h->ts && h->qe == h->te && h->rev) // PacBio-specific artifact (TODO: is this right when we skip target containment above?)
                    g->seq[qn].del = 1;
                continue;
            }
            p = asg_arc_pushp(g);
            *p = t;
        } else if (r == MA_HT_QCONT) g->seq[qn].del = 1;
    }
    asg_cleanup(g);
    std::cerr << "[M::" << __func__ << "] read " << g->n_arc << " arcs\n";
    return g;
}

void save_string_graph(const asg_t *g, const sdict_t *read_dict, const ma_sub_t *subreads, std::string graph_filename, const char *reads_filename)
{
    FILE *fp = fopen(graph_filename.c_str(), "w");

    // First we have to figure out which reads will be included in the graph. Include any which are
    // part of edges.
    set<size_t> used_read_indices;
    for (uint32_t i = 0; i < g->n_arc; ++i) {
        const asg_arc_t *p = &g->arc[i];
        size_t query_i = p->ul >> 33;
        size_t target_i = p->v >> 1;
        used_read_indices.insert(query_i);
        used_read_indices.insert(target_i);
    }

    // Also include any reads which are contigs.
    for (size_t i = 0; i < read_dict->n_seq; ++i) {
        if (is_read_illumina_contig(read_dict, i))
            used_read_indices.insert(i);
    }
    unordered_set<string> used_read_names;
    for (set<size_t>::iterator it = used_read_indices.begin(); it != used_read_indices.end(); ++it)
        used_read_names.insert(string(read_dict->seq[*it].name));

    // Now we can load in the sequences for the reads that we need.
    unordered_map<string,string> read_seqs;
    gzFile reads_file = reads_filename && strcmp(reads_filename, "-")? gzopen(reads_filename, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *ks = kseq_init(reads_file);
    while (kseq_read(ks) >= 0) {
        string read_name = ks->name.s;
        if (used_read_names.find(read_name) != used_read_names.end())  // if the read is used
            read_seqs[read_name] = ks->seq.s;
    }

    // Now we print the segment (S) lines.
    for (set<size_t>::iterator it = used_read_indices.begin(); it != used_read_indices.end(); ++it) {
        string read_name = read_dict->seq[*it].name;
        string read_seq = read_seqs[read_name];
        if (subreads) {
            const ma_sub_t *subread = &subreads[*it];
            read_name += ':';
            read_name += to_string(subread->s + 1);
            read_name += '-';
            read_name += to_string(subread->e);
            u_int32_t len = subread->e - subread->s;
            read_seq = read_seq.substr(subread->s, len);
        }
        fprintf(fp, "S\t%s\t%s\n", read_name.c_str(), read_seq.c_str());
    }

    // Then we print the link (L) lines.
    for (uint32_t i = 0; i < g->n_arc; ++i) {
        const asg_arc_t *p = &g->arc[i];

        size_t query_i = p->ul>>33;
        size_t target_i = p->v>>1;

        char *query_name = read_dict->seq[query_i].name;
        char *target_name = read_dict->seq[target_i].name;

        char query_strand = "+-"[p->ul>>32&1];
        char target_strand = "+-"[p->v&1];

        if (subreads) {
            const ma_sub_t *query_subread = &subreads[query_i];
            const ma_sub_t *target_subread = &subreads[target_i];
            fprintf(fp, "L\t%s:%d-%d\t%c\t%s:%d-%d\t%c",
                    query_name, query_subread->s + 1, query_subread->e, query_strand,
                    target_name, target_subread->s + 1, target_subread->e, target_strand);
        } else
            fprintf(fp, "L\t%s\t%c\t%s\t%c", query_name, query_strand, target_name, target_strand);
        fprintf(fp, "\t%dM\tSD:i:%d\tml:i:%d\tmr:f:%.4f\n", p->ol, (uint32_t)p->ul, p->ml, p->mr);
    }
    fclose(fp);
}

/*********************
 * Unitig generation *
 *********************/

#include "miniasm/kdq.h"
KDQ_INIT(uint64_t)

void destroy_unitig_graph(ma_ug_t *ug)
{
    uint32_t i;
    if (ug == 0) return;
    for (i = 0; i < ug->u.n; ++i) {
        free(ug->u.a[i].a);
        free(ug->u.a[i].s);
    }
    free(ug->u.a);
    destroy_string_graph(ug->g);
    free(ug);
}

void save_unitig_graph(const ma_ug_t *ug, const sdict_t *d, const ma_sub_t *sub, std::string graph_filename)
{
    uint32_t i, j, l;
    FILE *fp = fopen(graph_filename.c_str(), "w");
    char name[32];
    for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA
        ma_utg_t *p = &ug->u.a[i];
        sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
        fprintf(fp, "S\t%s\t%s\tLN:i:%d\n", name, p->s? p->s : "*", p->len);
        for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
            uint32_t x = p->a[j]>>33;
            if (sub) fprintf(fp, "a\t%s\t%d\t%s:%d-%d\t%c\t%d\n", name, l, d->seq[x].name, sub[x].s + 1, sub[x].e, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
            else fprintf(fp, "a\t%s\t%d\t%s\t%c\t%d\n", name, l, d->seq[x].name, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
        }
    }
    for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
        uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
        fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tSD:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
                (v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
    }
    for (i = 0; i < ug->u.n; ++i) { // summary of unitigs
        uint32_t cnt[2];
        ma_utg_t *u = &ug->u.a[i];
        if (u->start == UINT32_MAX) {
            fprintf(fp, "x\tutg%.6dc\t%d\t%d\n", i + 1, u->len, u->n);
        } else {
            for (j = 0; j < 2; ++j) cnt[j] = asg_arc_n(ug->g, i<<1|j);
            if (sub)
                fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s:%d-%d\t%c\t%s:%d-%d\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
                        d->seq[u->start>>1].name, sub[u->start>>1].s + 1, sub[u->start>>1].e, "+-"[u->start&1],
                        d->seq[u->end>>1].name, sub[u->end>>1].s + 1, sub[u->end>>1].e, "+-"[u->end&1]);
            else
                fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s\t%c\t%s\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
                        d->seq[u->start>>1].name, "+-"[u->start&1], d->seq[u->end>>1].name, "+-"[u->end&1]);
        }
    }
    fclose(fp);
}

#define arc_cnt(g, v) ((uint32_t)(g)->idx[(v)])
#define arc_first(g, v) ((g)->arc[(g)->idx[(v)]>>32])

ma_ug_t *make_unitig_graph(asg_t *g)
{
    int32_t *mark;
    uint32_t i, v, n_vtx = g->n_seq * 2;
    kdq_t(uint64_t) *q;
    ma_ug_t *ug;

    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = asg_init();
    mark = (int32_t*)calloc(n_vtx, 4);

    q = kdq_init(uint64_t);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t w, x, l, start, end, len;
        ma_utg_t *p;
        if (g->seq[v>>1].del || arc_cnt(g, v) == 0 || mark[v]) continue;
        mark[v] = 1;
        q->count = 0, start = v, end = v^1, len = 0;
        // forward
        w = v;
        while (1) {
            if (arc_cnt(g, w) != 1) break;
            x = arc_first(g, w).v; // w->x
            if (arc_cnt(g, x^1) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
            end = x^1, len += l;
            w = x;
            if (x == v) break;
        }
        if (start != (end^1) || kdq_size(q) == 0) { // linear unitig
            l = g->seq[end>>1].len;
            kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
            len += l;
        } else { // circular unitig
            start = end = UINT32_MAX;
            goto add_unitig; // then it is not necessary to do the backward
        }
        // backward
        x = v;
        while (1) { // similar to forward but not the same
            if (arc_cnt(g, x^1) != 1) break;
            w = arc_first(g, x^1).v ^ 1; // w->x
            if (arc_cnt(g, w) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
            start = w, len += l;
            x = w;
        }
add_unitig:
        if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        kv_pushp(ma_utg_t, ug->u, &p);
        p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
        p->m = p->n;
        kv_roundup32(p->m);
        p->a = (uint64_t*)malloc(8 * p->m);
        for (i = 0; i < kdq_size(q); ++i)
            p->a[i] = kdq_at(q, i);
    }
    kdq_destroy(uint64_t, q);

    // add arcs between unitigs; reusing mark for a different purpose
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;
    for (i = 0; i < ug->u.n; ++i) {
        if (ug->u.a[i].circ) continue;
        mark[ug->u.a[i].start] = i<<1 | 0;
        mark[ug->u.a[i].end] = i<<1 | 1;
    }
    for (i = 0; i < g->n_arc; ++i) {
        asg_arc_t *p = &g->arc[i];
        if (p->del) continue;
        if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
            asg_arc_t *q;
            uint32_t u = mark[p->ul>>32^1]^1;
            int l = ug->u.a[u>>1].len - p->ol;
            if (l < 0) l = 1;
            q = asg_arc_pushp(ug->g);
            q->ol = p->ol, q->del = 0;
            q->ul = (uint64_t)u<<32 | l;
            q->v = mark[p->v];
        }
    }
    for (i = 0; i < ug->u.n; ++i)
        asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
    asg_cleanup(ug->g);
    free(mark);
    return ug;
}

/*******************
 * Unitig sequence *
 *******************/


typedef struct {
    uint32_t utg:31, ori:1, start, len;
} utg_intv_t;

static char comp_tab[] = { // complement base
      0,   1,    2,     3,      4,   5,    6,     7,      8,   9,  10,    11,     12,  13,  14,    15,
     16,  17,  18,    19,     20,  21,  22,    23,     24,  25,  26,    27,     28,  29,  30,    31,
     32,  33,  34,    35,     36,  37,  38,    39,     40,  41,  42,    43,     44,  45,  46,    47,
     48,  49,  50,    51,     52,  53,  54,    55,     56,  57,  58,    59,     60,  61,  62,    63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',    91,     92,  93,  94,    95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};


int generate_unitig_seqs(ma_ug_t *g, const sdict_t *d, const ma_sub_t *sub, const char *reads_filename)
{
    gzFile fp = reads_filename && strcmp(reads_filename, "-")? gzopen(reads_filename, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) return -1;
    kseq_t *ks = kseq_init(fp);

    utg_intv_t *tmp = (utg_intv_t*)calloc(d->n_seq, sizeof(utg_intv_t));
    for (uint32_t i = 0; i < g->u.n; ++i) {
        ma_utg_t *u = &g->u.a[i];
        uint32_t l = 0;
        u->s = (char*)calloc(1, u->len + 1);
        memset(u->s, 'N', u->len);
        for (uint32_t j = 0; j < u->n; ++j) {
            utg_intv_t *t = &tmp[u->a[j]>>33];
            assert(t->len == 0);
            t->utg = i, t->ori = u->a[j]>>32&1;
            t->start = l, t->len = (uint32_t)u->a[j];
            l += t->len;
        }
    }

    while (kseq_read(ks) >= 0) {
        int32_t id;
        utg_intv_t *t;
        ma_utg_t *u;
        id = sd_get(d, ks->name.s);
        if (id < 0 || tmp[id].len == 0) continue;
        t = &tmp[id];
        u = &g->u.a[t->utg];
        if (sub) {
            assert(sub[id].e - sub[id].s <= ks->seq.l);
            memmove(ks->seq.s, ks->seq.s + sub[id].s, sub[id].e - sub[id].s);
            ks->seq.l = sub[id].e - sub[id].s;
        }
        if (!t->ori) { // forward strand
            for (uint32_t i = 0; i < t->len; ++i)
                u->s[t->start + i] = ks->seq.s[i];
        } else {
            for (uint32_t i = 0; i < t->len; ++i) {
                int c = (uint8_t)ks->seq.s[ks->seq.l - 1 - i];
                u->s[t->start + i] = c >= 128? 'N' : comp_tab[c];
            }
        }
    }
    free(tmp);

    kseq_destroy(ks);
    gzclose(fp);
    return 0;
}
