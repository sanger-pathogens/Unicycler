#ifndef SDICT_H
#define SDICT_H

#include <stdint.h>

typedef struct {
	char *name;
	uint32_t len, aux:31, del:1;
} sd_seq_t;

typedef struct {
	uint32_t n_seq, m_seq;
	sd_seq_t *seq;
	void *h;
} sdict_t;

sdict_t *init_seq_dict(void);
void destroy_seq_dict(sdict_t *d);
int32_t sd_put(sdict_t *d, const char *name, uint32_t len);
int32_t sd_get(const sdict_t *d, const char *name);
int32_t *sd_squeeze(sdict_t *d);

#endif
