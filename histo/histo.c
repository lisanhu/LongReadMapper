//
// Created by lisanhu on 9/30/18.
//

#include <stdbool.h>

#include "histo.h"

inline histo *histo_init(u32 cap) {
	histo *his = malloc(sizeof(histo));
	his->cap = cap;
	his->size = 0;
	his->entries = malloc(sizeof(entry) * cap);
	return his;
}

inline void histo_destroy(histo *h) {
	if (h->entries) {
		free(h->entries);
		h->entries = NULL;
	}
	h->cap = h->size = 0;
    free(h);
}

inline static u64 histo_key_hash(u64 key) {
	return key >> 4U;
//	return key;
}

inline static void histo_push(histo *h, u64 key, u64 val) {
	u32 cap = h->cap;
	u32 size = h->size;
	if (size + 1 == cap) {
		h->entries = realloc(h->entries, cap * 2 * sizeof(entry));
		h->cap = cap * 2;
	}

	entry e = {.key = key, .val = val, .bucket = histo_key_hash(key)};
	h->entries[h->size++] = e;
}

inline void histo_add(histo *h, u64 key) {
	u32 size = h->size;
	bool found = false;
	for (u32 i = 0; i < size; ++i) {
		if (h->entries[i].bucket == histo_key_hash(key)) {
			found = true;
			h->entries[i].val += 1;
			h->entries[i].key = key < h->entries[i].key ? key : h->entries[i].key;
		}
	}

	if (!found) {
		histo_push(h, key, 1);
	}
}

inline u64 histo_get(histo *h, u64 key) {
	u32 size = h->size, i;
	bool found = false;
	for (i = 0; i < size; ++i) {
		if (h->entries[i].key == histo_key_hash(key)) {
			found = true;
		}
	}

	if (!found) {
		return 0;
	} else {
		return h->entries[i].val;
	}
}

inline u64 histo_find_max(histo *h) {
	u64 r = 0, v = 0;
	for (u32 i = 0; i < h->size; ++i) {
		if (v < h->entries[i].val) {
			r = i; v = h->entries[i].val;
		}
	}
	return h->entries[r].key;
}

inline u64 histo_find_2_max(histo *h, entry *store) {
	store[0].bucket = store[0].val = store[0].key = 0;
    store[1].bucket = store[1].val = store[1].key = 0;

	for (u32 i = 0; i < h->size; ++i) {
		entry e = h->entries[i];
		if (store[1].val < e.val && store[0].val < e.val) {
		    store[1] = store[0];
		    store[0] = e;
		} else if (store[1].val < e.val && store[0].val >= e.val) store[1] = e;
	}
	return store[0].val + store[1].val;
}
