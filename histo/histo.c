//
// Created by lisanhu on 9/30/18.
//

#include <stdbool.h>

#include "histo.h"

// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})

// shorthand way to set value in hash with single line command.  Returns value
// returns 0=replaced existing item, 1=bucket empty (new key), 2-adding element previously deleted
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

inline histo *histo_init(u32 cap) {
//	KHASH_INIT(histo, u64, u64)
	histo *his = malloc(sizeof(histo));
	his->kh = kh_init(64);
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
	kh_destroy(64, h->kh);
	h->cap = h->size = 0;
}

inline static u64 histo_key_hash(u64 key) {
//	return key >> 4;
	return key;
}

inline static void histo_push(histo *h, u64 key, u64 val) {
	u32 cap = h->cap;
	u32 size = h->size;
	if (size + 1 == cap) {
		h->entries = realloc(h->entries, cap * 2 * sizeof(entry));
		h->cap = cap * 2;
	}

	entry e = {.key = histo_key_hash(key), .val = val};
	h->entries[h->size++] = e;
}

inline void histo_add(histo *h, u64 key) {
	u32 size = h->size;
	bool found = false;
	for (u32 i = 0; i < size; ++i) {
		if (h->entries[i].key == histo_key_hash(key)) {
			found = true;
			h->entries[i].val += 1;
		}
	}

	if (!found) {
		histo_push(h, key, 1);
	}
//	khiter_t k;
//	u64 ret = kh_get_val(64, h->kh, key, 0);
//	ret += 1;
//	kh_set(64, h->kh, key, ret);
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
	store[0].val = store[0].key = 0;
	store[1].val = store[1].key = 0;
//	khiter_t k;
//	for (k = kh_begin(h->kh); k != kh_end(h->kh); ++k) {
//		if (kh_exist(h->kh, k)) {
//			u64 key = kh_key(h->kh, k);
//			u64 val = kh_value(h->kh, k);
//			if (store[1].val <= val && store[0].val <= val) {
//				store[0].key = key;
//				store[0].val = val;
//			} else if (store[1].val <= val && store[0].val > val) {
//				store[1].key = key;
//				store[1].val = val;
//			}
//		}
//	}
	for (u32 i = 0; i < h->size; ++i) {
		entry e = h->entries[i];
		if (store[1].val <= e.val && store[0].val <= e.val) store[0] = e;
		if (store[1].val <= e.val && store[0].val > e.val) store[1] = e;
	}
	return store[0].val + store[1].val;
}
