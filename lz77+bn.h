#ifndef lz77_definition
#define lz77_definition

#include <errno.h>
#include <stdint.h>

// Naive LZ77 implementation inspired by CharGPT discussion
// and my personal passion to compressors in 198x

typedef struct lz77_s lz77_t;

enum {
    lz77_min_window = 10,
    lz77_max_window = 12,
    lz77_alphabet   = 1 << lz77_max_window
};

typedef struct lz77_binheap_s {
    int32_t  ns[lz77_alphabet]; // nodes
    int32_t  sx[lz77_alphabet]; // symbol -> ix in ns[]
    uint64_t fq[lz77_alphabet]; // frequency
    int32_t  nc;                // node count
} lz77_binheap_t;

typedef struct lz77_s {
    // `that` see: https://gist.github.com/leok7v/8d118985d3236b0069d419166f4111cf
    void*    that;  // caller supplied data
    errno_t  error; // sticky; for read()/write() compress() and decompress()
    // caller supplied read()/write() must error via .error field
    uint64_t (*read)(lz77_t*); //  reads 64 bits
    void     (*write)(lz77_t*, uint64_t b64); // writes 64 bits
    uint64_t written;
    lz77_binheap_t bh_txt;
    lz77_binheap_t bh_pos;
    lz77_binheap_t bh_len;
} lz77_t;

typedef struct lz77_if {
    // `window_bits` is a log2 of window size in bytes must be in range [10..20]
    void (*write_header)(lz77_t* lz77, size_t bytes, uint8_t window_bits);
    void (*compress)(lz77_t* lz77, const uint8_t* data, size_t bytes,
                     uint8_t window_bits);
    void (*read_header)(lz77_t* lz77, size_t *bytes, uint8_t *window_bits);
    void (*decompress)(lz77_t* lz77, uint8_t* data, size_t bytes,
                       uint8_t window_bits);
    // Writing and reading envelope of source data `bytes` and
    // `window_bits` is caller's responsibility.
} lz77_if;

extern lz77_if lz77;

#endif // lz77_definition

#if defined(lz77_implementation) && !defined(lz77_implemented)

#define lz77_implemented

#ifndef lz77_assert
#define lz77_assert(...) do { } while (0)
#endif

#ifndef lz77_println
#define lz77_println(...) do { } while (0)
#endif

#define lz77_if_error_return(lz);do {                   \
    lz77_assert(lz->error == 0);                        \
    if (lz->error) { return; }                          \
} while (0)

#define lz77_return_invalid(lz) do {                    \
    lz->error = EINVAL;                                 \
    return;                                             \
} while (0)

#ifdef lz77_historgram

static inline uint32_t lz77_bit_count(size_t v) {
    uint32_t count = 0;
    while (v) { count++; v >>= 1; }
    return count;
}

static size_t lz77_hist_len[64];
static size_t lz77_hist_pos[64];

#define lz77_init_histograms() do {                     \
    memset(lz77_hist_pos, 0x00, sizeof(lz77_hist_pos)); \
    memset(lz77_hist_len, 0x00, sizeof(lz77_hist_len)); \
} while (0)

#define  lz77_histogram_pos_len(pos, len) do {          \
    lz77_hist_pos[lz77_bit_count(pos)]++;               \
    lz77_hist_len[lz77_bit_count(len)]++;               \
} while (0)

#define lz77_dump_histograms() do {                                 \
    lz77_println("Histogram log2(len):");                           \
    for (int8_t i_ = 0; i_ < 64; i_++) {                            \
        if (lz77_hist_len[i_] > 0) {                                \
            lz77_println("len[%d]: %lld", i_, lz77_hist_len[i_]);   \
        }                                                           \
    }                                                               \
    lz77_println("Histogram log2(pos):");                           \
    for (int8_t i_ = 0; i_ < 64; i_++) {                            \
        if (lz77_hist_pos[i_] > 0) {                                \
            lz77_println("pos[%d]: %lld", i_, lz77_hist_pos[i_]);   \
        }                                                           \
    }                                                               \
} while (0)

#else

#define lz77_init_histograms()           do { } while (0)
#define lz77_histogram_pos_len(pos, len) do { } while (0)
#define lz77_dump_histograms()           do { } while (0)

#endif

// https://en.wikipedia.org/wiki/Binary_heap
// J. W. J. Williams 1964

static_assert(lz77_alphabet > 2 && (lz77_alphabet & (lz77_alphabet - 1)) == 0,
              "lz77_alphabet must be 2^n");

static inline void lz77_swap_int32(int32_t* a, int32_t* b) {
    int32_t t = *a; *a = *b; *b = t;
}

static inline void lz77_swap_int64(uint64_t* a, uint64_t* b) {
    uint64_t t = *a; *a = *b; *b = t;
}

static inline void lz77_binheap_swap(lz77_binheap_t* t, int32_t ix0, int32_t ix1) {
    lz77_assert(ix0 != ix1);
    const int32_t s0 = t->ns[ix0];
    const int32_t s1 = t->ns[ix1];
    lz77_assert(s0 != s1);
    lz77_assert(0 <= s0 && s0 < t->nc);
    lz77_assert(0 <= s1 && s1 < t->nc);
//  lz77_println("swap([%03d]: %02X %c, [%03d]: %02X %c) freq: %lld %lld",
//                ix0, s0, s0 >= 0x20 ? s0 : 0x20,
//                ix1, s1, s1 >= 0x20 ? s1 : 0x20,
//                t->fq[s0], t->fq[s1]);
    lz77_swap_int32(&t->ns[ix0], &t->ns[ix1]);
    lz77_swap_int32(&t->sx[s0], &t->sx[s1]);
    // do not need to swap fq[s0] fq[s1]
}

static inline int32_t lz77_binheap_up_heapify(lz77_binheap_t* t, int32_t ix) {
    lz77_assert(0 <= ix && ix < t->nc);
    while (ix > 0) {
        int32_t p = (ix - 1) / 2; // parent
        if (t->fq[t->ns[p]] >= t->fq[t->ns[ix]]) { break; }
        lz77_binheap_swap(t, ix, p);
        ix = p;
    }
    lz77_assert(t->sx[t->ns[ix]] == ix);
    return ix;
}

static int32_t lz77_binheap_add(lz77_binheap_t* t, int32_t sym) {
    lz77_assert(0 <= sym && sym < lz77_alphabet);
    lz77_assert(t->nc < lz77_alphabet);
    t->ns[t->nc] = sym;
    t->fq[sym] = 0;
    t->sx[sym] = t->nc;
    t->nc++;
    return lz77_binheap_up_heapify(t, t->nc - 1);
}

static inline int32_t lz77_binheap_inc_freq(lz77_binheap_t* t, int32_t sym) {
    lz77_assert(0 <= sym && sym < t->nc);
    int32_t ix = t->sx[sym];
    lz77_assert(0 <= ix && ix < t->nc, "ix: %d", ix);
    t->fq[sym]++;
//  lz77_println("%02X %c freq:%lld", sym, sym, t->fq[sym]);
    return lz77_binheap_up_heapify(t, ix);
}

static void lz77_binheap_init(lz77_binheap_t* t, int32_t nc) { // node count
    t->nc = 0; // no nodes
    memset(t->ns, 0xFF, sizeof(t->ns));
    memset(t->sx, 0xFF, sizeof(t->sx));
    memset(t->fq, 0, sizeof(t->fq));
    for (int32_t s = 0; s < nc; s++) {
        lz77_binheap_add(t, s);
    }
    lz77_assert(t->nc == nc);
}

static inline void lz77_write_bit(lz77_t* lz, uint64_t* b64,
        uint32_t* bp, uint64_t bit) {
    if (*bp == 64 && lz->error == 0) {
        lz->write(lz, *b64);
        *b64 = 0;
        *bp = 0;
        if (lz->error == 0) { lz->written += 8; }
    }
    *b64 |= bit << *bp;
    (*bp)++;
}

static inline void lz77_write_bits(lz77_t* lz, uint64_t* b64,
        uint32_t* bp, uint64_t bits, uint32_t n) {
    rt_assert(n <= 64);
    while (n > 0) {
        lz77_write_bit(lz, b64, bp, bits & 1);
        bits >>= 1;
        n--;
    }
}

static inline void lz77_write_number(lz77_t* lz, uint64_t* b64,
        uint32_t* bp, uint64_t bits, uint8_t base) {
    do {
        lz77_write_bits(lz, b64, bp, bits, base);
        bits >>= base;
        lz77_write_bit(lz, b64, bp, bits != 0); // continue bit
    } while (bits != 0);
}

static inline void lz77_flush(lz77_t* lz, uint64_t b64, uint32_t bp) {
    if (bp > 0 && lz->error == 0) {
        lz->write(lz, b64);
        if (lz->error == 0) { lz->written += 8; }
    }
}

static void lz77_write_header(lz77_t* lz, size_t bytes, uint8_t window_bits) {
    lz77_if_error_return(lz);
    if (window_bits < 10 || window_bits > 20) { lz77_return_invalid(lz); }
    lz->write(lz, (uint64_t)bytes);
    lz77_if_error_return(lz);
    lz->write(lz, (uint64_t)window_bits);
}

static void lz77_compress(lz77_t* lz, const uint8_t* data, size_t bytes,
        uint8_t window_bits) {
    lz77_if_error_return(lz);
    if (window_bits < 10 || window_bits > 20) { lz77_return_invalid(lz); }
    lz77_init_histograms();
    const size_t window = ((size_t)1U) << window_bits;
//  const uint8_t base = (window_bits - 4) / 2;
    const uint8_t base = 4;
    lz77_binheap_init(&lz->bh_txt, 0x80); // ascii text
    lz77_binheap_init(&lz->bh_pos, (int32_t)window);
    lz77_binheap_init(&lz->bh_len, (int32_t)window);
    uint64_t b64 = 0;
    uint32_t bp = 0;
    size_t i = 0;
    while (i < bytes) {
        // length and position of longest matching sequence
        size_t len = 0;
        size_t pos = 0;
        if (i >= 1) {
            size_t j = i - 1;
            size_t min_j = i > window ? i - window : 0;
            while (j > min_j) {
                rt_assert((i - j) < window);
                const size_t n = bytes - i;
                size_t k = 0;
                while (k < n && data[j + k] == data[i + k]) {
                    k++;
                }
                if (k > len) {
                    len = k;
                    pos = i - j;
                }
                j--;
            }
        }
        if (len > 2) {
//lz77_println("i: %lld pos: %lld len: %lld ->", i, pos, len);
//lz77_println("i: %lld pos: %lld len: %lld", i, lz->bh_pos.sx[pos], lz->bh_len.sx[len]);
            rt_assert(0 < pos && pos < window);
            rt_assert(0 < len);
            lz77_write_bits(lz, &b64, &bp, 0b11, 2); // flags
            lz77_if_error_return(lz);
            lz77_write_number(lz, &b64, &bp, lz->bh_pos.sx[pos], base);
            lz77_binheap_inc_freq(&lz->bh_pos, (int32_t)pos);
            lz77_if_error_return(lz);
            lz77_write_bit(lz, &b64, &bp, len >= window); // flag: long len
            lz77_if_error_return(lz);
            if (len >= window) {
                lz77_write_number(lz, &b64, &bp, len, base);
            } else {
                lz77_write_number(lz, &b64, &bp, lz->bh_len.sx[len], base);
                lz77_binheap_inc_freq(&lz->bh_len, (int32_t)len);
            }
            lz77_if_error_return(lz);
//          lz77_histogram_pos_len(pos, len);
            if (len < window) {
                lz77_histogram_pos_len(lz->bh_pos.sx[pos], lz->bh_len.sx[len]);
            }
            i += len;
        } else {
            const uint8_t b = data[i];
//lz77_println("i: %lld byte: %08X %c", i, b, b);
            // European texts are predominantly spaces and small ASCII letters:
            if (b < 0x80) {
                lz77_write_bit(lz, &b64, &bp, 0); // flags
                lz77_if_error_return(lz);
                // ASCII byte < 0x80 with 8th bit set to `0`
                const uint8_t bh = (uint8_t)lz->bh_txt.sx[b];
//              lz77_write_bits(lz, &b64, &bp, bh, 7);
                lz77_write_number(lz, &b64, &bp, bh, 2);
                lz77_binheap_inc_freq(&lz->bh_txt, b);
                lz77_if_error_return(lz);
            } else {
                lz77_write_bit(lz, &b64, &bp, 1); // flag: 1
                lz77_write_bit(lz, &b64, &bp, 0); // flag: 0
                lz77_if_error_return(lz);
                // only 7 bit because 8th bit is `1`
                const uint8_t bh = (uint8_t)lz->bh_txt.sx[b & 0x7F];
//              lz77_write_bits(lz, &b64, &bp, bh, 7);
                lz77_write_number(lz, &b64, &bp, bh, 2);
                lz77_binheap_inc_freq(&lz->bh_txt, b & 0x7F);
                lz77_if_error_return(lz);
            }
            i++;
        }
    }
    lz77_flush(lz, b64, bp);
    lz77_dump_histograms();
}

static inline uint64_t lz77_read_bit(lz77_t* lz, uint64_t* b64, uint32_t* bp) {
    if (*bp == 0) { *b64 = lz->read(lz); }
    uint64_t bit = (*b64 >> *bp) & 1;
    *bp = *bp == 63 ? 0 : *bp + 1;
    return bit;
}

static inline uint64_t lz77_read_bits(lz77_t* lz, uint64_t* b64,
        uint32_t* bp, uint32_t n) {
    rt_assert(n <= 64);
    uint64_t bits = 0;
    for (uint32_t i = 0; i < n && lz->error == 0; i++) {
        uint64_t bit = lz77_read_bit(lz, b64, bp);
        bits |= bit << i;
    }
    return bits;
}

static inline uint64_t lz77_read_number(lz77_t* lz, uint64_t* b64,
        uint32_t* bp, uint8_t base) {
    uint64_t bits = 0;
    uint64_t bit = 0;
    uint32_t shift = 0;
    do {
        bits |= (lz77_read_bits(lz, b64, bp, base) << shift);
        shift += base;
        bit = lz77_read_bit(lz, b64, bp);
    } while (bit && lz->error == 0);
    return bits;
}

static void lz77_read_header(lz77_t* lz, size_t *bytes, uint8_t *window_bits) {
    lz77_if_error_return(lz);
    *bytes = (size_t)lz->read(lz);
    *window_bits = (uint8_t)lz->read(lz);
    if (*window_bits < 10 || *window_bits > 20) { lz77_return_invalid(lz); }
}

static void lz77_decompress(lz77_t* lz, uint8_t* data, size_t bytes,
        uint8_t window_bits) {
    lz77_if_error_return(lz);
    uint64_t b64 = 0;
    uint32_t bp = 0;
    if (window_bits < 10 || window_bits > 20) { lz77_return_invalid(lz); }
    const size_t window = ((size_t)1U) << window_bits;
//  const uint8_t base = (window_bits - 4) / 2;
    const uint8_t base = 4;
    lz77_binheap_init(&lz->bh_txt, 0x80); // ascii text
    lz77_binheap_init(&lz->bh_pos, (int32_t)window);
    lz77_binheap_init(&lz->bh_len, (int32_t)window);
    size_t i = 0; // output data[i]
    while (i < bytes) {
        uint64_t bit0 = lz77_read_bit(lz, &b64, &bp);
        lz77_if_error_return(lz);
        if (bit0) {
            uint64_t bit1 = lz77_read_bit(lz, &b64, &bp);
            lz77_if_error_return(lz);
            if (bit1) {
                uint64_t pos = lz77_read_number(lz, &b64, &bp, base);
                lz77_if_error_return(lz);
                pos = lz->bh_pos.ns[pos];
                lz77_binheap_inc_freq(&lz->bh_pos, (int32_t)pos);
                uint64_t long_len = lz77_read_bit(lz, &b64, &bp);
                lz77_if_error_return(lz);
                uint64_t len = 0;
                if (long_len) {
                    len = lz77_read_number(lz, &b64, &bp, base);
                    lz77_if_error_return(lz);
                } else {
                    len = lz77_read_number(lz, &b64, &bp, base);
                    lz77_if_error_return(lz);
                    len = lz->bh_len.ns[len];
                    lz77_binheap_inc_freq(&lz->bh_len, (int32_t)len);
                }
//lz77_println("i: %lld pos: %lld len: %lld", i, pos, len);
                rt_assert(0 < pos && pos < window);
                if (!(0 < pos && pos < window)) { lz77_return_invalid(lz); }
                rt_assert(0 < len);
                if (len == 0) { lz77_return_invalid(lz); }
                // Cannot do memcpy() here because of possible overlap.
                // memcpy() may read more than one byte at a time.
                uint8_t* s = data - (size_t)pos;
                const size_t n = i + (size_t)len;
                while (i < n) { data[i] = s[i]; i++; }
            } else { // byte >= 0x80
//              uint8_t b = (uint8_t)lz77_read_bits(lz, &b64, &bp, 7);
                uint8_t b = (uint8_t)lz77_read_number(lz, &b64, &bp, 2);
                lz77_if_error_return(lz);
                uint8_t s = (uint8_t)lz->bh_txt.ns[b]; // symbol
                data[i] = 0x80 | s;
//lz77_println("i: %lld byte: %08X %c", i, s, s);
                lz77_binheap_inc_freq(&lz->bh_txt, s);
                i++;
            }
        } else { // ASCII byte < 0x80
//          uint8_t b = (uint8_t)lz77_read_bits(lz, &b64, &bp, 7);
            uint8_t b = (uint8_t)lz77_read_number(lz, &b64, &bp, 2);
            lz77_if_error_return(lz);
            uint8_t s = (uint8_t)lz->bh_txt.ns[b]; // symbol
//lz77_println("i: %lld byte: %08X %c", i, s, s);
            data[i] = s;
            lz77_binheap_inc_freq(&lz->bh_txt, s);
            i++;
        }
    }
}

lz77_if lz77 = {
    .write_header = lz77_write_header,
    .compress     = lz77_compress,
    .read_header  = lz77_read_header,
    .decompress   = lz77_decompress,
};

#endif // lz77_implementation

#if 0

Debug / base:3

compress    4096 ->       8   0.2%
compress    4096 ->       8   0.2%
compress      35 ->      24  68.6%
 verify   decompressed: Hello World Hello.World Hello World
 compress    7774 ->    2960  38.1% of "C:\Users\leo\github\leok7v\lz77\test.c"
 compress  314160 ->  114840  36.6% of "test/rt.h"
 compress  611835 ->  212888  34.8% of "test/ui.h"
 compress 8182289 -> 3519720  43.0% of "test/sqlite3.c"
 compress 1790976 ->  414384  23.1% of "C:\Users\leo\github\leok7v\lz77\bin\ARM64\debug\lz77.exe"


#endif