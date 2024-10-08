#ifndef lz77_definition
#define lz77_definition

#include <errno.h>
#include <math.h>
#include <stdint.h>

// Naive LZ77 implementation inspired by CharGPT discussion
// and my personal passion to compressors in 198x

typedef struct lz77_s lz77_t;

typedef struct lz77_s {
    // `that` see: https://gist.github.com/leok7v/8d118985d3236b0069d419166f4111cf
    void*    that;  // caller supplied data
    errno_t  error; // sticky; for read()/write() compress() and decompress()
    // caller supplied read()/write() must error via .error field
    uint64_t (*read)(lz77_t*); //  reads 64 bits
    void     (*write)(lz77_t*, uint64_t b64); // writes 64 bits
    uint64_t written;
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
#define lz77_histogram_pos_len(pos, bytes) do { } while (0)
#define lz77_dump_histograms()           do { } while (0)

#endif

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

typedef uint8_t map_entry_t[256]; // data[0] number of bytes [2..255]

typedef struct {
    map_entry_t entry[512 * 1024]; // 128MB
    int32_t entries;
    int32_t max_chain;
    int32_t max_bytes;
} map_t;

static uint32_t map_hash32(const uint8_t* data, int64_t bytes) {
    uint32_t hash  = 0x811c9dc5; // FNV_offset_basis for 32-bit
    uint32_t prime = 0x01000193; // FNV_prime for 32-bit
    if (bytes > 0) {
        for (int64_t i = 1; i < bytes; i++) {
            hash ^= (uint32_t)data[i];
            hash *= prime;
        }
    } else {
        for (int64_t i = 0; data[i] != 0; i++) {
            hash ^= (uint32_t)data[i];
            hash *= prime;
        }
    }
    return hash;
}

static uint64_t map_hash64(const uint8_t* data, int64_t bytes) {
    uint64_t hash  = 0xcbf29ce484222325; // FNV_offset_basis for 64-bit
    uint64_t prime = 0x100000001b3;      // FNV_prime for 64-bit
    if (bytes > 0) {
        for (int64_t i = 0; i < bytes; i++) {
            hash ^= (uint64_t)data[i];
            hash *= prime;
        }
    } else {
        for (int64_t i = 0; data[i] != 0; i++) {
            hash ^= (uint64_t)data[i];
            hash *= prime;
        }
    }
    return hash;
}

static void map_init(map_t* map) {
    for (int32_t i = 0; i < rt_countof(map->entry); i++) {
        map->entry[i][0] = 0;
    }
    map->entries = 0;
    map->max_chain = 0;
    map->max_bytes = 0;
}

static const uint8_t* map_get(const map_t* map, const uint8_t* data, uint8_t bytes) {
    uint64_t hash = map_hash64(data, bytes);
    size_t i = (size_t)hash % rt_countof(map->entry);
    while (map->entry[i][0] > 0) {
        if (map->entry[i][0] == bytes && memcmp(&map->entry[i][1], data, bytes) == 0) {
            return &map->entry[i][1];
        }
        i = (i + 1) % rt_countof(map->entry);
    }
    return null;
}

static void map_put(map_t* map, const uint8_t* data, uint8_t bytes) {
    rt_swear(2 <= bytes && bytes < rt_countof(map->entry[0]));
    rt_swear(map->entries < rt_countof(map->entry) * 3 / 4);
    uint64_t hash = map_hash64(data, bytes);
    size_t i = (size_t)hash % rt_countof(map->entry);
    int32_t rehash = 0;
    while (map->entry[i][0] > 0) {
        if (map->entry[i][0] == bytes &&
            memcmp(&map->entry[i][1], data, bytes) == 0) {
            return; // already exists
        }
        rehash++;
        i = (i + 1) % rt_countof(map->entry);
    }
    if (rehash > map->max_chain) { map->max_chain = rehash; }
    if (bytes  > map->max_bytes) { map->max_bytes = bytes; }
    map->entry[i][0] = bytes;
    memcpy(&map->entry[i][1], data, bytes);
    map->entries++;
}

static void map_clear(map_t *map) {
    for (int32_t i = 0; i < rt_countof(map->entry); i++) {
        map->entry[i][0] = 0;
    }
    map->entries = 0;
    map->max_chain = 0;
}

static uint64_t pos_freq[64 * 1024];
static uint64_t len_freq[64 * 1024];
static map_t map;  // word map
static map_t lens; // different bytes encountered
static map_t poss; // different pos encountered

static double lz77_entropy(uint64_t freq[], int32_t n) {
    double total = 0;
    double aha_entropy = 0.0;
    for (int32_t i = 0; i < n; i++) { total += (double)freq[i]; }
    for (int32_t i = 0; i < n; i++) {
        if (freq[i] > 0) {
            double p_i = (double)freq[i] / total;
            aha_entropy += p_i * log2(p_i);
        }
    }
    return -aha_entropy;
}

static void lz77_compress(lz77_t* lz, const uint8_t* data, size_t bytes,
        uint8_t window_bits) {
memset(pos_freq, 0x00, sizeof(pos_freq));
memset(len_freq, 0x00, sizeof(pos_freq));
map_init(&map);
map_init(&lens);
map_init(&poss);
    lz77_if_error_return(lz);
    if (window_bits < 10 || window_bits > 20) { lz77_return_invalid(lz); }
    lz77_init_histograms();
    const size_t window = ((size_t)1U) << window_bits;
    const uint8_t base = (window_bits - 4) / 2;
    uint64_t b64 = 0;
    uint32_t bp = 0;
    size_t i = 0;
    while (i < bytes) {
        // bytes and position of longest matching sequence
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
            rt_assert(0 < pos && pos < window);
            rt_assert(0 < len);
            lz77_write_bits(lz, &b64, &bp, 0b11, 2); // flags
            lz77_if_error_return(lz);
            lz77_write_number(lz, &b64, &bp, pos, base);
            lz77_if_error_return(lz);
            lz77_write_number(lz, &b64, &bp, len, base);
            lz77_if_error_return(lz);
            lz77_histogram_pos_len(pos, len);
if (len < window) { len_freq[len]++; }
if (pos < window) { pos_freq[pos]++; }
// lz77_println("\"%.*s\"", bytes, &data[i]);
if (len <= 255) {
    map_put(&map, &data[i], (uint8_t)len);
}
map_put(&lens, (uint8_t*)&len, (uint8_t)sizeof(len));
map_put(&poss, (uint8_t*)&pos, (uint8_t)sizeof(pos));
            i += len;
        } else {
            const uint8_t b = data[i];
            // European texts are predominantly spaces and small ASCII letters:
            if (b < 0x80) {
                lz77_write_bit(lz, &b64, &bp, 0); // flags
                lz77_if_error_return(lz);
                // ASCII byte < 0x80 with 8th bit set to `0`
                lz77_write_bits(lz, &b64, &bp, b, 7);
                lz77_if_error_return(lz);
            } else {
                lz77_write_bit(lz, &b64, &bp, 1); // flag: 1
                lz77_write_bit(lz, &b64, &bp, 0); // flag: 0
                lz77_if_error_return(lz);
                // only 7 bit because 8th bit is `1`
                lz77_write_bits(lz, &b64, &bp, b, 7);
                lz77_if_error_return(lz);
            }
            i++;
        }
    }
    lz77_flush(lz, b64, bp);
    lz77_dump_histograms();
    double len_bits = lz77_entropy(len_freq, (int32_t)window);
    double pos_bits = lz77_entropy(pos_freq, (int32_t)window);
    lz77_println("bits len: %.2f pos: %.2f words: %d "
                 "max chain: %d max bytes: %d #len: %d #pos: %d",
        len_bits, pos_bits, map.entries, map.max_chain, map.max_bytes,
        lens.entries, poss.entries);
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
    const uint8_t base = (window_bits - 4) / 2;
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
                uint64_t len = lz77_read_number(lz, &b64, &bp, base);
                lz77_if_error_return(lz);
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
                uint64_t b = lz77_read_bits(lz, &b64, &bp, 7);
                lz77_if_error_return(lz);
                data[i] = (uint8_t)b | 0x80;
                i++;
            }
        } else { // literal byte
            uint64_t b = lz77_read_bits(lz, &b64, &bp, 7); // ASCII byte < 0x80
            lz77_if_error_return(lz);
            data[i] = (uint8_t)b;
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

