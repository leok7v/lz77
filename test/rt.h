#ifndef rt_definition
#define rt_definition

// _________________________________ rt_std.h _________________________________

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <locale.h>
#include <malloc.h>
#include <signal.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define rt_stringify(x) #x
#define rt_tostring(x) rt_stringify(x)
#define rt_pragma(x) _Pragma(rt_tostring(x))

#if defined(__GNUC__) || defined(__clang__) // TODO: remove and fix code
#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#pragma GCC diagnostic ignored "-Wdeclaration-after-statement"
#pragma GCC diagnostic ignored "-Wfour-char-constants"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wunsafe-buffer-usage"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Waddress-of-packed-member"
#pragma GCC diagnostic ignored "-Wused-but-marked-unused" // because in debug only
#define rt_msvc_pragma(x)
#define rt_gcc_pragma(x) rt_pragma(x)
#else
#define rt_gcc_pragma(x)
#define rt_msvc_pragma(x) rt_pragma(x)
#endif

#ifdef _MSC_VER
    #define rt_suppress_constant_cond_exp _Pragma("warning(suppress: 4127)")
#else
    #define rt_suppress_constant_cond_exp
#endif

// Type aliases for floating-point types similar to <stdint.h>
typedef float  fp32_t;
typedef double fp64_t;
// "long fp64_t" is required by C standard but the bitness
// of it is not specified.

#ifdef __cplusplus
    #define rt_begin_c extern "C" {
    #define rt_end_c } // extern "C"
#else
    #define rt_begin_c // C headers compiled as C++
    #define rt_end_c
#endif

// rt_countof() and rt_countof() are suitable for
// small < 2^31 element arrays

#define rt_countof(a) ((int32_t)((int)(sizeof(a) / sizeof((a)[0]))))

#if defined(__GNUC__) || defined(__clang__)
    #define rt_force_inline __attribute__((always_inline))
#elif defined(_MSC_VER)
    #define rt_force_inline __forceinline
#endif

#ifndef __cplusplus
    #define null ((void*)0) // better than NULL which is zero
#else
    #define null nullptr
#endif

#if defined(_MSC_VER)
    #define rt_thread_local __declspec(thread)
#else
    #ifndef __cplusplus
        #define rt_thread_local _Thread_local // C99
    #else
        // C++ supports rt_thread_local keyword
    #endif
#endif

// rt_begin_packed rt_end_packed
// usage: typedef rt_begin_packed struct foo_s { ... } rt_end_packed foo_t;

#if defined(__GNUC__) || defined(__clang__)
#define rt_attribute_packed __attribute__((packed))
#define rt_begin_packed
#define rt_end_packed rt_attribute_packed
#else
#define rt_begin_packed rt_pragma( pack(push, 1) )
#define rt_end_packed rt_pragma( pack(pop) )
#define rt_attribute_packed
#endif

// usage: typedef struct rt_aligned_8 foo_s { ... } foo_t;

#if defined(__GNUC__) || defined(__clang__)
#define rt_aligned_8 __attribute__((aligned(8)))
#elif defined(_MSC_VER)
#define rt_aligned_8 __declspec(align(8))
#else
#define rt_aligned_8
#endif


// In callbacks the formal parameters are
// frequently unused. Also sometimes parameters
// are used in debug configuration only (e.g. rt_assert() checks)
// but not in release.
// C does not have anonymous parameters like C++
// Instead of:
//      void foo(param_type_t param) { (void)param; / *unused */ }
// use:
//      vod foo(param_type_t rt_unused(param)) { }

#if defined(__GNUC__) || defined(__clang__)
#define rt_unused(name) name __attribute__((unused))
#elif defined(_MSC_VER)
#define rt_unused(name) _Pragma("warning(suppress:  4100)") name
#else
#define rt_unused(name) name
#endif

// Because MS C compiler is unhappy about alloca() and
// does not implement (C99 optional) dynamic arrays on the stack:

#define rt_stackalloc(n) (_Pragma("warning(suppress: 6255 6263)") alloca(n))

// alloca() is messy and in general is a not a good idea.
// try to avoid if possible. Stack sizes vary from 64KB to 8MB in 2024.

// _________________________________ rt_str.h _________________________________

rt_begin_c

typedef struct rt_str64_t {
    char s[64];
} rt_str64_t;

typedef struct rt_str128_t {
    char s[128];
} rt_str128_t;

typedef struct rt_str1024_t {
    char s[1024];
} rt_str1024_t;

typedef struct rt_str32K_t {
    char s[32 * 1024];
} rt_str32K_t;

// truncating string printf:
// char s[100]; rt_str_printf(s, "Hello %s", "world");
// do not use with char* and char s[] parameters
// because rt_countof(s) will be sizeof(char*) not the size of the buffer.

#define rt_str_printf(s, ...) rt_str.format((s), rt_countof(s), "" __VA_ARGS__)

#define rt_strerr(r) (rt_str.error((r)).s) // use only as rt_str_printf() parameter

// The strings are expected to be UTF-8 encoded.
// Copy functions fatal fail if the destination buffer is too small.
// It is responsibility of the caller to make sure it won't happen.

typedef struct {
    char* (*drop_const)(const char* s); // because of strstr() and alike
    int32_t (*len)(const char* s);
    int32_t (*len16)(const uint16_t* utf16);
    int32_t (*utf8bytes)(const char* utf8, int32_t bytes); // 0 on error
    int32_t (*glyphs)(const char* utf8, int32_t bytes); // -1 on error
    bool (*starts)(const char* s1, const char* s2);  // s1 starts with s2
    bool (*ends)(const char* s1, const char* s2);    // s1 ends with s2
    bool (*istarts)(const char* s1, const char* s2); // ignore case
    bool (*iends)(const char* s1, const char* s2);   // ignore case
    // string truncation is fatal use strlen() to check at call site
    void (*lower)(char* d, int32_t capacity, const char* s); // ASCII only
    void (*upper)(char* d, int32_t capacity, const char* s); // ASCII only
    // utf8/utf16 conversion
    // If `chars` argument is -1, the function utf8_bytes includes the zero
    // terminating character in the conversion and the returned byte count.
    int32_t (*utf8_bytes)(const uint16_t* utf16, int32_t bytes); // bytes count
    // If `bytes` argument is -1, the function utf16_chars() includes the zero
    // terminating character in the conversion and the returned character count.
    int32_t (*utf16_chars)(const char* utf8, int32_t bytes); // chars count
    // utf8_bytes() and utf16_chars() return -1 on invalid UTF-8/UTF-16
    // utf8_bytes(L"", -1) returns 1 for zero termination
    // utf16_chars("", -1) returns 1 for zero termination
    // chars: -1 means both source and destination are zero terminated
    errno_t (*utf16to8)(char* utf8, int32_t capacity,
                        const uint16_t* utf16, int32_t chars);
    // bytes: -1 means both source and destination are zero terminated
    errno_t (*utf8to16)(uint16_t* utf16, int32_t capacity,
                        const char* utf8, int32_t bytes);
    // https://compart.com/en/unicode/U+1F41E
    // Lady Beetle: utf16 L"\xD83D\xDC1E" utf8 "\xF0\x9F\x90\x9E"
    //           surrogates:  high   low
    bool (*utf16_is_low_surrogate)(uint16_t utf16char);
    bool (*utf16_is_high_surrogate)(uint16_t utf16char);
    uint32_t (*utf32)(const char* utf8, int32_t bytes); // single codepoint
    // string formatting printf style:
    void (*format_va)(char* utf8, int32_t count, const char* format, va_list va);
    void (*format)(char* utf8, int32_t count, const char* format, ...);
    // format "dg" digit grouped; see below for known grouping separators:
    const char* (*grouping_separator)(void); // locale
    // Returned const char* pointer is short-living thread local and
    // intended to be used in the arguments list of .format() or .printf()
    // like functions, not stored or passed for prolonged call chains.
    // See implementation for details.
    rt_str64_t (*int64_dg)(int64_t v, bool uint, const char* gs);
    rt_str64_t (*int64)(int64_t v);      // with UTF-8 thin space
    rt_str64_t (*uint64)(uint64_t v);    // with UTF-8 thin space
    rt_str64_t (*int64_lc)(int64_t v);   // with locale separator
    rt_str64_t (*uint64_lc)(uint64_t v); // with locale separator
    rt_str128_t (*fp)(const char* format, fp64_t v); // respects locale
    // errors to strings
    rt_str1024_t (*error)(int32_t error);     // en-US
    rt_str1024_t (*error_nls)(int32_t error); // national locale string
    void (*test)(void);
} rt_str_if;

// Known grouping separators
// https://en.wikipedia.org/wiki/Decimal_separator#Digit_grouping
// coma "," separated decimal
// other commonly used separators:
// underscore "_" (Fortran, Kotlin)
// apostrophe "'" (C++14, Rebol, Red)
// backtick "`"
// space "\x20"
// thin_space "\xE2\x80\x89" Unicode: U+2009

extern rt_str_if rt_str;

rt_end_c
// ___________________________________ rt.h ___________________________________

// the rest is in alphabetical order (no inter dependencies)

// ________________________________ rt_args.h _________________________________

rt_begin_c

typedef struct {
    // On Unix it is responsibility of the main() to assign these values
    int32_t c;      // argc
    const char** v; // argv[argc]
    const char** env; // rt_args.env[] is null-terminated
    void    (*main)(int32_t argc, const char* argv[], const char** env);
    void    (*WinMain)(void); // windows specific
    int32_t (*option_index)(const char* option); // e.g. option: "--verbosity" or "-v"
    void    (*remove_at)(int32_t ix);
    /* argc=3 argv={"foo", "--verbose"} -> returns true; argc=1 argv={"foo"} */
    bool    (*option_bool)(const char* option);
    /* argc=3 argv={"foo", "--n", "153"} -> value==153, true; argc=1 argv={"foo"}
       also handles negative values (e.g. "-153") and hex (e.g. 0xBADF00D)
    */
    bool    (*option_int)(const char* option, int64_t *value);
    // for argc=3 argv={"foo", "--path", "bar"}
    //     option_str("--path", option)
    // returns option: "bar" and argc=1 argv={"foo"} */
    const char* (*option_str)(const char* option);
    // basename() for argc=3 argv={"/bin/foo.exe", ...} returns "foo":
    const char* (*basename)(void);
    void (*fini)(void);
    void (*test)(void);
} rt_args_if;

extern rt_args_if rt_args;

/* Usage:

    (both main() and WinMain() could be compiled at the same time on Windows):

    static int run(void);

    int main(int argc, char* argv[], char* envp[]) { // link.exe /SUBSYSTEM:CONSOLE
        rt_args.main(argc, argv, envp); // Initialize args with command-line parameters
        int r = run();
        rt_args.fini();
        return r;
    }

    #include "rt/rt_win32.h"

    int APIENTRY WinMain(HINSTANCE inst, HINSTANCE prev, char* command, int show) {
        // link.exe /SUBSYSTEM:WINDOWS
        rt_args.WinMain();
        int r = run();
        rt_args.fini();
        return 0;
    }

    static int run(void) {
        if (rt_args.option_bool("-v")) {
            rt_debug.verbosity.level = rt_debug.verbosity.verbose;
        }
        int64_t num = 0;
        if (rt_args.option_int("--number", &num)) {
            printf("--number: %ld\n", num);
        }
        const char* path = rt_args.option_str("--path");
        if (path != null) {
            printf("--path: %s\n", path);
        }
        printf("rt_args.basename(): %s\n", rt_args.basename());
        printf("rt_args.v[0]: %s\n", rt_args.v[0]);
        for (int i = 1; i < rt_args.c; i++) {
            printf("rt_args.v[%d]: %s\n", i, rt_args.v[i]);
        }
        return 0;
    }

    // Also see: https://github.com/leok7v/ut/blob/main/test/test1.c

*/

rt_end_c

// ______________________________ rt_backtrace.h ______________________________

// "bt" stands for Stack Back Trace (not British Telecom)

rt_begin_c

enum { rt_backtrace_max_depth = 32 };    // increase if not enough
enum { rt_backtrace_max_symbol = 1024 }; // MSFT symbol size limit

typedef struct thread_s* rt_thread_t;

typedef char rt_backtrace_symbol_t[rt_backtrace_max_symbol];
typedef char rt_backtrace_file_t[512];

typedef struct rt_backtrace_s {
    int32_t frames; // 0 if capture() failed
    uint32_t hash;
    errno_t  error; // last error set by capture() or symbolize()
    void* stack[rt_backtrace_max_depth];
    rt_backtrace_symbol_t symbol[rt_backtrace_max_depth];
    rt_backtrace_file_t file[rt_backtrace_max_depth];
    int32_t line[rt_backtrace_max_depth];
    bool symbolized;
} rt_backtrace_t;

//  calling .trace(bt, /*stop:*/"*")
//  stops backtrace dumping at any of the known Microsoft runtime
//  symbols:
//    "main",
//    "WinMain",
//    "BaseThreadInitThunk",
//    "RtlUserThreadStart",
//    "mainCRTStartup",
//    "WinMainCRTStartup",
//    "invoke_main"
// .trace(bt, null)
// provides complete backtrace to the bottom of stack

typedef struct {
    void (*capture)(rt_backtrace_t *bt, int32_t skip); // number of frames to skip
    void (*context)(rt_thread_t thread, const void* context, rt_backtrace_t *bt);
    void (*symbolize)(rt_backtrace_t *bt);
    // dump backtrace into rt_println():
    void (*trace)(const rt_backtrace_t* bt, const char* stop);
    void (*trace_self)(const char* stop);
    void (*trace_all_but_self)(void); // trace all threads
    const char* (*string)(const rt_backtrace_t* bt, char* text, int32_t count);
    void (*test)(void);
} rt_backtrace_if;

extern rt_backtrace_if rt_backtrace;

#define rt_backtrace_here() do {   \
    rt_backtrace_t bt_ = {0};      \
    rt_backtrace.capture(&bt_, 0); \
    rt_backtrace.symbolize(&bt_);  \
    rt_backtrace.trace(&bt_, "*"); \
} while (0)

rt_end_c

// _______________________________ rt_atomics.h _______________________________

// Will be deprecated soon after Microsoft fully supports <stdatomic.h>

rt_begin_c

typedef struct {
    void* (*exchange_ptr)(volatile void** a, void* v); // retuns previous value
    int32_t (*increment_int32)(volatile int32_t* a); // returns incremented
    int32_t (*decrement_int32)(volatile int32_t* a); // returns decremented
    int64_t (*increment_int64)(volatile int64_t* a); // returns incremented
    int64_t (*decrement_int64)(volatile int64_t* a); // returns decremented
    int32_t (*add_int32)(volatile int32_t* a, int32_t v); // returns result of add
    int64_t (*add_int64)(volatile int64_t* a, int64_t v); // returns result of add
    // returns the value held previously by "a" address:
    int32_t (*exchange_int32)(volatile int32_t* a, int32_t v);
    int64_t (*exchange_int64)(volatile int64_t* a, int64_t v);
    // compare_exchange functions compare the *a value with the comparand value.
    // If the *a is equal to the comparand value, the "v" value is stored in the address
    // specified by "a" otherwise, no operation is performed.
    // returns true if previous value *a was the same as "comparand"
    // false if *a was different from "comparand" and "a" was NOT modified.
    bool (*compare_exchange_int64)(volatile int64_t* a, int64_t comparand, int64_t v);
    bool (*compare_exchange_int32)(volatile int32_t* a, int32_t comparand, int32_t v);
    bool (*compare_exchange_ptr)(volatile void** a, void* comparand, void* v);
    void (*spinlock_acquire)(volatile int64_t* spinlock);
    void (*spinlock_release)(volatile int64_t* spinlock);
    int32_t (*load32)(volatile int32_t* a);
    int64_t (*load64)(volatile int64_t* a);
    void (*memory_fence)(void);
    void (*test)(void);
} rt_atomics_if;

extern rt_atomics_if rt_atomics;

rt_end_c

// ______________________________ rt_clipboard.h ______________________________

rt_begin_c

typedef struct ui_bitmap_s ui_bitmap_t;

typedef struct {
    errno_t (*put_text)(const char* s);
    errno_t (*get_text)(char* text, int32_t* bytes);
    errno_t (*put_image)(ui_bitmap_t* image); // only for Windows apps
    void (*test)(void);
} rt_clipboard_if;

extern rt_clipboard_if rt_clipboard;

rt_end_c

// ________________________________ rt_clock.h ________________________________

rt_begin_c

typedef struct {
    int32_t const nsec_in_usec; // nano in micro second
    int32_t const nsec_in_msec; // nano in milli
    int32_t const nsec_in_sec;
    int32_t const usec_in_msec; // micro in milli
    int32_t const msec_in_sec;  // milli in sec
    int32_t const usec_in_sec;  // micro in sec
    fp64_t   (*seconds)(void);      // since boot
    uint64_t (*nanoseconds)(void);  // since boot overflows in about 584.5 years
    uint64_t (*unix_microseconds)(void); // since January 1, 1970
    uint64_t (*unix_seconds)(void);      // since January 1, 1970
    uint64_t (*microseconds)(void); // NOT monotonic(!) UTC since epoch January 1, 1601
    uint64_t (*localtime)(void);    // local time microseconds since epoch
    void (*utc)(uint64_t microseconds, int32_t* year, int32_t* month,
        int32_t* day, int32_t* hh, int32_t* mm, int32_t* ss, int32_t* ms,
        int32_t* mc);
    void (*local)(uint64_t microseconds, int32_t* year, int32_t* month,
        int32_t* day, int32_t* hh, int32_t* mm, int32_t* ss, int32_t* ms,
        int32_t* mc);
    void (*test)(void);
} rt_clock_if;

extern rt_clock_if rt_clock;

rt_end_c


// _______________________________ rt_config.h ________________________________

rt_begin_c

// Persistent storage for configuration and other small data
// related to specific application.
// on Unix-like system ~/.name/key files are used.
// On Window User registry (could be .dot files/folders).
// "name" is customary basename of "rt_args.v[0]"

typedef struct {
    errno_t (*save)(const char* name, const char* key,
                    const void* data, int32_t bytes);
    int32_t (*size)(const char* name, const char* key);
    // load() returns number of actual loaded bytes:
    int32_t (*load)(const char* name, const char* key,
                    void* data, int32_t bytes);
    errno_t (*remove)(const char* name, const char* key);
    errno_t (*clean)(const char* name); // remove all subkeys
    void (*test)(void);
} rt_config_if;

extern rt_config_if rt_config;

rt_end_c


// ________________________________ rt_core.h _________________________________

rt_begin_c

typedef struct {
    errno_t (*err)(void); // errno or GetLastError()
    void (*set_err)(errno_t err); // errno = err or SetLastError()
    void (*abort)(void);
    void (*exit)(int32_t exit_code); // only 8 bits on posix
    void (*test)(void);
    struct {                                // posix
        errno_t const access_denied;        // EACCES
        errno_t const bad_file;             // EBADF
        errno_t const broken_pipe;          // EPIPE
        errno_t const device_not_ready;     // ENXIO
        errno_t const directory_not_empty;  // ENOTEMPTY
        errno_t const disk_full;            // ENOSPC
        errno_t const file_exists;          // EEXIST
        errno_t const file_not_found;       // ENOENT
        errno_t const insufficient_buffer;  // E2BIG
        errno_t const interrupted;          // EINTR
        errno_t const invalid_data;         // EINVAL
        errno_t const invalid_handle;       // EBADF
        errno_t const invalid_parameter;    // EINVAL
        errno_t const io_error;             // EIO
        errno_t const more_data;            // ENOBUFS
        errno_t const name_too_long;        // ENAMETOOLONG
        errno_t const no_child_process;     // ECHILD
        errno_t const not_a_directory;      // ENOTDIR
        errno_t const not_empty;            // ENOTEMPTY
        errno_t const out_of_memory;        // ENOMEM
        errno_t const path_not_found;       // ENOENT
        errno_t const pipe_not_connected;   // EPIPE
        errno_t const read_only_file;       // EROFS
        errno_t const resource_deadlock;    // EDEADLK
        errno_t const too_many_open_files;  // EMFILE
    } const error;
} rt_core_if;

extern rt_core_if rt_core;

rt_end_c

// ________________________________ rt_debug.h ________________________________

rt_begin_c

// debug interface essentially is:
// vfprintf(stderr, format, va)
// fprintf(stderr, format, ...)
// with the additional convience:
// 1. writing to system log (e.g. OutputDebugString() on Windows)
// 2. always appending \n at the end of the line and thus flushing buffer
// Warning: on Windows it is not real-time and subject to 30ms delays
//          that may or may not happen on some calls

typedef struct {
    int32_t level; // global verbosity (interpretation may vary)
    int32_t const quiet;    // 0
    int32_t const info;     // 1 basic information (errors and warnings)
    int32_t const verbose;  // 2 detailed diagnostic
    int32_t const debug;    // 3 including debug messages
    int32_t const trace;    // 4 everything (may include nested calls)
} rt_verbosity_if;

typedef struct {
    rt_verbosity_if verbosity;
    int32_t (*verbosity_from_string)(const char* s);
    // "T" connector for outside write; return false to proceed with default
    bool (*tee)(const char* s, int32_t count); // return true to intercept
    void (*output)(const char* s, int32_t count);
    void (*println_va)(const char* file, int32_t line, const char* func,
        const char* format, va_list va);
    void (*println)(const char* file, int32_t line, const char* func,
        const char* format, ...);
    void (*perrno)(const char* file, int32_t line,
        const char* func, int32_t err_no, const char* format, ...);
    void (*perror)(const char* file, int32_t line,
        const char* func, int32_t error, const char* format, ...);
    bool (*is_debugger_present)(void);
    void (*breakpoint)(void); // no-op if debugger is not present
    errno_t (*raise)(uint32_t exception);
    struct  {
        uint32_t const access_violation;
        uint32_t const datatype_misalignment;
        uint32_t const breakpoint;
        uint32_t const single_step;
        uint32_t const array_bounds;
        uint32_t const float_denormal_operand;
        uint32_t const float_divide_by_zero;
        uint32_t const float_inexact_result;
        uint32_t const float_invalid_operation;
        uint32_t const float_overflow;
        uint32_t const float_stack_check;
        uint32_t const float_underflow;
        uint32_t const int_divide_by_zero;
        uint32_t const int_overflow;
        uint32_t const priv_instruction;
        uint32_t const in_page_error;
        uint32_t const illegal_instruction;
        uint32_t const noncontinuable;
        uint32_t const stack_overflow;
        uint32_t const invalid_disposition;
        uint32_t const guard_page;
        uint32_t const invalid_handle;
        uint32_t const possible_deadlock;
    } exception;
    void (*test)(void);
} rt_debug_if;

#define rt_println(...) rt_debug.println(__FILE__, __LINE__, __func__, "" __VA_ARGS__)

extern rt_debug_if rt_debug;

rt_end_c

// ________________________________ rt_files.h ________________________________

rt_begin_c

// space for "short" 260 utf-16 characters filename in utf-8 format:
typedef struct rt_file_name_s { char s[1024]; } rt_file_name_t;

enum { rt_files_max_path = (int32_t)sizeof(rt_file_name_t) }; // *)

typedef struct rt_file_s rt_file_t;

typedef struct rt_files_stat_s {
    uint64_t created;
    uint64_t accessed;
    uint64_t updated;
    int64_t  size; // bytes
    int64_t  type; // device / folder / symlink
} rt_files_stat_t;

typedef struct rt_folder_s {
    uint8_t data[512]; // implementation specific
} rt_folder_t;

typedef struct {
    rt_file_t* const invalid; // (rt_file_t*)-1
    // rt_files_stat_t.type:
    int32_t const type_folder;
    int32_t const type_symlink;
    int32_t const type_device;
    // seek() methods:
    int32_t const seek_set;
    int32_t const seek_cur;
    int32_t const seek_end;
    // open() flags
    int32_t const o_rd; // read only
    int32_t const o_wr; // write only
    int32_t const o_rw; // != (o_rd | o_wr)
    int32_t const o_append;
    int32_t const o_create; // opens existing or creates new
    int32_t const o_excl;   // create fails if file exists
    int32_t const o_trunc;  // open always truncates to empty
    int32_t const o_sync;
    // known folders ids:
    struct { // known folders:
        int32_t const home;      // "c:\Users\<username>" or "/users/<username>"
        int32_t const desktop;
        int32_t const documents;
        int32_t const downloads;
        int32_t const music;
        int32_t const pictures;
        int32_t const videos;
        int32_t const shared;    // "c:\Users\Public"
        int32_t const bin;       // "c:\Program Files" aka "/bin" aka "/Applications"
        int32_t const data;      // "c:\ProgramData" aka "/var" aka "/data"
    } const folder;
    errno_t (*open)(rt_file_t* *file, const char* filename, int32_t flags);
    bool    (*is_valid)(rt_file_t* file); // checks both null and invalid
    errno_t (*seek)(rt_file_t* file, int64_t *position, int32_t method);
    errno_t (*stat)(rt_file_t* file, rt_files_stat_t* stat, bool follow_symlink);
    errno_t (*read)(rt_file_t* file, void* data, int64_t bytes, int64_t *transferred);
    errno_t (*write)(rt_file_t* file, const void* data, int64_t bytes, int64_t *transferred);
    errno_t (*flush)(rt_file_t* file);
    void    (*close)(rt_file_t* file);
    errno_t (*write_fully)(const char* filename, const void* data,
                           int64_t bytes, int64_t *transferred);
    bool (*exists)(const char* pathname); // does not guarantee any access writes
    bool (*is_folder)(const char* pathname);
    bool (*is_symlink)(const char* pathname);
    const char* (*basename)(const char* pathname); // c:\foo\bar.ext -> bar.ext
    errno_t (*mkdirs)(const char* pathname); // tries to deep create all folders in pathname
    errno_t (*rmdirs)(const char* pathname); // tries to remove folder and its subtree
    errno_t (*create_tmp)(char* file, int32_t count); // create temporary file
    errno_t (*chmod777)(const char* pathname); // and whole subtree new files and folders
    errno_t (*symlink)(const char* from, const char* to); // sym link "ln -s" **)
    errno_t (*link)(const char* from, const char* to); // hard link like "ln"
    errno_t (*unlink)(const char* pathname); // delete file or empty folder
    errno_t (*copy)(const char* from, const char* to); // allows overwriting
    errno_t (*move)(const char* from, const char* to); // allows overwriting
    errno_t (*cwd)(char* folder, int32_t count); // get current working dir
    errno_t (*chdir)(const char* folder); // set working directory
    const char* (*known_folder)(int32_t kf_id);
    const char* (*bin)(void);  // Windows: "c:\Program Files" Un*x: "/bin"
    const char* (*data)(void); // Windows: "c:\ProgramData" Un*x: /data or /var
    const char* (*tmp)(void);  // temporary folder (system or user)
    // There are better, native, higher performance ways to iterate thru
    // folders in Posix, Linux and Windows. The following is minimalistic
    // approach to folder content reading:
    errno_t (*opendir)(rt_folder_t* folder, const char* folder_name);
    const char* (*readdir)(rt_folder_t* folder, rt_files_stat_t* optional);
    void (*closedir)(rt_folder_t* folder);
    void (*test)(void);
} rt_files_if;

// *) rt_files_max_path is a compromise - way longer than Windows MAX_PATH of 260
// and somewhat shorter than 32 * 1024 Windows long path.
// Use with caution understanding that it is a limitation and where it is
// important heap may and should be used. Do not to rely on thread stack size
// (default: 1MB on Windows, Android Linux 64KB, 512 KB  on MacOS, Ubuntu 8MB)
//
// **) symlink on Win32 is only allowed in Admin elevated
//     processes and in Developer mode.

extern rt_files_if rt_files;

rt_end_c

// ______________________________ rt_generics.h _______________________________

rt_begin_c

// Most of ut/ui code is written the way of min(a,b) max(a,b)
// not having side effects on the arguments and thus evaluating
// them twice ain't a big deal. However, out of curiosity of
// usefulness of Generic() in C11 standard here is type-safe
// single evaluation of the arguments version of what could
// have been simple minimum and maximum macro definitions.
// Type safety comes with the cost of complexity in puritan
// or stoic language like C:

static inline int8_t   rt_max_int8(int8_t x, int8_t y)       { return x > y ? x : y; }
static inline int16_t  rt_max_int16(int16_t x, int16_t y)    { return x > y ? x : y; }
static inline int32_t  rt_max_int32(int32_t x, int32_t y)    { return x > y ? x : y; }
static inline int64_t  rt_max_int64(int64_t x, int64_t y)    { return x > y ? x : y; }
static inline uint8_t  rt_max_uint8(uint8_t x, uint8_t y)    { return x > y ? x : y; }
static inline uint16_t rt_max_uint16(uint16_t x, uint16_t y) { return x > y ? x : y; }
static inline uint32_t rt_max_uint32(uint32_t x, uint32_t y) { return x > y ? x : y; }
static inline uint64_t rt_max_uint64(uint64_t x, uint64_t y) { return x > y ? x : y; }
static inline fp32_t   rt_max_fp32(fp32_t x, fp32_t y)       { return x > y ? x : y; }
static inline fp64_t   rt_max_fp64(fp64_t x, fp64_t y)       { return x > y ? x : y; }

// MS cl.exe version 19.39.33523 has issues with "long":
// does not pick up int32_t/uint32_t types for "long" and "unsigned long"
// need to handle long / unsigned long separately:

static inline long          rt_max_long(long x, long y)                    { return x > y ? x : y; }
static inline unsigned long rt_max_ulong(unsigned long x, unsigned long y) { return x > y ? x : y; }

static inline int8_t   rt_min_int8(int8_t x, int8_t y)       { return x < y ? x : y; }
static inline int16_t  rt_min_int16(int16_t x, int16_t y)    { return x < y ? x : y; }
static inline int32_t  rt_min_int32(int32_t x, int32_t y)    { return x < y ? x : y; }
static inline int64_t  rt_min_int64(int64_t x, int64_t y)    { return x < y ? x : y; }
static inline uint8_t  rt_min_uint8(uint8_t x, uint8_t y)    { return x < y ? x : y; }
static inline uint16_t rt_min_uint16(uint16_t x, uint16_t y) { return x < y ? x : y; }
static inline uint32_t rt_min_uint32(uint32_t x, uint32_t y) { return x < y ? x : y; }
static inline uint64_t rt_min_uint64(uint64_t x, uint64_t y) { return x < y ? x : y; }
static inline fp32_t   rt_min_fp32(fp32_t x, fp32_t y)       { return x < y ? x : y; }
static inline fp64_t   rt_min_fp64(fp64_t x, fp64_t y)       { return x < y ? x : y; }

static inline long          rt_min_long(long x, long y)                    { return x < y ? x : y; }
static inline unsigned long rt_min_ulong(unsigned long x, unsigned long y) { return x < y ? x : y; }


static inline void rt_min_undefined(void) { }
static inline void rt_max_undefined(void) { }

#define rt_max(X, Y) _Generic((X) + (Y), \
    int8_t:   rt_max_int8,   \
    int16_t:  rt_max_int16,  \
    int32_t:  rt_max_int32,  \
    int64_t:  rt_max_int64,  \
    uint8_t:  rt_max_uint8,  \
    uint16_t: rt_max_uint16, \
    uint32_t: rt_max_uint32, \
    uint64_t: rt_max_uint64, \
    fp32_t:   rt_max_fp32,   \
    fp64_t:   rt_max_fp64,   \
    long:     rt_max_long,   \
    unsigned long: rt_max_ulong, \
    default:  rt_max_undefined)(X, Y)

#define rt_min(X, Y) _Generic((X) + (Y), \
    int8_t:   rt_min_int8,   \
    int16_t:  rt_min_int16,  \
    int32_t:  rt_min_int32,  \
    int64_t:  rt_min_int64,  \
    uint8_t:  rt_min_uint8,  \
    uint16_t: rt_min_uint16, \
    uint32_t: rt_min_uint32, \
    uint64_t: rt_min_uint64, \
    fp32_t:   rt_min_fp32,   \
    fp64_t:   rt_min_fp64,   \
    long:     rt_min_long,   \
    unsigned long: rt_min_ulong, \
    default:  rt_min_undefined)(X, Y)

// The expression (X) + (Y) is used in _Generic primarily for type promotion
// and compatibility between different types of the two operands. In C, when
// you perform an arithmetic operation like addition between two variables,
// the types of these variables undergo implicit conversions to a common type
// according to the usual arithmetic conversions defined in the C standard.
// This helps ensure that:
//
// Type Compatibility: The macro works correctly even if X and Y are of
// different types. By using (X) + (Y), both X and Y are promoted to a common
// type, which ensures that the macro selects the appropriate function
// that can handle this common type.
//
// Type Safety: It ensures that the selected function can handle the type
// resulting from the operation, thereby preventing type mismatches that
// could lead to undefined behavior or compilation errors.

typedef struct {
    void (*test)(void);
} rt_generics_if;

extern rt_generics_if rt_generics;

rt_end_c

// _______________________________ rt_glyphs.h ________________________________

// Square Four Corners https://www.compart.com/en/unicode/U+26F6
#define rt_glyph_square_four_corners                    "\xE2\x9B\xB6"

// Circled Cross Formee
// https://codepoints.net/U+1F902
#define rt_glyph_circled_cross_formee                   "\xF0\x9F\xA4\x82"

// White Large Square https://www.compart.com/en/unicode/U+2B1C
#define rt_glyph_white_large_square                     "\xE2\xAC\x9C"

// N-Ary Times Operator https://www.compart.com/en/unicode/U+2A09
#define rt_glyph_n_ary_times_operator                   "\xE2\xA8\x89"

// Heavy Minus Sign https://www.compart.com/en/unicode/U+2796
#define rt_glyph_heavy_minus_sign                       "\xE2\x9E\x96"

// Heavy Plus Sign https://www.compart.com/en/unicode/U+2795
#define rt_glyph_heavy_plus_sign                        "\xE2\x9E\x95"

// Clockwise Rightwards and Leftwards Open Circle Arrows with Circled One Overlay
// https://www.compart.com/en/unicode/U+1F502
// rt_glyph_clockwise_rightwards_and_leftwards_open_circle_arrows_with_circled_one_overlay...
#define rt_glyph_open_circle_arrows_one_overlay         "\xF0\x9F\x94\x82"

// Halfwidth Katakana-Hiragana Prolonged Sound Mark https://www.compart.com/en/unicode/U+FF70
#define rt_glyph_prolonged_sound_mark                   "\xEF\xBD\xB0"

// Fullwidth Plus Sign https://www.compart.com/en/unicode/U+FF0B
#define rt_glyph_fullwidth_plus_sign                    "\xEF\xBC\x8B"

// Fullwidth Hyphen-Minus https://www.compart.com/en/unicode/U+FF0D
#define rt_glyph_fullwidth_hyphen_minus                 "\xEF\xBC\x8D"


// Heavy Multiplication X https://www.compart.com/en/unicode/U+2716
#define rt_glyph_heavy_multiplication_x                 "\xE2\x9C\x96"

// Multiplication Sign https://www.compart.com/en/unicode/U+00D7
#define rt_glyph_multiplication_sign                    "\xC3\x97"

// Trigram For Heaven (caption menu button) https://www.compart.com/en/unicode/U+2630
#define rt_glyph_trigram_for_heaven                     "\xE2\x98\xB0"

// (tool bar drag handle like: msvc toolbars)
// Braille Pattern Dots-12345678  https://www.compart.com/en/unicode/U+28FF
#define rt_glyph_braille_pattern_dots_12345678          "\xE2\xA3\xBF"

// White Square Containing Black Medium Square
// https://compart.com/en/unicode/U+1F795
#define rt_glyph_white_square_containing_black_medium_square "\xF0\x9F\x9E\x95"

// White Medium Square
// https://compart.com/en/unicode/U+25FB
#define rt_glyph_white_medium_square                   "\xE2\x97\xBB"

// White Square with Upper Right Quadrant
// https://compart.com/en/unicode/U+25F3
#define rt_glyph_white_square_with_upper_right_quadrant "\xE2\x97\xB3"

// White Square with Upper Left Quadrant https://www.compart.com/en/unicode/U+25F0
#define rt_glyph_white_square_with_upper_left_quadrant "\xE2\x97\xB0"

// White Square with Lower Left Quadrant https://www.compart.com/en/unicode/U+25F1
#define rt_glyph_white_square_with_lower_left_quadrant "\xE2\x97\xB1"

// White Square with Lower Right Quadrant https://www.compart.com/en/unicode/U+25F2
#define rt_glyph_white_square_with_lower_right_quadrant "\xE2\x97\xB2"

// White Square with Upper Right Quadrant https://www.compart.com/en/unicode/U+25F3
#define rt_glyph_white_square_with_upper_right_quadrant "\xE2\x97\xB3"

// White Square with Vertical Bisecting Line https://www.compart.com/en/unicode/U+25EB
#define rt_glyph_white_square_with_vertical_bisecting_line "\xE2\x97\xAB"

// (White Square with Horizontal Bisecting Line)
// Squared Minus https://www.compart.com/en/unicode/U+229F
#define rt_glyph_squared_minus                          "\xE2\x8A\x9F"

// North East and South West Arrow https://www.compart.com/en/unicode/U+2922
#define rt_glyph_north_east_and_south_west_arrow        "\xE2\xA4\xA2"

// South East Arrow to Corner https://www.compart.com/en/unicode/U+21F2
#define rt_glyph_south_east_white_arrow_to_corner       "\xE2\x87\xB2"

// North West Arrow to Corner https://www.compart.com/en/unicode/U+21F1
#define rt_glyph_north_west_white_arrow_to_corner       "\xE2\x87\xB1"

// Leftwards Arrow to Bar https://www.compart.com/en/unicode/U+21E6
#define rt_glyph_leftwards_white_arrow_to_bar           "\xE2\x87\xA6"

// Rightwards Arrow to Bar https://www.compart.com/en/unicode/U+21E8
#define rt_glyph_rightwards_white_arrow_to_bar          "\xE2\x87\xA8"

// Upwards White Arrow https://www.compart.com/en/unicode/U+21E7
#define rt_glyph_upwards_white_arrow                    "\xE2\x87\xA7"

// Downwards White Arrow https://www.compart.com/en/unicode/U+21E9
#define rt_glyph_downwards_white_arrow                  "\xE2\x87\xA9"

// Leftwards White Arrow https://www.compart.com/en/unicode/U+21E4
#define rt_glyph_leftwards_white_arrow                  "\xE2\x87\xA4"

// Rightwards White Arrow https://www.compart.com/en/unicode/U+21E5
#define rt_glyph_rightwards_white_arrow                 "\xE2\x87\xA5"

// Upwards White Arrow on Pedestal https://www.compart.com/en/unicode/U+21EB
#define rt_glyph_upwards_white_arrow_on_pedestal        "\xE2\x87\xAB"

// Braille Pattern Dots-678 https://www.compart.com/en/unicode/U+28E0
#define rt_glyph_3_dots_tiny_right_bottom_triangle      "\xE2\xA3\xA0"

// Braille Pattern Dots-2345678 https://www.compart.com/en/unicode/U+28FE
// Combining the two into:
#define rt_glyph_dotted_right_bottom_triangle           "\xE2\xA3\xA0\xE2\xA3\xBE"

// Upper Right Drop-Shadowed White Square https://www.compart.com/en/unicode/U+2750
#define rt_glyph_upper_right_drop_shadowed_white_square "\xE2\x9D\x90"

// No-Break Space (NBSP)
// https://www.compart.com/en/unicode/U+00A0
#define rt_glyph_nbsp                                   "\xC2\xA0"

// Word Joiner (WJ)
// https://compart.com/en/unicode/U+2060
#define rt_glyph_word_joiner                            "\xE2\x81\xA0"

// Zero Width Space (ZWSP)
// https://compart.com/en/unicode/U+200B
#define rt_glyph_zwsp                                   "\xE2\x80\x8B"

// Zero Width Joiner (ZWJ)
// https://compart.com/en/unicode/u+200D
#define rt_glyph_zwj                                    "\xE2\x80\x8D"

// En Quad
// https://compart.com/en/unicode/U+2000
#define rt_glyph_en_quad "\xE2\x80\x80"

// Em Quad
// https://compart.com/en/unicode/U+2001
#define rt_glyph_em_quad "\xE2\x80\x81"

// En Space
// https://compart.com/en/unicode/U+2002
#define rt_glyph_en_space "\xE2\x80\x82"

// Em Space
// https://compart.com/en/unicode/U+2003
#define rt_glyph_em_space "\xE2\x80\x83"

// Hyphen https://www.compart.com/en/unicode/U+2010
#define rt_glyph_hyphen                                "\xE2\x80\x90"

// Non-Breaking Hyphen https://www.compart.com/en/unicode/U+2011
#define rt_glyph_non_breaking_hyphen                   "\xE2\x80\x91"

// Fullwidth Low Line https://www.compart.com/en/unicode/U+FF3F
#define rt_glyph_fullwidth_low_line                    "\xEF\xBC\xBF"

// #define rt_glyph_light_horizontal                     "\xE2\x94\x80"
// Light Horizontal https://www.compart.com/en/unicode/U+2500
#define rt_glyph_light_horizontal                     "\xE2\x94\x80"

// Three-Em Dash https://www.compart.com/en/unicode/U+2E3B
#define rt_glyph_three_em_dash                         "\xE2\xB8\xBB"

// Infinity https://www.compart.com/en/unicode/U+221E
#define rt_glyph_infinity                              "\xE2\x88\x9E"

// Black Large Circle https://www.compart.com/en/unicode/U+2B24
#define rt_glyph_black_large_circle                    "\xE2\xAC\xA4"

// Large Circle https://www.compart.com/en/unicode/U+25EF
#define rt_glyph_large_circle                          "\xE2\x97\xAF"

// Heavy Leftwards Arrow with Equilateral Arrowhead https://www.compart.com/en/unicode/U+1F818
#define rt_glyph_heavy_leftwards_arrow_with_equilateral_arrowhead           "\xF0\x9F\xA0\x98"

// Heavy Rightwards Arrow with Equilateral Arrowhead https://www.compart.com/en/unicode/U+1F81A
#define rt_glyph_heavy_rightwards_arrow_with_equilateral_arrowhead          "\xF0\x9F\xA0\x9A"

// Heavy Leftwards Arrow with Large Equilateral Arrowhead https://www.compart.com/en/unicode/U+1F81C
#define rt_glyph_heavy_leftwards_arrow_with_large_equilateral_arrowhead     "\xF0\x9F\xA0\x9C"

// Heavy Rightwards Arrow with Large Equilateral Arrowhead https://www.compart.com/en/unicode/U+1F81E
#define rt_glyph_heavy_rightwards_arrow_with_large_equilateral_arrowhead    "\xF0\x9F\xA0\x9E"

// CJK Unified Ideograph-5973: Kanji Onna "Female" https://www.compart.com/en/unicode/U+5973
#define rt_glyph_kanji_onna_female                                          "\xE2\xBC\xA5"

// Leftwards Arrow https://www.compart.com/en/unicode/U+2190
#define rt_glyph_leftward_arrow                                             "\xE2\x86\x90"

// Upwards Arrow https://www.compart.com/en/unicode/U+2191
#define rt_glyph_upwards_arrow                                              "\xE2\x86\x91"

// Rightwards Arrow
// https://www.compart.com/en/unicode/U+2192
#define rt_glyph_rightwards_arrow                                           "\xE2\x86\x92"

// Downwards Arrow https://www.compart.com/en/unicode/U+2193
#define rt_glyph_downwards_arrow                                            "\xE2\x86\x93"

// Thin Space https://www.compart.com/en/unicode/U+2009
#define rt_glyph_thin_space                                                 "\xE2\x80\x89"

// Medium Mathematical Space (MMSP) https://www.compart.com/en/unicode/U+205F
#define rt_glyph_mmsp                                                       "\xE2\x81\x9F"

// Three-Per-Em Space https://www.compart.com/en/unicode/U+2004
#define rt_glyph_three_per_em                                               "\xE2\x80\x84"

// Six-Per-Em Space https://www.compart.com/en/unicode/U+2006
#define rt_glyph_six_per_em                                                 "\xE2\x80\x86"

// Punctuation Space https://www.compart.com/en/unicode/U+2008
#define rt_glyph_punctuation                                                "\xE2\x80\x88"

// Hair Space https://www.compart.com/en/unicode/U+200A
#define rt_glyph_hair_space                                                 "\xE2\x80\x8A"

// Chinese "jin4" https://www.compart.com/en/unicode/U+58F9
#define rt_glyph_chinese_jin4                                               "\xE5\xA3\xB9"

// Chinese "gong" https://www.compart.com/en/unicode/U+8D70
#define rt_glyph_chinese_gong                                                "\xE8\xB5\xB0"

// https://www.compart.com/en/unicode/U+1F9F8
#define rt_glyph_teddy_bear                                                 "\xF0\x9F\xA7\xB8"

// https://www.compart.com/en/unicode/U+1F9CA
#define rt_glyph_ice_cube                                                   "\xF0\x9F\xA7\x8A"

// Speaker https://www.compart.com/en/unicode/U+1F508
#define rt_glyph_speaker                                                    "\xF0\x9F\x94\x88"

// Speaker with Cancellation Stroke https://www.compart.com/en/unicode/U+1F507
#define rt_glyph_mute                                                       "\xF0\x9F\x94\x87"

// TODO: this is used for Font Metric Visualization

// Full Block https://www.compart.com/en/unicode/U+2588
#define rt_glyph_full_block                             "\xE2\x96\x88"

// Black Square https://www.compart.com/en/unicode/U+25A0
#define rt_glyph_black_square                           "\xE2\x96\xA0"

// the appearance of a dragon walking
// CJK Unified Ideograph-9F98 https://www.compart.com/en/unicode/U+9F98
#define rt_glyph_walking_dragon                         "\xE9\xBE\x98"

// possibly highest "diacritical marks" character (Vietnamese)
// Latin Small Letter U with Horn and Hook Above https://www.compart.com/en/unicode/U+1EED
#define rt_glyph_u_with_horn_and_hook_above             "\xC7\xAD"

// possibly "long descender" character
// Latin Small Letter Qp Digraph https://www.compart.com/en/unicode/U+0239
#define rt_glyph_qp_digraph                             "\xC9\xB9"

// another possibly "long descender" character
// Cyrillic Small Letter Shha with Descender https://www.compart.com/en/unicode/U+0527
#define rt_glyph_shha_with_descender                    "\xD4\xA7"

// a"very long descender" character
// Tibetan Mark Caret Yig Mgo Phur Shad Ma https://www.compart.com/en/unicode/U+0F06
#define rt_glyph_caret_yig_mgo_phur_shad_ma             "\xE0\xBC\x86"

// Tibetan Vowel Sign Vocalic Ll https://www.compart.com/en/unicode/U+0F79
#define rt_glyph_vocalic_ll                             "\xE0\xBD\xB9"

// https://www.compart.com/en/unicode/U+1F4A3
#define rt_glyph_bomb "\xF0\x9F\x92\xA3"

// https://www.compart.com/en/unicode/U+1F4A1
#define rt_glyph_electric_light_bulb "\xF0\x9F\x92\xA1"

// https://www.compart.com/en/unicode/U+1F4E2
#define rt_glyph_public_address_loudspeaker "\xF0\x9F\x93\xA2"

// https://www.compart.com/en/unicode/U+1F517
#define rt_glyph_link_symbol "\xF0\x9F\x94\x97"

// https://www.compart.com/en/unicode/U+1F571
#define rt_glyph_black_skull_and_crossbones "\xF0\x9F\x95\xB1"

// https://www.compart.com/en/unicode/U+1F5B5
#define rt_glyph_screen "\xF0\x9F\x96\xB5"

// https://www.compart.com/en/unicode/U+1F5D7
#define rt_glyph_overlap "\xF0\x9F\x97\x97"

// https://www.compart.com/en/unicode/U+1F5D6
#define rt_glyph_maximize "\xF0\x9F\x97\x96"

// https://www.compart.com/en/unicode/U+1F5D5
#define rt_glyph_minimize "\xF0\x9F\x97\x95"

// Desktop Window
// https://compart.com/en/unicode/U+1F5D4
#define rt_glyph_desktop_window "\xF0\x9F\x97\x94"

// https://www.compart.com/en/unicode/U+1F5D9
#define rt_glyph_cancellation_x "\xF0\x9F\x97\x99"

// https://www.compart.com/en/unicode/U+1F5DF
#define rt_glyph_page_with_circled_text "\xF0\x9F\x97\x9F"

// https://www.compart.com/en/unicode/U+1F533
#define rt_glyph_white_square_button "\xF0\x9F\x94\xB3"

// https://www.compart.com/en/unicode/U+1F532
#define rt_glyph_black_square_button "\xF0\x9F\x94\xB2"

// https://www.compart.com/en/unicode/U+1F5F9
#define rt_glyph_ballot_box_with_bold_check "\xF0\x9F\x97\xB9"

// https://www.compart.com/en/unicode/U+1F5F8
#define rt_glyph_light_check_mark "\xF0\x9F\x97\xB8"

// https://compart.com/en/unicode/U+1F4BB
#define rt_glyph_personal_computer "\xF0\x9F\x92\xBB"

// https://compart.com/en/unicode/U+1F4DC
#define rt_glyph_desktop_computer "\xF0\x9F\x93\x9C"

// https://compart.com/en/unicode/U+1F4DD
#define rt_glyph_printer "\xF0\x9F\x93\x9D"

// https://compart.com/en/unicode/U+1F4F9
#define rt_glyph_video_camera "\xF0\x9F\x93\xB9"

// https://compart.com/en/unicode/U+1F4F8
#define rt_glyph_camera "\xF0\x9F\x93\xB8"

// https://compart.com/en/unicode/U+1F505
#define rt_glyph_high_brightness "\xF0\x9F\x94\x85"

// https://compart.com/en/unicode/U+1F506
#define rt_glyph_low_brightness "\xF0\x9F\x94\x86"

// https://compart.com/en/unicode/U+1F507
#define rt_glyph_speaker_with_cancellation_stroke "\xF0\x9F\x94\x87"

// https://compart.com/en/unicode/U+1F509
#define rt_glyph_speaker_with_one_sound_wave "\xF0\x9F\x94\x89"

// Right-Pointing Magnifying Glass
// https://compart.com/en/unicode/U+1F50E
#define rt_glyph_right_pointing_magnifying_glass "\xF0\x9F\x94\x8E"

// Radio Button
// https://compart.com/en/unicode/U+1F518
#define rt_glyph_radio_button "\xF0\x9F\x94\x98"

// https://compart.com/en/unicode/U+1F525
#define rt_glyph_fire "\xF0\x9F\x94\xA5"

// Gear
// https://compart.com/en/unicode/U+2699
#define rt_glyph_gear "\xE2\x9A\x99"

// Nut and Bolt
// https://compart.com/en/unicode/U+1F529
#define rt_glyph_nut_and_bolt "\xF0\x9F\x94\xA9"

// Hammer and Wrench
// https://compart.com/en/unicode/U+1F6E0
#define rt_glyph_hammer_and_wrench "\xF0\x9F\x9B\xA0"

// https://compart.com/en/unicode/U+1F53E
#define rt_glyph_upwards_button "\xF0\x9F\x94\xBE"

// https://compart.com/en/unicode/U+1F53F
#define rt_glyph_downwards_button "\xF0\x9F\x94\xBF"

// https://compart.com/en/unicode/U+1F5C7
#define rt_glyph_litter_in_bin_sign "\xF0\x9F\x97\x87"

// Checker Board
// https://compart.com/en/unicode/U+1F67E
#define rt_glyph_checker_board "\xF0\x9F\x9A\xBE"

// Reverse Checker Board
// https://compart.com/en/unicode/U+1F67F
#define rt_glyph_reverse_checker_board "\xF0\x9F\x9A\xBF"

// Clipboard
// https://compart.com/en/unicode/U+1F4CB
#define rt_glyph_clipboard "\xF0\x9F\x93\x8B"

// Two Joined Squares https://www.compart.com/en/unicode/U+29C9
#define rt_glyph_two_joined_squares "\xE2\xA7\x89"

// White Heavy Check Mark
// https://compart.com/en/unicode/U+2705
#define rt_glyph_white_heavy_check_mark "\xE2\x9C\x85"

// Negative Squared Cross Mark
// https://compart.com/en/unicode/U+274E
#define rt_glyph_negative_squared_cross_mark "\xE2\x9D\x8E"

// Lower Right Drop-Shadowed White Square
// https://compart.com/en/unicode/U+274F
#define rt_glyph_lower_right_drop_shadowed_white_square "\xE2\x9D\x8F"

// Upper Right Drop-Shadowed White Square
// https://compart.com/en/unicode/U+2750
#define rt_glyph_upper_right_drop_shadowed_white_square "\xE2\x9D\x90"

// Lower Right Shadowed White Square
// https://compart.com/en/unicode/U+2751
#define rt_glyph_lower_right_shadowed_white_square "\xE2\x9D\x91"

// Upper Right Shadowed White Square
// https://compart.com/en/unicode/U+2752
#define rt_glyph_upper_right_shadowed_white_square "\xE2\x9D\x92"

// Left Double Wiggly Fence
// https://compart.com/en/unicode/U+29DA
#define rt_glyph_left_double_wiggly_fence "\xE2\xA7\x9A"

// Right Double Wiggly Fence
// https://compart.com/en/unicode/U+29DB
#define rt_glyph_right_double_wiggly_fence "\xE2\xA7\x9B"

// Logical Or
// https://compart.com/en/unicode/U+2228
#define rt_glyph_logical_or "\xE2\x88\xA8"

// Logical And
// https://compart.com/en/unicode/U+2227
#define rt_glyph_logical_and "\xE2\x88\xA7"

// Double Vertical Bar (Pause)
// https://compart.com/en/unicode/U+23F8
#define rt_glyph_double_vertical_bar "\xE2\x8F\xB8"

// Black Square For Stop
// https://compart.com/en/unicode/U+23F9
#define rt_glyph_black_square_for_stop "\xE2\x8F\xB9"

// Black Circle For Record
// https://compart.com/en/unicode/U+23FA
#define rt_glyph_black_circle_for_record "\xE2\x8F\xBA"

// Negative Squared Latin Capital Letter "I"
// https://compart.com/en/unicode/U+1F158
#define rt_glyph_negative_squared_latin_capital_letter_i "\xF0\x9F\x85\x98"
#define rt_glyph_info rt_glyph_negative_squared_latin_capital_letter_i

// Circled Information Source
// https://compart.com/en/unicode/U+1F6C8
#define rt_glyph_circled_information_source "\xF0\x9F\x9B\x88"

// Information Source
// https://compart.com/en/unicode/U+2139
#define rt_glyph_information_source "\xE2\x84\xB9"

// Squared Cool
// https://compart.com/en/unicode/U+1F192
#define rt_glyph_squared_cool "\xF0\x9F\x86\x92"

// Squared OK
// https://compart.com/en/unicode/U+1F197
#define rt_glyph_squared_ok "\xF0\x9F\x86\x97"

// Squared Free
// https://compart.com/en/unicode/U+1F193
#define rt_glyph_squared_free "\xF0\x9F\x86\x93"

// Squared New
// https://compart.com/en/unicode/U+1F195
#define rt_glyph_squared_new "\xF0\x9F\x86\x95"

// Lady Beetle
// https://compart.com/en/unicode/U+1F41E
#define rt_glyph_lady_beetle "\xF0\x9F\x90\x9E"

// Brain
// https://compart.com/en/unicode/U+1F9E0
#define rt_glyph_brain "\xF0\x9F\xA7\xA0"

// South West Arrow with Hook
// https://www.compart.com/en/unicode/U+2926
#define rt_glyph_south_west_arrow_with_hook "\xE2\xA4\xA6"

// North West Arrow with Hook
// https://www.compart.com/en/unicode/U+2923
#define rt_glyph_north_west_arrow_with_hook "\xE2\xA4\xA3"

// White Sun with Rays
// https://www.compart.com/en/unicode/U+263C
#define rt_glyph_white_sun_with_rays "\xE2\x98\xBC"

// Black Sun with Rays
// https://www.compart.com/en/unicode/U+2600
#define rt_glyph_black_sun_with_rays "\xE2\x98\x80"

// Sun Behind Cloud
// https://www.compart.com/en/unicode/U+26C5
#define rt_glyph_sun_behind_cloud "\xE2\x9B\x85"

// White Sun
// https://www.compart.com/en/unicode/U+1F323
#define rt_glyph_white_sun "\xF0\x9F\x8C\xA3"

// Crescent Moon
// https://www.compart.com/en/unicode/U+1F319
#define rt_glyph_crescent_moon "\xF0\x9F\x8C\x99"

// Latin Capital Letter E with Cedilla and Breve
// https://compart.com/en/unicode/U+1E1C
#define rt_glyph_E_with_cedilla_and_breve "\xE1\xB8\x9C"

// Box Drawings Heavy Vertical and Horizontal
// https://compart.com/en/unicode/U+254B
#define rt_glyph_box_drawings_heavy_vertical_and_horizontal "\xE2\x95\x8B"

// Box Drawings Light Diagonal Cross
// https://compart.com/en/unicode/U+2573
#define rt_glyph_box_drawings_light_diagonal_cross "\xE2\x95\xB3"

// Combining Enclosing Square
// https://compart.com/en/unicode/U+20DE
#define rt_glyph_combining_enclosing_square "\xE2\x83\x9E"

// Combining Enclosing Screen
// https://compart.com/en/unicode/U+20E2
#define rt_glyph_combining_enclosing_screen "\xE2\x83\xA2"

// Combining Enclosing Keycap
// https://compart.com/en/unicode/U+20E3
#define rt_glyph_combining_enclosing_keycap "\xE2\x83\xA3"

// Combining Enclosing Circle
// https://compart.com/en/unicode/U+20DD
#define rt_glyph_combining_enclosing_circle "\xE2\x83\x9D"

// Frame with Picture
// https://compart.com/en/unicode/U+1F5BC
#define rt_glyph_frame_with_picture "\xF0\x9F\x96\xBC"
// with emoji variation selector: "\xF0\x9F\x96\xBC\xEF\xB8\x8F"

// Document with Picture
// https://compart.com/en/unicode/U+1F5BB
#define rt_glyph_document_with_picture "\xF0\x9F\x96\xBB"

// Frame with Tiles
// https://compart.com/en/unicode/U+1F5BD
#define rt_glyph_frame_with_tiles "\xF0\x9F\x96\xBD"

// Frame with an X
// https://compart.com/en/unicode/U+1F5BE
#define rt_glyph_frame_with_an_x "\xF0\x9F\x96\xBE"

// Left Right Arrow
// https://compart.com/en/unicode/U+2194
#define rt_glyph_left_right_arrow "\xE2\x86\x94"

// Up Down Arrow
// https://compart.com/en/unicode/U+2195
#define rt_glyph_up_down_arrow "\xE2\x86\x95"

// ________________________________ rt_heap.h _________________________________

rt_begin_c

// It is absolutely OK to use posix compliant
// malloc()/calloc()/realloc()/free() function calls with understanding
// that they introduce serialization points in multi-threaded applications
// and may be induce wait states that under pressure (all cores busy) may
// result in prolonged wait which may not be acceptable for real time
// processing pipelines.
//
// heap_if.functions may or may not be faster than malloc()/free() ...
//
// Some callers may find realloc() parameters more convenient to avoid
// anti-pattern
//      void* reallocated = realloc(p, new_size);
//      if (reallocated != null) { p = reallocated; }
// and avoid never ending discussion of legality and implementation
// compliance of the situation:
//      realloc(p /* when p == null */, ...)
//
// zero: true initializes allocated or reallocated tail memory to 0x00
// be careful with zeroing heap memory. It will result in virtual
// to physical memory mapping and may be expensive.

typedef struct rt_heap_s rt_heap_t;

typedef struct { // heap == null uses process serialized LFH
    errno_t (*alloc)(void* *a, int64_t bytes);
    errno_t (*alloc_zero)(void* *a, int64_t bytes);
    errno_t (*realloc)(void* *a, int64_t bytes);
    errno_t (*realloc_zero)(void* *a, int64_t bytes);
    void    (*free)(void* a);
    // heaps:
    rt_heap_t* (*create)(bool serialized);
    errno_t (*allocate)(rt_heap_t* heap, void* *a, int64_t bytes, bool zero);
    // reallocate may return ERROR_OUTOFMEMORY w/o changing 'a' *)
    errno_t (*reallocate)(rt_heap_t* heap, void* *a, int64_t bytes, bool zero);
    void    (*deallocate)(rt_heap_t* heap, void* a);
    int64_t (*bytes)(rt_heap_t* heap, void* a); // actual allocated size
    void    (*dispose)(rt_heap_t* heap);
    void    (*test)(void);
} rt_heap_if;

extern rt_heap_if rt_heap;

// *) zero in reallocate applies to the newly appended bytes

// On Windows rt_mem.heap is based on serialized LFH returned by GetProcessHeap()
// https://learn.microsoft.com/en-us/windows/win32/memory/low-fragmentation-heap
// threads can benefit from not serialized, not LFH if they allocate and free
// memory in time critical loops.

rt_end_c


// _______________________________ rt_loader.h ________________________________

rt_begin_c

// see:
// https://pubs.opengroup.org/onlinepubs/7908799/xsh/dlfcn.h.html

typedef struct {
    // mode:
    int32_t const local;
    int32_t const lazy;
    int32_t const now;
    int32_t const global;
    // "If the value of file is null, dlopen() provides a handle on a global
    //  symbol object." posix
    void* (*open)(const char* filename, int32_t mode);
    void* (*sym)(void* handle, const char* name);
    void  (*close)(void* handle);
    void (*test)(void);
} rt_loader_if;

extern rt_loader_if rt_loader;

rt_end_c

// _________________________________ rt_mem.h _________________________________

rt_begin_c

typedef struct {
    // whole file read only
    errno_t (*map_ro)(const char* filename, void** data, int64_t* bytes);
    // whole file read-write
    errno_t (*map_rw)(const char* filename, void** data, int64_t* bytes);
    void (*unmap)(void* data, int64_t bytes);
    // map_resource() maps data from resources, do NOT unmap!
    errno_t  (*map_resource)(const char* label, void** data, int64_t* bytes);
    int32_t (*page_size)(void); // 4KB or 64KB on Windows
    int32_t (*large_page_size)(void);  // 2MB on Windows
    // allocate() contiguous reserved virtual address range,
    // if possible committed to physical memory.
    // Memory guaranteed to be aligned to page boundary.
    // Memory is guaranteed to be initialized to zero on access.
    void* (*allocate)(int64_t bytes_multiple_of_page_size);
    void  (*deallocate)(void* a, int64_t bytes_multiple_of_page_size);
    void  (*test)(void);
} rt_mem_if;

extern rt_mem_if rt_mem;

rt_end_c


// _________________________________ rt_nls.h _________________________________

rt_begin_c

typedef struct { // i18n national language support
    void (*init)(void);
    const char* (*locale)(void);  // "en-US" "zh-CN" etc...
    // force locale for debugging and testing:
    errno_t (*set_locale)(const char* locale); // only for calling thread
    // nls(s) is same as string(strid(s), s)
    const char* (*str)(const char* defau1t); // returns localized string
    // strid("foo") returns -1 if there is no matching
    // ENGLISH NEUTRAL STRINGTABLE entry
    int32_t (*strid)(const char* s);
    // given strid > 0 returns localized string or default value
    const char* (*string)(int32_t strid, const char* defau1t);
} rt_nls_if;

extern rt_nls_if rt_nls;

rt_end_c

// _________________________________ rt_num.h _________________________________

rt_begin_c

typedef struct {
    uint64_t lo;
    uint64_t hi;
} rt_num128_t; // uint128_t may be supported by compiler

typedef struct {
    rt_num128_t (*add128)(const rt_num128_t a, const rt_num128_t b);
    rt_num128_t (*sub128)(const rt_num128_t a, const rt_num128_t b);
    rt_num128_t (*mul64x64)(uint64_t a, uint64_t b);
    uint64_t (*muldiv128)(uint64_t a, uint64_t b, uint64_t d);
    uint32_t (*gcd32)(uint32_t u, uint32_t v); // greatest common denominator
    // non-crypto strong pseudo-random number generators (thread safe)
    uint32_t (*random32)(uint32_t *state); // "Mulberry32"
    uint64_t (*random64)(uint64_t *state); // "Trust"
    // "FNV-1a" hash functions (if bytes == 0 expects zero terminated string)
    uint32_t (*hash32)(const char* s, int64_t bytes);
    uint64_t (*hash64)(const char* s, int64_t bytes);
    void     (*test)(void);
} rt_num_if;

extern rt_num_if rt_num;

rt_end_c


// _______________________________ rt_static.h ________________________________

rt_begin_c

// rt_static_init(unique_name) { code_to_execute_before_main }

#if defined(_MSC_VER)

#if defined(_WIN64) || defined(_M_X64)
#define _msvc_symbol_prefix_ ""
#else
#define _msvc_symbol_prefix_ "_"
#endif

#pragma comment(linker, "/include:rt_force_symbol_reference")

void* rt_force_symbol_reference(void* symbol);

#define _msvc_ctor_(sym_prefix, func)                                    \
  void func(void);                                                        \
  int32_t (* rt_array ## func)(void);                                     \
  int32_t func ## _wrapper(void);                                         \
  int32_t func ## _wrapper(void) { func();                                \
  rt_force_symbol_reference((void*)rt_array ## func);                     \
  rt_force_symbol_reference((void*)func ## _wrapper); return 0; }         \
  extern int32_t (* rt_array ## func)(void);                              \
  __pragma(comment(linker, "/include:" sym_prefix # func "_wrapper"))     \
  __pragma(section(".CRT$XCU", read))                                     \
  __declspec(allocate(".CRT$XCU"))                                        \
    int32_t (* rt_array ## func)(void) = func ## _wrapper;

#define rt_static_init2_(func, line) _msvc_ctor_(_msvc_symbol_prefix_, \
    func ## _constructor_##line)                                       \
    void func ## _constructor_##line(void)

#define rt_static_init1_(func, line) rt_static_init2_(func, line)

#define rt_static_init(func) rt_static_init1_(func, __LINE__)

#else
#define rt_static_init(n) __attribute__((constructor)) \
        static void _init_ ## n ## __LINE__ ## _ctor(void)
#endif

void rt_static_init_test(void);

rt_end_c

// _______________________________ rt_streams.h _______________________________

rt_begin_c

typedef struct rt_stream_if rt_stream_if;

typedef struct rt_stream_if {
    errno_t (*read)(rt_stream_if* s, void* data, int64_t bytes,
                    int64_t *transferred);
    errno_t (*write)(rt_stream_if* s, const void* data, int64_t bytes,
                     int64_t *transferred);
    void    (*close)(rt_stream_if* s); // optional
} rt_stream_if;

typedef struct {
    rt_stream_if   stream;
    const void* data_read;
    int64_t     bytes_read;
    int64_t     pos_read;
    void*       data_write;
    int64_t     bytes_write;
    int64_t     pos_write;
} rt_stream_memory_if;

typedef struct {
    void (*read_only)(rt_stream_memory_if* s,  const void* data, int64_t bytes);
    void (*write_only)(rt_stream_memory_if* s, void* data, int64_t bytes);
    void (*read_write)(rt_stream_memory_if* s, const void* read, int64_t read_bytes,
                                               void* write, int64_t write_bytes);
    void (*test)(void);
} rt_streams_if;

extern rt_streams_if rt_streams;

rt_end_c

// ______________________________ rt_processes.h ______________________________

rt_begin_c

typedef struct {
    const char* command;
    rt_stream_if* in;
    rt_stream_if* out;
    rt_stream_if* err;
    uint32_t exit_code;
    fp64_t   timeout; // seconds
} rt_processes_child_t;

// Process name may be an the executable filename with
// full, partial or absent pathname.
// Case insensitive on Windows.

typedef struct {
    const char* (*name)(void); // argv[0] like but full path
    uint64_t  (*pid)(const char* name); // 0 if process not found
    errno_t   (*pids)(const char* name, uint64_t* pids/*[size]*/, int32_t size,
                      int32_t *count); // return 0, error or ERROR_MORE_DATA
    errno_t   (*nameof)(uint64_t pid, char* name, int32_t count); // pathname
    bool      (*present)(uint64_t pid);
    errno_t   (*kill)(uint64_t pid, fp64_t timeout_seconds);
    errno_t   (*kill_all)(const char* name, fp64_t timeout_seconds);
    bool      (*is_elevated)(void); // Is process running as root/ Admin / System?
    errno_t   (*restart_elevated)(void); // retuns error or exits on success
    errno_t   (*run)(rt_processes_child_t* child);
    errno_t   (*popen)(const char* command, int32_t *xc, rt_stream_if* output,
                       fp64_t timeout_seconds); // <= 0 infinite
    // popen() does NOT guarantee stream zero termination on errors
    errno_t  (*spawn)(const char* command); // spawn fully detached process
    void (*test)(void);
} rt_processes_if;

extern rt_processes_if rt_processes;

rt_end_c

// _______________________________ rt_threads.h _______________________________

rt_begin_c

typedef struct rt_event_s* rt_event_t;

typedef struct {
    rt_event_t (*create)(void); // never returns null
    rt_event_t (*create_manual)(void); // never returns null
    void (*set)(rt_event_t e);
    void (*reset)(rt_event_t e);
    void (*wait)(rt_event_t e);
    // returns 0 on success or -1 on timeout
    int32_t (*wait_or_timeout)(rt_event_t e, fp64_t seconds); // seconds < 0 forever
    // returns event index or -1 on timeout or -2 on abandon
    int32_t (*wait_any)(int32_t n, rt_event_t events[]); // -1 on abandon
    int32_t (*wait_any_or_timeout)(int32_t n, rt_event_t e[], fp64_t seconds);
    void (*dispose)(rt_event_t e);
    void (*test)(void);
} rt_event_if;

extern rt_event_if rt_event;

typedef struct rt_aligned_8 rt_mutex_s { uint8_t content[40]; } rt_mutex_t;

typedef struct {
    void (*init)(rt_mutex_t* m);
    void (*lock)(rt_mutex_t* m);
    void (*unlock)(rt_mutex_t* m);
    void (*dispose)(rt_mutex_t* m);
    void (*test)(void);
} rt_mutex_if;

extern rt_mutex_if rt_mutex;

typedef struct thread_s* rt_thread_t;

typedef struct {
    rt_thread_t (*start)(void (*func)(void*), void* p); // never returns null
    errno_t     (*join)(rt_thread_t thread, fp64_t timeout_seconds); // < 0 forever
    void        (*detach)(rt_thread_t thread); // closes handle. thread is not joinable
    void        (*name)(const char* name); // names the thread
    void        (*realtime)(void); // bumps calling thread priority
    void        (*yield)(void);    // pthread_yield() / Win32: SwitchToThread()
    void        (*sleep_for)(fp64_t seconds);
    uint64_t    (*id_of)(rt_thread_t t);
    uint64_t    (*id)(void); // gettid()
    rt_thread_t (*self)(void); // Pseudo Handle may differ in access to .open(.id())
    errno_t     (*open)(rt_thread_t* t, uint64_t id);
    void        (*close)(rt_thread_t t);
    void        (*test)(void);
} rt_thread_if;

extern rt_thread_if rt_thread;

rt_end_c

// ________________________________ rt_vigil.h ________________________________

rt_begin_c

// better rt_assert() - augmented with printf format and parameters
// rt_swear() - release configuration rt_assert() in honor of:
// https://github.com/munificent/vigil

#define rt_static_assertion(condition) static_assert(condition, #condition)

typedef struct {
    int32_t (*failed_assertion)(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...);
    int32_t (*fatal_termination)(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...);
    int32_t (*fatal_if_error)(const char* file, int32_t line, const char* func,
        const char* condition, errno_t r, const char* format, ...);
    void (*test)(void);
} rt_vigil_if;

extern rt_vigil_if rt_vigil;

#if defined(DEBUG)
  #define rt_assert(b, ...) rt_suppress_constant_cond_exp           \
    /* const cond */                                                \
    (void)((!!(b)) || rt_vigil.failed_assertion(__FILE__, __LINE__, \
    __func__, #b, "" __VA_ARGS__))
#else
  #define rt_assert(b, ...) ((void)0)
#endif

// rt_swear() is runtime rt_assert() for both debug and release configurations

#define rt_swear(b, ...) rt_suppress_constant_cond_exp                 \
    /* const cond */                                                \
    (void)((!!(b)) || rt_vigil.failed_assertion(__FILE__, __LINE__, \
    __func__, #b, "" __VA_ARGS__))

#define rt_fatal(...) (void)(rt_vigil.fatal_termination(            \
    __FILE__, __LINE__,  __func__, "",  "" __VA_ARGS__))

#define rt_fatal_if(b, ...) rt_suppress_constant_cond_exp           \
    /* const cond */                                                \
    (void)((!(b)) || rt_vigil.fatal_termination(__FILE__, __LINE__, \
    __func__, #b, "" __VA_ARGS__))

#define rt_fatal_if_not(b, ...) rt_suppress_constant_cond_exp        \
    /* const cond */                                                 \
    (void)((!!(b)) || rt_vigil.fatal_termination(__FILE__, __LINE__, \
    __func__, #b, "" __VA_ARGS__))

#define rt_not_null(e, ...) rt_fatal_if((e) == null, "" __VA_ARGS__)

#define rt_fatal_if_error(r, ...) rt_suppress_constant_cond_exp      \
    /* const cond */                                                 \
    (void)(rt_vigil.fatal_if_error(__FILE__, __LINE__, __func__,     \
                                   #r, r, "" __VA_ARGS__))

#define rt_fatal_win32err(c, ...) rt_suppress_constant_cond_exp      \
    /* const cond */                                                 \
    (void)(rt_vigil.fatal_if_error(__FILE__, __LINE__, __func__,     \
                                   #c, rt_b2e(c), "" __VA_ARGS__))

rt_end_c
// ___________________________________ rt.h ___________________________________

// the rest is in alphabetical order (no inter dependencies)

// ________________________________ rt_work.h _________________________________

rt_begin_c

// Minimalistic "react"-like work_queue or work items and
// a thread based workers. See rt_worker_test() for usage.

typedef struct rt_event_s*     rt_event_t;
typedef struct rt_work_s       rt_work_t;
typedef struct rt_work_queue_s rt_work_queue_t;

typedef struct rt_work_s {
    rt_work_queue_t* queue; // queue where the call is or was last scheduled
    fp64_t when;       // proc() call will be made after or at this time
    void (*work)(rt_work_t* c);
    void*  data;       // extra data that will be passed to proc() call
    rt_event_t  done;  // if not null signalled after calling proc() or canceling
    rt_work_t*  next;  // next element in the queue (implementation detail)
    bool    canceled;  // set to true inside .cancel() call
} rt_work_t;

typedef struct rt_work_queue_s {
    rt_work_t* head;
    int64_t    lock; // spinlock
    rt_event_t changed; // if not null will be signaled when head changes
} rt_work_queue_t;

typedef struct rt_work_queue_if {
    void (*post)(rt_work_t* c);
    bool (*get)(rt_work_queue_t*, rt_work_t* *c);
    void (*call)(rt_work_t* c);
    void (*dispatch)(rt_work_queue_t* q); // all ready messages
    void (*cancel)(rt_work_t* c);
    void (*flush)(rt_work_queue_t* q); // cancel all requests in the queue
} rt_work_queue_if;

extern rt_work_queue_if  rt_work_queue;

typedef struct rt_worker_s {
    rt_work_queue_t queue;
    rt_thread_t     thread;
    rt_event_t      wake;
    volatile bool   quit;
} rt_worker_t;

typedef struct rt_worker_if {
    void    (*start)(rt_worker_t* tq);
    void    (*post)(rt_worker_t* tq, rt_work_t* w);
    errno_t (*join)(rt_worker_t* tq, fp64_t timeout);
    void    (*test)(void);
} rt_worker_if;

extern rt_worker_if rt_worker;

// worker thread waits for a queue's `wake` event with the timeout
// infinity or if queue is not empty delta time till the head
// item of the queue.
//
// Upon post() call the `wake` event is set and the worker thread
// wakes up and dispatches all the items with .when less then now
// calling function work() if it is not null and optionally signaling
// .done event if it is not null.
//
// When all ready items in the queue are processed worker thread locks
// the queue and if the head is present calculates next timeout based
// on .when time of the head or sets timeout to infinity if the queue
// is empty.
//
// Function .join() sets .quit to true signals .wake event and attempt
// to join the worker .thread with specified timeout.
// It is the responsibility of the caller to ensure that no other
// work is posted after calling .join() because it will be lost.

rt_end_c

/*
    Usage examples:

    // The dispatch_until() is just for testing purposes.
    // Usually rt_work_queue.dispatch(q) will be called inside each
    // iteration of message loop of a dispatch [UI] thread.

    static void dispatch_until(rt_work_queue_t* q, int32_t* i, const int32_t n);

    // simple way of passing a single pointer to call_later

    static void every_100ms(rt_work_t* w) {
        int32_t* i = (int32_t*)w->data;
        rt_println("i: %d", *i);
        (*i)++;
        w->when = rt_clock.seconds() + 0.100;
        rt_work_queue.post(w);
    }

    static void example_1(void) {
        rt_work_queue_t queue = {0};
        // if a single pointer will suffice
        int32_t i = 0;
        rt_work_t work = {
            .queue = &queue,
            .when  = rt_clock.seconds() + 0.100,
            .work  = every_100ms,
            .data  = &i
        };
        rt_work_queue.post(&work);
        dispatch_until(&queue, &i, 4);
    }

    // extending rt_work_t with extra data:

    typedef struct rt_work_ex_s {
        union {
            rt_work_t base;
            struct rt_work_s;
        };
        struct { int32_t a; int32_t b; } s;
        int32_t i;
    } rt_work_ex_t;

    static void every_200ms(rt_work_t* w) {
        rt_work_ex_t* ex = (rt_work_ex_t*)w;
        rt_println("ex { .i: %d, .s.a: %d .s.b: %d}", ex->i, ex->s.a, ex->s.b);
        ex->i++;
        const int32_t swap = ex->s.a; ex->s.a = ex->s.b; ex->s.b = swap;
        w->when = rt_clock.seconds() + 0.200;
        rt_work_queue.post(w);
    }

    static void example_2(void) {
        rt_work_queue_t queue = {0};
        rt_work_ex_t work = {
            .queue = &queue,
            .when  = rt_clock.seconds() + 0.200,
            .work  = every_200ms,
            .data  = null,
            .s = { .a = 1, .b = 2 },
            .i = 0
        };
        rt_work_queue.post(&work.base);
        dispatch_until(&queue, &work.i, 4);
    }

    static void dispatch_until(rt_work_queue_t* q, int32_t* i, const int32_t n) {
        while (q->head != null && *i < n) {
            rt_thread.sleep_for(0.0001); // 100 microseconds
            rt_work_queue.dispatch(q);
        }
        rt_work_queue.flush(q);
    }

    // worker:

    static void do_work(rt_work_t* w) {
        // TODO: something useful
    }

    static void worker_test(void) {
        rt_worker_t worker = { 0 };
        rt_worker.start(&worker);
        rt_work_t work = {
            .when  = rt_clock.seconds() + 0.010, // 10ms
            .done  = rt_event.create(),
            .work  = do_work
        };
        rt_worker.post(&worker, &work);
        rt_event.wait(work.done);    // await(work)
        rt_event.dispose(work.done); // responsibility of the caller
        rt_fatal_if_error(rt_worker.join(&worker, -1.0));
    }

    // Hint:
    // To monitor timing turn on MSVC Output / Show Timestamp (clock button)

*/

#endif // rt_definition

#ifdef rt_implementation
// ________________________________ rt_win32.h ________________________________

#ifdef WIN32

#pragma warning(push)
#pragma warning(disable: 4255) // no function prototype: '()' to '(void)'
#pragma warning(disable: 4459) // declaration of '...' hides global declaration

#pragma push_macro("UNICODE")
#define UNICODE // always because otherwise IME does not work

// ut:
#include <Windows.h>  // used by:
#include <Psapi.h>    // both rt_loader.c and rt_processes.c
#include <shellapi.h> // rt_processes.c
#include <winternl.h> // rt_processes.c
#include <initguid.h>     // for knownfolders
#include <KnownFolders.h> // rt_files.c
#include <AclAPI.h>       // rt_files.c
#include <ShlObj_core.h>  // rt_files.c
#include <Shlwapi.h>      // rt_files.c
// ui:
#include <commdlg.h>
#include <dbghelp.h>
#include <dwmapi.h>
#include <imm.h>
#include <ShellScalingApi.h>
#include <tlhelp32.h>
#include <VersionHelpers.h>
#include <windowsx.h>
#include <winnt.h>

#pragma pop_macro("UNICODE")

#pragma warning(pop)

#include <fcntl.h>

#define rt_export __declspec(dllexport)

// Win32 API BOOL -> errno_t translation

#define rt_b2e(call) ((errno_t)(call ? 0 : GetLastError()))

void rt_win32_close_handle(void* h);
/* translate ix to error */
errno_t rt_wait_ix2e(uint32_t r);


#endif // WIN32
// ___________________________________ rt.c ___________________________________

// #include "ut/macos.h" // TODO
// #include "ut/linux.h" // TODO


// ________________________________ rt_args.c _________________________________

static void* rt_args_memory;

static void rt_args_main(int32_t argc, const char* argv[], const char** env) {
    rt_swear(rt_args.c == 0 && rt_args.v == null && rt_args.env == null);
    rt_swear(rt_args_memory == null);
    rt_args.c = argc;
    rt_args.v = argv;
    rt_args.env = env;
}

static int32_t rt_args_option_index(const char* option) {
    for (int32_t i = 1; i < rt_args.c; i++) {
        if (strcmp(rt_args.v[i], "--") == 0) { break; } // no options after '--'
        if (strcmp(rt_args.v[i], option) == 0) { return i; }
    }
    return -1;
}

static void rt_args_remove_at(int32_t ix) {
    // returns new argc
    rt_assert(0 < rt_args.c);
    rt_assert(0 < ix && ix < rt_args.c); // cannot remove rt_args.v[0]
    for (int32_t i = ix; i < rt_args.c; i++) {
        rt_args.v[i] = rt_args.v[i + 1];
    }
    rt_args.v[rt_args.c - 1] = "";
    rt_args.c--;
}

static bool rt_args_option_bool(const char* option) {
    int32_t ix = rt_args_option_index(option);
    if (ix > 0) { rt_args_remove_at(ix); }
    return ix > 0;
}

static bool rt_args_option_int(const char* option, int64_t *value) {
    int32_t ix = rt_args_option_index(option);
    if (ix > 0 && ix < rt_args.c - 1) {
        const char* s = rt_args.v[ix + 1];
        int32_t base = (strstr(s, "0x") == s || strstr(s, "0X") == s) ? 16 : 10;
        const char* b = s + (base == 10 ? 0 : 2);
        char* e = null;
        errno = 0;
        int64_t v = strtoll(b, &e, base);
        if (errno == 0 && e > b && *e == 0) {
            *value = v;
        } else {
            ix = -1;
        }
    } else {
        ix = -1;
    }
    if (ix > 0) {
        rt_args_remove_at(ix); // remove option
        rt_args_remove_at(ix); // remove following number
    }
    return ix > 0;
}

static const char* rt_args_option_str(const char* option) {
    int32_t ix = rt_args_option_index(option);
    const char* s = null;
    if (ix > 0 && ix < rt_args.c - 1) {
        s = rt_args.v[ix + 1];
    } else {
        ix = -1;
    }
    if (ix > 0) {
        rt_args_remove_at(ix); // remove option
        rt_args_remove_at(ix); // remove following string
    }
    return ix > 0 ? s : null;
}

// Terminology: "quote" in the code and comments below
// actually refers to "fp64_t quote mark" and used for brevity.

// TODO: posix like systems
// Looks like all shells support quote marks but
// AFAIK MacOS bash and zsh also allow (and prefer) backslash escaped
// space character. Unclear what other escaping shell and posix compliant
// parser should support.
// Lengthy discussion here:
// https://stackoverflow.com/questions/1706551/parse-string-into-argv-argc

// Microsoft specific argument parsing:
// https://web.archive.org/web/20231115181633/http://learn.microsoft.com/en-us/cpp/c-language/parsing-c-command-line-arguments?view=msvc-170
// Alternative: just use CommandLineToArgvW()

typedef struct { const char* s; char* d; const char* e; } rt_args_pair_t;

static rt_args_pair_t rt_args_parse_backslashes(rt_args_pair_t p) {
    enum { quote = '"', backslash = '\\' };
    const char* s = p.s;
    char* d = p.d;
    rt_swear(*s == backslash);
    int32_t bsc = 0; // number of backslashes
    while (*s == backslash) { s++; bsc++; }
    if (*s == quote) {
        while (bsc > 1 && d < p.e) { *d++ = backslash; bsc -= 2; }
        if (bsc == 1 && d < p.e) { *d++ = *s++; }
    } else {
        // Backslashes are interpreted literally,
        // unless they immediately precede a quote:
        while (bsc > 0 && d < p.e) { *d++ = backslash; bsc--; }
    }
    return (rt_args_pair_t){ .s = s, .d = d, .e = p.e };
}

static rt_args_pair_t rt_args_parse_quoted(rt_args_pair_t p) {
    enum { quote = '"', backslash = '\\' };
    const char* s = p.s;
    char* d = p.d;
    rt_swear(*s == quote);
    s++; // opening quote (skip)
    while (*s != 0x00) {
        if (*s == backslash) {
            p = rt_args_parse_backslashes((rt_args_pair_t){
                        .s = s, .d = d, .e = p.e });
            s = p.s; d = p.d;
        } else if (*s == quote && s[1] == quote) {
            // Within a quoted string, a pair of quote is
            // interpreted as a single escaped quote.
            if (d < p.e) { *d++ = *s++; }
            s++; // 1 for 2 quotes
        } else if (*s == quote) {
            s++; // closing quote (skip)
            break;
        } else if (d < p.e) {
            *d++ = *s++;
        }
    }
    return (rt_args_pair_t){ .s = s, .d = d, .e = p.e };
}

static void rt_args_parse(const char* s) {
    rt_swear(s[0] != 0, "cannot parse empty string");
    rt_swear(rt_args.c == 0);
    rt_swear(rt_args.v == null);
    rt_swear(rt_args_memory == null);
    enum { quote = '"', backslash = '\\', tab = '\t', space = 0x20 };
    const int32_t len = (int32_t)strlen(s);
    // Worst-case scenario (possible to optimize with dry run of parse)
    // at least 2 characters per token in "a b c d e" plush null at the end:
    const int32_t k = ((len + 2) / 2 + 1) * (int32_t)sizeof(void*) + (int32_t)sizeof(void*);
    const int32_t n = k + (len + 2) * (int32_t)sizeof(char);
    rt_fatal_if_error(rt_heap.allocate(null, &rt_args_memory, n, true));
    rt_args.c = 0;
    rt_args.v = (const char**)rt_args_memory;
    char* d = (char*)(((char*)rt_args.v) + k);
    char* e = d + n; // end of memory
    // special rules for 1st argument:
    if (rt_args.c < n) { rt_args.v[rt_args.c++] = d; }
    if (*s == quote) {
        s++;
        while (*s != 0x00 && *s != quote && d < e) { *d++ = *s++; }
        if (*s == quote) { // // closing quote
            s++; // skip closing quote
            *d++ = 0x00;
        } else {
            while (*s != 0x00) { s++; }
        }
    } else {
        while (*s != 0x00 && *s != space && *s != tab && d < e) {
            *d++ = *s++;
        }
    }
    if (d < e) { *d++ = 0; }
    while (d < e) {
        while (*s == space || *s == tab) { s++; }
        if (*s == 0) { break; }
        if (*s == quote && s[1] == 0 && d < e) { // unbalanced single quote
            if (rt_args.c < n) { rt_args.v[rt_args.c++] = d; } // spec does not say what to do
            *d++ = *s++;
        } else if (*s == quote) { // quoted arg
            if (rt_args.c < n) { rt_args.v[rt_args.c++] = d; }
            rt_args_pair_t p = rt_args_parse_quoted(
                    (rt_args_pair_t){ .s = s, .d = d, .e = e });
            s = p.s; d = p.d;
        } else { // non-quoted arg (that can have quoted strings inside)
            if (rt_args.c < n) { rt_args.v[rt_args.c++] = d; }
            while (*s != 0) {
                if (*s == backslash) {
                    rt_args_pair_t p = rt_args_parse_backslashes(
                            (rt_args_pair_t){ .s = s, .d = d, .e = e });
                    s = p.s; d = p.d;
                } else if (*s == quote) {
                    rt_args_pair_t p = rt_args_parse_quoted(
                            (rt_args_pair_t){ .s = s, .d = d, .e = e });
                    s = p.s; d = p.d;
                } else if (*s == tab || *s == space) {
                    break;
                } else if (d < e) {
                    *d++ = *s++;
                }
            }
        }
        if (d < e) { *d++ = 0; }
    }
    if (rt_args.c < n) {
        rt_args.v[rt_args.c] = null;
    }
    rt_swear(rt_args.c < n, "not enough memory - adjust guestimates");
    rt_swear(d <= e, "not enough memory - adjust guestimates");
}

static const char* rt_args_basename(void) {
    static char basename[260];
    rt_swear(rt_args.c > 0);
    if (basename[0] == 0) {
        const char* s = rt_args.v[0];
        const char* b = s;
        while (*s != 0) {
            if (*s == '\\' || *s == '/') { b = s + 1; }
            s++;
        }
        int32_t n = rt_str.len(b);
        rt_swear(n < rt_countof(basename));
        strncpy(basename, b, rt_countof(basename) - 1);
        char* d = basename + n - 1;
        while (d > basename && *d != '.') { d--; }
        if (*d == '.') { *d = 0x00; }
    }
    return basename;
}

static void rt_args_fini(void) {
    rt_heap.deallocate(null, rt_args_memory); // can be null is parse() was not called
    rt_args_memory = null;
    rt_args.c = 0;
    rt_args.v = null;
}

static void rt_args_WinMain(void) {
    rt_swear(rt_args.c == 0 && rt_args.v == null && rt_args.env == null);
    rt_swear(rt_args_memory == null);
    const uint16_t* wcl = GetCommandLineW();
    int32_t n = (int32_t)rt_str.len16(wcl);
    char* cl = null;
    rt_fatal_if_error(rt_heap.allocate(null, (void**)&cl, n * 2 + 1, false));
    rt_str.utf16to8(cl, n * 2 + 1, wcl, -1);
    rt_args_parse(cl);
    rt_heap.deallocate(null, cl);
    rt_args.env = (const char**)(void*)_environ;
}

#ifdef RT_TESTS

// https://learn.microsoft.com/en-us/cpp/c-language/parsing-c-command-line-arguments
// Command-line input       argv[1]     argv[2]	    argv[3]
// "a b c" d e	            a b c       d           e
// "ab\"c" "\\" d           ab"c        \           d
// a\\\b d"e f"g h          a\\\b       de fg       h
// a\\\"b c d               a\"b        c           d
// a\\\\"b c" d e           a\\b c      d           e
// a"b"" c d                ab" c d

#ifndef __INTELLISENSE__ // confused data analysis

static void rt_args_test_verify(const char* cl, int32_t expected, ...) {
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
        rt_println("cl: `%s`", cl);
    }
    int32_t argc = rt_args.c;
    const char** argv = rt_args.v;
    void* memory = rt_args_memory;
    rt_args.c = 0;
    rt_args.v = null;
    rt_args_memory = null;
    rt_args_parse(cl);
    va_list va;
    va_start(va, expected);
    for (int32_t i = 0; i < expected; i++) {
        const char* s = va_arg(va, const char*);
//      if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
//          rt_println("rt_args.v[%d]: `%s` expected: `%s`", i, rt_args.v[i], s);
//      }
        // Warning 6385: reading data outside of array
        const char* ai = _Pragma("warning(suppress:  6385)")rt_args.v[i];
        rt_swear(strcmp(ai, s) == 0, "rt_args.v[%d]: `%s` expected: `%s`",
                 i, ai, s);
    }
    va_end(va);
    rt_args.fini();
    // restore command line arguments:
    rt_args.c = argc;
    rt_args.v = argv;
    rt_args_memory = memory;
}

#endif // __INTELLISENSE__

static void rt_args_test(void) {
    // The first argument (rt_args.v[0]) is treated specially.
    // It represents the program name. Because it must be a valid pathname,
    // parts surrounded by quote (") are allowed. The quote aren't included
    // in the rt_args.v[0] output. The parts surrounded by quote prevent
    // interpretation of a space or tab character as the end of the argument.
    // The escaping rules don't apply.
    rt_args_test_verify("\"c:\\foo\\bar\\snafu.exe\"", 1,
                     "c:\\foo\\bar\\snafu.exe");
    rt_args_test_verify("c:\\foo\\bar\\snafu.exe", 1,
                     "c:\\foo\\bar\\snafu.exe");
    rt_args_test_verify("foo.exe \"a b c\" d e", 4,
                     "foo.exe", "a b c", "d", "e");
    rt_args_test_verify("foo.exe \"ab\\\"c\" \"\\\\\" d", 4,
                     "foo.exe", "ab\"c", "\\", "d");
    rt_args_test_verify("foo.exe a\\\\\\b d\"e f\"g h", 4,
                     "foo.exe", "a\\\\\\b", "de fg", "h");
    rt_args_test_verify("foo.exe a\\\\\\b d\"e f\"g h", 4,
                     "foo.exe", "a\\\\\\b", "de fg", "h");
    rt_args_test_verify("foo.exe a\"b\"\" c d", 2, // unmatched quote
                     "foo.exe", "ab\" c d");
    // unbalanced quote and backslash:
    rt_args_test_verify("foo.exe \"",     2, "foo.exe", "\"");
    rt_args_test_verify("foo.exe \\",     2, "foo.exe", "\\");
    rt_args_test_verify("foo.exe \\\\",   2, "foo.exe", "\\\\");
    rt_args_test_verify("foo.exe \\\\\\", 2, "foo.exe", "\\\\\\");
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_args_test(void) {}

#endif

rt_args_if rt_args = {
    .main         = rt_args_main,
    .WinMain      = rt_args_WinMain,
    .option_index = rt_args_option_index,
    .remove_at    = rt_args_remove_at,
    .option_bool  = rt_args_option_bool,
    .option_int   = rt_args_option_int,
    .option_str   = rt_args_option_str,
    .basename     = rt_args_basename,
    .fini         = rt_args_fini,
    .test         = rt_args_test
};
// _______________________________ rt_atomics.c _______________________________

#include <stdatomic.h> // needs cl.exe /experimental:c11atomics command line

// see: https://developercommunity.visualstudio.com/t/C11--C17-include-stdatomich-issue/10620622

#pragma warning(push)
#pragma warning(disable: 4746) // volatile access of 'int32_var' is subject to /volatile:<iso|ms> setting; consider using __iso_volatile_load/store intrinsic functions

#ifndef UT_ATOMICS_HAS_STDATOMIC_H

static int32_t rt_atomics_increment_int32(volatile int32_t* a) {
    return InterlockedIncrement((volatile LONG*)a);
}

static int32_t rt_atomics_decrement_int32(volatile int32_t* a) {
    return InterlockedDecrement((volatile LONG*)a);
}

static int64_t rt_atomics_increment_int64(volatile int64_t* a) {
    return InterlockedIncrement64((__int64 volatile *)a);
}

static int64_t rt_atomics_decrement_int64(volatile int64_t* a) {
    return InterlockedDecrement64((__int64 volatile *)a);
}

static int32_t rt_atomics_add_int32(volatile int32_t* a, int32_t v) {
    return InterlockedAdd((LONG volatile *)a, v);
}

static int64_t rt_atomics_add_int64(volatile int64_t* a, int64_t v) {
    return InterlockedAdd64((__int64 volatile *)a, v);
}

static int64_t rt_atomics_exchange_int64(volatile int64_t* a, int64_t v) {
    return (int64_t)InterlockedExchange64((LONGLONG*)a, (LONGLONG)v);
}

static int32_t rt_atomics_exchange_int32(volatile int32_t* a, int32_t v) {
    rt_assert(sizeof(int32_t) == sizeof(unsigned long));
    return (int32_t)InterlockedExchange((volatile LONG*)a, (unsigned long)v);
}

static bool rt_atomics_compare_exchange_int64(volatile int64_t* a,
        int64_t comparand, int64_t v) {
    return (int64_t)InterlockedCompareExchange64((LONGLONG*)a,
        (LONGLONG)v, (LONGLONG)comparand) == comparand;
}

static bool rt_atomics_compare_exchange_int32(volatile int32_t* a,
        int32_t comparand, int32_t v) {
    return (int64_t)InterlockedCompareExchange((LONG*)a,
        (LONG)v, (LONG)comparand) == comparand;
}

static void memory_fence(void) {
#ifdef _M_ARM64
atomic_thread_fence(memory_order_seq_cst);
#else
_mm_mfence();
#endif
}

#else

// stdatomic.h version:

#ifndef __INTELLISENSE__ // IntelliSense chokes on _Atomic(_Type)
// __INTELLISENSE__ Defined as 1 during an IntelliSense compiler pass
// in the Visual Studio IDE. Otherwise, undefined. You can use this macro
// to guard code the IntelliSense compiler doesn't understand,
// or use it to toggle between the build and IntelliSense compiler.


// _strong() operations are the same as _explicit(..., memory_order_seq_cst)
// memory_order_seq_cst stands for Sequentially Consistent Ordering
//
// This is the strongest memory order, providing the guarantee that
// all sequentially consistent operations appear to be executed in
// the same order on all threads (cores)
//
// int_fast32_t: Fastest integer type with at least 32 bits.
// int_least32_t: Smallest integer type with at least 32 bits.

rt_static_assertion(sizeof(int32_t) == sizeof(int_fast32_t));
rt_static_assertion(sizeof(int32_t) == sizeof(int_least32_t));

static int32_t rt_atomics_increment_int32(volatile int32_t* a) {
    return atomic_fetch_add((volatile atomic_int_fast32_t*)a, 1) + 1;
}

static int32_t rt_atomics_decrement_int32(volatile int32_t* a) {
    return atomic_fetch_sub((volatile atomic_int_fast32_t*)a, 1) - 1;
}

static int64_t rt_atomics_increment_int64(volatile int64_t* a) {
    return atomic_fetch_add((volatile atomic_int_fast64_t*)a, 1) + 1;
}

static int64_t rt_atomics_decrement_int64(volatile int64_t* a) {
    return atomic_fetch_sub((volatile atomic_int_fast64_t*)a, 1) - 1;
}

static int32_t rt_atomics_add_int32(volatile int32_t* a, int32_t v) {
    return atomic_fetch_add((volatile atomic_int_fast32_t*)a, v) + v;
}

static int64_t rt_atomics_add_int64(volatile int64_t* a, int64_t v) {
    return atomic_fetch_add((volatile atomic_int_fast64_t*)a, v) + v;
}

static int64_t rt_atomics_exchange_int64(volatile int64_t* a, int64_t v) {
    return atomic_exchange((volatile atomic_int_fast64_t*)a, v);
}

static int32_t rt_atomics_exchange_int32(volatile int32_t* a, int32_t v) {
    return atomic_exchange((volatile atomic_int_fast32_t*)a, v);
}

static bool rt_atomics_compare_exchange_int64(volatile int64_t* a,
    int64_t comparand, int64_t v) {
    return atomic_compare_exchange_strong((volatile atomic_int_fast64_t*)a,
        &comparand, v);
}

// Code here is not "seen" by IntelliSense but is compiled normally.
static bool rt_atomics_compare_exchange_int32(volatile int32_t* a,
    int32_t comparand, int32_t v) {
    return atomic_compare_exchange_strong((volatile atomic_int_fast32_t*)a,
        &comparand, v);
}

static void memory_fence(void) { atomic_thread_fence(memory_order_seq_cst); }

#endif // __INTELLISENSE__

#endif // UT_ATOMICS_HAS_STDATOMIC_H

static int32_t rt_atomics_load_int32(volatile int32_t* a) {
    return rt_atomics.add_int32(a, 0);
}

static int64_t rt_atomics_load_int64(volatile int64_t* a) {
    return rt_atomics.add_int64(a, 0);
}

static void* rt_atomics_exchange_ptr(volatile void* *a, void* v) {
    rt_static_assertion(sizeof(void*) == sizeof(uint64_t));
    return (void*)(intptr_t)rt_atomics.exchange_int64((int64_t*)a, (int64_t)v);
}

static bool rt_atomics_compare_exchange_ptr(volatile void* *a, void* comparand, void* v) {
    rt_static_assertion(sizeof(void*) == sizeof(int64_t));
    return rt_atomics.compare_exchange_int64((int64_t*)a,
        (int64_t)comparand, (int64_t)v);
}

#pragma push_macro("rt_sync_bool_compare_and_swap")
#pragma push_macro("rt_builtin_cpu_pause")

// https://en.wikipedia.org/wiki/Spinlock

#define rt_sync_bool_compare_and_swap(p, old_val, new_val)          \
    (_InterlockedCompareExchange64(p, new_val, old_val) == old_val)

// https://stackoverflow.com/questions/37063700/mm-pause-usage-in-gcc-on-intel
#define rt_builtin_cpu_pause() do { YieldProcessor(); } while (0)

static void spinlock_acquire(volatile int64_t* spinlock) {
    // Very basic implementation of a spinlock. This is currently
    // only used to guarantee thread-safety during context initialization
    // and shutdown (which are both executed very infrequently and
    // have minimal thread contention).
    // Not a performance champion (because of mem_fence()) but serves
    // the purpose. mem_fence() can be reduced to mem_sfence()... sigh
    while (!rt_sync_bool_compare_and_swap(spinlock, 0, 1)) {
        while (*spinlock) { rt_builtin_cpu_pause(); }
    }
    rt_atomics.memory_fence();
    // not strictly necessary on strong mem model Intel/AMD but
    // see: https://cfsamsonbooks.gitbook.io/explaining-atomics-in-rust/
    //      Fig 2 Inconsistent C11 execution of SB and 2+2W
    rt_assert(*spinlock == 1);
}

#pragma pop_macro("rt_builtin_cpu_pause")
#pragma pop_macro("rt_sync_bool_compare_and_swap")

static void spinlock_release(volatile int64_t* spinlock) {
    rt_assert(*spinlock == 1);
    *spinlock = 0;
    // tribute to lengthy Linus discussion going since 2006:
    rt_atomics.memory_fence();
}

static void rt_atomics_test(void) {
    #ifdef RT_TESTS
    volatile int32_t int32_var = 0;
    volatile int64_t int64_var = 0;
    volatile void* ptr_var = null;
    int64_t spinlock = 0;
    void* old_ptr = rt_atomics.exchange_ptr(&ptr_var, (void*)123);
    rt_swear(old_ptr == null);
    rt_swear(ptr_var == (void*)123);
    int32_t incremented_int32 = rt_atomics.increment_int32(&int32_var);
    rt_swear(incremented_int32 == 1);
    rt_swear(int32_var == 1);
    int32_t decremented_int32 = rt_atomics.decrement_int32(&int32_var);
    rt_swear(decremented_int32 == 0);
    rt_swear(int32_var == 0);
    int64_t incremented_int64 = rt_atomics.increment_int64(&int64_var);
    rt_swear(incremented_int64 == 1);
    rt_swear(int64_var == 1);
    int64_t decremented_int64 = rt_atomics.decrement_int64(&int64_var);
    rt_swear(decremented_int64 == 0);
    rt_swear(int64_var == 0);
    int32_t added_int32 = rt_atomics.add_int32(&int32_var, 5);
    rt_swear(added_int32 == 5);
    rt_swear(int32_var == 5);
    int64_t added_int64 = rt_atomics.add_int64(&int64_var, 10);
    rt_swear(added_int64 == 10);
    rt_swear(int64_var == 10);
    int32_t old_int32 = rt_atomics.exchange_int32(&int32_var, 3);
    rt_swear(old_int32 == 5);
    rt_swear(int32_var == 3);
    int64_t old_int64 = rt_atomics.exchange_int64(&int64_var, 6);
    rt_swear(old_int64 == 10);
    rt_swear(int64_var == 6);
    bool int32_exchanged = rt_atomics.compare_exchange_int32(&int32_var, 3, 4);
    rt_swear(int32_exchanged);
    rt_swear(int32_var == 4);
    bool int64_exchanged = rt_atomics.compare_exchange_int64(&int64_var, 6, 7);
    rt_swear(int64_exchanged);
    rt_swear(int64_var == 7);
    ptr_var = (void*)0x123;
    bool ptr_exchanged = rt_atomics.compare_exchange_ptr(&ptr_var,
        (void*)0x123, (void*)0x456);
    rt_swear(ptr_exchanged);
    rt_swear(ptr_var == (void*)0x456);
    rt_atomics.spinlock_acquire(&spinlock);
    rt_swear(spinlock == 1);
    rt_atomics.spinlock_release(&spinlock);
    rt_swear(spinlock == 0);
    int32_t loaded_int32 = rt_atomics.load32(&int32_var);
    rt_swear(loaded_int32 == int32_var);
    int64_t loaded_int64 = rt_atomics.load64(&int64_var);
    rt_swear(loaded_int64 == int64_var);
    rt_atomics.memory_fence();
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
    #endif
}

#ifndef __INTELLISENSE__ // IntelliSense chokes on _Atomic(_Type)

rt_static_assertion(sizeof(void*) == sizeof(int64_t));
rt_static_assertion(sizeof(void*) == sizeof(uintptr_t));

rt_atomics_if rt_atomics = {
    .exchange_ptr    = rt_atomics_exchange_ptr,
    .increment_int32 = rt_atomics_increment_int32,
    .decrement_int32 = rt_atomics_decrement_int32,
    .increment_int64 = rt_atomics_increment_int64,
    .decrement_int64 = rt_atomics_decrement_int64,
    .add_int32 = rt_atomics_add_int32,
    .add_int64 = rt_atomics_add_int64,
    .exchange_int32  = rt_atomics_exchange_int32,
    .exchange_int64  = rt_atomics_exchange_int64,
    .compare_exchange_int64 = rt_atomics_compare_exchange_int64,
    .compare_exchange_int32 = rt_atomics_compare_exchange_int32,
    .compare_exchange_ptr = rt_atomics_compare_exchange_ptr,
    .load32 = rt_atomics_load_int32,
    .load64 = rt_atomics_load_int64,
    .spinlock_acquire = spinlock_acquire,
    .spinlock_release = spinlock_release,
    .memory_fence = memory_fence,
    .test = rt_atomics_test
};

#endif // __INTELLISENSE__

// 2024-03-20 latest windows runtime and toolchain cl.exe
// ... VC\Tools\MSVC\14.39.33519\include
// see:
//     vcruntime_c11_atomic_support.h
//     vcruntime_c11_stdatomic.h
//     stdatomic.h
// https://developercommunity.visualstudio.com/t/C11--C17-include--issue/10620622
// cl.exe /std:c11 /experimental:c11atomics
// command line option are required
// even in C17 mode in spring of 2024

#pragma warning(pop)

// ______________________________ rt_backtrace.c ______________________________

static void* rt_backtrace_process;
static DWORD rt_backtrace_pid;

typedef rt_begin_packed struct symbol_info_s {
    SYMBOL_INFO info; char name[rt_backtrace_max_symbol];
} rt_end_packed symbol_info_t;

#pragma push_macro("rt_backtrace_load_dll")

#define rt_backtrace_load_dll(fn) do {              \
    if (GetModuleHandleA(fn) == null) {      \
        rt_fatal_win32err(LoadLibraryA(fn)); \
    }                                        \
} while (0)

static void rt_backtrace_init(void) {
    if (rt_backtrace_process == null) {
        rt_backtrace_load_dll("dbghelp.dll");
        rt_backtrace_load_dll("imagehlp.dll");
        DWORD options = SymGetOptions();
//      options |= SYMOPT_DEBUG;
        options |= SYMOPT_NO_PROMPTS;
        options |= SYMOPT_LOAD_LINES;
        options |= SYMOPT_UNDNAME;
        options |= SYMOPT_LOAD_ANYTHING;
        rt_swear(SymSetOptions(options));
        rt_backtrace_pid = GetProcessId(GetCurrentProcess());
        rt_swear(rt_backtrace_pid != 0);
        rt_backtrace_process = OpenProcess(PROCESS_ALL_ACCESS, false,
                                           rt_backtrace_pid);
        rt_swear(rt_backtrace_process != null);
        rt_swear(SymInitialize(rt_backtrace_process, null, true), "%s",
                            rt_str.error(rt_core.err()));
    }
}

#pragma pop_macro("rt_backtrace_load_dll")

static void rt_backtrace_capture(rt_backtrace_t* bt, int32_t skip) {
    rt_backtrace_init();
    SetLastError(0);
    bt->frames = CaptureStackBackTrace(1 + skip, rt_countof(bt->stack),
        bt->stack, (DWORD*)&bt->hash);
    bt->error = rt_core.err();
}

static bool rt_backtrace_function(DWORD64 pc, SYMBOL_INFO* si) {
    // find DLL exported function
    bool found = false;
    const DWORD64 module_base = SymGetModuleBase64(rt_backtrace_process, pc);
    if (module_base != 0) {
        const DWORD flags = GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT |
                            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS;
        HMODULE module_handle = null;
        if (GetModuleHandleExA(flags, (const char*)pc, &module_handle)) {
            DWORD bytes = 0;
            IMAGE_EXPORT_DIRECTORY* dir = (IMAGE_EXPORT_DIRECTORY*)
                    ImageDirectoryEntryToDataEx(module_handle, true,
                            IMAGE_DIRECTORY_ENTRY_EXPORT, &bytes, null);
            if (dir) {
                uint8_t* m = (uint8_t*)module_handle;
                DWORD* functions = (DWORD*)(m + dir->AddressOfFunctions);
                DWORD* names = (DWORD*)(m + dir->AddressOfNames);
                WORD* ordinals = (WORD*)(m + dir->AddressOfNameOrdinals);
                DWORD64 address = 0; // closest address
                DWORD64 min_distance = (DWORD64)-1;
                const char* function = NULL; // closest function name
                for (DWORD i = 0; i < dir->NumberOfNames; i++) {
                    // function address
                    DWORD64 fa = (DWORD64)(m + functions[ordinals[i]]);
                    if (fa <= pc) {
                        DWORD64 distance = pc - fa;
                        if (distance < min_distance) {
                            min_distance = distance;
                            address = fa;
                            function = (const char*)(m + names[i]);
                        }
                    }
                }
                if (function != null) {
                    si->ModBase = (uint64_t)m;
                    snprintf(si->Name, si->MaxNameLen - 1, "%s", function);
                    si->Name[si->MaxNameLen - 1] = 0x00;
                    si->NameLen = (DWORD)strlen(si->Name);
                    si->Address = address;
                    found = true;
                }
            }
        }
    }
    return found;
}

// SimpleStackWalker::showVariablesAt() can be implemented if needed like this:
// https://accu.org/journals/overload/29/165/orr/
// https://github.com/rogerorr/articles/tree/main/Debugging_Optimised_Code
// https://github.com/rogerorr/articles/blob/main/Debugging_Optimised_Code/SimpleStackWalker.cpp#L301

static const void rt_backtrace_symbolize_inline_frame(rt_backtrace_t* bt,
        int32_t i, DWORD64 pc, DWORD inline_context, symbol_info_t* si) {
    si->info.Name[0] = 0;
    si->info.NameLen = 0;
    bt->file[i][0] = 0;
    bt->line[i] = 0;
    bt->symbol[i][0] = 0;
    DWORD64 displacement = 0;
    if (SymFromInlineContext(rt_backtrace_process, pc, inline_context,
                            &displacement, &si->info)) {
        rt_str_printf(bt->symbol[i], "%s", si->info.Name);
    } else {
        bt->error = rt_core.err();
    }
    IMAGEHLP_LINE64 li = { .SizeOfStruct = sizeof(IMAGEHLP_LINE64) };
    DWORD offset = 0;
    if (SymGetLineFromInlineContext(rt_backtrace_process,
                                    pc, inline_context, 0,
                                    &offset, &li)) {
        rt_str_printf(bt->file[i], "%s", li.FileName);
        bt->line[i] = li.LineNumber;
    }
}

// Too see kernel addresses in Stack Back Traces:
//
// Windows Registry Editor Version 5.00
// [HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Memory Management]
// "DisablePagingExecutive"=dword:00000001
//
// https://learn.microsoft.com/en-us/previous-versions/windows/it-pro/windows-server-2003/cc757875(v=ws.10)

static int32_t rt_backtrace_symbolize_frame(rt_backtrace_t* bt, int32_t i) {
    const DWORD64 pc = (DWORD64)bt->stack[i];
    symbol_info_t si = {
        .info = { .SizeOfStruct = sizeof(SYMBOL_INFO),
                  .MaxNameLen = rt_countof(si.name) }
    };
    bt->file[i][0] = 0;
    bt->line[i] = 0;
    bt->symbol[i][0] = 0;
    DWORD64 offsetFromSymbol = 0;
    const DWORD inline_count =
        SymAddrIncludeInlineTrace(rt_backtrace_process, pc);
    if (inline_count > 0) {
        DWORD ic = 0; // inline context
        DWORD fi = 0; // frame index
        if (SymQueryInlineTrace(rt_backtrace_process,
                                pc, 0, pc, pc, &ic, &fi)) {
            for (DWORD k = 0; k < inline_count; k++, ic++) {
                rt_backtrace_symbolize_inline_frame(bt, i, pc, ic, &si);
                i++;
            }
        }
    } else {
        if (SymFromAddr(rt_backtrace_process, pc, &offsetFromSymbol, &si.info)) {
            rt_str_printf(bt->symbol[i], "%s", si.info.Name);
            DWORD d = 0; // displacement
            IMAGEHLP_LINE64 ln = { .SizeOfStruct = sizeof(IMAGEHLP_LINE64) };
            if (SymGetLineFromAddr64(rt_backtrace_process, pc, &d, &ln)) {
                bt->line[i] = ln.LineNumber;
                rt_str_printf(bt->file[i], "%s", ln.FileName);
            } else {
                bt->error = rt_core.err();
                if (rt_backtrace_function(pc, &si.info)) {
                    GetModuleFileNameA((HANDLE)si.info.ModBase, bt->file[i],
                        rt_countof(bt->file[i]) - 1);
                    bt->file[i][rt_countof(bt->file[i]) - 1] = 0;
                    bt->line[i]    = 0;
                } else  {
                    bt->file[i][0] = 0x00;
                    bt->line[i]    = 0;
                }
            }
            i++;
        } else {
            bt->error = rt_core.err();
            if (rt_backtrace_function(pc, &si.info)) {
                rt_str_printf(bt->symbol[i], "%s", si.info.Name);
                GetModuleFileNameA((HANDLE)si.info.ModBase, bt->file[i],
                    rt_countof(bt->file[i]) - 1);
                bt->file[i][rt_countof(bt->file[i]) - 1] = 0;
                bt->error = 0;
                i++;
            } else {
                // will not do i++
            }
        }
    }
    return i;
}

static void rt_backtrace_symbolize_backtrace(rt_backtrace_t* bt) {
    rt_assert(!bt->symbolized);
    bt->error = 0;
    rt_backtrace_init();
    // rt_backtrace_symbolize_frame() may produce zero, one or many frames
    int32_t n = bt->frames;
    void* stack[rt_countof(bt->stack)];
    memcpy(stack, bt->stack, n * sizeof(stack[0]));
    bt->frames = 0;
    for (int32_t i = 0; i < n && bt->frames < rt_countof(bt->stack); i++) {
        bt->stack[bt->frames] = stack[i];
        bt->frames = rt_backtrace_symbolize_frame(bt, i);
    }
    bt->symbolized = true;
}

static void rt_backtrace_symbolize(rt_backtrace_t* bt) {
    if (!bt->symbolized) { rt_backtrace_symbolize_backtrace(bt); }
}

static const char* rt_backtrace_stops[] = {
    "main",
    "WinMain",
    "BaseThreadInitThunk",
    "RtlUserThreadStart",
    "mainCRTStartup",
    "WinMainCRTStartup",
    "invoke_main",
    "NdrInterfacePointerMemorySize",
    null
};

static void rt_backtrace_trace(const rt_backtrace_t* bt, const char* stop) {
    #pragma push_macro("rt_backtrace_glyph_called_from")
    #define rt_backtrace_glyph_called_from rt_glyph_north_west_arrow_with_hook
    rt_assert(bt->symbolized, "need rt_backtrace.symbolize(bt)");
    const char** alt = stop != null && strcmp(stop, "*") == 0 ?
                       rt_backtrace_stops : null;
    for (int32_t i = 0; i < bt->frames; i++) {
        rt_debug.println(bt->file[i], bt->line[i], bt->symbol[i],
            rt_backtrace_glyph_called_from "%s",
            i == i < bt->frames - 1 ? "\n" : ""); // extra \n for last line
        if (stop != null && strcmp(bt->symbol[i], stop) == 0) { break; }
        const char** s = alt;
        while (s != null && *s != null && strcmp(bt->symbol[i], *s) != 0) { s++; }
        if (s != null && *s != null)  { break; }
    }
    #pragma pop_macro("rt_backtrace_glyph_called_from")
}


static const char* rt_backtrace_string(const rt_backtrace_t* bt,
        char* text, int32_t count) {
    rt_assert(bt->symbolized, "need rt_backtrace.symbolize(bt)");
    char s[1024];
    char* p = text;
    int32_t n = count;
    for (int32_t i = 0; i < bt->frames && n > 128; i++) {
        int32_t line = bt->line[i];
        const char* file = bt->file[i];
        const char* name = bt->symbol[i];
        if (file[0] != 0 && name[0] != 0) {
            rt_str_printf(s, "%s(%d): %s\n", file, line, name);
        } else if (file[0] == 0 && name[0] != 0) {
            rt_str_printf(s, "%s\n", name);
        }
        s[rt_countof(s) - 1] = 0;
        int32_t k = (int32_t)strlen(s);
        if (k < n) {
            memcpy(p, s, (size_t)k + 1);
            p += k;
            n -= k;
        }
    }
    return text;
}

typedef struct { char name[32]; } rt_backtrace_thread_name_t;

static rt_backtrace_thread_name_t rt_backtrace_thread_name(HANDLE thread) {
    rt_backtrace_thread_name_t tn;
    tn.name[0] = 0;
    wchar_t* thread_name = null;
    if (SUCCEEDED(GetThreadDescription(thread, &thread_name))) {
        rt_str.utf16to8(tn.name, rt_countof(tn.name), thread_name, -1);
        LocalFree(thread_name);
    }
    return tn;
}

static void rt_backtrace_context(rt_thread_t thread, const void* ctx,
        rt_backtrace_t* bt) {
    CONTEXT* context = (CONTEXT*)ctx;
    STACKFRAME64 stack_frame = { 0 };
    int machine_type = IMAGE_FILE_MACHINE_UNKNOWN;
    #if defined(_M_IX86)
        #error "Unsupported platform"
    #elif defined(_M_ARM64)
        machine_type = IMAGE_FILE_MACHINE_ARM64;
        stack_frame = (STACKFRAME64){
            .AddrPC    = {.Offset = context->Pc, .Mode = AddrModeFlat},
            .AddrFrame = {.Offset = context->Fp, .Mode = AddrModeFlat},
            .AddrStack = {.Offset = context->Sp, .Mode = AddrModeFlat}
        };
    #elif defined(_M_X64)
        machine_type = IMAGE_FILE_MACHINE_AMD64;
        stack_frame = (STACKFRAME64){
            .AddrPC    = {.Offset = context->Rip, .Mode = AddrModeFlat},
            .AddrFrame = {.Offset = context->Rbp, .Mode = AddrModeFlat},
            .AddrStack = {.Offset = context->Rsp, .Mode = AddrModeFlat}
        };
    #elif defined(_M_IA64)
        int machine_type = IMAGE_FILE_MACHINE_IA64;
        stack_frame = (STACKFRAME64){
            .AddrPC     = {.Offset = context->StIIP, .Mode = AddrModeFlat},
            .AddrFrame  = {.Offset = context->IntSp, .Mode = AddrModeFlat},
            .AddrBStore = {.Offset = context->RsBSP, .Mode = AddrModeFlat},
            .AddrStack  = {.Offset = context->IntSp, .Mode = AddrModeFlat}
        }
    #elif defined(_M_ARM64)
        machine_type = IMAGE_FILE_MACHINE_ARM64;
        stack_frame = (STACKFRAME64){
            .AddrPC    = {.Offset = context->Pc, .Mode = AddrModeFlat},
            .AddrFrame = {.Offset = context->Fp, .Mode = AddrModeFlat},
            .AddrStack = {.Offset = context->Sp, .Mode = AddrModeFlat}
        };
    #else
        #error "Unsupported platform"
    #endif
    rt_backtrace_init();
    while (StackWalk64(machine_type, rt_backtrace_process,
            (HANDLE)thread, &stack_frame, context, null,
            SymFunctionTableAccess64, SymGetModuleBase64, null)) {
        DWORD64 pc = stack_frame.AddrPC.Offset;
        if (pc == 0) { break; }
        if (bt->frames < rt_countof(bt->stack)) {
            bt->stack[bt->frames] = (void*)pc;
            bt->frames = rt_backtrace_symbolize_frame(bt, bt->frames);
        }
    }
    bt->symbolized = true;
}

static void rt_backtrace_thread(HANDLE thread, rt_backtrace_t* bt) {
    bt->frames = 0;
    // cannot suspend callers thread
    rt_swear(rt_thread.id_of(thread) != rt_thread.id());
    if (SuspendThread(thread) != (DWORD)-1) {
        CONTEXT context = { .ContextFlags = CONTEXT_FULL };
        GetThreadContext(thread, &context);
        rt_backtrace.context(thread, &context, bt);
        if (ResumeThread(thread) == (DWORD)-1) {
            rt_println("ResumeThread() failed %s", rt_str.error(rt_core.err()));
            ExitProcess(0xBD);
        }
    }
}

static void rt_backtrace_trace_self(const char* stop) {
    rt_backtrace_t bt = {{0}};
    rt_backtrace.capture(&bt, 2);
    rt_backtrace.symbolize(&bt);
    rt_backtrace.trace(&bt, stop);
}

static void rt_backtrace_trace_all_but_self(void) {
    rt_backtrace_init();
    rt_assert(rt_backtrace_process != null && rt_backtrace_pid != 0);
    HANDLE snapshot = CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0);
    if (snapshot == INVALID_HANDLE_VALUE) {
        rt_println("CreateToolhelp32Snapshot failed %s",
                rt_str.error(rt_core.err()));
    } else {
        THREADENTRY32 te = { .dwSize = sizeof(THREADENTRY32) };
        if (!Thread32First(snapshot, &te)) {
            rt_println("Thread32First failed %s", rt_str.error(rt_core.err()));
        } else {
            do {
                if (te.th32OwnerProcessID == rt_backtrace_pid) {
                    static const DWORD flags = THREAD_ALL_ACCESS |
                       THREAD_SUSPEND_RESUME | THREAD_GET_CONTEXT;
                    uint32_t tid = te.th32ThreadID;
                    if (tid != (uint32_t)rt_thread.id()) {
                        HANDLE thread = OpenThread(flags, false, tid);
                        if (thread != null) {
                            rt_backtrace_t bt = {0};
                            rt_backtrace_thread(thread, &bt);
                            rt_backtrace_thread_name_t tn = rt_backtrace_thread_name(thread);
                            rt_debug.println(">Thread", tid, tn.name,
                                "id 0x%08X (%d)", tid, tid);
                            if (bt.frames > 0) {
                                rt_backtrace.trace(&bt, "*");
                            }
                            rt_debug.println("<Thread", tid, tn.name, "");
                            rt_win32_close_handle(thread);
                        }
                    }
                }
            } while (Thread32Next(snapshot, &te));
        }
        rt_win32_close_handle(snapshot);
    }
}

#ifdef RT_TESTS

static bool (*rt_backtrace_debug_tee)(const char* s, int32_t count);

static char  rt_backtrace_test_output[16 * 1024];
static char* rt_backtrace_test_output_p;

static bool rt_backtrace_tee(const char* s, int32_t count) {
    if (count > 0 && s[count - 1] == 0) { // zero terminated
        int32_t k = (int32_t)(uintptr_t)(
            rt_backtrace_test_output_p - rt_backtrace_test_output);
        int32_t space = rt_countof(rt_backtrace_test_output) - k;
        if (count < space) {
            memcpy(rt_backtrace_test_output_p, s, count);
            rt_backtrace_test_output_p += count - 1; // w/o 0x00
        }
    } else {
        rt_debug.breakpoint(); // incorrect output() cannot append
    }
    return true; // intercepted, do not do OutputDebugString()
}

static void rt_backtrace_test_thread(void* e) {
    rt_event.wait(*(rt_event_t*)e);
}

static void rt_backtrace_test(void) {
    rt_backtrace_debug_tee = rt_debug.tee;
    rt_backtrace_test_output_p = rt_backtrace_test_output;
    rt_backtrace_test_output[0] = 0x00;
    rt_debug.tee = rt_backtrace_tee;
    rt_backtrace_t bt = {{0}};
    rt_backtrace.capture(&bt, 0);
    // rt_backtrace_test <- rt_core_test <- run <- main
    rt_swear(bt.frames >= 3);
    rt_backtrace.symbolize(&bt);
    rt_backtrace.trace(&bt, null);
    rt_backtrace.trace(&bt, "main");
    rt_backtrace.trace(&bt, null);
    rt_backtrace.trace(&bt, "main");
    rt_event_t e = rt_event.create();
    rt_thread_t thread = rt_thread.start(rt_backtrace_test_thread, &e);
    rt_backtrace.trace_all_but_self();
    rt_event.set(e);
    rt_thread.join(thread, -1.0);
    rt_event.dispose(e);
    rt_debug.tee = rt_backtrace_debug_tee;
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
        rt_debug.output(rt_backtrace_test_output,
            (int32_t)strlen(rt_backtrace_test_output) + 1);
    }
    rt_swear(strstr(rt_backtrace_test_output, "rt_backtrace_test") != null,
          "%s", rt_backtrace_test_output);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_backtrace_test(void) { }

#endif

rt_backtrace_if rt_backtrace = {
    .capture            = rt_backtrace_capture,
    .context            = rt_backtrace_context,
    .symbolize          = rt_backtrace_symbolize,
    .trace              = rt_backtrace_trace,
    .trace_self         = rt_backtrace_trace_self,
    .trace_all_but_self = rt_backtrace_trace_all_but_self,
    .string             = rt_backtrace_string,
    .test               = rt_backtrace_test
};

// ______________________________ rt_clipboard.c ______________________________

static errno_t rt_clipboard_put_text(const char* utf8) {
    int32_t chars = rt_str.utf16_chars(utf8, -1);
    int32_t bytes = (chars + 1) * 2;
    uint16_t* utf16 = null;
    errno_t r = rt_heap.alloc((void**)&utf16, (size_t)bytes);
    if (utf16 != null) {
        rt_str.utf8to16(utf16, bytes, utf8, -1);
        rt_assert(utf16[chars - 1] == 0);
        const int32_t n = (int32_t)rt_str.len16(utf16) + 1;
        r = OpenClipboard(GetDesktopWindow()) ? 0 : rt_core.err();
        if (r != 0) { rt_println("OpenClipboard() failed %s", rt_strerr(r)); }
        if (r == 0) {
            r = EmptyClipboard() ? 0 : rt_core.err();
            if (r != 0) { rt_println("EmptyClipboard() failed %s", rt_strerr(r)); }
        }
        void* global = null;
        if (r == 0) {
            global = GlobalAlloc(GMEM_MOVEABLE, (size_t)n * 2);
            r = global != null ? 0 : rt_core.err();
            if (r != 0) { rt_println("GlobalAlloc() failed %s", rt_strerr(r)); }
        }
        if (r == 0) {
            char* d = (char*)GlobalLock(global);
            rt_not_null(d);
            memcpy(d, utf16, (size_t)n * 2);
            r = rt_b2e(SetClipboardData(CF_UNICODETEXT, global));
            GlobalUnlock(global);
            if (r != 0) {
                rt_println("SetClipboardData() failed %s", rt_strerr(r));
                GlobalFree(global);
            } else {
                // do not free global memory. It's owned by system clipboard now
            }
        }
        if (r == 0) {
            r = rt_b2e(CloseClipboard());
            if (r != 0) {
                rt_println("CloseClipboard() failed %s", rt_strerr(r));
            }
        }
        rt_heap.free(utf16);
    }
    return r;
}

static errno_t rt_clipboard_get_text(char* utf8, int32_t* bytes) {
    rt_not_null(bytes);
    errno_t r = rt_b2e(OpenClipboard(GetDesktopWindow()));
    if (r != 0) { rt_println("OpenClipboard() failed %s", rt_strerr(r)); }
    if (r == 0) {
        HANDLE global = GetClipboardData(CF_UNICODETEXT);
        if (global == null) {
            r = rt_core.err();
        } else {
            uint16_t* utf16 = (uint16_t*)GlobalLock(global);
            if (utf16 != null) {
                int32_t utf8_bytes = rt_str.utf8_bytes(utf16, -1);
                if (utf8 != null) {
                    char* decoded = (char*)malloc((size_t)utf8_bytes);
                    if (decoded == null) {
                        r = ERROR_OUTOFMEMORY;
                    } else {
                        rt_str.utf16to8(decoded, utf8_bytes, utf16, -1);
                        int32_t n = rt_min(*bytes, utf8_bytes);
                        memcpy(utf8, decoded, (size_t)n);
                        free(decoded);
                        if (n < utf8_bytes) {
                            r = ERROR_INSUFFICIENT_BUFFER;
                        }
                    }
                }
                *bytes = utf8_bytes;
                GlobalUnlock(global);
            }
        }
        r = rt_b2e(CloseClipboard());
    }
    return r;
}

#ifdef RT_TESTS

static void rt_clipboard_test(void) {
    rt_fatal_if_error(rt_clipboard.put_text("Hello Clipboard"));
    char text[256];
    int32_t bytes = rt_countof(text);
    rt_fatal_if_error(rt_clipboard.get_text(text, &bytes));
    rt_swear(strcmp(text, "Hello Clipboard") == 0);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_clipboard_test(void) {
}

#endif

rt_clipboard_if rt_clipboard = {
    .put_text   = rt_clipboard_put_text,
    .get_text   = rt_clipboard_get_text,
    .put_image  = null, // implemented in ui.app
    .test       = rt_clipboard_test
};

// ________________________________ rt_clock.c ________________________________

enum {
    rt_clock_nsec_in_usec = 1000, // nano in micro
    rt_clock_nsec_in_msec = rt_clock_nsec_in_usec * 1000, // nano in milli
    rt_clock_nsec_in_sec  = rt_clock_nsec_in_msec * 1000,
    rt_clock_usec_in_msec = 1000, // micro in mill
    rt_clock_msec_in_sec  = 1000, // milli in sec
    rt_clock_usec_in_sec  = rt_clock_usec_in_msec * rt_clock_msec_in_sec // micro in sec
};

static uint64_t rt_clock_microseconds_since_epoch(void) { // NOT monotonic
    FILETIME ft; // time in 100ns interval (tenth of microsecond)
    // since 12:00 A.M. January 1, 1601 Coordinated Universal Time (UTC)
    GetSystemTimePreciseAsFileTime(&ft);
    uint64_t microseconds =
        (((uint64_t)ft.dwHighDateTime) << 32 | ft.dwLowDateTime) / 10;
    rt_assert(microseconds > 0);
    return microseconds;
}

static uint64_t rt_clock_localtime(void) {
    TIME_ZONE_INFORMATION tzi; // UTC = local time + bias
    GetTimeZoneInformation(&tzi);
    uint64_t bias = (uint64_t)tzi.Bias * 60LL * 1000 * 1000; // in microseconds
    return rt_clock_microseconds_since_epoch() - bias;
}

static void rt_clock_utc(uint64_t microseconds,
        int32_t* year, int32_t* month, int32_t* day,
        int32_t* hh, int32_t* mm, int32_t* ss, int32_t* ms, int32_t* mc) {
    uint64_t time_in_100ns = microseconds * 10;
    FILETIME mst = { (DWORD)(time_in_100ns & 0xFFFFFFFF),
                     (DWORD)(time_in_100ns >> 32) };
    SYSTEMTIME utc;
    FileTimeToSystemTime(&mst, &utc);
    *year = utc.wYear;
    *month = utc.wMonth;
    *day = utc.wDay;
    *hh = utc.wHour;
    *mm = utc.wMinute;
    *ss = utc.wSecond;
    *ms = utc.wMilliseconds;
    *mc = microseconds % 1000;
}

static void rt_clock_local(uint64_t microseconds,
        int32_t* year, int32_t* month, int32_t* day,
        int32_t* hh, int32_t* mm, int32_t* ss, int32_t* ms, int32_t* mc) {
    uint64_t time_in_100ns = microseconds * 10;
    FILETIME mst = { (DWORD)(time_in_100ns & 0xFFFFFFFF), (DWORD)(time_in_100ns >> 32) };
    SYSTEMTIME utc;
    FileTimeToSystemTime(&mst, &utc);
    DYNAMIC_TIME_ZONE_INFORMATION tzi;
    GetDynamicTimeZoneInformation(&tzi);
    SYSTEMTIME lt = {0};
    SystemTimeToTzSpecificLocalTimeEx(&tzi, &utc, &lt);
    *year = lt.wYear;
    *month = lt.wMonth;
    *day = lt.wDay;
    *hh = lt.wHour;
    *mm = lt.wMinute;
    *ss = lt.wSecond;
    *ms = lt.wMilliseconds;
    *mc = microseconds % 1000;
}

static fp64_t rt_clock_seconds(void) { // since_boot
    LARGE_INTEGER qpc;
    QueryPerformanceCounter(&qpc);
    static fp64_t one_over_freq;
    if (one_over_freq == 0) {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        one_over_freq = 1.0 / (fp64_t)frequency.QuadPart;
    }
    return (fp64_t)qpc.QuadPart * one_over_freq;
}

// Max duration in nanoseconds=2^64 - 1 nanoseconds
//                          2^64 - 1 ns        1 sec          1 min
// Max Duration in Hours =  ----------- x  ------------ x -------------
//                          10^9 ns / s    60 sec / min   60 min / hour
//
//                              1 hour
// Max Duration in Days =  ---------------
//                          24 hours / day
//
// it would take approximately 213,503 days (or about 584.5 years)
// for rt_clock.nanoseconds() to overflow
//
// for divider = rt_num.gcd32(nsec_in_sec, freq) below and 10MHz timer
// the actual duration is shorter because of (mul == 100)
//    (uint64_t)qpc.QuadPart * mul
// 64 bit overflow and is about 5.8 years.
//
// In a long running code like services is advisable to use
// rt_clock.nanoseconds() to measure only deltas and pay close attention
// to the wrap around despite of 5 years monotony

static uint64_t rt_clock_nanoseconds(void) {
    LARGE_INTEGER qpc;
    QueryPerformanceCounter(&qpc);
    static uint32_t freq;
    static uint32_t mul = rt_clock_nsec_in_sec;
    if (freq == 0) {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        rt_assert(frequency.HighPart == 0);
        // even 1GHz frequency should fit into 32 bit unsigned
        rt_assert(frequency.HighPart == 0, "%08lX%%08lX",
               frequency.HighPart, frequency.LowPart);
        // known values: 10,000,000 and 3,000,000 10MHz, 3MHz
        rt_assert(frequency.LowPart % (1000 * 1000) == 0);
        // if we start getting weird frequencies not
        // multiples of MHz rt_num.gcd() approach may need
        // to be revised in favor of rt_num.muldiv64x64()
        freq = frequency.LowPart;
        rt_assert(freq != 0 && freq < (uint32_t)rt_clock.nsec_in_sec);
        // to avoid rt_num.muldiv128:
        uint32_t divider = rt_num.gcd32((uint32_t)rt_clock.nsec_in_sec, freq);
        freq /= divider;
        mul  /= divider;
    }
    uint64_t ns_mul_freq = (uint64_t)qpc.QuadPart * mul;
    return freq == 1 ? ns_mul_freq : ns_mul_freq / freq;
}

// Difference between 1601 and 1970 in microseconds:

static const uint64_t rt_clock_epoch_diff_usec = 11644473600000000ULL;

static uint64_t rt_clock_unix_microseconds(void) {
    return rt_clock.microseconds() - rt_clock_epoch_diff_usec;
}

static uint64_t rt_clock_unix_seconds(void) {
    return rt_clock.unix_microseconds() / (uint64_t)rt_clock.usec_in_sec;
}

static void rt_clock_test(void) {
    #ifdef RT_TESTS
    // TODO: implement more tests
    uint64_t t0 = rt_clock.nanoseconds();
    uint64_t t1 = rt_clock.nanoseconds();
    int32_t count = 0;
    while (t0 == t1 && count < 1024) {
        t1 = rt_clock.nanoseconds();
        count++;
    }
    rt_swear(t0 != t1, "count: %d t0: %lld t1: %lld", count, t0, t1);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
    #endif
}

rt_clock_if rt_clock = {
    .nsec_in_usec      = rt_clock_nsec_in_usec,
    .nsec_in_msec      = rt_clock_nsec_in_msec,
    .nsec_in_sec       = rt_clock_nsec_in_sec,
    .usec_in_msec      = rt_clock_usec_in_msec,
    .msec_in_sec       = rt_clock_msec_in_sec,
    .usec_in_sec       = rt_clock_usec_in_sec,
    .seconds           = rt_clock_seconds,
    .nanoseconds       = rt_clock_nanoseconds,
    .unix_microseconds = rt_clock_unix_microseconds,
    .unix_seconds      = rt_clock_unix_seconds,
    .microseconds      = rt_clock_microseconds_since_epoch,
    .localtime         = rt_clock_localtime,
    .utc               = rt_clock_utc,
    .local             = rt_clock_local,
    .test              = rt_clock_test
};

// _______________________________ rt_config.c ________________________________

// On Unix the implementation should keep KV pairs in
// key-named files inside .name/ folder

static const char* rt_config_apps = "Software\\ui\\apps";

static const DWORD rt_config_access =
    KEY_READ|KEY_WRITE|KEY_SET_VALUE|KEY_QUERY_VALUE|
    KEY_ENUMERATE_SUB_KEYS|DELETE;

static errno_t rt_config_get_reg_key(const char* name, HKEY *key) {
    char path[256] = {0};
    rt_str_printf(path, "%s\\%s", rt_config_apps, name);
    errno_t r = RegOpenKeyExA(HKEY_CURRENT_USER, path, 0, rt_config_access, key);
    if (r != 0) {
        const DWORD option = REG_OPTION_NON_VOLATILE;
        r = RegCreateKeyExA(HKEY_CURRENT_USER, path, 0, null, option,
                            rt_config_access, null, key, null);
    }
    return r;
}

static errno_t rt_config_save(const char* name,
        const char* key, const void* data, int32_t bytes) {
    errno_t r = 0;
    HKEY k = null;
    r = rt_config_get_reg_key(name, &k);
    if (k != null) {
        r = RegSetValueExA(k, key, 0, REG_BINARY,
            (const uint8_t*)data, (DWORD)bytes);
        rt_fatal_if_error(RegCloseKey(k));
    }
    return r;
}

static errno_t rt_config_remove(const char* name, const char* key) {
    errno_t r = 0;
    HKEY k = null;
    r = rt_config_get_reg_key(name, &k);
    if (k != null) {
        r = RegDeleteValueA(k, key);
        rt_fatal_if_error(RegCloseKey(k));
    }
    return r;
}

static errno_t rt_config_clean(const char* name) {
    errno_t r = 0;
    HKEY k = null;
    if (RegOpenKeyExA(HKEY_CURRENT_USER, rt_config_apps,
                                      0, rt_config_access, &k) == 0) {
       r = RegDeleteTreeA(k, name);
       rt_fatal_if_error(RegCloseKey(k));
    }
    return r;
}

static int32_t rt_config_size(const char* name, const char* key) {
    int32_t bytes = -1;
    HKEY k = null;
    errno_t r = rt_config_get_reg_key(name, &k);
    if (k != null) {
        DWORD type = REG_BINARY;
        DWORD cb = 0;
        r = RegQueryValueExA(k, key, null, &type, null, &cb);
        if (r == ERROR_FILE_NOT_FOUND) {
            bytes = 0; // do not report data_size() often used this way
        } else if (r != 0) {
            rt_println("%s.RegQueryValueExA(\"%s\") failed %s",
                name, key, rt_strerr(r));
            bytes = 0; // on any error behave as empty data
        } else {
            bytes = (int32_t)cb;
        }
        rt_fatal_if_error(RegCloseKey(k));
    }
    return bytes;
}

static int32_t rt_config_load(const char* name,
        const char* key, void* data, int32_t bytes) {
    int32_t read = -1;
    HKEY k = null;
    errno_t r = rt_config_get_reg_key(name, &k);
    if (k != null) {
        DWORD type = REG_BINARY;
        DWORD cb = (DWORD)bytes;
        r = RegQueryValueExA(k, key, null, &type, (uint8_t*)data, &cb);
        if (r == ERROR_MORE_DATA) {
            // returns -1 ui_app.data_size() should be used
        } else if (r != 0) {
            if (r != ERROR_FILE_NOT_FOUND) {
                rt_println("%s.RegQueryValueExA(\"%s\") failed %s",
                    name, key, rt_strerr(r));
            }
            read = 0; // on any error behave as empty data
        } else {
            read = (int32_t)cb;
        }
        rt_fatal_if_error(RegCloseKey(k));
    }
    return read;
}

#ifdef RT_TESTS

static void rt_config_test(void) {
    const char* name = strrchr(rt_args.v[0], '\\');
    if (name == null) { name = strrchr(rt_args.v[0], '/'); }
    name = name != null ? name + 1 : rt_args.v[0];
    rt_swear(name != null);
    const char* key = "test";
    const char data[] = "data";
    int32_t bytes = sizeof(data);
    rt_swear(rt_config.save(name, key, data, bytes) == 0);
    char read[256];
    rt_swear(rt_config.load(name, key, read, bytes) == bytes);
    int32_t size = rt_config.size(name, key);
    rt_swear(size == bytes);
    rt_swear(rt_config.remove(name, key) == 0);
    rt_swear(rt_config.clean(name) == 0);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_config_test(void) { }

#endif

rt_config_if rt_config = {
    .save   = rt_config_save,
    .size   = rt_config_size,
    .load   = rt_config_load,
    .remove = rt_config_remove,
    .clean  = rt_config_clean,
    .test   = rt_config_test
};

// ________________________________ rt_core.c _________________________________

// abort does NOT call atexit() functions and
// does NOT flush rt_streams. Also Win32 runtime
// abort() attempt to show Abort/Retry/Ignore
// MessageBox - thus ExitProcess()

static void rt_core_abort(void) { ExitProcess(ERROR_FATAL_APP_EXIT); }

static void rt_core_exit(int32_t exit_code) { exit(exit_code); }

// TODO: consider r = HRESULT_FROM_WIN32() and r = HRESULT_CODE(hr);
// this separates posix error codes from win32 error codes


static errno_t rt_core_err(void) { return (errno_t)GetLastError(); }

static void rt_core_seterr(errno_t err) { SetLastError((DWORD)err); }

rt_static_init(runtime) {
    SetErrorMode(
        // The system does not display the critical-error-handler message box.
        // Instead, the system sends the error to the calling process:
        SEM_FAILCRITICALERRORS|
        // The system automatically fixes memory alignment faults and
        // makes them invisible to the application.
        SEM_NOALIGNMENTFAULTEXCEPT|
        // The system does not display the Windows Error Reporting dialog.
        SEM_NOGPFAULTERRORBOX|
        // The OpenFile function does not display a message box when it fails
        // to find a file. Instead, the error is returned to the caller.
        // This error mode overrides the OF_PROMPT flag.
        SEM_NOOPENFILEERRORBOX);
}

#ifdef RT_TESTS

static void rt_core_test(void) { // in alphabetical order
    rt_args.test();
    rt_atomics.test();
    rt_backtrace.test();
    rt_clipboard.test();
    rt_clock.test();
    rt_config.test();
    rt_debug.test();
    rt_event.test();
    rt_files.test();
    rt_generics.test();
    rt_heap.test();
    rt_loader.test();
    rt_mem.test();
    rt_mutex.test();
    rt_num.test();
    rt_processes.test();
    rt_static_init_test();
    rt_str.test();
    rt_streams.test();
    rt_thread.test();
    rt_vigil.test();
    rt_worker.test();
}

#else

static void rt_core_test(void) { }

#endif

rt_core_if rt_core = {
    .err     = rt_core_err,
    .set_err = rt_core_seterr,
    .abort   = rt_core_abort,
    .exit    = rt_core_exit,
    .test    = rt_core_test,
    .error   = {                                              // posix
        .access_denied          = ERROR_ACCESS_DENIED,        // EACCES
        .bad_file               = ERROR_BAD_FILE_TYPE,        // EBADF
        .broken_pipe            = ERROR_BROKEN_PIPE,          // EPIPE
        .device_not_ready       = ERROR_NOT_READY,            // ENXIO
        .directory_not_empty    = ERROR_DIR_NOT_EMPTY,        // ENOTEMPTY
        .disk_full              = ERROR_DISK_FULL,            // ENOSPC
        .file_exists            = ERROR_FILE_EXISTS,          // EEXIST
        .file_not_found         = ERROR_FILE_NOT_FOUND,       // ENOENT
        .insufficient_buffer    = ERROR_INSUFFICIENT_BUFFER,  // E2BIG
        .interrupted            = ERROR_OPERATION_ABORTED,    // EINTR
        .invalid_data           = ERROR_INVALID_DATA,         // EINVAL
        .invalid_handle         = ERROR_INVALID_HANDLE,       // EBADF
        .invalid_parameter      = ERROR_INVALID_PARAMETER,    // EINVAL
        .io_error               = ERROR_IO_DEVICE,            // EIO
        .more_data              = ERROR_MORE_DATA,            // ENOBUFS
        .name_too_long          = ERROR_FILENAME_EXCED_RANGE, // ENAMETOOLONG
        .no_child_process       = ERROR_NO_PROC_SLOTS,        // ECHILD
        .not_a_directory        = ERROR_DIRECTORY,            // ENOTDIR
        .not_empty              = ERROR_DIR_NOT_EMPTY,        // ENOTEMPTY
        .out_of_memory          = ERROR_OUTOFMEMORY,          // ENOMEM
        .path_not_found         = ERROR_PATH_NOT_FOUND,       // ENOENT
        .pipe_not_connected     = ERROR_PIPE_NOT_CONNECTED,   // EPIPE
        .read_only_file         = ERROR_WRITE_PROTECT,        // EROFS
        .resource_deadlock      = ERROR_LOCK_VIOLATION,       // EDEADLK
        .too_many_open_files    = ERROR_TOO_MANY_OPEN_FILES,  // EMFILE
    }
};

#pragma comment(lib, "advapi32")
#pragma comment(lib, "ntdll")
#pragma comment(lib, "psapi")
#pragma comment(lib, "shell32")
#pragma comment(lib, "shlwapi")
#pragma comment(lib, "kernel32")
#pragma comment(lib, "user32") // clipboard
#pragma comment(lib, "imm32")  // Internationalization input method
#pragma comment(lib, "ole32")  // rt_files.known_folder CoMemFree
#pragma comment(lib, "dbghelp")
#pragma comment(lib, "imagehlp")



// ________________________________ rt_debug.c ________________________________

static const char* rt_debug_abbreviate(const char* file) {
    const char* fn = strrchr(file, '\\');
    if (fn == null) { fn = strrchr(file, '/'); }
    return fn != null ? fn + 1 : file;
}

#ifdef WINDOWS

static int32_t rt_debug_max_file_line;
static int32_t rt_debug_max_function;

static void rt_debug_output(const char* s, int32_t count) {
    bool intercepted = false;
    if (rt_debug.tee != null) { intercepted = rt_debug.tee(s, count); }
    if (!intercepted) {
        // For link.exe /Subsystem:Windows code stdout/stderr are often closed
        if (stderr != null && fileno(stderr) >= 0) {
            fprintf(stderr, "%s", s);
        }
        // SetConsoleCP(CP_UTF8) is not guaranteed to be called
        uint16_t* wide = rt_stackalloc((count + 1) * sizeof(uint16_t));
        rt_str.utf8to16(wide, count, s, -1);
        OutputDebugStringW(wide);
    }
}

static void rt_debug_println_va(const char* file, int32_t line, const char* func,
        const char* format, va_list va) {
    if (func == null) { func = ""; }
    char file_line[1024];
    if (line == 0 && file == null || file[0] == 0x00) {
        file_line[0] = 0x00;
    } else {
        if (file == null) { file = ""; } // backtrace can have null files
        // full path is useful in MSVC debugger output pane (clickable)
        // for all other scenarios short filename without path is preferable:
        const char* name = IsDebuggerPresent() ? file : rt_files.basename(file);
        snprintf(file_line, rt_countof(file_line) - 1, "%s(%d):", name, line);
    }
    file_line[rt_countof(file_line) - 1] = 0x00; // always zero terminated'
    rt_debug_max_file_line = rt_max(rt_debug_max_file_line,
                                    (int32_t)strlen(file_line));
    rt_debug_max_function  = rt_max(rt_debug_max_function,
                                    (int32_t)strlen(func));
    char prefix[2 * 1024];
    // snprintf() does not guarantee zero termination on truncation
    snprintf(prefix, rt_countof(prefix) - 1, "%-*s %-*s",
            rt_debug_max_file_line, file_line,
            rt_debug_max_function,  func);
    prefix[rt_countof(prefix) - 1] = 0; // zero terminated
    char text[2 * 1024];
    if (format != null && format[0] != 0) {
        #if defined(__GNUC__) || defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wformat-nonliteral"
        #endif
        vsnprintf(text, rt_countof(text) - 1, format, va);
        text[rt_countof(text) - 1] = 0;
        #if defined(__GNUC__) || defined(__clang__)
        #pragma GCC diagnostic pop
        #endif
    } else {
        text[0] = 0;
    }
    char output[4 * 1024];
    snprintf(output, rt_countof(output) - 1, "%s %s", prefix, text);
    output[rt_countof(output) - 2] = 0;
    // strip trailing \n which can be remnant of fprintf("...\n")
    int32_t n = (int32_t)strlen(output);
    while (n > 0 && (output[n - 1] == '\n' || output[n - 1] == '\r')) {
        output[n - 1] = 0;
        n--;
    }
    rt_assert(n + 1 < rt_countof(output));
    // Win32 OutputDebugString() needs \n
    output[n + 0] = '\n';
    output[n + 1] = 0;
    rt_debug.output(output, n + 2); // including 0x00
}

#else // posix version:

static void rt_debug_vprintf(const char* file, int32_t line, const char* func,
        const char* format, va_list va) {
    fprintf(stderr, "%s(%d): %s ", file, line, func);
    vfprintf(stderr, format, va);
    fprintf(stderr, "\n");
}

#endif

static void rt_debug_perrno(const char* file, int32_t line,
    const char* func, int32_t err_no, const char* format, ...) {
    if (err_no != 0) {
        if (format != null && format[0] != 0) {
            va_list va;
            va_start(va, format);
            rt_debug.println_va(file, line, func, format, va);
            va_end(va);
        }
        rt_debug.println(file, line, func, "errno: %d %s", err_no, strerror(err_no));
    }
}

static void rt_debug_perror(const char* file, int32_t line,
    const char* func, int32_t error, const char* format, ...) {
    if (error != 0) {
        if (format != null && format[0] != 0) {
            va_list va;
            va_start(va, format);
            rt_debug.println_va(file, line, func, format, va);
            va_end(va);
        }
        rt_debug.println(file, line, func, "error: %s", rt_strerr(error));
    }
}

static void rt_debug_println(const char* file, int32_t line, const char* func,
        const char* format, ...) {
    va_list va;
    va_start(va, format);
    rt_debug.println_va(file, line, func, format, va);
    va_end(va);
}

static bool rt_debug_is_debugger_present(void) { return IsDebuggerPresent(); }

static void rt_debug_breakpoint(void) {
    if (rt_debug.is_debugger_present()) { DebugBreak(); }
}

static errno_t rt_debug_raise(uint32_t exception) {
    rt_core.set_err(0);
    RaiseException(exception, EXCEPTION_NONCONTINUABLE, 0, null);
    return rt_core.err();
}

static int32_t rt_debug_verbosity_from_string(const char* s) {
    char* n = null;
    long v = strtol(s, &n, 10);
    if (stricmp(s, "quiet") == 0) {
        return rt_debug.verbosity.quiet;
    } else if (stricmp(s, "info") == 0) {
        return rt_debug.verbosity.info;
    } else if (stricmp(s, "verbose") == 0) {
        return rt_debug.verbosity.verbose;
    } else if (stricmp(s, "debug") == 0) {
        return rt_debug.verbosity.debug;
    } else if (stricmp(s, "trace") == 0) {
        return rt_debug.verbosity.trace;
    } else if (n > s && rt_debug.verbosity.quiet <= v &&
               v <= rt_debug.verbosity.trace) {
        return v;
    } else {
        rt_fatal("invalid verbosity: %s", s);
        return rt_debug.verbosity.quiet;
    }
}

static void rt_debug_test(void) {
    #ifdef RT_TESTS
    // not clear what can be tested here
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
    #endif
}

#ifndef STATUS_POSSIBLE_DEADLOCK
#define STATUS_POSSIBLE_DEADLOCK 0xC0000194uL
#endif

rt_debug_if rt_debug = {
    .verbosity = {
        .level   =  0,
        .quiet   =  0,
        .info    =  1,
        .verbose =  2,
        .debug   =  3,
        .trace   =  4,
    },
    .verbosity_from_string = rt_debug_verbosity_from_string,
    .tee                   = null,
    .output                = rt_debug_output,
    .println               = rt_debug_println,
    .println_va            = rt_debug_println_va,
    .perrno                = rt_debug_perrno,
    .perror                = rt_debug_perror,
    .is_debugger_present   = rt_debug_is_debugger_present,
    .breakpoint            = rt_debug_breakpoint,
    .raise                 = rt_debug_raise,
    .exception             = {
        .access_violation        = EXCEPTION_ACCESS_VIOLATION,
        .datatype_misalignment   = EXCEPTION_DATATYPE_MISALIGNMENT,
        .breakpoint              = EXCEPTION_BREAKPOINT,
        .single_step             = EXCEPTION_SINGLE_STEP,
        .array_bounds            = EXCEPTION_ARRAY_BOUNDS_EXCEEDED,
        .float_denormal_operand  = EXCEPTION_FLT_DENORMAL_OPERAND,
        .float_divide_by_zero    = EXCEPTION_FLT_DIVIDE_BY_ZERO,
        .float_inexact_result    = EXCEPTION_FLT_INEXACT_RESULT,
        .float_invalid_operation = EXCEPTION_FLT_INVALID_OPERATION,
        .float_overflow          = EXCEPTION_FLT_OVERFLOW,
        .float_stack_check       = EXCEPTION_FLT_STACK_CHECK,
        .float_underflow         = EXCEPTION_FLT_UNDERFLOW,
        .int_divide_by_zero      = EXCEPTION_INT_DIVIDE_BY_ZERO,
        .int_overflow            = EXCEPTION_INT_OVERFLOW,
        .priv_instruction        = EXCEPTION_PRIV_INSTRUCTION,
        .in_page_error           = EXCEPTION_IN_PAGE_ERROR,
        .illegal_instruction     = EXCEPTION_ILLEGAL_INSTRUCTION,
        .noncontinuable          = EXCEPTION_NONCONTINUABLE_EXCEPTION,
        .stack_overflow          = EXCEPTION_STACK_OVERFLOW,
        .invalid_disposition     = EXCEPTION_INVALID_DISPOSITION,
        .guard_page              = EXCEPTION_GUARD_PAGE,
        .invalid_handle          = EXCEPTION_INVALID_HANDLE,
        .possible_deadlock       = EXCEPTION_POSSIBLE_DEADLOCK
    },
    .test                  = rt_debug_test
};

// ________________________________ rt_files.c ________________________________

// TODO: test FILE_APPEND_DATA
// https://learn.microsoft.com/en-us/windows/win32/fileio/appending-one-file-to-another-file?redirectedfrom=MSDN

// are posix and Win32 seek in agreement?
rt_static_assertion(SEEK_SET == FILE_BEGIN);
rt_static_assertion(SEEK_CUR == FILE_CURRENT);
rt_static_assertion(SEEK_END == FILE_END);

#ifndef O_SYNC
#define O_SYNC (0x10000)
#endif

static errno_t rt_files_open(rt_file_t* *file, const char* fn, int32_t f) {
    DWORD access = (f & rt_files.o_wr) ? GENERIC_WRITE :
                   (f & rt_files.o_rw) ? GENERIC_READ | GENERIC_WRITE :
                                      GENERIC_READ;
    access |= (f & rt_files.o_append) ? FILE_APPEND_DATA : 0;
    DWORD disposition =
        (f & rt_files.o_create) ? ((f & rt_files.o_excl)  ? CREATE_NEW :
                                (f & rt_files.o_trunc) ? CREATE_ALWAYS :
                                                      OPEN_ALWAYS) :
            (f & rt_files.o_trunc) ? TRUNCATE_EXISTING : OPEN_EXISTING;
    const DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
    DWORD attr = FILE_ATTRIBUTE_NORMAL;
    attr |= (f & O_SYNC) ? FILE_FLAG_WRITE_THROUGH : 0;
    *file = CreateFileA(fn, access, share, null, disposition, attr, null);
    return *file != INVALID_HANDLE_VALUE ? 0 : rt_core.err();
}

static bool rt_files_is_valid(rt_file_t* file) { // both null and rt_files.invalid
    return file != rt_files.invalid && file != null;
}

static errno_t rt_files_seek(rt_file_t* file, int64_t *position, int32_t method) {
    LARGE_INTEGER distance_to_move = { .QuadPart = *position };
    LARGE_INTEGER p = { 0 }; // pointer
    errno_t r = rt_b2e(SetFilePointerEx(file, distance_to_move, &p, (DWORD)method));
    if (r == 0) { *position = p.QuadPart; }
    return r;
}

static inline uint64_t rt_files_ft_to_us(FILETIME ft) { // us (microseconds)
    return (ft.dwLowDateTime | (((uint64_t)ft.dwHighDateTime) << 32)) / 10;
}

static int64_t rt_files_a2t(DWORD a) {
    int64_t type = 0;
    if (a & FILE_ATTRIBUTE_REPARSE_POINT) {
        type |= rt_files.type_symlink;
    }
    if (a & FILE_ATTRIBUTE_DIRECTORY) {
        type |= rt_files.type_folder;
    }
    if (a & FILE_ATTRIBUTE_DEVICE) {
        type |= rt_files.type_device;
    }
    return type;
}

#ifdef FILES_LINUX_PATH_BY_FD

static int get_final_path_name_by_fd(int fd, char *buffer, int32_t bytes) {
    swear(bytes >= 0);
    char fd_path[16 * 1024];
    // /proc/self/fd/* is a symbolic link
    snprintf(fd_path, sizeof(fd_path), "/proc/self/fd/%d", fd);
    size_t len = readlink(fd_path, buffer, bytes - 1);
    if (len != -1) { buffer[len] = 0x00; } // Null-terminate the result
    return len == -1 ? errno : 0;
}

#endif

static errno_t rt_files_stat(rt_file_t* file, rt_files_stat_t* s,
                             bool follow_symlink) {
    errno_t r = 0;
    BY_HANDLE_FILE_INFORMATION fi;
    rt_fatal_win32err(GetFileInformationByHandle(file, &fi));
    const bool symlink =
        (fi.dwFileAttributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0;
    if (follow_symlink && symlink) {
        const DWORD flags = FILE_NAME_NORMALIZED | VOLUME_NAME_DOS;
        DWORD n = GetFinalPathNameByHandleA(file, null, 0, flags);
        if (n == 0) {
            r = rt_core.err();
        } else {
            char* name = null;
            r = rt_heap.allocate(null, (void**)&name, (int64_t)n + 2, false);
            if (r == 0) {
                n = GetFinalPathNameByHandleA(file, name, n + 1, flags);
                if (n == 0) {
                    r = rt_core.err();
                } else {
                    rt_file_t* f = rt_files.invalid;
                    r = rt_files.open(&f, name, rt_files.o_rd);
                    if (r == 0) { // keep following:
                        r = rt_files.stat(f, s, follow_symlink);
                        rt_files.close(f);
                    }
                }
                rt_heap.deallocate(null, name);
            }
        }
    } else {
        s->size = (int64_t)((uint64_t)fi.nFileSizeLow |
                          (((uint64_t)fi.nFileSizeHigh) << 32));
        s->created  = rt_files_ft_to_us(fi.ftCreationTime); // since epoch
        s->accessed = rt_files_ft_to_us(fi.ftLastAccessTime);
        s->updated  = rt_files_ft_to_us(fi.ftLastWriteTime);
        s->type = rt_files_a2t(fi.dwFileAttributes);
    }
    return r;
}

static errno_t rt_files_read(rt_file_t* file, void* data, int64_t bytes, int64_t *transferred) {
    errno_t r = 0;
    *transferred = 0;
    while (bytes > 0 && r == 0) {
        DWORD chunk_size = (DWORD)(bytes > UINT32_MAX ? UINT32_MAX : bytes);
        DWORD bytes_read = 0;
        r = rt_b2e(ReadFile(file, data, chunk_size, &bytes_read, null));
        if (r == 0) {
            *transferred += bytes_read;
            bytes -= bytes_read;
            data = (uint8_t*)data + bytes_read;
        }
    }
    return r;
}

static errno_t rt_files_write(rt_file_t* file, const void* data, int64_t bytes, int64_t *transferred) {
    errno_t r = 0;
    *transferred = 0;
    while (bytes > 0 && r == 0) {
        DWORD chunk_size = (DWORD)(bytes > UINT32_MAX ? UINT32_MAX : bytes);
        DWORD bytes_read = 0;
        r = rt_b2e(WriteFile(file, data, chunk_size, &bytes_read, null));
        if (r == 0) {
            *transferred += bytes_read;
            bytes -= bytes_read;
            data = (const uint8_t*)data + bytes_read;
        }
    }
    return r;
}

static errno_t rt_files_flush(rt_file_t* file) {
    return rt_b2e(FlushFileBuffers(file));
}

static void rt_files_close(rt_file_t* file) {
    rt_win32_close_handle(file);
}

static errno_t rt_files_write_fully(const char* filename, const void* data,
                                 int64_t bytes, int64_t *transferred) {
    if (transferred != null) { *transferred = 0; }
    errno_t r = 0;
    const DWORD access = GENERIC_WRITE;
    const DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
    const DWORD flags = FILE_ATTRIBUTE_NORMAL | FILE_FLAG_WRITE_THROUGH;
    HANDLE file = CreateFileA(filename, access, share, null, CREATE_ALWAYS,
                              flags, null);
    if (file == INVALID_HANDLE_VALUE) {
        r = rt_core.err();
    } else {
        int64_t written = 0;
        const uint8_t* p = (const uint8_t*)data;
        while (r == 0 && bytes > 0) {
            uint64_t write = bytes >= UINT32_MAX ?
                (uint64_t)(UINT32_MAX) - 0xFFFFuLL : (uint64_t)bytes;
            rt_assert(0 < write && write < (uint64_t)UINT32_MAX);
            DWORD chunk = 0;
            r = rt_b2e(WriteFile(file, p, (DWORD)write, &chunk, null));
            written += chunk;
            bytes -= chunk;
        }
        if (transferred != null) { *transferred = written; }
        errno_t rc = rt_b2e(FlushFileBuffers(file));
        if (r == 0) { r = rc; }
        rt_win32_close_handle(file);
    }
    return r;
}

static errno_t rt_files_unlink(const char* pathname) {
    if (rt_files.is_folder(pathname)) {
        return rt_b2e(RemoveDirectoryA(pathname));
    } else {
        return rt_b2e(DeleteFileA(pathname));
    }
}

static errno_t rt_files_create_tmp(char* fn, int32_t count) {
    // create temporary file (not folder!) see folders_test() about racing
    rt_swear(fn != null && count > 0);
    const char* tmp = rt_files.tmp();
    errno_t r = 0;
    if (count < (int32_t)strlen(tmp) + 8) {
        r = ERROR_BUFFER_OVERFLOW;
    } else {
        rt_assert(count > (int32_t)strlen(tmp) + 8);
        // If GetTempFileNameA() succeeds, the return value is the length,
        // in chars, of the string copied to lpBuffer, not including the
        // terminating null character.If the function fails,
        // the return value is zero.
        if (count > (int32_t)strlen(tmp) + 8) {
            char prefix[4] = { 0 };
            r = GetTempFileNameA(tmp, prefix, 0, fn) == 0 ? rt_core.err() : 0;
            if (r == 0) {
                rt_assert(rt_files.exists(fn) && !rt_files.is_folder(fn));
            } else {
                rt_println("GetTempFileNameA() failed %s", rt_strerr(r));
            }
        } else {
            r = ERROR_BUFFER_OVERFLOW;
        }
    }
    return r;
}

#pragma push_macro("files_acl_args")
#pragma push_macro("files_get_acl")
#pragma push_macro("files_set_acl")

#define rt_files_acl_args(acl) DACL_SECURITY_INFORMATION, null, null, acl, null

#define rt_files_get_acl(obj, type, acl, sd) (errno_t)(         \
    (type == SE_FILE_OBJECT ? GetNamedSecurityInfoA((char*)obj, \
             SE_FILE_OBJECT, rt_files_acl_args(acl), &sd) :     \
    (type == SE_KERNEL_OBJECT) ? GetSecurityInfo((HANDLE)obj,   \
             SE_KERNEL_OBJECT, rt_files_acl_args(acl), &sd) :   \
    ERROR_INVALID_PARAMETER))

#define rt_files_set_acl(obj, type, acl) (errno_t)(             \
    (type == SE_FILE_OBJECT ? SetNamedSecurityInfoA((char*)obj, \
             SE_FILE_OBJECT, rt_files_acl_args(acl)) :          \
    (type == SE_KERNEL_OBJECT) ? SetSecurityInfo((HANDLE)obj,   \
             SE_KERNEL_OBJECT, rt_files_acl_args(acl)) :        \
    ERROR_INVALID_PARAMETER))

static errno_t rt_files_acl_add_ace(ACL* acl, SID* sid, uint32_t mask,
                                 ACL** free_me, byte flags) {
    ACL_SIZE_INFORMATION info = {0};
    ACL* bigger = null;
    uint32_t bytes_needed = sizeof(ACCESS_ALLOWED_ACE) + GetLengthSid(sid)
                          - sizeof(DWORD);
    errno_t r = rt_b2e(GetAclInformation(acl, &info, sizeof(ACL_SIZE_INFORMATION),
        AclSizeInformation));
    if (r == 0 && info.AclBytesFree < bytes_needed) {
        const int64_t bytes = (int64_t)(info.AclBytesInUse + bytes_needed);
        r = rt_heap.allocate(null, (void**)&bigger, bytes, true);
        if (r == 0) {
            r = rt_b2e(InitializeAcl((ACL*)bigger,
                    info.AclBytesInUse + bytes_needed, ACL_REVISION));
        }
    }
    if (r == 0 && bigger != null) {
        for (int32_t i = 0; i < (int32_t)info.AceCount; i++) {
            ACCESS_ALLOWED_ACE* ace = null;
            r = rt_b2e(GetAce(acl, (DWORD)i, (void**)&ace));
            if (r != 0) { break; }
            r = rt_b2e(AddAce(bigger, ACL_REVISION, MAXDWORD, ace,
                           ace->Header.AceSize));
            if (r != 0) { break; }
        }
    }
    if (r == 0) {
        ACCESS_ALLOWED_ACE* ace = null;
        r = rt_heap.allocate(null, (void**)&ace, bytes_needed, true);
        if (r == 0) {
            ace->Header.AceFlags = flags;
            ace->Header.AceType = ACCESS_ALLOWED_ACE_TYPE;
            ace->Header.AceSize = (WORD)bytes_needed;
            ace->Mask = mask;
            ace->SidStart = sizeof(ACCESS_ALLOWED_ACE);
            memcpy(&ace->SidStart, sid, GetLengthSid(sid));
            r = rt_b2e(AddAce(bigger != null ? bigger : acl, ACL_REVISION, MAXDWORD,
                           ace, bytes_needed));
            rt_heap.deallocate(null, ace);
        }
    }
    *free_me = bigger;
    return r;
}

static errno_t rt_files_lookup_sid(ACCESS_ALLOWED_ACE* ace) {
    // handy for debugging
    SID* sid = (SID*)&ace->SidStart;
    DWORD l1 = 128, l2 = 128;
    char account[128];
    char group[128];
    SID_NAME_USE use;
    errno_t r = rt_b2e(LookupAccountSidA(null, sid, account,
                                     &l1, group, &l2, &use));
    if (r == 0) {
        rt_println("%s/%s: type: %d, mask: 0x%X, flags:%d",
                group, account,
                ace->Header.AceType, ace->Mask, ace->Header.AceFlags);
    } else {
        rt_println("LookupAccountSidA() failed %s", rt_strerr(r));
    }
    return r;
}

static errno_t rt_files_add_acl_ace(void* obj, int32_t obj_type,
                                 int32_t sid_type, uint32_t mask) {
    uint8_t stack[SECURITY_MAX_SID_SIZE] = {0};
    DWORD n = rt_countof(stack);
    SID* sid = (SID*)stack;
    errno_t r = rt_b2e(CreateWellKnownSid((WELL_KNOWN_SID_TYPE)sid_type,
                                       null, sid, &n));
    if (r != 0) {
        return ERROR_INVALID_PARAMETER;
    }
    ACL* acl = null;
    void* sd = null;
    r = rt_files_get_acl(obj, obj_type, &acl, sd);
    if (r == 0) {
        ACCESS_ALLOWED_ACE* found = null;
        for (int32_t i = 0; i < acl->AceCount; i++) {
            ACCESS_ALLOWED_ACE* ace = null;
            r = rt_b2e(GetAce(acl, (DWORD)i, (void**)&ace));
            if (r != 0) { break; }
            if (EqualSid((SID*)&ace->SidStart, sid)) {
                if (ace->Header.AceType == ACCESS_ALLOWED_ACE_TYPE &&
                   (ace->Header.AceFlags & INHERITED_ACE) == 0) {
                    found = ace;
                } else if (ace->Header.AceType !=
                           ACCESS_ALLOWED_ACE_TYPE) {
                    rt_println("%d ACE_TYPE is not supported.",
                             ace->Header.AceType);
                    r = ERROR_INVALID_PARAMETER;
                }
                break;
            }
        }
        if (r == 0 && found) {
            if ((found->Mask & mask) != mask) {
//              rt_println("updating existing ace");
                found->Mask |= mask;
                r = rt_files_set_acl(obj, obj_type, acl);
            } else {
//              rt_println("desired access is already allowed by ace");
            }
        } else if (r == 0) {
//          rt_println("inserting new ace");
            ACL* new_acl = null;
            byte flags = obj_type == SE_FILE_OBJECT ?
                CONTAINER_INHERIT_ACE | OBJECT_INHERIT_ACE : 0;
            r = rt_files_acl_add_ace(acl, sid, mask, &new_acl, flags);
            if (r == 0) {
                r = rt_files_set_acl(obj, obj_type, (new_acl != null ? new_acl : acl));
            }
            if (new_acl != null) { rt_heap.deallocate(null, new_acl); }
        }
    }
    if (sd != null) { LocalFree(sd); }
    return r;
}

#pragma pop_macro("files_set_acl")
#pragma pop_macro("files_get_acl")
#pragma pop_macro("files_acl_args")

static errno_t rt_files_chmod777(const char* pathname) {
    SID_IDENTIFIER_AUTHORITY SIDAuthWorld = SECURITY_WORLD_SID_AUTHORITY;
    PSID everyone = null; // Create a well-known SID for the Everyone group.
    rt_fatal_win32err(AllocateAndInitializeSid(&SIDAuthWorld, 1,
             SECURITY_WORLD_RID, 0, 0, 0, 0, 0, 0, 0, &everyone));
    EXPLICIT_ACCESSA ea[1] = { { 0 } };
    // Initialize an EXPLICIT_ACCESS structure for an ACE.
    ea[0].grfAccessPermissions = 0xFFFFFFFF;
    ea[0].grfAccessMode  = GRANT_ACCESS; // The ACE will allow everyone all access.
    ea[0].grfInheritance = SUB_CONTAINERS_AND_OBJECTS_INHERIT;
    ea[0].Trustee.TrusteeForm = TRUSTEE_IS_SID;
    ea[0].Trustee.TrusteeType = TRUSTEE_IS_WELL_KNOWN_GROUP;
    ea[0].Trustee.ptstrName  = (LPSTR)everyone;
    // Create a new ACL that contains the new ACEs.
    ACL* acl = null;
    rt_fatal_if_error(SetEntriesInAclA(1, ea, null, &acl));
    // Initialize a security descriptor.
    uint8_t stack[SECURITY_DESCRIPTOR_MIN_LENGTH] = {0};
    SECURITY_DESCRIPTOR* sd = (SECURITY_DESCRIPTOR*)stack;
    rt_fatal_win32err(InitializeSecurityDescriptor(sd,
        SECURITY_DESCRIPTOR_REVISION));
    // Add the ACL to the security descriptor.
    rt_fatal_win32err(SetSecurityDescriptorDacl(sd,
        /* present flag: */ true, acl, /* not a default DACL: */  false));
    // Change the security attributes
    errno_t r = rt_b2e(SetFileSecurityA(pathname, DACL_SECURITY_INFORMATION, sd));
    if (r != 0) {
        rt_println("chmod777(%s) failed %s", pathname, rt_strerr(r));
    }
    if (everyone != null) { FreeSid(everyone); }
    if (acl != null) { LocalFree(acl); }
    return r;
}

// https://docs.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createdirectorya
// "If lpSecurityAttributes is null, the directory gets a default security
//  descriptor. The ACLs in the default security descriptor for a directory
//  are inherited from its parent directory."

static errno_t rt_files_mkdirs(const char* dir) {
    const int32_t n = (int32_t)strlen(dir) + 1;
    char* s = null;
    errno_t r = rt_heap.allocate(null, (void**)&s, n, true);
    const char* next = strchr(dir, '\\');
    if (next == null) { next = strchr(dir, '/'); }
    while (r == 0 && next != null) {
        if (next > dir && *(next - 1) != ':') {
            memcpy(s, dir, (size_t)(next - dir));
            r = rt_b2e(CreateDirectoryA(s, null));
            if (r == ERROR_ALREADY_EXISTS) { r = 0; }
        }
        if (r == 0) {
            const char* prev = ++next;
            next = strchr(prev, '\\');
            if (next == null) { next = strchr(prev, '/'); }
        }
    }
    if (r == 0) {
        r = rt_b2e(CreateDirectoryA(dir, null));
    }
    rt_heap.deallocate(null, s);
    return r == ERROR_ALREADY_EXISTS ? 0 : r;
}

#pragma push_macro("rt_files_realloc_path")
#pragma push_macro("rt_files_append_name")

#define rt_files_realloc_path(r, pn, pnc, fn, name) do {                \
    const int32_t bytes = (int32_t)(strlen(fn) + strlen(name) + 3);     \
    if (bytes > pnc) {                                                  \
        r = rt_heap.reallocate(null, (void**)&pn, bytes, false);        \
        if (r != 0) {                                                   \
            pnc = bytes;                                                \
        } else {                                                        \
            rt_heap.deallocate(null, pn);                               \
            pn = null;                                                  \
        }                                                               \
    }                                                                   \
} while (0)

#define rt_files_append_name(pn, pnc, fn, name) do {     \
    if (strcmp(fn, "\\") == 0 || strcmp(fn, "/") == 0) { \
        rt_str.format(pn, pnc, "\\%s", name);            \
    } else {                                             \
        rt_str.format(pn, pnc, "%.*s\\%s", k, fn, name); \
    }                                                    \
} while (0)

static errno_t rt_files_rmdirs(const char* fn) {
    rt_files_stat_t st;
    rt_folder_t folder;
    errno_t r = rt_files.opendir(&folder, fn);
    if (r == 0) {
        int32_t k = (int32_t)strlen(fn);
        // remove trailing backslash (except if it is root: "/" or "\\")
        if (k > 1 && (fn[k - 1] == '/' || fn[k - 1] == '\\')) {
            k--;
        }
        int32_t pnc = 64 * 1024; // pathname "pn" capacity in bytes
        char* pn = null;
        r = rt_heap.allocate(null, (void**)&pn, pnc, false);
        while (r == 0) {
            // recurse into sub folders and remove them first
            // do NOT follow symlinks - it could be disastrous
            const char* name = rt_files.readdir(&folder, &st);
            if (name == null) { break; }
            if (strcmp(name, ".") != 0 && strcmp(name, "..") != 0 &&
                (st.type & rt_files.type_symlink) == 0 &&
                (st.type & rt_files.type_folder) != 0) {
                rt_files_realloc_path(r, pn, pnc, fn, name);
                if (r == 0) {
                    rt_files_append_name(pn, pnc, fn, name);
                    r = rt_files.rmdirs(pn);
                }
            }
        }
        rt_files.closedir(&folder);
        r = rt_files.opendir(&folder, fn);
        while (r == 0) {
            const char* name = rt_files.readdir(&folder, &st);
            if (name == null) { break; }
            // symlinks are already removed as normal files
            if (strcmp(name, ".") != 0 && strcmp(name, "..") != 0 &&
                (st.type & rt_files.type_folder) == 0) {
                rt_files_realloc_path(r, pn, pnc, fn, name);
                if (r == 0) {
                    rt_files_append_name(pn, pnc, fn, name);
                    r = rt_files.unlink(pn);
                    if (r != 0) {
                        rt_println("remove(%s) failed %s", pn, rt_strerr(r));
                    }
                }
            }
        }
        rt_heap.deallocate(null, pn);
        rt_files.closedir(&folder);
    }
    if (r == 0) { r = rt_files.unlink(fn); }
    return r;
}

#pragma pop_macro("rt_files_append_name")
#pragma pop_macro("rt_files_realloc_path")

static bool rt_files_exists(const char* path) {
    return PathFileExistsA(path);
}

static bool rt_files_is_folder(const char* path) {
    return PathIsDirectoryA(path);
}

static bool rt_files_is_symlink(const char* filename) {
    DWORD attributes = GetFileAttributesA(filename);
    return attributes != INVALID_FILE_ATTRIBUTES &&
          (attributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0;
}

static const char* rt_files_basename(const char* pathname) {
    const char* bn = strrchr(pathname, '\\');
    if (bn == null) { bn = strrchr(pathname, '/'); }
    return bn != null ? bn + 1 : pathname;
}

static errno_t rt_files_copy(const char* s, const char* d) {
    return rt_b2e(CopyFileA(s, d, false));
}

static errno_t rt_files_move(const char* s, const char* d) {
    static const DWORD flags =
        MOVEFILE_REPLACE_EXISTING |
        MOVEFILE_COPY_ALLOWED |
        MOVEFILE_WRITE_THROUGH;
    return rt_b2e(MoveFileExA(s, d, flags));
}

static errno_t rt_files_link(const char* from, const char* to) {
    // note reverse order of parameters:
    return rt_b2e(CreateHardLinkA(to, from, null));
}

static errno_t rt_files_symlink(const char* from, const char* to) {
    // The correct order of parameters for CreateSymbolicLinkA is:
    // CreateSymbolicLinkA(symlink_to_create, existing_file, flags);
    DWORD flags = rt_files.is_folder(from) ? SYMBOLIC_LINK_FLAG_DIRECTORY : 0;
    return rt_b2e(CreateSymbolicLinkA(to, from, flags));
}

static const char* rt_files_known_folder(int32_t kf) {
    // known folder ids order must match enum see:
    static const GUID* kf_ids[] = {
        &FOLDERID_Profile,
        &FOLDERID_Desktop,
        &FOLDERID_Documents,
        &FOLDERID_Downloads,
        &FOLDERID_Music,
        &FOLDERID_Pictures,
        &FOLDERID_Videos,
        &FOLDERID_Public,
        &FOLDERID_ProgramFiles,
        &FOLDERID_ProgramData
    };
    static rt_file_name_t known_folders[rt_countof(kf_ids)];
    rt_fatal_if(!(0 <= kf && kf < rt_countof(kf_ids)), "invalid kf=%d", kf);
    if (known_folders[kf].s[0] == 0) {
        uint16_t* path = null;
        rt_fatal_if_error(SHGetKnownFolderPath(kf_ids[kf], 0, null, &path));
        const int32_t n = rt_countof(known_folders[kf].s);
        rt_str.utf16to8(known_folders[kf].s, n, path, -1);
        CoTaskMemFree(path);
	}
    return known_folders[kf].s;
}

static const char* rt_files_bin(void) {
    return rt_files_known_folder(rt_files.folder.bin);
}

static const char* rt_files_data(void) {
    return rt_files_known_folder(rt_files.folder.data);
}

static const char* rt_files_tmp(void) {
    static char tmp[rt_files_max_path];
    if (tmp[0] == 0) {
        // If GetTempPathA() succeeds, the return value is the length,
        // in chars, of the string copied to lpBuffer, not including
        // the terminating null character. If the function fails, the
        // return value is zero.
        errno_t r = GetTempPathA(rt_countof(tmp), tmp) == 0 ? rt_core.err() : 0;
        rt_fatal_if(r != 0, "GetTempPathA() failed %s", rt_strerr(r));
    }
    return tmp;
}

static errno_t rt_files_cwd(char* fn, int32_t count) {
    rt_swear(count > 1);
    DWORD bytes = (DWORD)(count - 1);
    errno_t r = rt_b2e(GetCurrentDirectoryA(bytes, fn));
    fn[count - 1] = 0; // always
    return r;
}

static errno_t rt_files_chdir(const char* fn) {
    return rt_b2e(SetCurrentDirectoryA(fn));
}

typedef struct rt_files_dir_s {
    HANDLE handle;
    WIN32_FIND_DATAA find; // On Win64: 320 bytes
} rt_files_dir_t;

rt_static_assertion(sizeof(rt_files_dir_t) <= sizeof(rt_folder_t));

static errno_t rt_files_opendir(rt_folder_t* folder, const char* folder_name) {
    rt_files_dir_t* d = (rt_files_dir_t*)(void*)folder;
    int32_t n = (int32_t)strlen(folder_name);
    char* fn = null;
    // extra room for "\*" suffix
    errno_t r = rt_heap.allocate(null, (void**)&fn, (int64_t)n + 3, false);
    if (r == 0) {
        rt_str.format(fn, n + 3, "%s\\*", folder_name);
        fn[n + 2] = 0;
        d->handle = FindFirstFileA(fn, &d->find);
        if (d->handle == INVALID_HANDLE_VALUE) { r = rt_core.err(); }
        rt_heap.deallocate(null, fn);
    }
    return r;
}

static uint64_t rt_files_ft2us(FILETIME* ft) { // 100ns units to microseconds:
    return (((uint64_t)ft->dwHighDateTime) << 32 | ft->dwLowDateTime) / 10;
}

static const char* rt_files_readdir(rt_folder_t* folder, rt_files_stat_t* s) {
    const char* fn = null;
    rt_files_dir_t* d = (rt_files_dir_t*)(void*)folder;
    if (FindNextFileA(d->handle, &d->find)) {
        fn = d->find.cFileName;
        // Ensure zero termination
        d->find.cFileName[rt_countof(d->find.cFileName) - 1] = 0x00;
        if (s != null) {
            s->accessed = rt_files_ft2us(&d->find.ftLastAccessTime);
            s->created = rt_files_ft2us(&d->find.ftCreationTime);
            s->updated = rt_files_ft2us(&d->find.ftLastWriteTime);
            s->type = rt_files_a2t(d->find.dwFileAttributes);
            s->size = (int64_t)((((uint64_t)d->find.nFileSizeHigh) << 32) |
                                  (uint64_t)d->find.nFileSizeLow);
        }
    }
    return fn;
}

static void rt_files_closedir(rt_folder_t* folder) {
    rt_files_dir_t* d = (rt_files_dir_t*)(void*)folder;
    rt_fatal_win32err(FindClose(d->handle));
}

#pragma push_macro("files_test_failed")

#ifdef RT_TESTS

// TODO: change rt_fatal_if() to swear()

#define rt_files_test_failed " failed %s", rt_strerr(rt_core.err())

#pragma push_macro("verbose") // --verbosity trace

#define verbose(...) do {                                       \
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) { \
        rt_println(__VA_ARGS__);                                   \
    }                                                           \
} while (0)

static void folders_dump_time(const char* label, uint64_t us) {
    int32_t year = 0;
    int32_t month = 0;
    int32_t day = 0;
    int32_t hh = 0;
    int32_t mm = 0;
    int32_t ss = 0;
    int32_t ms = 0;
    int32_t mc = 0;
    rt_clock.local(us, &year, &month, &day, &hh, &mm, &ss, &ms, &mc);
    rt_println("%-7s: %04d-%02d-%02d %02d:%02d:%02d.%03d:%03d",
            label, year, month, day, hh, mm, ss, ms, mc);
}

static void folders_test(void) {
    uint64_t now = rt_clock.microseconds(); // microseconds since epoch
    uint64_t before = now - 1 * (uint64_t)rt_clock.usec_in_sec; // one second earlier
    uint64_t after  = now + 2 * (uint64_t)rt_clock.usec_in_sec; // two seconds later
    int32_t year = 0;
    int32_t month = 0;
    int32_t day = 0;
    int32_t hh = 0;
    int32_t mm = 0;
    int32_t ss = 0;
    int32_t ms = 0;
    int32_t mc = 0;
    rt_clock.local(now, &year, &month, &day, &hh, &mm, &ss, &ms, &mc);
    verbose("now: %04d-%02d-%02d %02d:%02d:%02d.%03d:%03d",
             year, month, day, hh, mm, ss, ms, mc);
    // Test cwd, setcwd
    const char* tmp = rt_files.tmp();
    char cwd[256] = { 0 };
    rt_fatal_if(rt_files.cwd(cwd, sizeof(cwd)) != 0, "rt_files.cwd() failed");
    rt_fatal_if(rt_files.chdir(tmp) != 0, "rt_files.chdir(\"%s\") failed %s",
                tmp, rt_strerr(rt_core.err()));
    // there is no racing free way to create temporary folder
    // without having a temporary file for the duration of folder usage:
    char tmp_file[rt_files_max_path]; // create_tmp() is thread safe race free:
    errno_t r = rt_files.create_tmp(tmp_file, rt_countof(tmp_file));
    rt_fatal_if(r != 0, "rt_files.create_tmp() failed %s", rt_strerr(r));
    char tmp_dir[rt_files_max_path];
    rt_str_printf(tmp_dir, "%s.dir", tmp_file);
    r = rt_files.mkdirs(tmp_dir);
    rt_fatal_if(r != 0, "rt_files.mkdirs(%s) failed %s", tmp_dir, rt_strerr(r));
    verbose("%s", tmp_dir);
    rt_folder_t folder;
    char pn[rt_files_max_path] = { 0 };
    rt_str_printf(pn, "%s/file", tmp_dir);
    // cannot test symlinks because they are only
    // available to Administrators and in Developer mode
//  char sym[rt_files_max_path] = { 0 };
    char hard[rt_files_max_path] = { 0 };
    char sub[rt_files_max_path] = { 0 };
    rt_str_printf(hard, "%s/hard", tmp_dir);
    rt_str_printf(sub, "%s/subd", tmp_dir);
    const char* content = "content";
    int64_t transferred = 0;
    r = rt_files.write_fully(pn, content, (int64_t)strlen(content), &transferred);
    rt_fatal_if(r != 0, "rt_files.write_fully(\"%s\") failed %s", pn, rt_strerr(r));
    rt_swear(transferred == (int64_t)strlen(content));
    r = rt_files.link(pn, hard);
    rt_fatal_if(r != 0, "rt_files.link(\"%s\", \"%s\") failed %s",
                      pn, hard, rt_strerr(r));
    r = rt_files.mkdirs(sub);
    rt_fatal_if(r != 0, "rt_files.mkdirs(\"%s\") failed %s", sub, rt_strerr(r));
    r = rt_files.opendir(&folder, tmp_dir);
    rt_fatal_if(r != 0, "rt_files.opendir(\"%s\") failed %s", tmp_dir, rt_strerr(r));
    for (;;) {
        rt_files_stat_t st = { 0 };
        const char* name = rt_files.readdir(&folder, &st);
        if (name == null) { break; }
        uint64_t at = st.accessed;
        uint64_t ct = st.created;
        uint64_t ut = st.updated;
        rt_swear(ct <= at && ct <= ut);
        rt_clock.local(ct, &year, &month, &day, &hh, &mm, &ss, &ms, &mc);
        bool is_folder = st.type & rt_files.type_folder;
        bool is_symlink = st.type & rt_files.type_symlink;
        int64_t bytes = st.size;
        verbose("%s: %04d-%02d-%02d %02d:%02d:%02d.%03d:%03d %lld bytes %s%s",
                name, year, month, day, hh, mm, ss, ms, mc,
                bytes, is_folder ? "[folder]" : "", is_symlink ? "[symlink]" : "");
        if (strcmp(name, "file") == 0 || strcmp(name, "hard") == 0) {
            rt_swear(bytes == (int64_t)strlen(content),
                    "size of \"%s\": %lld is incorrect expected: %d",
                    name, bytes, transferred);
        }
        if (strcmp(name, ".") == 0 || strcmp(name, "..") == 0) {
            rt_swear(is_folder, "\"%s\" is_folder: %d", name, is_folder);
        } else {
            rt_swear((strcmp(name, "subd") == 0) == is_folder,
                  "\"%s\" is_folder: %d", name, is_folder);
            // empirically timestamps are imprecise on NTFS
            rt_swear(at >= before, "access: %lld  >= %lld", at, before);
            if (ct < before || ut < before || at >= after || ct >= after || ut >= after) {
                rt_println("file: %s", name);
                folders_dump_time("before", before);
                folders_dump_time("create", ct);
                folders_dump_time("update", ut);
                folders_dump_time("access", at);
            }
            rt_swear(ct >= before, "create: %lld  >= %lld", ct, before);
            rt_swear(ut >= before, "update: %lld  >= %lld", ut, before);
            // and no later than 2 seconds since folders_test()
            rt_swear(at < after, "access: %lld  < %lld", at, after);
            rt_swear(ct < after, "create: %lld  < %lld", ct, after);
            rt_swear(at < after, "update: %lld  < %lld", ut, after);
        }
    }
    rt_files.closedir(&folder);
    r = rt_files.rmdirs(tmp_dir);
    rt_fatal_if(r != 0, "rt_files.rmdirs(\"%s\") failed %s",
                     tmp_dir, rt_strerr(r));
    r = rt_files.unlink(tmp_file);
    rt_fatal_if(r != 0, "rt_files.unlink(\"%s\") failed %s",
                     tmp_file, rt_strerr(r));
    rt_fatal_if(rt_files.chdir(cwd) != 0, "rt_files.chdir(\"%s\") failed %s",
             cwd, rt_strerr(rt_core.err()));
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#pragma pop_macro("verbose")

static void rt_files_test_append_thread(void* p) {
    rt_file_t* f = (rt_file_t*)p;
    uint8_t data[256] = {0};
    for (int i = 0; i < 256; i++) { data[i] = (uint8_t)i; }
    int64_t transferred = 0;
    rt_fatal_if(rt_files.write(f, data, rt_countof(data), &transferred) != 0 ||
             transferred != rt_countof(data), "rt_files.write()" rt_files_test_failed);
}

static void rt_files_test(void) {
    folders_test();
    uint64_t now = rt_clock.microseconds(); // epoch time
    char tf[256]; // temporary file
    rt_fatal_if(rt_files.create_tmp(tf, rt_countof(tf)) != 0,
            "rt_files.create_tmp()" rt_files_test_failed);
    uint8_t data[256] = {0};
    int64_t transferred = 0;
    for (int i = 0; i < 256; i++) { data[i] = (uint8_t)i; }
    {
        rt_file_t* f = rt_files.invalid;
        rt_fatal_if(rt_files.open(&f, tf,
                 rt_files.o_wr | rt_files.o_create | rt_files.o_trunc) != 0 ||
                !rt_files.is_valid(f), "rt_files.open()" rt_files_test_failed);
        rt_fatal_if(rt_files.write_fully(tf, data, rt_countof(data), &transferred) != 0 ||
                 transferred != rt_countof(data),
                "rt_files.write_fully()" rt_files_test_failed);
        rt_fatal_if(rt_files.open(&f, tf, rt_files.o_rd) != 0 ||
                !rt_files.is_valid(f), "rt_files.open()" rt_files_test_failed);
        for (int32_t i = 0; i < 256; i++) {
            for (int32_t j = 1; j < 256 - i; j++) {
                uint8_t test[rt_countof(data)] = { 0 };
                int64_t position = i;
                rt_fatal_if(rt_files.seek(f, &position, rt_files.seek_set) != 0 ||
                         position != i,
                        "rt_files.seek(position: %lld) failed %s",
                         position, rt_strerr(rt_core.err()));
                rt_fatal_if(rt_files.read(f, test, j, &transferred) != 0 ||
                         transferred != j,
                        "rt_files.read() transferred: %lld failed %s",
                        transferred, rt_strerr(rt_core.err()));
                for (int32_t k = 0; k < j; k++) {
                    rt_swear(test[k] == data[i + k],
                         "Data mismatch at position: %d, length %d"
                         "test[%d]: 0x%02X != data[%d + %d]: 0x%02X ",
                          i, j,
                          k, test[k], i, k, data[i + k]);
                }
            }
        }
        rt_swear((rt_files.o_rd | rt_files.o_wr) != rt_files.o_rw);
        rt_fatal_if(rt_files.open(&f, tf, rt_files.o_rw) != 0 || !rt_files.is_valid(f),
                "rt_files.open()" rt_files_test_failed);
        for (int32_t i = 0; i < 256; i++) {
            uint8_t val = ~data[i];
            int64_t pos = i;
            rt_fatal_if(rt_files.seek(f, &pos, rt_files.seek_set) != 0 || pos != i,
                    "rt_files.seek() failed %s", rt_core.err());
            rt_fatal_if(rt_files.write(f, &val, 1, &transferred) != 0 ||
                     transferred != 1, "rt_files.write()" rt_files_test_failed);
            pos = i;
            rt_fatal_if(rt_files.seek(f, &pos, rt_files.seek_set) != 0 || pos != i,
                    "rt_files.seek(pos: %lld i: %d) failed %s", pos, i, rt_core.err());
            uint8_t read_val = 0;
            rt_fatal_if(rt_files.read(f, &read_val, 1, &transferred) != 0 ||
                     transferred != 1, "rt_files.read()" rt_files_test_failed);
            rt_swear(read_val == val, "Data mismatch at position %d", i);
        }
        rt_files_stat_t s = { 0 };
        rt_files.stat(f, &s, false);
        uint64_t before = now - 1 * (uint64_t)rt_clock.usec_in_sec; // one second before now
        uint64_t after  = now + 2 * (uint64_t)rt_clock.usec_in_sec; // two seconds after
        rt_swear(before <= s.created  && s.created  <= after,
             "before: %lld created: %lld after: %lld", before, s.created, after);
        rt_swear(before <= s.accessed && s.accessed <= after,
             "before: %lld created: %lld accessed: %lld", before, s.accessed, after);
        rt_swear(before <= s.updated  && s.updated  <= after,
             "before: %lld created: %lld updated: %lld", before, s.updated, after);
        rt_files.close(f);
        rt_fatal_if(rt_files.open(&f, tf, rt_files.o_wr | rt_files.o_create | rt_files.o_trunc) != 0 ||
                !rt_files.is_valid(f), "rt_files.open()" rt_files_test_failed);
        rt_files.stat(f, &s, false);
        rt_swear(s.size == 0, "File is not empty after truncation. .size: %lld", s.size);
        rt_files.close(f);
    }
    {  // Append test with threads
        rt_file_t* f = rt_files.invalid;
        rt_fatal_if(rt_files.open(&f, tf, rt_files.o_rw | rt_files.o_append) != 0 ||
                !rt_files.is_valid(f), "rt_files.open()" rt_files_test_failed);
        rt_thread_t thread1 = rt_thread.start(rt_files_test_append_thread, f);
        rt_thread_t thread2 = rt_thread.start(rt_files_test_append_thread, f);
        rt_thread.join(thread1, -1);
        rt_thread.join(thread2, -1);
        rt_files.close(f);
    }
    {   // write_fully, exists, is_folder, mkdirs, rmdirs, create_tmp, chmod777
        rt_fatal_if(rt_files.write_fully(tf, data, rt_countof(data), &transferred) != 0 ||
                 transferred != rt_countof(data),
                "rt_files.write_fully() failed %s", rt_core.err());
        rt_fatal_if(!rt_files.exists(tf), "file \"%s\" does not exist", tf);
        rt_fatal_if(rt_files.is_folder(tf), "%s is a folder", tf);
        rt_fatal_if(rt_files.chmod777(tf) != 0, "rt_files.chmod777(\"%s\") failed %s",
                 tf, rt_strerr(rt_core.err()));
        char folder[256] = { 0 };
        rt_str_printf(folder, "%s.folder\\subfolder", tf);
        rt_fatal_if(rt_files.mkdirs(folder) != 0, "rt_files.mkdirs(\"%s\") failed %s",
            folder, rt_strerr(rt_core.err()));
        rt_fatal_if(!rt_files.is_folder(folder), "\"%s\" is not a folder", folder);
        rt_fatal_if(rt_files.chmod777(folder) != 0, "rt_files.chmod777(\"%s\") failed %s",
                 folder, rt_strerr(rt_core.err()));
        rt_fatal_if(rt_files.rmdirs(folder) != 0, "rt_files.rmdirs(\"%s\") failed %s",
                 folder, rt_strerr(rt_core.err()));
        rt_fatal_if(rt_files.exists(folder), "folder \"%s\" still exists", folder);
    }
    {   // getcwd, chdir
        const char* tmp = rt_files.tmp();
        char cwd[256] = { 0 };
        rt_fatal_if(rt_files.cwd(cwd, sizeof(cwd)) != 0, "rt_files.cwd() failed");
        rt_fatal_if(rt_files.chdir(tmp) != 0, "rt_files.chdir(\"%s\") failed %s",
                 tmp, rt_strerr(rt_core.err()));
        // symlink
        if (rt_processes.is_elevated()) {
            char sym_link[rt_files_max_path];
            rt_str_printf(sym_link, "%s.sym_link", tf);
            rt_fatal_if(rt_files.symlink(tf, sym_link) != 0,
                "rt_files.symlink(\"%s\", \"%s\") failed %s",
                tf, sym_link, rt_strerr(rt_core.err()));
            rt_fatal_if(!rt_files.is_symlink(sym_link), "\"%s\" is not a sym_link", sym_link);
            rt_fatal_if(rt_files.unlink(sym_link) != 0, "rt_files.unlink(\"%s\") failed %s",
                    sym_link, rt_strerr(rt_core.err()));
        } else {
            rt_println("Skipping rt_files.symlink test: process is not elevated");
        }
        // hard link
        char hard_link[rt_files_max_path];
        rt_str_printf(hard_link, "%s.hard_link", tf);
        rt_fatal_if(rt_files.link(tf, hard_link) != 0,
            "rt_files.link(\"%s\", \"%s\") failed %s",
            tf, hard_link, rt_strerr(rt_core.err()));
        rt_fatal_if(!rt_files.exists(hard_link), "\"%s\" does not exist", hard_link);
        rt_fatal_if(rt_files.unlink(hard_link) != 0, "rt_files.unlink(\"%s\") failed %s",
                 hard_link, rt_strerr(rt_core.err()));
        rt_fatal_if(rt_files.exists(hard_link), "\"%s\" still exists", hard_link);
        // copy, move:
        rt_fatal_if(rt_files.copy(tf, "copied_file") != 0,
            "rt_files.copy(\"%s\", 'copied_file') failed %s",
            tf, rt_strerr(rt_core.err()));
        rt_fatal_if(!rt_files.exists("copied_file"), "'copied_file' does not exist");
        rt_fatal_if(rt_files.move("copied_file", "moved_file") != 0,
            "rt_files.move('copied_file', 'moved_file') failed %s",
            rt_strerr(rt_core.err()));
        rt_fatal_if(rt_files.exists("copied_file"), "'copied_file' still exists");
        rt_fatal_if(!rt_files.exists("moved_file"), "'moved_file' does not exist");
        rt_fatal_if(rt_files.unlink("moved_file") != 0,
                "rt_files.unlink('moved_file') failed %s",
                 rt_strerr(rt_core.err()));
        rt_fatal_if(rt_files.chdir(cwd) != 0, "rt_files.chdir(\"%s\") failed %s",
                    cwd, rt_strerr(rt_core.err()));
    }
    rt_fatal_if(rt_files.unlink(tf) != 0);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_files_test(void) {}

#endif // RT_TESTS

#pragma pop_macro("files_test_failed")

rt_files_if rt_files = {
    .invalid  = (rt_file_t*)INVALID_HANDLE_VALUE,
    // rt_files_stat_t.type:
    .type_folder  = 0x00000010, // FILE_ATTRIBUTE_DIRECTORY
    .type_symlink = 0x00000400, // FILE_ATTRIBUTE_REPARSE_POINT
    .type_device  = 0x00000040, // FILE_ATTRIBUTE_DEVICE
    // seek() methods:
    .seek_set = SEEK_SET,
    .seek_cur = SEEK_CUR,
    .seek_end = SEEK_END,
    // open() flags: missing O_RSYNC, O_DSYNC, O_NONBLOCK, O_NOCTTY
    .o_rd     = O_RDONLY,
    .o_wr     = O_WRONLY,
    .o_rw     = O_RDWR,
    .o_append = O_APPEND,
    .o_create = O_CREAT,
    .o_excl   = O_EXCL,
    .o_trunc  = O_TRUNC,
    .o_sync   = O_SYNC,
    // known folders ids:
    .folder = {
        .home      = 0, // c:\Users\<username>
        .desktop   = 1,
        .documents = 2,
        .downloads = 3,
        .music     = 4,
        .pictures  = 5,
        .videos    = 6,
        .shared    = 7, // c:\Users\Public
        .bin       = 8, // c:\Program Files
        .data      = 9  // c:\ProgramData
    },
    // methods:
    .open         = rt_files_open,
    .is_valid     = rt_files_is_valid,
    .seek         = rt_files_seek,
    .stat         = rt_files_stat,
    .read         = rt_files_read,
    .write        = rt_files_write,
    .flush        = rt_files_flush,
    .close        = rt_files_close,
    .write_fully  = rt_files_write_fully,
    .exists       = rt_files_exists,
    .is_folder    = rt_files_is_folder,
    .is_symlink   = rt_files_is_symlink,
    .mkdirs       = rt_files_mkdirs,
    .rmdirs       = rt_files_rmdirs,
    .create_tmp   = rt_files_create_tmp,
    .chmod777     = rt_files_chmod777,
    .unlink       = rt_files_unlink,
    .link         = rt_files_link,
    .symlink      = rt_files_symlink,
    .basename     = rt_files_basename,
    .copy         = rt_files_copy,
    .move         = rt_files_move,
    .cwd          = rt_files_cwd,
    .chdir        = rt_files_chdir,
    .known_folder = rt_files_known_folder,
    .bin          = rt_files_bin,
    .data         = rt_files_data,
    .tmp          = rt_files_tmp,
    .opendir      = rt_files_opendir,
    .readdir      = rt_files_readdir,
    .closedir     = rt_files_closedir,
    .test         = rt_files_test
};

// ______________________________ rt_generics.c _______________________________

#ifdef RT_TESTS

static void rt_generics_test(void) {
    {
        int8_t a = 10, b = 20;
        rt_swear(rt_max(a++, b++) == 20);
        rt_swear(rt_min(a++, b++) == 11);
    }
    {
        int32_t a = 10, b = 20;
        rt_swear(rt_max(a++, b++) == 20);
        rt_swear(rt_min(a++, b++) == 11);
    }
    {
        fp32_t a = 1.1f, b = 2.2f;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        fp64_t a = 1.1, b = 2.2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        fp32_t a = 1.1f, b = 2.2f;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        fp64_t a = 1.1, b = 2.2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        char a = 1, b = 2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        unsigned char a = 1, b = 2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    // MS cl.exe version 19.39.33523 has issues with "long":
    // does not pick up int32_t/uint32_t types for "long" and "unsigned long"
    {
        long int a = 1, b = 2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        unsigned long a = 1, b = 2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    {
        long long a = 1, b = 2;
        rt_swear(rt_max(a, b) == b);
        rt_swear(rt_min(a, b) == a);
    }
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_generics_test(void) { }

#endif

rt_generics_if rt_generics = {
    .test = rt_generics_test
};


// ________________________________ rt_heap.c _________________________________

static errno_t rt_heap_alloc(void* *a, int64_t bytes) {
    return rt_heap.allocate(null, a, bytes, false);
}

static errno_t rt_heap_alloc_zero(void* *a, int64_t bytes) {
    return rt_heap.allocate(null, a, bytes, true);
}

static errno_t rt_heap_realloc(void* *a, int64_t bytes) {
    return rt_heap.reallocate(null, a, bytes, false);
}

static errno_t rt_heap_realloc_zero(void* *a, int64_t bytes) {
    return rt_heap.reallocate(null, a, bytes, true);
}

static void rt_heap_free(void* a) {
    rt_heap.deallocate(null, a);
}

static rt_heap_t* rt_heap_create(bool serialized) {
    const DWORD options = serialized ? 0 : HEAP_NO_SERIALIZE;
    return (rt_heap_t*)HeapCreate(options, 0, 0);
}

static void rt_heap_dispose(rt_heap_t* h) {
    rt_fatal_win32err(HeapDestroy((HANDLE)h));
}

static inline HANDLE rt_heap_or_process_heap(rt_heap_t* h) {
    static HANDLE process_heap;
    if (process_heap == null) { process_heap = GetProcessHeap(); }
    return h != null ? (HANDLE)h : process_heap;
}

static errno_t rt_heap_allocate(rt_heap_t* h, void* *p, int64_t bytes, bool zero) {
    rt_swear(bytes > 0);
    #ifdef DEBUG
        static bool enabled;
        if (!enabled) {
            enabled = true;
            HeapSetInformation(null, HeapEnableTerminationOnCorruption, null, 0);
        }
    #endif
    const DWORD flags = zero ? HEAP_ZERO_MEMORY : 0;
    *p = HeapAlloc(rt_heap_or_process_heap(h), flags, (SIZE_T)bytes);
    return *p == null ? ERROR_OUTOFMEMORY : 0;
}

static errno_t rt_heap_reallocate(rt_heap_t* h, void* *p, int64_t bytes,
        bool zero) {
    rt_swear(bytes > 0);
    const DWORD flags = zero ? HEAP_ZERO_MEMORY : 0;
    void* a = *p == null ? // HeapReAlloc(..., null, bytes) may not work
        HeapAlloc(rt_heap_or_process_heap(h), flags, (SIZE_T)bytes) :
        HeapReAlloc(rt_heap_or_process_heap(h), flags, *p, (SIZE_T)bytes);
    if (a != null) { *p = a; }
    return a == null ? ERROR_OUTOFMEMORY : 0;
}

static void rt_heap_deallocate(rt_heap_t* h, void* a) {
    rt_fatal_win32err(HeapFree(rt_heap_or_process_heap(h), 0, a));
}

static int64_t rt_heap_bytes(rt_heap_t* h, void* a) {
    SIZE_T bytes = HeapSize(rt_heap_or_process_heap(h), 0, a);
    rt_fatal_if(bytes == (SIZE_T)-1);
    return (int64_t)bytes;
}

#ifdef RT_TESTS

static void rt_heap_test(void) {
    // TODO: allocate, reallocate deallocate, create, dispose
    void*   a[1024]; // addresses
    int32_t b[1024]; // bytes
    uint32_t seed = 0x1;
    for (int i = 0; i < 1024; i++) {
        b[i] = (int32_t)(rt_num.random32(&seed) % 1024) + 1;
        errno_t r = rt_heap.alloc(&a[i], b[i]);
        rt_swear(r == 0);
    }
    for (int i = 0; i < 1024; i++) {
        rt_heap.free(a[i]);
    }
    HeapCompact(rt_heap_or_process_heap(null), 0);
    // "There is no extended error information for HeapValidate;
    //  do not call GetLastError."
    rt_swear(HeapValidate(rt_heap_or_process_heap(null), 0, null));
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_heap_test(void) { }

#endif

rt_heap_if rt_heap = {
    .alloc        = rt_heap_alloc,
    .alloc_zero   = rt_heap_alloc_zero,
    .realloc      = rt_heap_realloc,
    .realloc_zero = rt_heap_realloc_zero,
    .free         = rt_heap_free,
    .create       = rt_heap_create,
    .allocate     = rt_heap_allocate,
    .reallocate   = rt_heap_reallocate,
    .deallocate   = rt_heap_deallocate,
    .bytes        = rt_heap_bytes,
    .dispose      = rt_heap_dispose,
    .test         = rt_heap_test
};

// _______________________________ rt_loader.c ________________________________

// This is oversimplified Win32 version completely ignoring mode.

// I bit more Posix compliant version is here:
// https://github.com/dlfcn-win32/dlfcn-win32/blob/master/src/dlfcn.c
// POSIX says that if the value of file is NULL, a handle on a global
// symbol object must be provided. That object must be able to access
// all symbols from the original program file, and any objects loaded
// with the RTLD_GLOBAL flag.
// The return value from GetModuleHandle( ) allows us to retrieve
// symbols only from the original program file. EnumProcessModules() is
// used to access symbols from other libraries. For objects loaded
// with the RTLD_LOCAL flag, we create our own list later on. They are
// excluded from EnumProcessModules() iteration.

static void* rt_loader_all;

static void* rt_loader_sym_all(const char* name) {
    void* sym = null;
    DWORD bytes = 0;
    rt_fatal_win32err(EnumProcessModules(GetCurrentProcess(),
                                         null, 0, &bytes));
    rt_assert(bytes % sizeof(HMODULE) == 0);
    rt_assert(bytes / sizeof(HMODULE) < 1024); // OK to allocate 8KB on stack
    HMODULE* modules = null;
    rt_fatal_if_error(rt_heap.allocate(null, (void**)&modules, bytes, false));
    rt_fatal_win32err(EnumProcessModules(GetCurrentProcess(),
                                         modules, bytes, &bytes));
    const int32_t n = bytes / (int32_t)sizeof(HMODULE);
    for (int32_t i = 0; i < n && sym != null; i++) {
        sym = rt_loader.sym(modules[i], name);
    }
    if (sym == null) {
        sym = rt_loader.sym(GetModuleHandleA(null), name);
    }
    rt_heap.deallocate(null, modules);
    return sym;
}

static void* rt_loader_open(const char* filename, int32_t rt_unused(mode)) {
    return filename == null ? &rt_loader_all : (void*)LoadLibraryA(filename);
}

static void* rt_loader_sym(void* handle, const char* name) {
    return handle == &rt_loader_all ?
            (void*)rt_loader_sym_all(name) :
            (void*)GetProcAddress((HMODULE)handle, name);
}

static void rt_loader_close(void* handle) {
    if (handle != &rt_loader_all) {
        rt_fatal_win32err(FreeLibrary(handle));
    }
}

#ifdef RT_TESTS

// manually test exported function once and comment out because of
// creating .lib out of each .exe is annoying

#undef RT_LOADER_TEST_EXPORTED_FUNCTION

#ifdef RT_LOADER_TEST_EXPORTED_FUNCTION


static int32_t rt_loader_test_calls_count;

rt_export void rt_loader_test_exported_function(void);

void rt_loader_test_exported_function(void) { rt_loader_test_calls_count++; }

#endif

static void rt_loader_test(void) {
    void* global = rt_loader.open(null, rt_loader.local);
    rt_loader.close(global);
    // NtQueryTimerResolution - http://undocumented.ntinternals.net/
    typedef long (__stdcall *query_timer_resolution_t)(
        long* minimum_resolution,
        long* maximum_resolution,
        long* current_resolution);
    void* nt_dll = rt_loader.open("ntdll", rt_loader.local);
    query_timer_resolution_t query_timer_resolution =
        (query_timer_resolution_t)rt_loader.sym(nt_dll, "NtQueryTimerResolution");
    // in 100ns = 0.1us units
    long min_resolution = 0;
    long max_resolution = 0; // lowest possible delay between timer events
    long cur_resolution = 0;
    rt_fatal_if(query_timer_resolution(
        &min_resolution, &max_resolution, &cur_resolution) != 0);
//  if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
//      rt_println("timer resolution min: %.3f max: %.3f cur: %.3f millisecond",
//          min_resolution / 10.0 / 1000.0,
//          max_resolution / 10.0 / 1000.0,
//          cur_resolution / 10.0 / 1000.0);
//      // Interesting observation cur_resolution sometimes 15.625ms or 1.0ms
//  }
    rt_loader.close(nt_dll);
#ifdef RT_LOADER_TEST_EXPORTED_FUNCTION
    rt_loader_test_calls_count = 0;
    rt_loader_test_exported_function(); // to make sure it is linked in
    rt_swear(rt_loader_test_calls_count == 1);
    typedef void (*foo_t)(void);
    foo_t foo = (foo_t)rt_loader.sym(global, "rt_loader_test_exported_function");
    foo();
    rt_swear(rt_loader_test_calls_count == 2);
#endif
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_loader_test(void) {}

#endif

enum {
    rt_loader_local  = 0,       // RTLD_LOCAL  All symbols are not made available for relocation processing by other modules.
    rt_loader_lazy   = 1,       // RTLD_LAZY   Relocations are performed at an implementation-dependent time.
    rt_loader_now    = 2,       // RTLD_NOW    Relocations are performed when the object is loaded.
    rt_loader_global = 0x00100, // RTLD_GLOBAL All symbols are available for relocation processing of other modules.
};

rt_loader_if rt_loader = {
    .local  = rt_loader_local,
    .lazy   = rt_loader_lazy,
    .now    = rt_loader_now,
    .global = rt_loader_global,
    .open   = rt_loader_open,
    .sym    = rt_loader_sym,
    .close  = rt_loader_close,
    .test   = rt_loader_test
};

// _________________________________ rt_mem.c _________________________________

static errno_t rt_mem_map_view_of_file(HANDLE file,
        void* *data, int64_t *bytes, bool rw) {
    errno_t r = 0;
    void* address = null;
    HANDLE mapping = CreateFileMapping(file, null,
        rw ? PAGE_READWRITE : PAGE_READONLY,
        (uint32_t)(*bytes >> 32), (uint32_t)*bytes, null);
    if (mapping == null) {
        r = rt_core.err();
    } else {
        DWORD access = rw ? FILE_MAP_ALL_ACCESS : FILE_MAP_READ;
        address = MapViewOfFile(mapping, access, 0, 0, (SIZE_T)*bytes);
        if (address == null) { r = rt_core.err(); }
        rt_win32_close_handle(mapping);
    }
    if (r == 0) {
        *data = address;
    } else {
        *data = null;
        *bytes = 0;
    }
    return r;
}

// see: https://learn.microsoft.com/en-us/windows/win32/secauthz/enabling-and-disabling-privileges-in-c--

static errno_t rt_mem_set_token_privilege(void* token,
            const char* name, bool e) {
    TOKEN_PRIVILEGES tp = { .PrivilegeCount = 1 };
    tp.Privileges[0].Attributes = e ? SE_PRIVILEGE_ENABLED : 0;
    rt_fatal_win32err(LookupPrivilegeValueA(null, name, &tp.Privileges[0].Luid));
    return rt_b2e(AdjustTokenPrivileges(token, false, &tp,
               sizeof(TOKEN_PRIVILEGES), null, null));
}

static errno_t rt_mem_adjust_process_privilege_manage_volume_name(void) {
    // see: https://devblogs.microsoft.com/oldnewthing/20160603-00/?p=93565
    const uint32_t access = TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY;
    const HANDLE process = GetCurrentProcess();
    HANDLE token = null;
    errno_t r = rt_b2e(OpenProcessToken(process, access, &token));
    if (r == 0) {
        const char* se_manage_volume_name = "SeManageVolumePrivilege";
        r = rt_mem_set_token_privilege(token, se_manage_volume_name, true);
        rt_win32_close_handle(token);
    }
    return r;
}

static errno_t rt_mem_map_file(const char* filename, void* *data,
        int64_t *bytes, bool rw) {
    if (rw) { // for SetFileValidData() call:
        (void)rt_mem_adjust_process_privilege_manage_volume_name();
    }
    errno_t r = 0;
    const DWORD access = GENERIC_READ | (rw ? GENERIC_WRITE : 0);
    const DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
    const DWORD disposition = rw ? OPEN_ALWAYS : OPEN_EXISTING;
    const DWORD flags = FILE_ATTRIBUTE_NORMAL;
    HANDLE file = CreateFileA(filename, access, share, null, disposition,
                              flags, null);
    if (file == INVALID_HANDLE_VALUE) {
        r = rt_core.err();
    } else {
        LARGE_INTEGER eof = { .QuadPart = 0 };
        rt_fatal_win32err(GetFileSizeEx(file, &eof));
        if (rw && *bytes > eof.QuadPart) { // increase file size
            const LARGE_INTEGER size = { .QuadPart = *bytes };
            r = r != 0 ? r : (rt_b2e(SetFilePointerEx(file, size, null, FILE_BEGIN)));
            r = r != 0 ? r : (rt_b2e(SetEndOfFile(file)));
            // the following not guaranteed to work but helps with sparse files
            r = r != 0 ? r : (rt_b2e(SetFileValidData(file, *bytes)));
            // SetFileValidData() only works for Admin (verified) or System accounts
            if (r == ERROR_PRIVILEGE_NOT_HELD) { r = 0; } // ignore
            // SetFileValidData() is also semi-security hole because it allows to read
            // previously not zeroed disk content of other files
            const LARGE_INTEGER zero = { .QuadPart = 0 }; // rewind stream:
            r = r != 0 ? r : (rt_b2e(SetFilePointerEx(file, zero, null, FILE_BEGIN)));
        } else {
            *bytes = eof.QuadPart;
        }
        r = r != 0 ? r : rt_mem_map_view_of_file(file, data, bytes, rw);
        rt_win32_close_handle(file);
    }
    return r;
}

static errno_t rt_mem_map_ro(const char* filename, void* *data, int64_t *bytes) {
    return rt_mem_map_file(filename, data, bytes, false);
}

static errno_t rt_mem_map_rw(const char* filename, void* *data, int64_t *bytes) {
    return rt_mem_map_file(filename, data, bytes, true);
}

static void rt_mem_unmap(void* data, int64_t bytes) {
    rt_assert(data != null && bytes > 0);
    (void)bytes; /* unused only need for posix version */
    if (data != null && bytes > 0) {
        rt_fatal_win32err(UnmapViewOfFile(data));
    }
}

static errno_t rt_mem_map_resource(const char* label, void* *data, int64_t *bytes) {
    HRSRC res = FindResourceA(null, label, (const char*)RT_RCDATA);
    // "LockResource does not actually lock memory; it is just used to
    // obtain a pointer to the memory containing the resource data.
    // The name of the function comes from versions prior to Windows XP,
    // when it was used to lock a global memory block allocated by LoadResource."
    if (res != null) { *bytes = SizeofResource(null, res); }
    HGLOBAL g = res != null ? LoadResource(null, res) : null;
    *data = g != null ? LockResource(g) : null;
    return *data != null ? 0 : rt_core.err();
}

static int32_t rt_mem_page_size(void) {
    static SYSTEM_INFO system_info;
    if (system_info.dwPageSize == 0) {
        GetSystemInfo(&system_info);
    }
    return (int32_t)system_info.dwPageSize;
}

static int rt_mem_large_page_size(void) {
    static SIZE_T large_page_minimum = 0;
    if (large_page_minimum == 0) {
        large_page_minimum = GetLargePageMinimum();
    }
    return (int32_t)large_page_minimum;
}

static void* rt_mem_allocate(int64_t bytes_multiple_of_page_size) {
    rt_assert(bytes_multiple_of_page_size > 0);
    SIZE_T bytes = (SIZE_T)bytes_multiple_of_page_size;
    SIZE_T page_size = (SIZE_T)rt_mem_page_size();
    rt_assert(bytes % page_size == 0);
    errno_t r = 0;
    void* a = null;
    if (bytes_multiple_of_page_size < 0 || bytes % page_size != 0) {
        SetLastError(ERROR_INVALID_PARAMETER);
        r = EINVAL;
    } else {
        const DWORD type = MEM_COMMIT | MEM_RESERVE;
        const DWORD physical = type | MEM_PHYSICAL;
        a = VirtualAlloc(null, bytes, physical, PAGE_READWRITE);
        if (a == null) {
            a = VirtualAlloc(null, bytes, type, PAGE_READWRITE);
        }
        if (a == null) {
            r = rt_core.err();
            if (r != 0) {
                rt_println("VirtualAlloc(%lld) failed %s", bytes, rt_strerr(r));
            }
        } else {
            r = VirtualLock(a, bytes) ? 0 : rt_core.err();
            if (r == ERROR_WORKING_SET_QUOTA) {
                // The default size is 345 pages (for example,
                // this is 1,413,120 bytes on systems with a 4K page size).
                SIZE_T min_mem = 0, max_mem = 0;
                r = rt_b2e(GetProcessWorkingSetSize(GetCurrentProcess(), &min_mem, &max_mem));
                if (r != 0) {
                    rt_println("GetProcessWorkingSetSize() failed %s", rt_strerr(r));
                } else {
                    max_mem =  max_mem + bytes * 2LL;
                    max_mem = (max_mem + page_size - 1) / page_size * page_size +
                               page_size * 16;
                    if (min_mem < max_mem) { min_mem = max_mem; }
                    r = rt_b2e(SetProcessWorkingSetSize(GetCurrentProcess(),
                            min_mem, max_mem));
                    if (r != 0) {
                        rt_println("SetProcessWorkingSetSize(%lld, %lld) failed %s",
                            (uint64_t)min_mem, (uint64_t)max_mem, rt_strerr(r));
                    } else {
                        r = rt_b2e(VirtualLock(a, bytes));
                    }
                }
            }
            if (r != 0) {
                rt_println("VirtualLock(%lld) failed %s", bytes, rt_strerr(r));
            }
        }
    }
    if (r != 0) {
        rt_println("mem_alloc_pages(%lld) failed %s", bytes, rt_strerr(r));
        rt_assert(a == null);
    }
    return a;
}

static void rt_mem_deallocate(void* a, int64_t bytes_multiple_of_page_size) {
    rt_assert(bytes_multiple_of_page_size > 0);
    SIZE_T bytes = (SIZE_T)bytes_multiple_of_page_size;
    errno_t r = 0;
    SIZE_T page_size = (SIZE_T)rt_mem_page_size();
    if (bytes_multiple_of_page_size < 0 || bytes % page_size != 0) {
        r = EINVAL;
        rt_println("failed %s", rt_strerr(r));
    } else {
        if (a != null) {
            // in case it was successfully locked
            r = rt_b2e(VirtualUnlock(a, bytes));
            if (r != 0) {
                rt_println("VirtualUnlock() failed %s", rt_strerr(r));
            }
            // If the "dwFreeType" parameter is MEM_RELEASE, "dwSize" parameter
            // must be the base address returned by the VirtualAlloc function when
            // the region of pages is reserved.
            r = rt_b2e(VirtualFree(a, 0, MEM_RELEASE));
            if (r != 0) { rt_println("VirtuaFree() failed %s", rt_strerr(r)); }
        }
    }
}

static void rt_mem_test(void) {
    #ifdef RT_TESTS
    rt_swear(rt_args.c > 0);
    void* data = null;
    int64_t bytes = 0;
    rt_swear(rt_mem.map_ro(rt_args.v[0], &data, &bytes) == 0);
    rt_swear(data != null && bytes != 0);
    rt_mem.unmap(data, bytes);
    // TODO: page_size large_page_size allocate deallocate
    // TODO: test heap functions
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
    #endif
}

rt_mem_if rt_mem = {
    .map_ro          = rt_mem_map_ro,
    .map_rw          = rt_mem_map_rw,
    .unmap           = rt_mem_unmap,
    .map_resource    = rt_mem_map_resource,
    .page_size       = rt_mem_page_size,
    .large_page_size = rt_mem_large_page_size,
    .allocate        = rt_mem_allocate,
    .deallocate      = rt_mem_deallocate,
    .test            = rt_mem_test
};

// _________________________________ rt_nls.c _________________________________

// Simplistic Win32 implementation of national language support.
// Windows NLS family of functions is very complicated and has
// difficult history of LANGID vs LCID etc... See:
// ResolveLocaleName()
// GetThreadLocale()
// SetThreadLocale()
// GetUserDefaultLocaleName()
// WM_SETTINGCHANGE lParam="intl"
// and many others...

enum {
    rt_nls_str_count_max = 1024,
    rt_nls_str_mem_max = 64 * rt_nls_str_count_max
};

static char  rt_nls_strings_memory[rt_nls_str_mem_max]; // increase if overflows
static char* rt_nls_strings_free = rt_nls_strings_memory;

static int32_t rt_nls_strings_count;

static const char* rt_nls_ls[rt_nls_str_count_max]; // localized strings
static const char* rt_nls_ns[rt_nls_str_count_max]; // neutral language strings

static uint16_t* rt_nls_load_string(int32_t strid, LANGID lang_id) {
    rt_assert(0 <= strid && strid < rt_countof(rt_nls_ns));
    uint16_t* r = null;
    int32_t block = strid / 16 + 1;
    int32_t index  = strid % 16;
    HRSRC res = FindResourceExW(((HMODULE)null), RT_STRING,
        MAKEINTRESOURCEW(block), lang_id);
//  rt_println("FindResourceExA(block=%d lang_id=%04X)=%p", block, lang_id, res);
    uint8_t* memory = res == null ? null : (uint8_t*)LoadResource(null, res);
    uint16_t* ws = memory == null ? null : (uint16_t*)LockResource(memory);
//  rt_println("LockResource(block=%d lang_id=%04X)=%p", block, lang_id, ws);
    if (ws != null) {
        for (int32_t i = 0; i < 16 && r == null; i++) {
            if (ws[0] != 0) {
                int32_t count = (int32_t)ws[0];  // String size in characters.
                ws++;
                rt_assert(ws[count - 1] == 0, "use rc.exe /n command line option");
                if (i == index) { // the string has been found
//                  rt_println("%04X found %s", lang_id, utf16to8(ws));
                    r = ws;
                }
                ws += count;
            } else {
                ws++;
            }
        }
    }
    return r;
}

static const char* rt_nls_save_string(uint16_t* utf16) {
    const int32_t bytes = rt_str.utf8_bytes(utf16, -1);
    rt_swear(bytes > 1);
    char* s = rt_nls_strings_free;
    uintptr_t left = (uintptr_t)rt_countof(rt_nls_strings_memory) -
        (uintptr_t)(rt_nls_strings_free - rt_nls_strings_memory);
    rt_fatal_if(left < (uintptr_t)bytes, "string_memory[] overflow");
    rt_str.utf16to8(s, (int32_t)left, utf16, -1);
    rt_assert((int32_t)strlen(s) == bytes - 1, "utf16to8() does not truncate");
    rt_nls_strings_free += bytes;
    return s;
}

static const char* rt_nls_localized_string(int32_t strid) {
    rt_swear(0 < strid && strid < rt_countof(rt_nls_ns));
    const char* s = null;
    if (0 < strid && strid < rt_countof(rt_nls_ns)) {
        if (rt_nls_ls[strid] != null) {
            s = rt_nls_ls[strid];
        } else {
            LCID lc_id = GetThreadLocale();
            LANGID lang_id = LANGIDFROMLCID(lc_id);
            uint16_t* utf16 = rt_nls_load_string(strid, lang_id);
            if (utf16 == null) { // try default dialect:
                LANGID primary = PRIMARYLANGID(lang_id);
                lang_id = MAKELANGID(primary, SUBLANG_NEUTRAL);
                utf16 = rt_nls_load_string(strid, lang_id);
            }
            if (utf16 != null && utf16[0] != 0x0000) {
                s = rt_nls_save_string(utf16);
                rt_nls_ls[strid] = s;
            }
        }
    }
    return s;
}

static int32_t rt_nls_strid(const char* s) {
    int32_t strid = -1;
    for (int32_t i = 1; i < rt_nls_strings_count && strid == -1; i++) {
        if (rt_nls_ns[i] != null && strcmp(s, rt_nls_ns[i]) == 0) {
            strid = i;
            rt_nls_localized_string(strid); // to save it, ignore result
        }
    }
    return strid;
}

static const char* rt_nls_string(int32_t strid, const char* defau1t) {
    const char* r = rt_nls_localized_string(strid);
    return r == null ? defau1t : r;
}

static const char* rt_nls_str(const char* s) {
    int32_t id = rt_nls_strid(s);
    return id < 0 ? s : rt_nls_string(id, s);
}

static const char* rt_nls_locale(void) {
    uint16_t utf16[LOCALE_NAME_MAX_LENGTH + 1];
    LCID lc_id = GetThreadLocale();
    int32_t n = LCIDToLocaleName(lc_id, utf16, rt_countof(utf16),
        LOCALE_ALLOW_NEUTRAL_NAMES);
    static char ln[LOCALE_NAME_MAX_LENGTH * 4 + 1];
    ln[0] = 0;
    if (n == 0) {
        errno_t r = rt_core.err();
        rt_println("LCIDToLocaleName(0x%04X) failed %s", lc_id, rt_str.error(r));
    } else {
        rt_str.utf16to8(ln, rt_countof(ln), utf16, -1);
    }
    return ln;
}

static errno_t rt_nls_set_locale(const char* locale) {
    errno_t r = 0;
    uint16_t utf16[LOCALE_NAME_MAX_LENGTH + 1];
    rt_str.utf8to16(utf16, rt_countof(utf16), locale, -1);
    uint16_t rln[LOCALE_NAME_MAX_LENGTH + 1]; // resolved locale name
    int32_t n = (int32_t)ResolveLocaleName(utf16, rln, (DWORD)rt_countof(rln));
    if (n == 0) {
        r = rt_core.err();
        rt_println("ResolveLocaleName(\"%s\") failed %s", locale, rt_str.error(r));
    } else {
        LCID lc_id = LocaleNameToLCID(rln, LOCALE_ALLOW_NEUTRAL_NAMES);
        if (lc_id == 0) {
            r = rt_core.err();
            rt_println("LocaleNameToLCID(\"%s\") failed %s", locale, rt_str.error(r));
        } else {
            rt_fatal_win32err(SetThreadLocale(lc_id));
            memset((void*)rt_nls_ls, 0, sizeof(rt_nls_ls)); // start all over
        }
    }
    return r;
}

static void rt_nls_init(void) {
    static_assert(rt_countof(rt_nls_ns) % 16 == 0, 
                 "rt_countof(ns) must be multiple of 16");
    LANGID lang_id = MAKELANGID(LANG_ENGLISH, SUBLANG_NEUTRAL);
    for (int32_t strid = 0; strid < rt_countof(rt_nls_ns); strid += 16) {
        int32_t block = strid / 16 + 1;
        HRSRC res = FindResourceExW(((HMODULE)null), RT_STRING,
            MAKEINTRESOURCEW(block), lang_id);
        uint8_t* memory = res == null ? null : (uint8_t*)LoadResource(null, res);
        uint16_t* ws = memory == null ? null : (uint16_t*)LockResource(memory);
        if (ws == null) { break; }
        for (int32_t i = 0; i < 16; i++) {
            int32_t ix = strid + i;
            uint16_t count = ws[0];
            if (count > 0) {
                ws++;
                rt_fatal_if(ws[count - 1] != 0, "use rc.exe /n");
                rt_nls_ns[ix] = rt_nls_save_string(ws);
                rt_nls_strings_count = ix + 1;
//              rt_println("ns[%d] := %d \"%s\"", ix, strlen(rt_nls_ns[ix]), rt_nls_ns[ix]);
                ws += count;
            } else {
                ws++;
            }
        }
    }
}

rt_nls_if rt_nls = {
    .init       = rt_nls_init,
    .strid      = rt_nls_strid,
    .str        = rt_nls_str,
    .string     = rt_nls_string,
    .locale     = rt_nls_locale,
    .set_locale = rt_nls_set_locale,
};
// _________________________________ rt_num.c _________________________________

#include <intrin.h>
//#include <immintrin.h> // _tzcnt_u32

static inline rt_num128_t rt_num_add128_inline(const rt_num128_t a, const rt_num128_t b) {
    rt_num128_t r = a;
    r.hi += b.hi;
    r.lo += b.lo;
    if (r.lo < b.lo) { r.hi++; } // carry
    return r;
}

static inline rt_num128_t rt_num_sub128_inline(const rt_num128_t a, const rt_num128_t b) {
    rt_num128_t r = a;
    r.hi -= b.hi;
    if (r.lo < b.lo) { r.hi--; } // borrow
    r.lo -= b.lo;
    return r;
}

static rt_num128_t rt_num_add128(const rt_num128_t a, const rt_num128_t b) {
    return rt_num_add128_inline(a, b);
}

static rt_num128_t rt_num_sub128(const rt_num128_t a, const rt_num128_t b) {
    return rt_num_sub128_inline(a, b);
}

static rt_num128_t rt_num_mul64x64(uint64_t a, uint64_t b) {
    uint64_t a_lo = (uint32_t)a;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b;
    uint64_t b_hi = b >> 32;
    uint64_t low = a_lo * b_lo;
    uint64_t cross1 = a_hi * b_lo;
    uint64_t cross2 = a_lo * b_hi;
    uint64_t high = a_hi * b_hi;
    // this cannot overflow as (2^32-1)^2 + 2^32-1 < 2^64-1
    cross1 += low >> 32;
    // this one can overflow
    cross1 += cross2;
    // propagate the carry if any
    high += ((uint64_t)(cross1 < cross2 != 0)) << 32;
    high = high + (cross1 >> 32);
    low = ((cross1 & 0xFFFFFFFF) << 32) + (low & 0xFFFFFFFF);
    return (rt_num128_t){.lo = low, .hi = high };
}

static inline void rt_num_shift128_left_inline(rt_num128_t* n) {
    const uint64_t top = (1ULL << 63);
    n->hi = (n->hi << 1) | ((n->lo & top) ? 1 : 0);
    n->lo = (n->lo << 1);
}

static inline void rt_num_shift128_right_inline(rt_num128_t* n) {
    const uint64_t top = (1ULL << 63);
    n->lo = (n->lo >> 1) | ((n->hi & 0x1) ? top : 0);
    n->hi = (n->hi >> 1);
}

static inline bool rt_num_less128_inline(const rt_num128_t a, const rt_num128_t b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

static inline bool rt_num_uint128_high_bit(const rt_num128_t a) {
    return (int64_t)a.hi < 0;
}

static uint64_t rt_num_muldiv128(uint64_t a, uint64_t b, uint64_t divisor) {
    rt_swear(divisor > 0, "divisor: %lld", divisor);
    rt_num128_t r = rt_num.mul64x64(a, b); // reminder: a * b
    uint64_t q = 0; // quotient
    if (r.hi >= divisor) {
        q = UINT64_MAX; // overflow
    } else {
        int32_t  shift = 0;
        rt_num128_t d = { .hi = 0, .lo = divisor };
        while (!rt_num_uint128_high_bit(d) && rt_num_less128_inline(d, r)) {
            rt_num_shift128_left_inline(&d);
            shift++;
        }
        rt_assert(shift <= 64);
        while (shift >= 0 && (d.hi != 0 || d.lo != 0)) {
            if (!rt_num_less128_inline(r, d)) {
                r = rt_num_sub128_inline(r, d);
                rt_assert(shift < 64);
                q |= (1ULL << shift);
            }
            rt_num_shift128_right_inline(&d);
            shift--;
        }
    }
    return q;
}

static uint32_t rt_num_gcd32(uint32_t u, uint32_t v) {
    #pragma push_macro("rt_trailing_zeros")
    #ifdef _M_ARM64
    #define rt_trailing_zeros(x) (_CountTrailingZeros(x))
    #else
    #define rt_trailing_zeros(x) ((int32_t)_tzcnt_u32(x))
    #endif
    if (u == 0) {
        return v;
    } else if (v == 0) {
        return u;
    }
    uint32_t i = rt_trailing_zeros(u);  u >>= i;
    uint32_t j = rt_trailing_zeros(v);  v >>= j;
    uint32_t k = rt_min(i, j);
    for (;;) {
        rt_assert(u % 2 == 1, "u = %d should be odd", u);
        rt_assert(v % 2 == 1, "v = %d should be odd", v);
        if (u > v) { uint32_t swap = u; u = v; v = swap; }
        v -= u;
        if (v == 0) { return u << k; }
        v >>= rt_trailing_zeros(v);
    }
    #pragma pop_macro("rt_trailing_zeros")
}

static uint32_t rt_num_random32(uint32_t* state) {
    // https://gist.github.com/tommyettinger/46a874533244883189143505d203312c
    static rt_thread_local bool started; // first seed must be odd
    if (!started) { started = true; *state |= 1; }
    uint32_t z = (*state += 0x6D2B79F5UL);
    z = (z ^ (z >> 15)) * (z | 1UL);
    z ^= z + (z ^ (z >> 7)) * (z | 61UL);
    return z ^ (z >> 14);
}

static uint64_t rt_num_random64(uint64_t *state) {
    // https://gist.github.com/tommyettinger/e6d3e8816da79b45bfe582384c2fe14a
    static rt_thread_local bool started; // first seed must be odd
    if (!started) { started = true; *state |= 1; }
	const uint64_t s = *state;
	const uint64_t z = (s ^ s >> 25) * (*state += 0x6A5D39EAE12657AAULL);
	return z ^ (z >> 22);
}

// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function

static uint32_t rt_num_hash32(const char *data, int64_t len) {
    uint32_t hash  = 0x811c9dc5;  // FNV_offset_basis for 32-bit
    uint32_t prime = 0x01000193; // FNV_prime for 32-bit
    if (len > 0) {
        for (int64_t i = 1; i < len; i++) {
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

static uint64_t rt_num_hash64(const char *data, int64_t len) {
    uint64_t hash  = 0xcbf29ce484222325; // FNV_offset_basis for 64-bit
    uint64_t prime = 0x100000001b3;      // FNV_prime for 64-bit
    if (len > 0) {
        for (int64_t i = 0; i < len; i++) {
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

static uint32_t ctz_2(uint32_t x) {
    if (x == 0) return 32;
    unsigned n = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        n++;
    }
    return n;
}

static void rt_num_test(void) {
    #ifdef RT_TESTS
    {
        rt_swear(rt_num.gcd32(1000000000, 24000000) == 8000000);
        // https://asecuritysite.com/encryption/nprimes?y=64
        // https://www.rapidtables.com/convert/number/decimal-to-hex.html
        uint64_t p = 15843490434539008357u; // prime
        uint64_t q = 16304766625841520833u; // prime
        // pq: 258324414073910997987910483408576601381
        //     0xC25778F20853A9A1EC0C27C467C45D25
        rt_num128_t pq = {.hi = 0xC25778F20853A9A1uLL,
                       .lo = 0xEC0C27C467C45D25uLL };
        rt_num128_t p_q = rt_num.mul64x64(p, q);
        rt_swear(p_q.hi == pq.hi && pq.lo == pq.lo);
        uint64_t p1 = rt_num.muldiv128(p, q, q);
        uint64_t q1 = rt_num.muldiv128(p, q, p);
        rt_swear(p1 == p);
        rt_swear(q1 == q);
    }
    #ifdef DEBUG
    enum { n = 100 };
    #else
    enum { n = 10000 };
    #endif
    uint64_t seed64 = 1;
    for (int32_t i = 0; i < n; i++) {
        uint64_t p = rt_num.random64(&seed64);
        uint64_t q = rt_num.random64(&seed64);
        uint64_t p1 = rt_num.muldiv128(p, q, q);
        uint64_t q1 = rt_num.muldiv128(p, q, p);
        rt_swear(p == p1, "0%16llx (0%16llu) != 0%16llx (0%16llu)", p, p1);
        rt_swear(q == q1, "0%16llx (0%16llu) != 0%16llx (0%16llu)", p, p1);
    }
    uint32_t seed32 = 1;
    for (int32_t i = 0; i < n; i++) {
        uint64_t p = rt_num.random32(&seed32);
        uint64_t q = rt_num.random32(&seed32);
        uint64_t r = rt_num.muldiv128(p, q, 1);
        rt_swear(r == p * q);
        // division by the maximum uint64_t value:
        r = rt_num.muldiv128(p, q, UINT64_MAX);
        rt_swear(r == 0);
    }
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
    #endif
}

rt_num_if rt_num = {
    .add128    = rt_num_add128,
    .sub128    = rt_num_sub128,
    .mul64x64  = rt_num_mul64x64,
    .muldiv128 = rt_num_muldiv128,
    .gcd32     = rt_num_gcd32,
    .random32  = rt_num_random32,
    .random64  = rt_num_random64,
    .hash32    = rt_num_hash32,
    .hash64    = rt_num_hash64,
    .test      = rt_num_test
};

// ______________________________ rt_processes.c ______________________________

typedef struct rt_processes_pidof_lambda_s rt_processes_pidof_lambda_t;

typedef struct rt_processes_pidof_lambda_s {
    bool (*each)(rt_processes_pidof_lambda_t* p, uint64_t pid); // returns true to continue
    uint64_t* pids;
    size_t size;  // pids[size]
    size_t count; // number of valid pids in the pids
    fp64_t timeout;
    errno_t error;
} rt_processes_pidof_lambda_t;

static int32_t rt_processes_for_each_pidof(const char* pname, rt_processes_pidof_lambda_t* la) {
    char stack[1024]; // avoid alloca()
    int32_t n = rt_str.len(pname);
    rt_fatal_if(n + 5 >= rt_countof(stack), "name is too long: %s", pname);
    const char* name = pname;
    // append ".exe" if not present:
    if (!rt_str.iends(pname, ".exe")) {
        int32_t k = (int32_t)strlen(pname) + 5;
        char* exe = stack;
        rt_str.format(exe, k, "%s.exe", pname);
        name = exe;
    }
    const char* base = strrchr(name, '\\');
    if (base != null) {
        base++; // advance past "\\"
    } else {
        base = name;
    }
    uint16_t wn[1024];
    rt_fatal_if(strlen(base) >= rt_countof(wn), "name too long: %s", base);
    rt_str.utf8to16(wn, rt_countof(wn), base, -1);
    size_t count = 0;
    uint64_t pid = 0;
    uint8_t* data = null;
    ULONG bytes = 0;
    errno_t r = NtQuerySystemInformation(SystemProcessInformation, data, 0, &bytes);
    #pragma push_macro("STATUS_INFO_LENGTH_MISMATCH")
    #define STATUS_INFO_LENGTH_MISMATCH      0xC0000004
    while (r == (errno_t)STATUS_INFO_LENGTH_MISMATCH) {
        // bytes == 420768 on Windows 11 which may be a bit
        // too much for stack alloca()
        // add little extra if new process is spawned in between calls.
        bytes += sizeof(SYSTEM_PROCESS_INFORMATION) * 32;
        r = rt_heap.reallocate(null, (void**)&data, bytes, false);
        if (r == 0) {
            r = NtQuerySystemInformation(SystemProcessInformation, data, bytes, &bytes);
        } else {
            rt_assert(r == (errno_t)ERROR_NOT_ENOUGH_MEMORY);
        }
    }
    #pragma pop_macro("STATUS_INFO_LENGTH_MISMATCH")
    if (r == 0 && data != null) {
        SYSTEM_PROCESS_INFORMATION* proc = (SYSTEM_PROCESS_INFORMATION*)data;
        while (proc != null) {
            uint16_t* img = proc->ImageName.Buffer; // last name only, not a pathname!
            bool match = img != null && wcsicmp(img, wn) == 0;
            if (match) {
                pid = (uint64_t)proc->UniqueProcessId; // HANDLE .UniqueProcessId
                if (base != name) {
                    char path[rt_files_max_path];
                    match = rt_processes.nameof(pid, path, rt_countof(path)) == 0 &&
                            rt_str.iends(path, name);
//                  rt_println("\"%s\" -> \"%s\" match: %d", name, path, match);
                }
            }
            if (match) {
                if (la != null && count < la->size && la->pids != null) {
                    la->pids[count] = pid;
                }
                count++;
                if (la != null && la->each != null && !la->each(la, pid)) {
                    break;
                }
            }
            proc = proc->NextEntryOffset != 0 ? (SYSTEM_PROCESS_INFORMATION*)
                ((uint8_t*)proc + proc->NextEntryOffset) : null;
        }
    }
    if (data != null) { rt_heap.deallocate(null, data); }
    rt_assert(count <= (uint64_t)INT32_MAX);
    return (int32_t)count;
}

static errno_t rt_processes_nameof(uint64_t pid, char* name, int32_t count) {
    rt_assert(name != null && count > 0);
    errno_t r = 0;
    name[0] = 0;
    HANDLE p = OpenProcess(PROCESS_ALL_ACCESS, false, (DWORD)pid);
    if (p != null) {
        r = rt_b2e(GetModuleFileNameExA(p, null, name, count));
        name[count - 1] = 0; // ensure zero termination
        rt_win32_close_handle(p);
    } else {
        r = ERROR_NOT_FOUND;
    }
    return r;
}

static bool rt_processes_present(uint64_t pid) {
    void* h = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, 0, (DWORD)pid);
    bool b = h != null;
    if (h != null) { rt_win32_close_handle(h); }
    return b;
}

static bool rt_processes_first_pid(rt_processes_pidof_lambda_t* lambda, uint64_t pid) {
    lambda->pids[0] = pid;
    return false;
}

static uint64_t rt_processes_pid(const char* pname) {
    uint64_t first[1] = {0};
    rt_processes_pidof_lambda_t lambda = {
        .each = rt_processes_first_pid,
        .pids = first,
        .size  = 1,
        .count = 0,
        .timeout = 0,
        .error = 0
    };
    rt_processes_for_each_pidof(pname, &lambda);
    return first[0];
}

static bool rt_processes_store_pid(rt_processes_pidof_lambda_t* lambda, uint64_t pid) {
    if (lambda->pids != null && lambda->count < lambda->size) {
        lambda->pids[lambda->count++] = pid;
    }
    return true; // always - need to count all
}

static errno_t rt_processes_pids(const char* pname, uint64_t* pids/*[size]*/,
        int32_t size, int32_t *count) {
    *count = 0;
    rt_processes_pidof_lambda_t lambda = {
        .each = rt_processes_store_pid,
        .pids = pids,
        .size = (size_t)size,
        .count = 0,
        .timeout = 0,
        .error = 0
    };
    *count = rt_processes_for_each_pidof(pname, &lambda);
    return (int32_t)lambda.count == *count ? 0 : ERROR_MORE_DATA;
}

static errno_t rt_processes_kill(uint64_t pid, fp64_t timeout) {
    DWORD milliseconds = timeout < 0 ? INFINITE : (DWORD)(timeout * 1000);
    enum { access = PROCESS_QUERY_LIMITED_INFORMATION |
                    PROCESS_TERMINATE | SYNCHRONIZE };
    rt_assert((DWORD)pid == pid); // Windows... HANDLE vs DWORD in different APIs
    errno_t r = ERROR_NOT_FOUND;
    HANDLE h = OpenProcess(access, 0, (DWORD)pid);
    if (h != null) {
        char path[rt_files_max_path];
        path[0] = 0;
        r = rt_b2e(TerminateProcess(h, ERROR_PROCESS_ABORTED));
        if (r == 0) {
            DWORD ix = WaitForSingleObject(h, milliseconds);
            r = rt_wait_ix2e(ix);
        } else {
            DWORD bytes = rt_countof(path);
            errno_t rq = rt_b2e(QueryFullProcessImageNameA(h, 0, path, &bytes));
            if (rq != 0) {
                rt_println("QueryFullProcessImageNameA(pid=%d, h=%p) "
                        "failed %s", pid, h, rt_strerr(rq));
            }
        }
        rt_win32_close_handle(h);
        if (r == ERROR_ACCESS_DENIED) { // special case
            rt_thread.sleep_for(0.015); // need to wait a bit
            HANDLE retry = OpenProcess(access, 0, (DWORD)pid);
            // process may have died before we have chance to terminate it:
            if (retry == null) {
                rt_println("TerminateProcess(pid=%d, h=%p, im=%s) "
                        "failed but zombie died after: %s",
                        pid, h, path, rt_strerr(r));
                r = 0;
            } else {
                rt_win32_close_handle(retry);
            }
        }
        if (r != 0) {
            rt_println("TerminateProcess(pid=%d, h=%p, im=%s) failed %s",
                pid, h, path, rt_strerr(r));
        }
    }
    if (r != 0) { errno = r; }
    return r;
}

static bool rt_processes_kill_one(rt_processes_pidof_lambda_t* lambda, uint64_t pid) {
    errno_t r = rt_processes_kill(pid, lambda->timeout);
    if (r != 0) { lambda->error = r; }
    return true; // keep going
}

static errno_t rt_processes_kill_all(const char* name, fp64_t timeout) {
    rt_processes_pidof_lambda_t lambda = {
        .each = rt_processes_kill_one,
        .pids = null,
        .size  = 0,
        .count = 0,
        .timeout = timeout,
        .error = 0
    };
    int32_t c = rt_processes_for_each_pidof(name, &lambda);
    return c == 0 ? ERROR_NOT_FOUND : lambda.error;
}

static bool rt_processes_is_elevated(void) { // Is process running as Admin / System ?
    BOOL elevated = false;
    PSID administrators_group = null;
    // Allocate and initialize a SID of the administrators group.
    SID_IDENTIFIER_AUTHORITY administrators_group_authority = SECURITY_NT_AUTHORITY;
    errno_t r = rt_b2e(AllocateAndInitializeSid(&administrators_group_authority, 2,
                SECURITY_BUILTIN_DOMAIN_RID, DOMAIN_ALIAS_RID_ADMINS,
                0, 0, 0, 0, 0, 0, &administrators_group));
    if (r != 0) {
        rt_println("AllocateAndInitializeSid() failed %s", rt_strerr(r));
    }
    PSID system_ops = null;
    SID_IDENTIFIER_AUTHORITY system_ops_authority = SECURITY_NT_AUTHORITY;
    r = rt_b2e(AllocateAndInitializeSid(&system_ops_authority, 2,
            SECURITY_BUILTIN_DOMAIN_RID, DOMAIN_ALIAS_RID_SYSTEM_OPS,
            0, 0, 0, 0, 0, 0, &system_ops));
    if (r != 0) {
        rt_println("AllocateAndInitializeSid() failed %s", rt_strerr(r));
    }
    if (administrators_group != null) {
        r = rt_b2e(CheckTokenMembership(null, administrators_group, &elevated));
    }
    if (system_ops != null && !elevated) {
        r = rt_b2e(CheckTokenMembership(null, administrators_group, &elevated));
    }
    if (administrators_group != null) { FreeSid(administrators_group); }
    if (system_ops != null) { FreeSid(system_ops); }
    if (r != 0) {
        rt_println("failed %s", rt_strerr(r));
    }
    return elevated;
}

static errno_t rt_processes_restart_elevated(void) {
    errno_t r = 0;
    if (!rt_processes.is_elevated()) {
        const char* path = rt_processes.name();
        SHELLEXECUTEINFOA sei = { sizeof(sei) };
        sei.lpVerb = "runas";
        sei.lpFile = path;
        sei.hwnd = null;
        sei.nShow = SW_NORMAL;
        r = rt_b2e(ShellExecuteExA(&sei));
        if (r == ERROR_CANCELLED) {
            rt_println("The user unable or refused to allow privileges elevation");
        } else if (r == 0) {
            rt_core.exit(0); // second copy of the app is running now
        }
    }
    return r;
}

static void rt_processes_close_pipes(STARTUPINFOA* si,
        HANDLE *read_out,
        HANDLE *read_err,
        HANDLE *write_in) {
    if (si->hStdOutput != INVALID_HANDLE_VALUE) { rt_win32_close_handle(si->hStdOutput); }
    if (si->hStdError  != INVALID_HANDLE_VALUE) { rt_win32_close_handle(si->hStdError);  }
    if (si->hStdInput  != INVALID_HANDLE_VALUE) { rt_win32_close_handle(si->hStdInput);  }
    if (*read_out != INVALID_HANDLE_VALUE) { rt_win32_close_handle(*read_out); }
    if (*read_err != INVALID_HANDLE_VALUE) { rt_win32_close_handle(*read_err); }
    if (*write_in != INVALID_HANDLE_VALUE) { rt_win32_close_handle(*write_in); }
}

static errno_t rt_processes_child_read(rt_stream_if* out, HANDLE pipe) {
    char data[32 * 1024]; // Temporary buffer for reading
    DWORD available = 0;
    errno_t r = rt_b2e(PeekNamedPipe(pipe, null, sizeof(data), null,
                                 &available, null));
    if (r != 0) {
        if (r != ERROR_BROKEN_PIPE) { // unexpected!
//          rt_println("PeekNamedPipe() failed %s", rt_strerr(r));
        }
        // process has exited and closed the pipe
        rt_assert(r == ERROR_BROKEN_PIPE);
    } else if (available > 0) {
        DWORD bytes_read = 0;
        r = rt_b2e(ReadFile(pipe, data, sizeof(data), &bytes_read, null));
//      rt_println("r: %d bytes_read: %d", r, bytes_read);
        if (out != null) {
            if (r == 0) {
                r = out->write(out, data, bytes_read, null);
            }
        } else {
            // no one interested - drop on the floor
        }
    }
    return r;
}

static errno_t rt_processes_child_write(rt_stream_if* in, HANDLE pipe) {
    errno_t r = 0;
    if (in != null) {
        uint8_t  memory[32 * 1024]; // Temporary buffer for reading
        uint8_t* data = memory;
        int64_t bytes_read = 0;
        in->read(in, data, sizeof(data), &bytes_read);
        while (r == 0 && bytes_read > 0) {
            DWORD bytes_written = 0;
            r = rt_b2e(WriteFile(pipe, data, (DWORD)bytes_read,
                             &bytes_written, null));
            rt_println("r: %d bytes_written: %d", r, bytes_written);
            rt_assert((int32_t)bytes_written <= bytes_read);
            data += bytes_written;
            bytes_read -= bytes_written;
        }
    }
    return r;
}

static errno_t rt_processes_run(rt_processes_child_t* child) {
    const fp64_t deadline = rt_clock.seconds() + child->timeout;
    errno_t r = 0;
    STARTUPINFOA si = {
        .cb = sizeof(STARTUPINFOA),
        .hStdInput  = INVALID_HANDLE_VALUE,
        .hStdOutput = INVALID_HANDLE_VALUE,
        .hStdError  = INVALID_HANDLE_VALUE,
        .dwFlags     = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES,
        .wShowWindow = SW_HIDE
    };
    SECURITY_ATTRIBUTES sa = { sizeof(sa), null, true };  // Inheritable handles
    PROCESS_INFORMATION pi = {0};
    HANDLE read_out = INVALID_HANDLE_VALUE;
    HANDLE read_err = INVALID_HANDLE_VALUE;
    HANDLE write_in = INVALID_HANDLE_VALUE;
    errno_t ro = rt_b2e(CreatePipe(&read_out, &si.hStdOutput, &sa, 0));
    errno_t re = rt_b2e(CreatePipe(&read_err, &si.hStdError,  &sa, 0));
    errno_t ri = rt_b2e(CreatePipe(&si.hStdInput, &write_in,  &sa, 0));
    if (ro != 0 || re != 0 || ri != 0) {
        rt_processes_close_pipes(&si, &read_out, &read_err, &write_in);
        if (ro != 0) { rt_println("CreatePipe() failed %s", rt_strerr(ro)); r = ro; }
        if (re != 0) { rt_println("CreatePipe() failed %s", rt_strerr(re)); r = re; }
        if (ri != 0) { rt_println("CreatePipe() failed %s", rt_strerr(ri)); r = ri; }
    }
    if (r == 0) {
        r = rt_b2e(CreateProcessA(null, rt_str.drop_const(child->command),
                null, null, true, CREATE_NO_WINDOW, null, null, &si, &pi));
        if (r != 0) {
            rt_println("CreateProcess() failed %s", rt_strerr(r));
            rt_processes_close_pipes(&si, &read_out, &read_err, &write_in);
        }
    }
    if (r == 0) {
        // not relevant: stdout can be written in other threads
        rt_win32_close_handle(pi.hThread);
        pi.hThread = null;
        // need to close si.hStdO* handles on caller side so,
        // when the process closes handles of the pipes, EOF happens
        // on caller side with io result ERROR_BROKEN_PIPE
        // indicating no more data can be read or written
        rt_win32_close_handle(si.hStdOutput);
        rt_win32_close_handle(si.hStdError);
        rt_win32_close_handle(si.hStdInput);
        si.hStdOutput = INVALID_HANDLE_VALUE;
        si.hStdError  = INVALID_HANDLE_VALUE;
        si.hStdInput  = INVALID_HANDLE_VALUE;
        bool done = false;
        while (!done && r == 0) {
            if (child->timeout > 0 && rt_clock.seconds() > deadline) {
                r = rt_b2e(TerminateProcess(pi.hProcess, ERROR_SEM_TIMEOUT));
                if (r != 0) {
                    rt_println("TerminateProcess() failed %s", rt_strerr(r));
                } else {
                    done = true;
                }
            }
            if (r == 0) { r = rt_processes_child_write(child->in, write_in); }
            if (r == 0) { r = rt_processes_child_read(child->out, read_out); }
            if (r == 0) { r = rt_processes_child_read(child->err, read_err); }
            if (!done) {
                DWORD ix = WaitForSingleObject(pi.hProcess, 0);
                // ix == 0 means process has exited (or terminated)
                // r == ERROR_BROKEN_PIPE process closed one of the handles
                done = ix == WAIT_OBJECT_0 || r == ERROR_BROKEN_PIPE;
            }
            // to avoid tight loop 100% cpu utilization:
            if (!done) { rt_thread.yield(); }
        }
        // broken pipe actually signifies EOF on the pipe
        if (r == ERROR_BROKEN_PIPE) { r = 0; } // not an error
//      if (r != 0) { rt_println("pipe loop failed %s", rt_strerr(r));}
        DWORD xc = 0;
        errno_t rx = rt_b2e(GetExitCodeProcess(pi.hProcess, &xc));
        if (rx == 0) {
            child->exit_code = xc;
        } else {
            rt_println("GetExitCodeProcess() failed %s", rt_strerr(rx));
            if (r != 0) { r = rx; } // report earliest error
        }
        rt_processes_close_pipes(&si, &read_out, &read_err, &write_in);
        // expected never to fail
        rt_win32_close_handle(pi.hProcess);
    }
    return r;
}

typedef struct {
    rt_stream_if stream;
    rt_stream_if* output;
    errno_t error;
} rt_processes_io_merge_out_and_err_if;

static errno_t rt_processes_merge_write(rt_stream_if* stream, const void* data,
        int64_t bytes, int64_t* transferred) {
    if (transferred != null) { *transferred = 0; }
    rt_processes_io_merge_out_and_err_if* s =
        (rt_processes_io_merge_out_and_err_if*)stream;
    if (s->output != null && bytes > 0) {
        s->error = s->output->write(s->output, data, bytes, transferred);
    }
    return s->error;
}

static errno_t rt_processes_open(const char* command, int32_t *exit_code,
        rt_stream_if* output,  fp64_t timeout) {
    rt_not_null(output);
    rt_processes_io_merge_out_and_err_if merge_out_and_err = {
        .stream ={ .write = rt_processes_merge_write },
        .output = output,
        .error = 0
    };
    rt_processes_child_t child = {
        .command = command,
        .in = null,
        .out = &merge_out_and_err.stream,
        .err = &merge_out_and_err.stream,
        .exit_code = 0,
        .timeout = timeout
    };
    errno_t r = rt_processes.run(&child);
    if (exit_code != null) { *exit_code = (int32_t)child.exit_code; }
    uint8_t zero = 0; // zero termination
    merge_out_and_err.stream.write(&merge_out_and_err.stream, &zero, 1, null);
    if (r == 0 && merge_out_and_err.error != 0) {
        r = merge_out_and_err.error; // zero termination is not guaranteed
    }
    return r;
}

static errno_t rt_processes_spawn(const char* command) {
    errno_t r = 0;
    STARTUPINFOA si = {
        .cb = sizeof(STARTUPINFOA),
        .dwFlags     = STARTF_USESHOWWINDOW
                     | CREATE_NEW_PROCESS_GROUP
                     | DETACHED_PROCESS,
        .wShowWindow = SW_HIDE,
        .hStdInput  = INVALID_HANDLE_VALUE,
        .hStdOutput = INVALID_HANDLE_VALUE,
        .hStdError  = INVALID_HANDLE_VALUE
    };
    const DWORD flags = CREATE_BREAKAWAY_FROM_JOB
                | CREATE_NO_WINDOW
                | CREATE_NEW_PROCESS_GROUP
                | DETACHED_PROCESS;
    PROCESS_INFORMATION pi = { .hProcess = null, .hThread = null };
    r = rt_b2e(CreateProcessA(null, rt_str.drop_const(command), null, null,
            /*bInheritHandles:*/false, flags, null, null, &si, &pi));
    if (r == 0) { // Close handles immediately
        rt_win32_close_handle(pi.hProcess);
        rt_win32_close_handle(pi.hThread);
    } else {
        rt_println("CreateProcess() failed %s", rt_strerr(r));
    }
    return r;
}

static const char* rt_processes_name(void) {
    static char mn[rt_files_max_path];
    if (mn[0] == 0) {
        rt_fatal_win32err(GetModuleFileNameA(null, mn, rt_countof(mn)));
    }
    return mn;
}

#ifdef RT_TESTS

#pragma push_macro("verbose") // --verbosity trace

#define verbose(...) do {                                       \
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) { \
        rt_println(__VA_ARGS__);                                   \
    }                                                           \
} while (0)

static void rt_processes_test(void) {
    #ifdef RT_TESTS // in alphabetical order
    const char* names[] = { "svchost", "RuntimeBroker", "conhost" };
    for (int32_t j = 0; j < rt_countof(names); j++) {
        int32_t size  = 0;
        int32_t count = 0;
        uint64_t* pids = null;
        errno_t r = rt_processes.pids(names[j], null, size, &count);
        while (r == ERROR_MORE_DATA && count > 0) {
            size = count * 2; // set of processes may change rapidly
            r = rt_heap.reallocate(null, (void**)&pids,
                                  (int64_t)sizeof(uint64_t) * (int64_t)size,
                                  false);
            if (r == 0) {
                r = rt_processes.pids(names[j], pids, size, &count);
            }
        }
        if (r == 0 && count > 0) {
            for (int32_t i = 0; i < count; i++) {
                char path[256] = {0};
                #pragma warning(suppress: 6011) // dereferencing null
                r = rt_processes.nameof(pids[i], path, rt_countof(path));
                if (r != ERROR_NOT_FOUND) {
                    rt_assert(r == 0 && path[0] != 0);
                    verbose("%6d %s %s", pids[i], path, rt_strerr(r));
                }
            }
        }
        rt_heap.deallocate(null, pids);
    }
    // test popen()
    int32_t xc = 0;
    char data[32 * 1024];
    rt_stream_memory_if output;
    rt_streams.write_only(&output, data, rt_countof(data));
    const char* cmd = "cmd /c dir 2>nul >nul";
    errno_t r = rt_processes.popen(cmd, &xc, &output.stream, 99999.0);
    verbose("r: %d xc: %d output:\n%s", r, xc, data);
    rt_streams.write_only(&output, data, rt_countof(data));
    cmd = "cmd /c dir \"folder that does not exist\\\"";
    r = rt_processes.popen(cmd, &xc, &output.stream, 99999.0);
    verbose("r: %d xc: %d output:\n%s", r, xc, data);
    rt_streams.write_only(&output, data, rt_countof(data));
    cmd = "cmd /c dir";
    r = rt_processes.popen(cmd, &xc, &output.stream, 99999.0);
    verbose("r: %d xc: %d output:\n%s", r, xc, data);
    rt_streams.write_only(&output, data, rt_countof(data));
    cmd = "cmd /c timeout 1";
    r = rt_processes.popen(cmd, &xc, &output.stream, 1.0E-9);
    verbose("r: %d xc: %d output:\n%s", r, xc, data);
    #endif
}

#pragma pop_macro("verbose")

#else

static void rt_processes_test(void) { }

#endif

rt_processes_if rt_processes = {
    .pid                 = rt_processes_pid,
    .pids                = rt_processes_pids,
    .nameof              = rt_processes_nameof,
    .present             = rt_processes_present,
    .kill                = rt_processes_kill,
    .kill_all            = rt_processes_kill_all,
    .is_elevated         = rt_processes_is_elevated,
    .restart_elevated    = rt_processes_restart_elevated,
    .run                 = rt_processes_run,
    .popen               = rt_processes_open,
    .spawn               = rt_processes_spawn,
    .name                = rt_processes_name,
    .test                = rt_processes_test
};

// _______________________________ rt_static.c ________________________________

static void*   rt_static_symbol_reference[1024];
static int32_t rt_static_symbol_reference_count;

void* rt_force_symbol_reference(void* symbol) {
    rt_assert(rt_static_symbol_reference_count <= rt_countof(rt_static_symbol_reference),
        "increase size of rt_static_symbol_reference[%d] to at least %d",
        rt_countof(rt_static_symbol_reference), rt_static_symbol_reference);
    if (rt_static_symbol_reference_count < rt_countof(rt_static_symbol_reference)) {
        rt_static_symbol_reference[rt_static_symbol_reference_count] = symbol;
//      rt_println("rt_static_symbol_reference[%d] = %p", rt_static_symbol_reference_count,
//               rt_static_symbol_reference[symbol_reference_count]);
        rt_static_symbol_reference_count++;
    }
    return symbol;
}

// test rt_static_init() { code } that will be executed in random
// order but before main()

#ifdef RT_TESTS

static int32_t rt_static_init_function_called;

static void rt_force_inline rt_static_init_function(void) {
    rt_static_init_function_called = 1;
}

rt_static_init(static_init_test) { rt_static_init_function(); }

void rt_static_init_test(void) {
    rt_fatal_if(rt_static_init_function_called != 1,
        "static_init_function() expected to be called before main()");
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

void rt_static_init_test(void) {}

#endif

// _________________________________ rt_str.c _________________________________

static char* rt_str_drop_const(const char* s) {
    #if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wcast-qual"
    #endif
    return (char*)s;
    #if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic pop
    #endif
}

static int32_t rt_str_len(const char* s) { return (int32_t)strlen(s); }

static int32_t rt_str_utf16len(const uint16_t* utf16) {
    return (int32_t)wcslen(utf16);
}

static int32_t rt_str_utf8bytes(const char* s, int32_t b) {
    rt_assert(b >= 1, "should not be called with bytes < 1");
    const uint8_t* const u = (const uint8_t*)s;
    // based on:
    // https://stackoverflow.com/questions/66715611/check-for-valid-utf-8-encoding-in-c
    if (b >= 1 && (u[0] & 0x80u) == 0x00u) {
        return 1;
    } else if (b > 1) {
        uint32_t c = (u[0] << 8) | u[1];
        // TODO: 0xC080 is a hack - consider removing
        if (c == 0xC080) { return 2; } // 0xC080 as not zero terminating '\0'
        if (0xC280 <= c && c <= 0xDFBF && (c & 0xE0C0) == 0xC080) { return 2; }
        if (b > 2) {
            c = (c << 8) | u[2];
            // reject utf16 surrogates:
            if (0xEDA080 <= c && c <= 0xEDBFBF) { return 0; }
            if (0xE0A080 <= c && c <= 0xEFBFBF && (c & 0xF0C0C0) == 0xE08080) {
                return 3;
            }
            if (b > 3) {
                c = (c << 8) | u[3];
                if (0xF0908080 <= c && c <= 0xF48FBFBF &&
                    (c & 0xF8C0C0C0) == 0xF0808080) {
                    return 4;
                }
            }
        }
    }
    return 0; // invalid utf8 sequence
}

static int32_t rt_str_glyphs(const char* utf8, int32_t bytes) {
    rt_swear(bytes >= 0);
    bool ok = true;
    int32_t i = 0;
    int32_t k = 1;
    while (i < bytes && ok) {
        const int32_t b = rt_str.utf8bytes(utf8 + i, bytes - i);
        ok = 0 < b && i + b <= bytes;
        if (ok) { i += b; k++; }
    }
    return ok ? k - 1 : -1;
}

static void rt_str_lower(char* d, int32_t capacity, const char* s) {
    int32_t n = rt_str.len(s);
    rt_swear(capacity > n);
    for (int32_t i = 0; i < n; i++) { d[i] = (char)tolower(s[i]); }
    d[n] = 0;
}

static void rt_str_upper(char* d, int32_t capacity, const char* s) {
    int32_t n = rt_str.len(s);
    rt_swear(capacity > n);
    for (int32_t i = 0; i < n; i++) { d[i] = (char)toupper(s[i]); }
    d[n] = 0;
}

static bool rt_str_starts(const char* s1, const char* s2) {
    int32_t n1 = (int32_t)strlen(s1);
    int32_t n2 = (int32_t)strlen(s2);
    return n1 >= n2 && memcmp(s1, s2, n2) == 0;
}

static bool rt_str_ends(const char* s1, const char* s2) {
    int32_t n1 = (int32_t)strlen(s1);
    int32_t n2 = (int32_t)strlen(s2);
    return n1 >= n2 && memcmp(s1 + n1 - n2, s2, n2) == 0;
}

static bool rt_str_i_starts(const char* s1, const char* s2) {
    int32_t n1 = (int32_t)strlen(s1);
    int32_t n2 = (int32_t)strlen(s2);
    return n1 >= n2 && strnicmp(s1, s2, n2) == 0;

}

static bool rt_str_i_ends(const char* s1, const char* s2) {
    int32_t n1 = (int32_t)strlen(s1);
    int32_t n2 = (int32_t)strlen(s2);
    return n1 >= n2 && strnicmp(s1 + n1 - n2, s2, n2) == 0;
}

static int32_t rt_str_utf8_bytes(const uint16_t* utf16, int32_t chars) {
    // If `chars` argument is -1, the function utf8_bytes includes the zero
    // terminating character in the conversion and the returned byte count.
    // Function will fail (return 0) on incomplete surrogate pairs like
    // 0xD83D without following 0xDC1E https://compart.com/en/unicode/U+1F41E
    if (chars == 0) { return 0; }
    if (chars < 0 && utf16[0] == 0x0000) { return 1; }
    const int32_t required_bytes_count =
        WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
        utf16, chars, null, 0, null, null);
    if (required_bytes_count == 0) {
        errno_t r = rt_core.err();
        rt_println("WideCharToMultiByte() failed %s", rt_strerr(r));
        rt_core.set_err(r);
    }
    return required_bytes_count == 0 ? -1 : required_bytes_count;
}

static int32_t rt_str_utf16_chars(const char* utf8, int32_t bytes) {
    // If `bytes` argument is -1, the function utf16_chars() includes the zero
    // terminating character in the conversion and the returned character count.
    if (bytes == 0) { return 0; }
    if (bytes < 0 && utf8[0] == 0x00) { return 1; }
    const int32_t required_wide_chars_count =
        MultiByteToWideChar(CP_UTF8, 0, utf8, bytes, null, 0);
    if (required_wide_chars_count == 0) {
        errno_t r = rt_core.err();
        rt_println("MultiByteToWideChar() failed %s", rt_strerr(r));
        rt_core.set_err(r);
    }
    return required_wide_chars_count == 0 ? -1 : required_wide_chars_count;
}

static errno_t rt_str_utf16to8(char* utf8, int32_t capacity,
        const uint16_t* utf16, int32_t chars) {
    if (chars == 0) { return 0; }
    if (chars < 0 && utf16[0] == 0x0000) {
        rt_swear(capacity >= 1);
        utf8[0] = 0x00;
        return 0;
    }
    const int32_t required = rt_str.utf8_bytes(utf16, chars);
    errno_t r = required < 0 ? rt_core.err() : 0;
    if (r == 0) {
        rt_swear(required > 0 && capacity >= required);
        int32_t bytes = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                            utf16, chars, utf8, capacity, null, null);
        rt_swear(required == bytes);
    }
    return r;
}

static errno_t rt_str_utf8to16(uint16_t* utf16, int32_t capacity,
        const char* utf8, int32_t bytes) {
    const int32_t required = rt_str.utf16_chars(utf8, bytes);
    errno_t r = required < 0 ? rt_core.err() : 0;
    if (r == 0) {
        rt_swear(required >= 0 && capacity >= required);
        int32_t count = MultiByteToWideChar(CP_UTF8, 0, utf8, bytes,
                                            utf16, capacity);
        rt_swear(required == count);
#if 0 // TODO: incorrect need output != input
        if (count > 0 && !IsNormalizedString(NormalizationC, utf16, count)) {
            rt_core.set_err(0);
            int32_t n = NormalizeString(NormalizationC, utf16, count, utf16, count);
            if (n <= 0) {
                r = rt_core.err();
                rt_println("NormalizeString() failed %s", rt_strerr(r));
            }
        }
#endif 
    }
    return r;
}

static bool rt_str_utf16_is_low_surrogate(uint16_t utf16char) {
    return 0xDC00 <= utf16char && utf16char <= 0xDFFF;
}

static bool rt_str_utf16_is_high_surrogate(uint16_t utf16char) {
    return 0xD800 <= utf16char && utf16char <= 0xDBFF;
}

static uint32_t rt_str_utf32(const char* utf8, int32_t bytes) {
    uint32_t utf32 = 0;
    if ((utf8[0] & 0x80) == 0) {
        utf32 = utf8[0];
        rt_swear(bytes == 1);
    } else if ((utf8[0] & 0xE0) == 0xC0) {
        utf32  = (utf8[0] & 0x1F) << 6;
        utf32 |= (utf8[1] & 0x3F);
        rt_swear(bytes == 2);
    } else if ((utf8[0] & 0xF0) == 0xE0) {
        utf32  = (utf8[0] & 0x0F) << 12;
        utf32 |= (utf8[1] & 0x3F) <<  6;
        utf32 |= (utf8[2] & 0x3F);
        rt_swear(bytes == 3);
    } else if ((utf8[0] & 0xF8) == 0xF0) {
        utf32  = (utf8[0] & 0x07) << 18;
        utf32 |= (utf8[1] & 0x3F) << 12;
        utf32 |= (utf8[2] & 0x3F) <<  6;
        utf32 |= (utf8[3] & 0x3F);
        rt_swear(bytes == 4);
    } else {
        rt_swear(false);
    }
    return utf32;
}

static void rt_str_format_va(char* utf8, int32_t count, const char* format,
        va_list va) {
    #if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wformat-nonliteral"
    #endif
    vsnprintf(utf8, (size_t)count, format, va);
    utf8[count - 1] = 0;
    #if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic pop
    #endif
}

static void rt_str_format(char* utf8, int32_t count, const char* format, ...) {
    va_list va;
    va_start(va, format);
    rt_str.format_va(utf8, count, format, va);
    va_end(va);
}

static rt_str1024_t rt_str_error_for_language(int32_t error, LANGID language) {
    DWORD flags = FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS;
    HMODULE module = null;
    HRESULT hr = 0 <= error && error <= 0xFFFF ?
        HRESULT_FROM_WIN32((uint32_t)error) : (HRESULT)error;
    if ((error & 0xC0000000U) == 0xC0000000U) {
        // https://stackoverflow.com/questions/25566234/how-to-convert-specific-ntstatus-value-to-the-hresult
        static HMODULE ntdll; // RtlNtStatusToDosError implies linking to ntdll
        if (ntdll == null) { ntdll = GetModuleHandleA("ntdll.dll"); }
        if (ntdll == null) { ntdll = LoadLibraryA("ntdll.dll"); }
        module = ntdll;
        hr = HRESULT_FROM_WIN32(RtlNtStatusToDosError((NTSTATUS)error));
        flags |= FORMAT_MESSAGE_FROM_HMODULE;
    }
    rt_str1024_t text;
    uint16_t utf16[rt_countof(text.s)];
    DWORD count = FormatMessageW(flags, module, hr, language,
            utf16, rt_countof(utf16) - 1, (va_list*)null);
    utf16[rt_countof(utf16) - 1] = 0; // always
    // If FormatMessageW() succeeds, the return value is the number of utf16
    // characters stored in the output buffer, excluding the terminating zero.
    if (count > 0) {
        rt_swear(count < rt_countof(utf16));
        utf16[count] = 0;
        // remove trailing '\r\n'
        int32_t k = count;
        if (k > 0 && utf16[k - 1] == '\n') { utf16[k - 1] = 0; }
        k = (int32_t)rt_str.len16(utf16);
        if (k > 0 && utf16[k - 1] == '\r') { utf16[k - 1] = 0; }
        char message[rt_countof(text.s)];
        const int32_t bytes = rt_str.utf8_bytes(utf16, -1);
        if (bytes >= rt_countof(message)) {
            rt_str_printf(message, "error message is too long: %d bytes", bytes);
        } else {
            rt_str.utf16to8(message, rt_countof(message), utf16, -1);
        }
        // truncating printf to string:
        rt_str_printf(text.s, "0x%08X(%d) \"%s\"", error, error, message);
    } else {
        rt_str_printf(text.s, "0x%08X(%d)", error, error);
    }
    return text;
}

static rt_str1024_t rt_str_error(int32_t error) {
    const LANGID language = MAKELANGID(LANG_ENGLISH, SUBLANG_DEFAULT);
    return rt_str_error_for_language(error, language);
}

static rt_str1024_t rt_str_error_nls(int32_t error) {
    const LANGID language = MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT);
    return rt_str_error_for_language(error, language);
}

static const char* rt_str_grouping_separator(void) {
    #ifdef WINDOWS
        // en-US Windows 10/11:
        // grouping_separator == ","
        // decimal_separator  == "."
        static char grouping_separator[8];
        if (grouping_separator[0] == 0x00) {
            errno_t r = rt_b2e(GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_STHOUSAND,
                grouping_separator, sizeof(grouping_separator)));
            rt_swear(r == 0 && grouping_separator[0] != 0);
        }
        return grouping_separator;
    #else
        // en-US Windows 10/11:
        // grouping_separator == ""
        // decimal_separator  == "."
        struct lconv *locale_info = localeconv();
        const char* grouping_separator = null;
        if (grouping_separator == null) {
            grouping_separator = locale_info->thousands_sep;
            swear(grouping_separator != null);
        }
        return grouping_separator;
    #endif
}

// Posix and Win32 C runtime:
//   #include <locale.h>
//   struct lconv *locale_info = localeconv();
//   const char* grouping_separator = locale_info->thousands_sep;
//   const char* decimal_separator = locale_info->decimal_point;
// en-US Windows 1x:
// grouping_separator == ""
// decimal_separator  == "."
//
// Win32 API:
//   rt_b2e(GetLocaleInfo(LOCALE_USER_DEFAULT, LOCALE_STHOUSAND,
//       grouping_separator, sizeof(grouping_separator)));
//   rt_b2e(GetLocaleInfo(LOCALE_USER_DEFAULT, LOCALE_SDECIMAL,
//       decimal_separator, sizeof(decimal_separator)));
// en-US Windows 1x:
// grouping_separator == ","
// decimal_separator  == "."

static rt_str64_t rt_str_int64_dg(int64_t v, // digit_grouped
        bool uint, const char* gs) { // grouping separator: gs
    // sprintf format %`lld may not be implemented or
    // does not respect locale or UI separators...
    // Do it hard way:
    const int32_t m = (int32_t)strlen(gs);
    rt_swear(m < 5); // utf-8 4 bytes max
    // 64 calls per thread 32 or less bytes each because:
    // "18446744073709551615" 21 characters + 6x4 groups:
    // "18'446'744'073'709'551'615" 27 characters
    rt_str64_t text;
    enum { max_text_bytes = rt_countof(text.s) };
    int64_t abs64 = v < 0 ? -v : v; // incorrect for INT64_MIN
    uint64_t n = uint ? (uint64_t)v :
        (v != INT64_MIN ? (uint64_t)abs64 : (uint64_t)INT64_MIN);
    int32_t i = 0;
    int32_t groups[8]; // 2^63 - 1 ~= 9 x 10^19 upto 7 groups of 3 digits
    do {
        groups[i] = n % 1000;
        n = n / 1000;
        i++;
    } while (n > 0);
    const int32_t gc = i - 1; // group count
    char* s = text.s;
    if (v < 0 && !uint) { *s++ = '-'; } // sign
    int32_t r = max_text_bytes - 1;
    while (i > 0) {
        i--;
        rt_assert(r > 3 + m);
        if (i == gc) {
            rt_str.format(s, r, "%d%s", groups[i], gc > 0 ? gs : "");
        } else {
            rt_str.format(s, r, "%03d%s", groups[i], i > 0 ? gs : "");
        }
        int32_t k = (int32_t)strlen(s);
        r -= k;
        s += k;
    }
    *s = 0;
    return text;
}

static rt_str64_t rt_str_int64(int64_t v) {
    return rt_str_int64_dg(v, false, rt_glyph_hair_space);
}

static rt_str64_t rt_str_uint64(uint64_t v) {
    return rt_str_int64_dg(v, true, rt_glyph_hair_space);
}

static rt_str64_t rt_str_int64_lc(int64_t v) {
    return rt_str_int64_dg(v, false, rt_str_grouping_separator());
}

static rt_str64_t rt_str_uint64_lc(uint64_t v) {
    return rt_str_int64_dg(v, true, rt_str_grouping_separator());
}

static rt_str128_t rt_str_fp(const char* format, fp64_t v) {
    static char decimal_separator[8];
    if (decimal_separator[0] == 0) {
        errno_t r = rt_b2e(GetLocaleInfoA(LOCALE_USER_DEFAULT, LOCALE_SDECIMAL,
            decimal_separator, sizeof(decimal_separator)));
        rt_swear(r == 0 && decimal_separator[0] != 0);
    }
    rt_swear(strlen(decimal_separator) <= 4);
    rt_str128_t f; // formatted float point
    // snprintf format does not handle thousands separators on all know runtimes
    // and respects setlocale() on Un*x systems but in MS runtime only when
    // _snprintf_l() is used.
    f.s[0] = 0x00;
    rt_str.format(f.s, rt_countof(f.s), format, v);
    f.s[rt_countof(f.s) - 1] = 0x00;
    rt_str128_t text;
    char* s = f.s;
    char* d = text.s;
    while (*s != 0x00) {
        if (*s == '.') {
            const char* sep = decimal_separator;
            while (*sep != 0x00) { *d++ = *sep++; }
            s++;
        } else {
            *d++ = *s++;
        }
    }
    *d = 0x00;
    // TODO: It's possible to handle mantissa grouping but...
    // Not clear if human expects it in 5 digits or 3 digits chunks
    // and unfortunately locale does not specify how
    return text;
}

#ifdef RT_TESTS

static void rt_str_test(void) {
    rt_swear(rt_str.len("hello") == 5);
    rt_swear(rt_str.starts("hello world", "hello"));
    rt_swear(rt_str.ends("hello world", "world"));
    rt_swear(rt_str.istarts("hello world", "HeLlO"));
    rt_swear(rt_str.iends("hello world", "WoRlD"));
    char ls[20] = {0};
    rt_str.lower(ls, rt_countof(ls), "HeLlO WoRlD");
    rt_swear(strcmp(ls, "hello world") == 0);
    char upper[11] = {0};
    rt_str.upper(upper, rt_countof(upper), "hello12345");
    rt_swear(strcmp(upper,  "HELLO12345") == 0);
    #pragma push_macro("glyph_chinese_one")
    #pragma push_macro("glyph_chinese_two")
    #pragma push_macro("glyph_teddy_bear")
    #pragma push_macro("glyph_ice_cube")
    #define glyph_chinese_one "\xE5\xA3\xB9"
    #define glyph_chinese_two "\xE8\xB4\xB0"
    #define glyph_teddy_bear  "\xF0\x9F\xA7\xB8"
    #define glyph_ice_cube    "\xF0\x9F\xA7\x8A"
    const char* utf8_str =
            glyph_teddy_bear
            "0"
            rt_glyph_chinese_jin4 rt_glyph_chinese_gong
            "3456789 "
            glyph_ice_cube;
    rt_swear(rt_str.utf8bytes("\x01", 1) == 1);
    rt_swear(rt_str.utf8bytes("\x7F", 1) == 1);
    rt_swear(rt_str.utf8bytes("\x80", 1) == 0);
//  swear(rt_str.utf8bytes(glyph_chinese_one, 0) == 0);
    rt_swear(rt_str.utf8bytes(glyph_chinese_one, 1) == 0);
    rt_swear(rt_str.utf8bytes(glyph_chinese_one, 2) == 0);
    rt_swear(rt_str.utf8bytes(glyph_chinese_one, 3) == 3);
    rt_swear(rt_str.utf8bytes(glyph_teddy_bear,  4) == 4);
    #pragma pop_macro("glyph_ice_cube")
    #pragma pop_macro("glyph_teddy_bear")
    #pragma pop_macro("glyph_chinese_two")
    #pragma pop_macro("glyph_chinese_one")
    uint16_t wide_str[100] = {0};
    rt_str.utf8to16(wide_str, rt_countof(wide_str), utf8_str, -1);
    char utf8[100] = {0};
    rt_str.utf16to8(utf8, rt_countof(utf8), wide_str, -1);
    uint16_t utf16[100];
    rt_str.utf8to16(utf16, rt_countof(utf16), utf8, -1);
    char narrow_str[100] = {0};
    rt_str.utf16to8(narrow_str, rt_countof(narrow_str), utf16, -1);
    rt_swear(strcmp(narrow_str, utf8_str) == 0);
    char formatted[100];
    rt_str.format(formatted, rt_countof(formatted), "n: %d, s: %s", 42, "test");
    rt_swear(strcmp(formatted, "n: 42, s: test") == 0);
    // numeric values digit grouping format:
    rt_swear(strcmp("0", rt_str.int64_dg(0, true, ",").s) == 0);
    rt_swear(strcmp("-1", rt_str.int64_dg(-1, false, ",").s) == 0);
    rt_swear(strcmp("999", rt_str.int64_dg(999, true, ",").s) == 0);
    rt_swear(strcmp("-999", rt_str.int64_dg(-999, false, ",").s) == 0);
    rt_swear(strcmp("1,001", rt_str.int64_dg(1001, true, ",").s) == 0);
    rt_swear(strcmp("-1,001", rt_str.int64_dg(-1001, false, ",").s) == 0);
    rt_swear(strcmp("18,446,744,073,709,551,615",
        rt_str.int64_dg(UINT64_MAX, true, ",").s) == 0
    );
    rt_swear(strcmp("9,223,372,036,854,775,807",
        rt_str.int64_dg(INT64_MAX, false, ",").s) == 0
    );
    rt_swear(strcmp("-9,223,372,036,854,775,808",
        rt_str.int64_dg(INT64_MIN, false, ",").s) == 0
    );
    //  see:
    // https://en.wikipedia.org/wiki/Single-precision_floating-point_format
    uint32_t pi_fp32 = 0x40490FDBULL; // 3.14159274101257324
    rt_swear(strcmp("3.141592741",
                rt_str.fp("%.9f", *(fp32_t*)&pi_fp32).s) == 0,
          "%s", rt_str.fp("%.9f", *(fp32_t*)&pi_fp32).s
    );
    //  3.141592741
    //  ********^ (*** true digits ^ first rounded digit)
    //    123456 (%.6f)
    //
    //  https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    uint64_t pi_fp64 = 0x400921FB54442D18ULL;
    rt_swear(strcmp("3.141592653589793116",
                rt_str.fp("%.18f", *(fp64_t*)&pi_fp64).s) == 0,
          "%s", rt_str.fp("%.18f", *(fp64_t*)&pi_fp64).s
    );
    //  3.141592653589793116
    //  *****************^ (*** true digits ^ first rounded digit)
    //    123456789012345 (%.15f)
    //  https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    //
    //  actual "pi" first 64 digits:
    //  3.1415926535897932384626433832795028841971693993751058209749445923
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_str_test(void) {}

#endif

rt_str_if rt_str = {
    .drop_const              = rt_str_drop_const,
    .len                     = rt_str_len,
    .len16                   = rt_str_utf16len,
    .utf8bytes               = rt_str_utf8bytes,
    .glyphs                  = rt_str_glyphs,
    .lower                   = rt_str_lower,
    .upper                   = rt_str_upper,
    .starts                  = rt_str_starts,
    .ends                    = rt_str_ends,
    .istarts                 = rt_str_i_starts,
    .iends                   = rt_str_i_ends,
    .utf8_bytes              = rt_str_utf8_bytes,
    .utf16_chars             = rt_str_utf16_chars,
    .utf16to8                = rt_str_utf16to8,
    .utf8to16                = rt_str_utf8to16,
    .utf16_is_low_surrogate  = rt_str_utf16_is_low_surrogate,
    .utf16_is_high_surrogate = rt_str_utf16_is_high_surrogate,
    .utf32                   = rt_str_utf32,
    .format                  = rt_str_format,
    .format_va               = rt_str_format_va,
    .error                   = rt_str_error,
    .error_nls               = rt_str_error_nls,
    .grouping_separator      = rt_str_grouping_separator,
    .int64_dg                = rt_str_int64_dg,
    .int64                   = rt_str_int64,
    .uint64                  = rt_str_uint64,
    .int64_lc                = rt_str_int64,
    .uint64_lc               = rt_str_uint64,
    .fp                      = rt_str_fp,
    .test                    = rt_str_test
};

// _______________________________ rt_streams.c _______________________________

static errno_t rt_streams_memory_read(rt_stream_if* stream, void* data, int64_t bytes,
        int64_t *transferred) {
    rt_swear(bytes > 0);
    rt_stream_memory_if* s = (rt_stream_memory_if*)stream;
    rt_swear(0 <= s->pos_read && s->pos_read <= s->bytes_read,
          "bytes: %lld stream .pos: %lld .bytes: %lld",
          bytes, s->pos_read, s->bytes_read);
    int64_t transfer = rt_min(bytes, s->bytes_read - s->pos_read);
    memcpy(data, (const uint8_t*)s->data_read + s->pos_read, (size_t)transfer);
    s->pos_read += transfer;
    if (transferred != null) { *transferred = transfer; }
    return 0;
}

static errno_t rt_streams_memory_write(rt_stream_if* stream, const void* data, int64_t bytes,
        int64_t *transferred) {
    rt_swear(bytes > 0);
    rt_stream_memory_if* s = (rt_stream_memory_if*)stream;
    rt_swear(0 <= s->pos_write && s->pos_write <= s->bytes_write,
          "bytes: %lld stream .pos: %lld .bytes: %lld",
          bytes, s->pos_write, s->bytes_write);
    bool overflow = s->bytes_write - s->pos_write <= 0;
    int64_t transfer = rt_min(bytes, s->bytes_write - s->pos_write);
    memcpy((uint8_t*)s->data_write + s->pos_write, data, (size_t)transfer);
    s->pos_write += transfer;
    if (transferred != null) { *transferred = transfer; }
    return overflow ? ERROR_INSUFFICIENT_BUFFER : 0;
}

static void rt_streams_read_only(rt_stream_memory_if* s,
        const void* data, int64_t bytes) {
    s->stream.read = rt_streams_memory_read;
    s->stream.write = null;
    s->data_read = data;
    s->bytes_read = bytes;
    s->pos_read = 0;
    s->data_write = null;
    s->bytes_write = 0;
    s->pos_write = 0;
}

static void rt_streams_write_only(rt_stream_memory_if* s,
        void* data, int64_t bytes) {
    s->stream.read = null;
    s->stream.write = rt_streams_memory_write;
    s->data_read = null;
    s->bytes_read = 0;
    s->pos_read = 0;
    s->data_write = data;
    s->bytes_write = bytes;
    s->pos_write = 0;
}

static void rt_streams_read_write(rt_stream_memory_if* s,
        const void* read, int64_t read_bytes,
        void* write, int64_t write_bytes) {
    s->stream.read = rt_streams_memory_read;
    s->stream.write = rt_streams_memory_write;
    s->data_read = read;
    s->bytes_read = read_bytes;
    s->pos_read = 0;
    s->pos_read = 0;
    s->data_write = write;
    s->bytes_write = write_bytes;
    s->pos_write = 0;
}

#ifdef RT_TESTS

static void rt_streams_test(void) {
    {   // read test
        uint8_t memory[256];
        for (int32_t i = 0; i < rt_countof(memory); i++) { memory[i] = (uint8_t)i; }
        for (int32_t i = 1; i < rt_countof(memory) - 1; i++) {
            rt_stream_memory_if ms; // memory stream
            rt_streams.read_only(&ms, memory, sizeof(memory));
            uint8_t data[256];
            for (int32_t j = 0; j < rt_countof(data); j++) { data[j] = 0xFF; }
            int64_t transferred = 0;
            errno_t r = ms.stream.read(&ms.stream, data, i, &transferred);
            rt_swear(r == 0 && transferred == i);
            for (int32_t j = 0; j < i; j++) { rt_swear(data[j] == memory[j]); }
            for (int32_t j = i; j < rt_countof(data); j++) { rt_swear(data[j] == 0xFF); }
        }
    }
    {   // write test
        // TODO: implement
    }
    {   // read/write test
        // TODO: implement
    }
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_streams_test(void) { }

#endif

rt_streams_if rt_streams = {
    .read_only  = rt_streams_read_only,
    .write_only = rt_streams_write_only,
    .read_write = rt_streams_read_write,
    .test       = rt_streams_test
};

// _______________________________ rt_threads.c _______________________________

// events:

static rt_event_t rt_event_create(void) {
    HANDLE e = CreateEvent(null, false, false, null);
    rt_not_null(e);
    return (rt_event_t)e;
}

static rt_event_t rt_event_create_manual(void) {
    HANDLE e = CreateEvent(null, true, false, null);
    rt_not_null(e);
    return (rt_event_t)e;
}

static void rt_event_set(rt_event_t e) {
    rt_fatal_win32err(SetEvent((HANDLE)e));
}

static void rt_event_reset(rt_event_t e) {
    rt_fatal_win32err(ResetEvent((HANDLE)e));
}

static int32_t rt_event_wait_or_timeout(rt_event_t e, fp64_t seconds) {
    uint32_t ms = seconds < 0 ? INFINITE : (uint32_t)(seconds * 1000.0 + 0.5);
    DWORD i = WaitForSingleObject(e, ms);
    rt_swear(i != WAIT_FAILED, "i: %d", i);
    errno_t r = rt_wait_ix2e(i);
    if (r != 0) { rt_swear(i == WAIT_TIMEOUT || i == WAIT_ABANDONED); }
    return i == WAIT_TIMEOUT ? -1 : (i == WAIT_ABANDONED ? -2 : i);
}

static void rt_event_wait(rt_event_t e) { rt_event_wait_or_timeout(e, -1); }

static int32_t rt_event_wait_any_or_timeout(int32_t n,
        rt_event_t events[], fp64_t s) {
    rt_swear(n < 64); // Win32 API limit
    const uint32_t ms = s < 0 ? INFINITE : (uint32_t)(s * 1000.0 + 0.5);
    const HANDLE* es = (const HANDLE*)events;
    DWORD i = WaitForMultipleObjects((DWORD)n, es, false, ms);
    rt_swear(i != WAIT_FAILED, "i: %d", i);
    errno_t r = rt_wait_ix2e(i);
    if (r != 0) { rt_swear(i == WAIT_TIMEOUT || i == WAIT_ABANDONED); }
    return i == WAIT_TIMEOUT ? -1 : (i == WAIT_ABANDONED ? -2 : i);
}

static int32_t rt_event_wait_any(int32_t n, rt_event_t e[]) {
    return rt_event_wait_any_or_timeout(n, e, -1);
}

static void rt_event_dispose(rt_event_t h) {
    rt_win32_close_handle(h);
}

#ifdef RT_TESTS

// test:

// check if the elapsed time is within the expected range
static void rt_event_test_check_time(fp64_t start, fp64_t expected) {
    fp64_t elapsed = rt_clock.seconds() - start;
    // Old Windows scheduler is prone to 2x16.6ms ~ 33ms delays (observed)
    rt_swear(elapsed >= expected - 0.04 && elapsed <= expected + 0.250,
          "expected: %f elapsed %f seconds", expected, elapsed);
}

static void rt_event_test(void) {
    rt_event_t event = rt_event.create();
    fp64_t start = rt_clock.seconds();
    rt_event.set(event);
    rt_event.wait(event);
    rt_event_test_check_time(start, 0); // Event should be immediate
    rt_event.reset(event);
    start = rt_clock.seconds();
    const fp64_t timeout_seconds = 1.0 / 8.0;
    int32_t result = rt_event.wait_or_timeout(event, timeout_seconds);
    rt_event_test_check_time(start, timeout_seconds);
    rt_swear(result == -1); // Timeout expected
    enum { count = 5 };
    rt_event_t events[count];
    for (int32_t i = 0; i < rt_countof(events); i++) {
        events[i] = rt_event.create_manual();
    }
    start = rt_clock.seconds();
    rt_event.set(events[2]); // Set the third event
    int32_t index = rt_event.wait_any(rt_countof(events), events);
    rt_swear(index == 2);
    rt_event_test_check_time(start, 0);
    rt_swear(index == 2); // Third event should be triggered
    rt_event.reset(events[2]); // Reset the third event
    start = rt_clock.seconds();
    result = rt_event.wait_any_or_timeout(rt_countof(events), events, timeout_seconds);
    rt_swear(result == -1);
    rt_event_test_check_time(start, timeout_seconds);
    rt_swear(result == -1); // Timeout expected
    // Clean up
    rt_event.dispose(event);
    for (int32_t i = 0; i < rt_countof(events); i++) {
        rt_event.dispose(events[i]);
    }
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_event_test(void) { }

#endif

rt_event_if rt_event = {
    .create              = rt_event_create,
    .create_manual       = rt_event_create_manual,
    .set                 = rt_event_set,
    .reset               = rt_event_reset,
    .wait                = rt_event_wait,
    .wait_or_timeout     = rt_event_wait_or_timeout,
    .wait_any            = rt_event_wait_any,
    .wait_any_or_timeout = rt_event_wait_any_or_timeout,
    .dispose             = rt_event_dispose,
    .test                = rt_event_test
};

// mutexes:

rt_static_assertion(sizeof(CRITICAL_SECTION) == sizeof(rt_mutex_t));

static void rt_mutex_init(rt_mutex_t* m) {
    CRITICAL_SECTION* cs = (CRITICAL_SECTION*)m;
    rt_fatal_win32err(InitializeCriticalSectionAndSpinCount(cs, 4096));
}

static void rt_mutex_lock(rt_mutex_t* m) {
    EnterCriticalSection((CRITICAL_SECTION*)m);
}

static void rt_mutex_unlock(rt_mutex_t* m) {
    LeaveCriticalSection((CRITICAL_SECTION*)m);
}

static void rt_mutex_dispose(rt_mutex_t* m) {
    DeleteCriticalSection((CRITICAL_SECTION*)m);
}

// test:

// check if the elapsed time is within the expected range
static void rt_mutex_test_check_time(fp64_t start, fp64_t expected) {
    fp64_t elapsed = rt_clock.seconds() - start;
    // Old Windows scheduler is prone to 2x16.6ms ~ 33ms delays
    rt_swear(elapsed >= expected - 0.04 && elapsed <= expected + 0.04,
          "expected: %f elapsed %f seconds", expected, elapsed);
}

static void rt_mutex_test_lock_unlock(void* arg) {
    rt_mutex_t* mutex = (rt_mutex_t*)arg;
    rt_mutex.lock(mutex);
    rt_thread.sleep_for(0.01); // Hold the mutex for 10ms
    rt_mutex.unlock(mutex);
}

static void rt_mutex_test(void) {
    rt_mutex_t mutex;
    rt_mutex.init(&mutex);
    fp64_t start = rt_clock.seconds();
    rt_mutex.lock(&mutex);
    rt_mutex.unlock(&mutex);
    // Lock and unlock should be immediate
    rt_mutex_test_check_time(start, 0);
    enum { count = 5 };
    rt_thread_t ts[count];
    for (int32_t i = 0; i < rt_countof(ts); i++) {
        ts[i] = rt_thread.start(rt_mutex_test_lock_unlock, &mutex);
    }
    // Wait for all threads to finish
    for (int32_t i = 0; i < rt_countof(ts); i++) {
        rt_thread.join(ts[i], -1);
    }
    rt_mutex.dispose(&mutex);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

rt_mutex_if rt_mutex = {
    .init    = rt_mutex_init,
    .lock    = rt_mutex_lock,
    .unlock  = rt_mutex_unlock,
    .dispose = rt_mutex_dispose,
    .test    = rt_mutex_test
};

// threads:

static void* rt_thread_ntdll(void) {
    static HMODULE ntdll;
    if (ntdll == null) {
        ntdll = (void*)GetModuleHandleA("ntdll.dll");
    }
    if (ntdll == null) {
        ntdll = rt_loader.open("ntdll.dll", 0);
    }
    rt_not_null(ntdll);
    return ntdll;
}

static fp64_t rt_thread_ns2ms(int64_t ns) {
    return (fp64_t)ns / (fp64_t)rt_clock.nsec_in_msec;
}

static void rt_thread_set_timer_resolution(uint64_t nanoseconds) {
    typedef int32_t (*query_timer_resolution_t)(ULONG* minimum_resolution,
        ULONG* maximum_resolution, ULONG* actual_resolution);
    typedef int32_t (*set_timer_resolution_t)(ULONG requested_resolution,
        BOOLEAN set, ULONG* actual_resolution); // ntdll.dll
    void* nt_dll = rt_thread_ntdll();
    query_timer_resolution_t query_timer_resolution =  (query_timer_resolution_t)
        rt_loader.sym(nt_dll, "NtQueryTimerResolution");
    set_timer_resolution_t set_timer_resolution = (set_timer_resolution_t)
        rt_loader.sym(nt_dll, "NtSetTimerResolution");
    unsigned long min100ns = 16 * 10 * 1000;
    unsigned long max100ns =  1 * 10 * 1000;
    unsigned long cur100ns =  0;
    rt_fatal_if(query_timer_resolution(&min100ns, &max100ns, &cur100ns) != 0);
    uint64_t max_ns = max100ns * 100uLL;
//  uint64_t min_ns = min100ns * 100uLL;
//  uint64_t cur_ns = cur100ns * 100uLL;
    // max resolution is lowest possible delay between timer events
//  if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
//      rt_println("timer resolution min: %.3f max: %.3f cur: %.3f"
//          " ms (milliseconds)",
//          rt_thread_ns2ms(min_ns),
//          rt_thread_ns2ms(max_ns),
//          rt_thread_ns2ms(cur_ns));
//  }
    // note that maximum resolution is actually < minimum
    nanoseconds = rt_max(max_ns, nanoseconds);
    unsigned long ns = (unsigned long)((nanoseconds + 99) / 100);
    rt_fatal_if(set_timer_resolution(ns, true, &cur100ns) != 0);
    rt_fatal_if(query_timer_resolution(&min100ns, &max100ns, &cur100ns) != 0);
//  if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
//      min_ns = min100ns * 100uLL;
//      max_ns = max100ns * 100uLL; // the smallest interval
//      cur_ns = cur100ns * 100uLL;
//      rt_println("timer resolution min: %.3f max: %.3f cur: %.3f ms (milliseconds)",
//          rt_thread_ns2ms(min_ns),
//          rt_thread_ns2ms(max_ns),
//          rt_thread_ns2ms(cur_ns));
//  }
}

static void rt_thread_power_throttling_disable_for_process(void) {
    static bool disabled_for_the_process;
    if (!disabled_for_the_process) {
        PROCESS_POWER_THROTTLING_STATE pt = { 0 };
        pt.Version = PROCESS_POWER_THROTTLING_CURRENT_VERSION;
        pt.ControlMask = PROCESS_POWER_THROTTLING_EXECUTION_SPEED;
        pt.StateMask = 0;
        rt_fatal_win32err(SetProcessInformation(GetCurrentProcess(),
            ProcessPowerThrottling, &pt, sizeof(pt)));
        // PROCESS_POWER_THROTTLING_IGNORE_TIMER_RESOLUTION
        // does not work on Win10. There is no easy way to
        // distinguish Windows 11 from 10 (Microsoft great engineering)
        pt.ControlMask = PROCESS_POWER_THROTTLING_IGNORE_TIMER_RESOLUTION;
        pt.StateMask = 0;
        // ignore error on Windows 10:
        (void)SetProcessInformation(GetCurrentProcess(),
            ProcessPowerThrottling, &pt, sizeof(pt));
        disabled_for_the_process = true;
    }
}

static void rt_thread_power_throttling_disable_for_thread(HANDLE thread) {
    THREAD_POWER_THROTTLING_STATE pt = { 0 };
    pt.Version = THREAD_POWER_THROTTLING_CURRENT_VERSION;
    pt.ControlMask = THREAD_POWER_THROTTLING_EXECUTION_SPEED;
    pt.StateMask = 0;
    rt_fatal_win32err(SetThreadInformation(thread,
        ThreadPowerThrottling, &pt, sizeof(pt)));
}

static void rt_thread_disable_power_throttling(void) {
    rt_thread_power_throttling_disable_for_process();
    rt_thread_power_throttling_disable_for_thread(GetCurrentThread());
}

static const char* rt_thread_rel2str(int32_t rel) {
    switch (rel) {
        case RelationProcessorCore   : return "ProcessorCore   ";
        case RelationNumaNode        : return "NumaNode        ";
        case RelationCache           : return "Cache           ";
        case RelationProcessorPackage: return "ProcessorPackage";
        case RelationGroup           : return "Group           ";
        case RelationProcessorDie    : return "ProcessorDie    ";
        case RelationNumaNodeEx      : return "NumaNodeEx      ";
        case RelationProcessorModule : return "ProcessorModule ";
        default: rt_assert(false, "fix me"); return "???";
    }
}

static uint64_t rt_thread_next_physical_processor_affinity_mask(void) {
    static volatile int32_t initialized;
    static int32_t init;
    static int32_t next = 1; // next physical core to use
    static int32_t cores = 0; // number of physical processors (cores)
    static uint64_t any;
    static uint64_t affinity[64]; // mask for each physical processor
    bool set_to_true = rt_atomics.compare_exchange_int32(&init, false, true);
    if (set_to_true) {
        // Concept D: 6 cores, 12 logical processors: 27 lpi entries
        static SYSTEM_LOGICAL_PROCESSOR_INFORMATION lpi[64];
        DWORD bytes = 0;
        GetLogicalProcessorInformation(null, &bytes);
        rt_assert(bytes % sizeof(lpi[0]) == 0);
        // number of lpi entries == 27 on 6 core / 12 logical processors system
        int32_t n = bytes / sizeof(lpi[0]);
        rt_assert(bytes <= sizeof(lpi), "increase lpi[%d]", n);
        rt_fatal_win32err(GetLogicalProcessorInformation(&lpi[0], &bytes));
        for (int32_t i = 0; i < n; i++) {
//          if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
//              rt_println("[%2d] affinity mask 0x%016llX relationship=%d %s", i,
//                  lpi[i].ProcessorMask, lpi[i].Relationship,
//                  rt_thread_rel2str(lpi[i].Relationship));
//          }
            if (lpi[i].Relationship == RelationProcessorCore) {
                rt_assert(cores < rt_countof(affinity), "increase affinity[%d]", cores);
                if (cores < rt_countof(affinity)) {
                    any |= lpi[i].ProcessorMask;
                    affinity[cores] = lpi[i].ProcessorMask;
                    cores++;
                }
            }
        }
        initialized = true;
    } else {
        while (initialized == 0) { rt_thread.sleep_for(1 / 1024.0); }
        rt_assert(any != 0); // should not ever happen
        if (any == 0) { any = (uint64_t)(-1LL); }
    }
    uint64_t mask = next < cores ? affinity[next] : any;
    rt_assert(mask != 0);
    // assume last physical core is least popular
    if (next < cores) { next++; } // not circular
    return mask;
}

static void rt_thread_realtime(void) {
    rt_fatal_win32err(SetPriorityClass(GetCurrentProcess(),
        REALTIME_PRIORITY_CLASS));
    rt_fatal_win32err(SetThreadPriority(GetCurrentThread(),
        THREAD_PRIORITY_TIME_CRITICAL));
    rt_fatal_win32err(SetThreadPriorityBoost(GetCurrentThread(),
        /* bDisablePriorityBoost = */ false));
    // desired: 0.5ms = 500us (microsecond) = 50,000ns
    rt_thread_set_timer_resolution((uint64_t)rt_clock.nsec_in_usec * 500ULL);
    rt_fatal_win32err(SetThreadAffinityMask(GetCurrentThread(),
        rt_thread_next_physical_processor_affinity_mask()));
    rt_thread_disable_power_throttling();
}

static void rt_thread_yield(void) { SwitchToThread(); }

static rt_thread_t rt_thread_start(void (*func)(void*), void* p) {
    rt_thread_t t = (rt_thread_t)CreateThread(null, 0,
        (LPTHREAD_START_ROUTINE)(void*)func, p, 0, null);
    rt_not_null(t);
    return t;
}

static bool is_handle_valid(void* h) {
    DWORD flags = 0;
    return GetHandleInformation(h, &flags);
}

static errno_t rt_thread_join(rt_thread_t t, fp64_t timeout) {
    rt_not_null(t);
    rt_fatal_if(!is_handle_valid(t));
    const uint32_t ms = timeout < 0 ? INFINITE : (uint32_t)(timeout * 1000.0 + 0.5);
    DWORD ix = WaitForSingleObject(t, (DWORD)ms);
    errno_t r = rt_wait_ix2e(ix);
    rt_assert(r != ERROR_REQUEST_ABORTED, "AFAIK thread can`t be ABANDONED");
    if (r == 0) {
        rt_win32_close_handle(t);
    } else {
        rt_println("failed to join thread %p %s", t, rt_strerr(r));
    }
    return r;
}

static void rt_thread_detach(rt_thread_t t) {
    rt_not_null(t);
    rt_fatal_if(!is_handle_valid(t));
    rt_win32_close_handle(t);
}

static void rt_thread_name(const char* name) {
    uint16_t stack[128];
    rt_fatal_if(rt_str.len(name) >= rt_countof(stack), "name too long: %s", name);
    rt_str.utf8to16(stack, rt_countof(stack), name, -1);
    HRESULT r = SetThreadDescription(GetCurrentThread(), stack);
    // notoriously returns 0x10000000 for no good reason whatsoever
    rt_fatal_if(!SUCCEEDED(r));
}

static void rt_thread_sleep_for(fp64_t seconds) {
    rt_assert(seconds >= 0);
    if (seconds < 0) { seconds = 0; }
    int64_t ns100 = (int64_t)(seconds * 1.0e+7); // in 0.1 us aka 100ns
    typedef int32_t (__stdcall *nt_delay_execution_t)(BOOLEAN alertable,
        PLARGE_INTEGER DelayInterval);
    static nt_delay_execution_t NtDelayExecution;
    // delay in 100-ns units. negative value means delay relative to current.
    LARGE_INTEGER delay = {0}; // delay in 100-ns units.
    delay.QuadPart = -ns100; // negative value means delay relative to current.
    if (NtDelayExecution == null) {
        void* ntdll = rt_thread_ntdll();
        NtDelayExecution = (nt_delay_execution_t)
            rt_loader.sym(ntdll, "NtDelayExecution");
        rt_not_null(NtDelayExecution);
    }
    // If "alertable" is set, sleep_for() can break earlier
    // as a result of NtAlertThread call.
    NtDelayExecution(false, &delay);
}

static uint64_t rt_thread_id_of(rt_thread_t t) {
    return (uint64_t)GetThreadId((HANDLE)t);
}

static uint64_t rt_thread_id(void) {
    return (uint64_t)GetThreadId(GetCurrentThread());
}

static rt_thread_t rt_thread_self(void) {
    // GetCurrentThread() returns pseudo-handle, not a real handle
    // if real handle is ever needed may do
    // rt_thread_t t = rt_thread.open(rt_thread.id()) and
    // rt_thread.close(t) instead.
    return (rt_thread_t)GetCurrentThread();
}

static errno_t rt_thread_open(rt_thread_t *t, uint64_t id) {
    // GetCurrentThread() returns pseudo-handle, not a real handle.
    // if real handle is ever needed do rt_thread_id_of() instead
    // but don't forget to do rt_thread.close() after that.
    *t = (rt_thread_t)OpenThread(THREAD_ALL_ACCESS, false, (DWORD)id);
    return *t == null ? rt_core.err() : 0;
}

static void rt_thread_close(rt_thread_t t) {
    rt_not_null(t);
    rt_win32_close_handle((HANDLE)t);
}

#ifdef RT_TESTS

// test: https://en.wikipedia.org/wiki/Dining_philosophers_problem

typedef struct rt_thread_philosophers_s rt_thread_philosophers_t;

typedef struct {
    rt_thread_philosophers_t* ps;
    rt_mutex_t  fork;
    rt_mutex_t* left_fork;
    rt_mutex_t* right_fork;
    rt_thread_t thread;
    uint64_t    id;
} rt_thread_philosopher_t;

typedef struct rt_thread_philosophers_s {
    rt_thread_philosopher_t philosopher[3];
    rt_event_t fed_up[3];
    uint32_t seed;
    volatile bool enough;
} rt_thread_philosophers_t;

#pragma push_macro("verbose") // --verbosity trace

#define verbose(...) do {                                 \
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) { \
        rt_println(__VA_ARGS__);                             \
    }                                                     \
} while (0)

static void rt_thread_philosopher_think(rt_thread_philosopher_t* p) {
    verbose("philosopher %d is thinking.", p->id);
    // Random think time between .1 and .3 seconds
    fp64_t seconds = (rt_num.random32(&p->ps->seed) % 30 + 1) / 100.0;
    rt_thread.sleep_for(seconds);
}

static void rt_thread_philosopher_eat(rt_thread_philosopher_t* p) {
    verbose("philosopher %d is eating.", p->id);
    // Random eat time between .1 and .2 seconds
    fp64_t seconds = (rt_num.random32(&p->ps->seed) % 20 + 1) / 100.0;
    rt_thread.sleep_for(seconds);
}

// To avoid deadlocks in the Three Philosophers problem, we can implement
// the Tanenbaum's solution, which ensures that one of the philosophers
// (e.g., the last one) tries to pick up the right fork first, while the
// others pick up the left fork first. This breaks the circular wait
// condition and prevents deadlock.

// If the philosopher is the last one (p->id == n - 1) they will try to pick
// up the right fork first and then the left fork. All other philosophers will
// pick up the left fork first and then the right fork, as before. This change
// ensures that at least one philosopher will be able to eat, breaking the
// circular wait condition and preventing deadlock.

static void rt_thread_philosopher_routine(void* arg) {
    rt_thread_philosopher_t* p = (rt_thread_philosopher_t*)arg;
    enum { n = rt_countof(p->ps->philosopher) };
    rt_thread.name("philosopher");
    rt_thread.realtime();
    while (!p->ps->enough) {
        rt_thread_philosopher_think(p);
        if (p->id == n - 1) { // Last philosopher picks up the right fork first
            rt_mutex.lock(p->right_fork);
            verbose("philosopher %d picked up right fork.", p->id);
            rt_mutex.lock(p->left_fork);
            verbose("philosopher %d picked up left fork.", p->id);
        } else { // Other philosophers pick up the left fork first
            rt_mutex.lock(p->left_fork);
            verbose("philosopher %d picked up left fork.", p->id);
            rt_mutex.lock(p->right_fork);
            verbose("philosopher %d picked up right fork.", p->id);
        }
        rt_thread_philosopher_eat(p);
        rt_mutex.unlock(p->right_fork);
        verbose("philosopher %d put down right fork.", p->id);
        rt_mutex.unlock(p->left_fork);
        verbose("philosopher %d put down left fork.", p->id);
        rt_event.set(p->ps->fed_up[p->id]);
    }
}

static void rt_thread_detached_sleep(void* rt_unused(p)) {
    rt_thread.sleep_for(1000.0); // seconds
}

static void rt_thread_detached_loop(void* rt_unused(p)) {
    uint64_t sum = 0;
    for (uint64_t i = 0; i < UINT64_MAX; i++) { sum += i; }
    // make sure that compiler won't get rid of the loop:
    rt_swear(sum == 0x8000000000000001ULL, "sum: %llu 0x%16llX", sum, sum);
}

static void rt_thread_test(void) {
    rt_thread_philosophers_t ps = { .seed = 1 };
    enum { n = rt_countof(ps.philosopher) };
    // Initialize mutexes for forks
    for (int32_t i = 0; i < n; i++) {
        rt_thread_philosopher_t* p = &ps.philosopher[i];
        p->id = i;
        p->ps = &ps;
        rt_mutex.init(&p->fork);
        p->left_fork = &p->fork;
        ps.fed_up[i] = rt_event.create();
    }
    // Create and start philosopher threads
    for (int32_t i = 0; i < n; i++) {
        rt_thread_philosopher_t* p = &ps.philosopher[i];
        rt_thread_philosopher_t* r = &ps.philosopher[(i + 1) % n];
        p->right_fork = r->left_fork;
        p->thread = rt_thread.start(rt_thread_philosopher_routine, p);
    }
    // wait for all philosophers being fed up:
    for (int32_t i = 0; i < n; i++) { rt_event.wait(ps.fed_up[i]); }
    ps.enough = true;
    // join all philosopher threads
    for (int32_t i = 0; i < n; i++) {
        rt_thread_philosopher_t* p = &ps.philosopher[i];
        rt_thread.join(p->thread, -1);
    }
    // Dispose of mutexes and events
    for (int32_t i = 0; i < n; ++i) {
        rt_thread_philosopher_t* p = &ps.philosopher[i];
        rt_mutex.dispose(&p->fork);
        rt_event.dispose(ps.fed_up[i]);
    }
    // detached threads are hacky and not that swell of an idea
    // but sometimes can be useful for 1. quick hacks 2. threads
    // that execute blocking calls that e.g. write logs to the
    // internet service that hangs.
    // test detached threads
    rt_thread_t detached_sleep = rt_thread.start(rt_thread_detached_sleep, null);
    rt_thread.detach(detached_sleep);
    rt_thread_t detached_loop = rt_thread.start(rt_thread_detached_loop, null);
    rt_thread.detach(detached_loop);
    // leave detached threads sleeping and running till ExitProcess(0)
    // that should NOT hang.
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#pragma pop_macro("verbose")

#else
static void rt_thread_test(void) { }
#endif

rt_thread_if rt_thread = {
    .start     = rt_thread_start,
    .join      = rt_thread_join,
    .detach    = rt_thread_detach,
    .name      = rt_thread_name,
    .realtime  = rt_thread_realtime,
    .yield     = rt_thread_yield,
    .sleep_for = rt_thread_sleep_for,
    .id_of     = rt_thread_id_of,
    .id        = rt_thread_id,
    .self      = rt_thread_self,
    .open      = rt_thread_open,
    .close     = rt_thread_close,
    .test      = rt_thread_test
};
// ________________________________ rt_vigil.c ________________________________

#include <stdio.h>
#include <string.h>

static void rt_vigil_breakpoint_and_abort(void) {
    rt_debug.breakpoint(); // only if debugger is present
    rt_debug.raise(rt_debug.exception.noncontinuable);
    rt_core.abort();
}

static int32_t rt_vigil_failed_assertion(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...) {
    va_list va;
    va_start(va, format);
    rt_debug.println_va(file, line, func, format, va);
    va_end(va);
    rt_debug.println(file, line, func, "assertion failed: %s\n", condition);
    // avoid warnings: conditional expression always true and unreachable code
    const bool always_true = rt_core.abort != null;
    if (always_true) { rt_vigil_breakpoint_and_abort(); }
    return 0;
}

static int32_t rt_vigil_fatal_termination_va(const char* file, int32_t line,
        const char* func, const char* condition, errno_t r,
        const char* format, va_list va) {
    const int32_t er = rt_core.err();
    const int32_t en = errno;
    rt_debug.println_va(file, line, func, format, va);
    if (r != er && r != 0) {
        rt_debug.perror(file, line, func, r, "");
    }
    // report last errors:
    if (er != 0) { rt_debug.perror(file, line, func, er, ""); }
    if (en != 0) { rt_debug.perrno(file, line, func, en, ""); }
    if (condition != null && condition[0] != 0) {
        rt_debug.println(file, line, func, "FATAL: %s\n", condition);
    } else {
        rt_debug.println(file, line, func, "FATAL\n");
    }
    const bool always_true = rt_core.abort != null;
    if (always_true) { rt_vigil_breakpoint_and_abort(); }
    return 0;
}

static int32_t rt_vigil_fatal_termination(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...) {
    va_list va;
    va_start(va, format);
    rt_vigil_fatal_termination_va(file, line, func, condition, 0, format, va);
    va_end(va);
    return 0;
}

static int32_t rt_vigil_fatal_if_error(const char* file, int32_t line,
    const char* func, const char* condition, errno_t r,
    const char* format, ...) {
    if (r != 0) {
        va_list va;
        va_start(va, format);
        rt_vigil_fatal_termination_va(file, line, func, condition, r, format, va);
        va_end(va);
    }
    return 0;
}

#ifdef RT_TESTS

static rt_vigil_if  rt_vigil_test_saved;
static int32_t      rt_vigil_test_failed_assertion_count;

#pragma push_macro("rt_vigil")
// intimate knowledge of vigil.*() functions used in macro definitions
#define rt_vigil rt_vigil_test_saved

static int32_t rt_vigil_test_failed_assertion(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...) {
    rt_fatal_if_not(strcmp(file,  __FILE__) == 0, "file: %s", file);
    rt_fatal_if_not(line > __LINE__, "line: %s", line);
    rt_assert(strcmp(func, "rt_vigil_test") == 0, "func: %s", func);
    rt_fatal_if(condition == null || condition[0] == 0);
    rt_fatal_if(format == null || format[0] == 0);
    rt_vigil_test_failed_assertion_count++;
    if (rt_debug.verbosity.level >= rt_debug.verbosity.trace) {
        va_list va;
        va_start(va, format);
        rt_debug.println_va(file, line, func, format, va);
        va_end(va);
        rt_debug.println(file, line, func, "assertion failed: %s (expected)\n",
                     condition);
    }
    return 0;
}

static int32_t rt_vigil_test_fatal_calls_count;

static int32_t rt_vigil_test_fatal_termination(const char* file, int32_t line,
        const char* func, const char* condition, const char* format, ...) {
    const int32_t er = rt_core.err();
    const int32_t en = errno;
    rt_assert(er == 2, "rt_core.err: %d expected 2", er);
    rt_assert(en == 2, "errno: %d expected 2", en);
    rt_fatal_if_not(strcmp(file,  __FILE__) == 0, "file: %s", file);
    rt_fatal_if_not(line > __LINE__, "line: %s", line);
    rt_assert(strcmp(func, "rt_vigil_test") == 0, "func: %s", func);
    rt_assert(strcmp(condition, "") == 0); // not yet used expected to be ""
    rt_assert(format != null && format[0] != 0);
    rt_vigil_test_fatal_calls_count++;
    if (rt_debug.verbosity.level > rt_debug.verbosity.trace) {
        va_list va;
        va_start(va, format);
        rt_debug.println_va(file, line, func, format, va);
        va_end(va);
        if (er != 0) { rt_debug.perror(file, line, func, er, ""); }
        if (en != 0) { rt_debug.perrno(file, line, func, en, ""); }
        if (condition != null && condition[0] != 0) {
            rt_debug.println(file, line, func, "FATAL: %s (testing)\n", condition);
        } else {
            rt_debug.println(file, line, func, "FATAL (testing)\n");
        }
    }
    return 0;
}

#pragma pop_macro("rt_vigil")

static void rt_vigil_test(void) {
    rt_vigil_test_saved = rt_vigil;
    int32_t en = errno;
    int32_t er = rt_core.err();
    errno = 2; // ENOENT
    rt_core.set_err(2); // ERROR_FILE_NOT_FOUND
    rt_vigil.failed_assertion  = rt_vigil_test_failed_assertion;
    rt_vigil.fatal_termination = rt_vigil_test_fatal_termination;
    int32_t count = rt_vigil_test_fatal_calls_count;
    rt_fatal("testing: %s call", "fatal()");
    rt_assert(rt_vigil_test_fatal_calls_count == count + 1);
    count = rt_vigil_test_failed_assertion_count;
    rt_assert(false, "testing: rt_assert(%s)", "false");
    #ifdef DEBUG // verify that rt_assert() is only compiled in DEBUG:
        rt_fatal_if_not(rt_vigil_test_failed_assertion_count == count + 1);
    #else // not RELEASE buid:
        rt_fatal_if_not(rt_vigil_test_failed_assertion_count == count);
    #endif
    count = rt_vigil_test_failed_assertion_count;
    rt_swear(false, "testing: swear(%s)", "false");
    // swear() is triggered in both debug and release configurations:
    rt_fatal_if_not(rt_vigil_test_failed_assertion_count == count + 1);
    errno = en;
    rt_core.set_err(er);
    rt_vigil = rt_vigil_test_saved;
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

#else

static void rt_vigil_test(void) { }

#endif

rt_vigil_if rt_vigil = {
    .failed_assertion  = rt_vigil_failed_assertion,
    .fatal_termination = rt_vigil_fatal_termination,
    .fatal_if_error    = rt_vigil_fatal_if_error,
    .test              = rt_vigil_test
};

// ________________________________ rt_win32.c ________________________________

void rt_win32_close_handle(void* h) {
    #pragma warning(suppress: 6001) // shut up overzealous IntelliSense
    rt_fatal_win32err(CloseHandle((HANDLE)h));
}

// WAIT_ABANDONED only reported for mutexes not events
// WAIT_FAILED means event was invalid handle or was disposed
// by another thread while the calling thread was waiting for it.

/* translate ix to error */
errno_t rt_wait_ix2e(uint32_t r) {
    const int32_t ix = (int32_t)r;
    return (errno_t)(
          (int32_t)WAIT_OBJECT_0 <= ix && ix <= WAIT_OBJECT_0 + 63 ? 0 :
          (ix == WAIT_ABANDONED ? ERROR_REQUEST_ABORTED :
            (ix == WAIT_TIMEOUT ? ERROR_TIMEOUT :
              (ix == WAIT_FAILED) ? rt_core.err() : ERROR_INVALID_HANDLE
            )
          )
    );
}

// ________________________________ rt_work.c _________________________________

static void rt_work_queue_no_duplicates(rt_work_t* w) {
    rt_work_t* e = w->queue->head;
    bool found = false;
    while (e != null && !found) {
        found = e == w;
        if (!found) { e = e->next; }
    }
    rt_swear(!found);
}

static void rt_work_queue_post(rt_work_t* w) {
    rt_assert(w->queue != null && w != null && w->when >= 0.0);
    rt_work_queue_t* q = w->queue;
    rt_atomics.spinlock_acquire(&q->lock);
    rt_work_queue_no_duplicates(w); // under lock
    //  Enqueue in time sorted order least ->time first to save
    //  time searching in fetching from queue which is more frequent.
    rt_work_t* p = null;
    rt_work_t* e = q->head;
    while (e != null && e->when <= w->when) {
        p = e;
        e = e->next;
    }
    w->next = e;
    bool head = p == null;
    if (head) {
        q->head = w;
    } else {
        p->next = w;
    }
    rt_atomics.spinlock_release(&q->lock);
    if (head && q->changed != null) { rt_event.set(q->changed); }
}

static void rt_work_queue_cancel(rt_work_t* w) {
    rt_swear(!w->canceled && w->queue != null && w->queue->head != null);
    rt_work_queue_t* q = w->queue;
    rt_atomics.spinlock_acquire(&q->lock);
    rt_work_t* p = null;
    rt_work_t* e = q->head;
    bool changed = false; // head changed
    while (e != null && !w->canceled) {
        if (e == w) {
            changed = p == null;
            if (changed) {
                q->head = e->next;
        } else {
                p->next = e->next;
            }
            e->next = null;
            e->canceled = true;
        } else {
            p = e;
            e = e->next;
        }
    }
    rt_atomics.spinlock_release(&q->lock);
    rt_swear(w->canceled);
    if (w->done != null) { rt_event.set(w->done); }
    if (changed && q->changed != null) { rt_event.set(q->changed); }
}

static void rt_work_queue_flush(rt_work_queue_t* q) {
    while (q->head != null) { rt_work_queue.cancel(q->head); }
}

static bool rt_work_queue_get(rt_work_queue_t* q, rt_work_t* *r) {
    rt_work_t* w = null;
    rt_atomics.spinlock_acquire(&q->lock);
    bool changed = q->head != null && q->head->when <= rt_clock.seconds();
    if (changed) {
        w = q->head;
        q->head = w->next;
        w->next = null;
    }
    rt_atomics.spinlock_release(&q->lock);
    *r = w;
    if (changed && q->changed != null) { rt_event.set(q->changed); }
    return w != null;
}

static void rt_work_queue_call(rt_work_t* w) {
    if (w->work != null) { w->work(w); }
    if (w->done != null) { rt_event.set(w->done); }
}

static void rt_work_queue_dispatch(rt_work_queue_t* q) {
    rt_work_t* w = null;
    while (rt_work_queue.get(q, &w)) { rt_work_queue.call(w); }
}

rt_work_queue_if rt_work_queue = {
    .post     = rt_work_queue_post,
    .get      = rt_work_queue_get,
    .call     = rt_work_queue_call,
    .dispatch = rt_work_queue_dispatch,
    .cancel   = rt_work_queue_cancel,
    .flush    = rt_work_queue_flush
};

static void rt_worker_thread(void* p) {
    rt_thread.name("worker");
    rt_worker_t* worker = (rt_worker_t*)p;
    rt_work_queue_t* q = &worker->queue;
    while (!worker->quit) {
        rt_work_queue.dispatch(q);
        fp64_t timeout = -1.0; // forever
        rt_atomics.spinlock_acquire(&q->lock);
        if (q->head != null) {
            timeout = rt_max(0, q->head->when - rt_clock.seconds());
        }
        rt_atomics.spinlock_release(&q->lock);
        // if another item is inserted into head after unlocking
        // the `wake` event guaranteed to be signalled
        if (!worker->quit && timeout != 0) {
            rt_event.wait_or_timeout(worker->wake, timeout);
        }
    }
    rt_work_queue.dispatch(q);
}

static void rt_worker_start(rt_worker_t* worker) {
    rt_assert(worker->wake == null && !worker->quit);
    worker->wake  = rt_event.create();
    worker->queue = (rt_work_queue_t){
        .head = null, .lock = 0, .changed = worker->wake
    };
    worker->thread = rt_thread.start(rt_worker_thread, worker);
}

static errno_t rt_worker_join(rt_worker_t* worker, fp64_t to) {
    worker->quit = true;
    rt_event.set(worker->wake);
    errno_t r = rt_thread.join(worker->thread, to);
    if (r == 0) {
        rt_event.dispose(worker->wake);
        worker->wake = null;
        worker->thread = null;
        worker->quit = false;
        rt_swear(worker->queue.head == null);
    }
    return r;
}

static void rt_worker_post(rt_worker_t* worker, rt_work_t* w) {
    rt_assert(!worker->quit && worker->wake != null && worker->thread != null);
    w->queue = &worker->queue;
    rt_work_queue.post(w);
}

static void rt_worker_test(void);

rt_worker_if rt_worker = {
    .start = rt_worker_start,
    .post  = rt_worker_post,
    .join  = rt_worker_join,
    .test  = rt_worker_test
};

#ifdef RT_TESTS

// tests:

// keep in mind that rt_println() may be blocking and is a subject
// of "astronomical" wait state times in order of dozens of ms.

static int32_t rt_test_called;

static void rt_never_called(rt_work_t* rt_unused(w)) {
    rt_test_called++;
}

static void rt_work_queue_test_1(void) {
    rt_test_called = 0;
    // testing insertion time ordering of two events into queue
    const fp64_t now = rt_clock.seconds();
    rt_work_queue_t q = {0};
    rt_work_t c1 = {
        .queue = &q,
        .work = rt_never_called,
        .when = now + 1.0
    };
    rt_work_t c2 = {
        .queue = &q,
        .work = rt_never_called,
        .when = now + 0.5
    };
    rt_work_queue.post(&c1);
    rt_swear(q.head == &c1 && q.head->next == null);
    rt_work_queue.post(&c2);
    rt_swear(q.head == &c2 && q.head->next == &c1);
    rt_work_queue.flush(&q);
    // test that canceled events are not dispatched
    rt_swear(rt_test_called == 0 && c1.canceled && c2.canceled && q.head == null);
    c1.canceled = false;
    c2.canceled = false;
    // test the rt_work_queue.cancel() function
    rt_work_queue.post(&c1);
    rt_work_queue.post(&c2);
    rt_swear(q.head == &c2 && q.head->next == &c1);
    rt_work_queue.cancel(&c2);
    rt_swear(c2.canceled && q.head == &c1 && q.head->next == null);
    c2.canceled = false;
    rt_work_queue.post(&c2);
    rt_work_queue.cancel(&c1);
    rt_swear(c1.canceled && q.head == &c2 && q.head->next == null);
    rt_work_queue.flush(&q);
    rt_swear(rt_test_called == 0 && c1.canceled && c2.canceled && q.head == null);
}

// simple way of passing a single pointer to call_later

static fp64_t rt_test_work_start; // makes timing debug traces easier to read

static void rt_every_millisecond(rt_work_t* w) {
    int32_t* i = (int32_t*)w->data;
    fp64_t now = rt_clock.seconds();
    if (rt_debug.verbosity.level > rt_debug.verbosity.info) {
        const fp64_t since_start = now - rt_test_work_start;
        const fp64_t dt = w->when - rt_test_work_start;
        rt_println("%d now: %.6f time: %.6f", *i, since_start, dt);
    }
    (*i)++;
    // read rt_clock.seconds() again because rt_println() above could block
    w->when = rt_clock.seconds() + 0.001;
    rt_work_queue.post(w);
}

static void rt_work_queue_test_2(void) {
    rt_thread.realtime();
    rt_test_work_start = rt_clock.seconds();
    rt_work_queue_t q = {0};
    // if a single pointer will suffice
    int32_t i = 0;
    rt_work_t c = {
        .queue = &q,
        .work = rt_every_millisecond,
        .when = rt_test_work_start + 0.001,
        .data = &i
    };
    rt_work_queue.post(&c);
    while (q.head != null && i < 8) {
        rt_thread.sleep_for(0.0001); // 100 microseconds
        rt_work_queue.dispatch(&q);
    }
    rt_work_queue.flush(&q);
    rt_swear(q.head == null);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) {
        rt_println("called: %d times", i);
    }
}

// extending rt_work_t with extra data:

typedef struct rt_work_ex_s {
    // nameless union opens up base fields into rt_work_ex_t
    // it is not necessary at all
    union {
        rt_work_t base;
        struct rt_work_s;
    };
    struct { int32_t a; int32_t b; } s;
    int32_t i;
} rt_work_ex_t;

static void rt_every_other_millisecond(rt_work_t* w) {
    rt_work_ex_t* ex = (rt_work_ex_t*)w;
    fp64_t now = rt_clock.seconds();
    if (rt_debug.verbosity.level > rt_debug.verbosity.info) {
        const fp64_t since_start = now - rt_test_work_start;
        const fp64_t dt  = w->when - rt_test_work_start;
        rt_println(".i: %d .extra: {.a: %d .b: %d} now: %.6f time: %.6f",
                ex->i, ex->s.a, ex->s.b, since_start, dt);
    }
    ex->i++;
    const int32_t swap = ex->s.a; ex->s.a = ex->s.b; ex->s.b = swap;
    // read rt_clock.seconds() again because rt_println() above could block
    w->when = rt_clock.seconds() + 0.002;
    rt_work_queue.post(w);
}

static void rt_work_queue_test_3(void) {
    rt_thread.realtime();
    rt_static_assertion(offsetof(rt_work_ex_t, base) == 0);
    const fp64_t now = rt_clock.seconds();
    rt_work_queue_t q = {0};
    rt_work_ex_t ex = {
        .queue = &q,
        .work = rt_every_other_millisecond,
        .when = now + 0.002,
        .s = { .a = 1, .b = 2 },
        .i = 0
    };
    rt_work_queue.post(&ex.base);
    while (q.head != null && ex.i < 8) {
        rt_thread.sleep_for(0.0001); // 100 microseconds
        rt_work_queue.dispatch(&q);
    }
    rt_work_queue.flush(&q);
    rt_swear(q.head == null);
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) {
        rt_println("called: %d times", ex.i);
    }
}

static void rt_work_queue_test(void) {
    rt_work_queue_test_1();
    rt_work_queue_test_2();
    rt_work_queue_test_3();
    if (rt_debug.verbosity.level > rt_debug.verbosity.quiet) { rt_println("done"); }
}

static int32_t rt_test_do_work_called;

static void rt_test_do_work(rt_work_t* rt_unused(w)) {
    rt_test_do_work_called++;
}

static void rt_worker_test(void) {
//  uncomment one of the following lines to see the output
//  rt_debug.verbosity.level = rt_debug.verbosity.info;
//  rt_debug.verbosity.level = rt_debug.verbosity.verbose;
    rt_work_queue_test(); // first test rt_work_queue
    rt_worker_t worker = { 0 };
    rt_worker.start(&worker);
    rt_work_t asap = {
        .when = 0, // A.S.A.P.
        .done = rt_event.create(),
        .work = rt_test_do_work
    };
    rt_work_t later = {
        .when = rt_clock.seconds() + 0.010, // 10ms
        .done = rt_event.create(),
        .work = rt_test_do_work
    };
    rt_worker.post(&worker, &asap);
    rt_worker.post(&worker, &later);
    // because `asap` and `later` are local variables
    // code needs to wait for them to be processed inside
    // this function before they goes out of scope
    rt_event.wait(asap.done); // await(asap)
    rt_event.dispose(asap.done); // responsibility of the caller
    // wait for later:
    rt_event.wait(later.done); // await(later)
    rt_event.dispose(later.done); // responsibility of the caller
    // quit the worker thread:
    rt_fatal_if_error(rt_worker.join(&worker, -1.0));
    // does worker respect .when dispatch time?
    rt_swear(rt_clock.seconds() >= later.when);
}

#else

static void rt_work_queue_test(void) {}
static void rt_worker_test(void) {}

#endif

#endif // rt_implementation

