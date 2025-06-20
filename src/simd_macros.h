#ifndef SIMD_MACROS_H_
#define SIMD_MACROS_H_


#if !defined(USE_SCALAR) && !defined(USE_SSE) && !defined(USE_AVX2) && !defined(USE_AVX512)
  #if defined(__AVX512F__)
    #define USE_AVX512
  #elif defined(__AVX2__)
    #define USE_AVX2
  #elif defined(__SSE4_2__)
    #define USE_SSE
  #else
    #define USE_SCALAR
  #endif
#endif

#ifdef USE_SCALAR
typedef int SIMD_INT;

#define SET1_INT(x) (x)
#define LOAD_INT(p) (*(p))
#define STORE_INT(p, v) (*(p) = (v))
#define LOADU_INT(p) (*(p))
#define STOREU_INT(p, v) (*(p) = (v))
#define ADD_INT(a, b) ((a) + (b))
#define MAX_INT(a, b) std::max((a), (b))
#define CMP_GT_INT32(a, b) ((a) > (b) ? 1 : 0)
#define COUNT_EQUAL_BYTES(a, b) \
    ((a & 0xFF000000) == (b & 0xFF000000)) + \
    ((a & 0x00FF0000) == (b & 0x00FF0000)) + \
    ((a & 0x0000FF00) == (b & 0x0000FF00)) + \
    ((a & 0x000000FF) == (b & 0x000000FF))
#define INT_PER_BLOCK 1
#define BYTES_PER_BLOCK 4

#elif defined(USE_SSE)
#include <immintrin.h>

typedef __m128i SIMD_INT;

#define SET1_INT _mm_set1_epi32
#define LOAD_INT _mm_load_si128
#define STORE_INT _mm_store_si128
#define LOADU_INT _mm_loadu_si128
#define STOREU_INT _mm_storeu_si128
#define ADD_INT _mm_add_epi32
#define MAX_INT _mm_max_epi32
#define CMP_GT_INT32(a, b) _mm_movemask_epi8(_mm_cmpgt_epi32(a, b))
#define COUNT_EQUAL_BYTES(a, b) _mm_popcnt_u32(_mm_movemask_epi8(_mm_cmpeq_epi8(a, b)))
#define INT_PER_BLOCK 4
#define BYTES_PER_BLOCK 16

#elif defined(USE_AVX2)
#include <immintrin.h>

typedef __m256i SIMD_INT;

#define SET1_INT _mm256_set1_epi32
#define LOAD_INT _mm256_load_si256
#define STORE_INT _mm256_store_si256
#define LOADU_INT _mm256_loadu_si256
#define STOREU_INT _mm256_storeu_si256
#define ADD_INT _mm256_add_epi32
#define MAX_INT _mm256_max_epi32
#define CMP_GT_INT32 _mm256_cmpgt_epi32_mask
#define COUNT_EQUAL_BYTES(a, b) _mm_popcnt_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)))
#define INT_PER_BLOCK 8
#define BYTES_PER_BLOCK 32

#elif defined(USE_AVX512)
#include <immintrin.h>

typedef __m512i SIMD_INT;

#define SET1_INT _mm512_set1_epi32
#define LOAD_INT _mm512_load_si512
#define STORE_INT _mm512_store_si512
#define LOADU_INT _mm512_loadu_si512
#define STOREU_INT _mm512_storeu_si512
#define ADD_INT _mm512_add_epi32
#define MAX_INT _mm512_max_epi32
#define CMP_GT_INT32 _mm512_cmpgt_epi32_mask
#define COUNT_EQUAL_BYTES(a, b) _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(a, b))
#define INT_PER_BLOCK 16
#define BYTES_PER_BLOCK 64

#else
#error "Either USE_SCALAR, USE_SSE, USE_AVX2, or USE_AVX512 must be defined"
#endif


#endif /* SIMD_MACROS_H_ */
