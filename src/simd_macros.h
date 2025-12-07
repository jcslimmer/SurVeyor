#ifndef SIMD_MACROS_H_
#define SIMD_MACROS_H_

#include <cstdint>
#include <cstring>
#include <algorithm>

#if !defined(USE_SCALAR) && !defined(USE_SSE) && !defined(USE_AVX2) //&& !defined(USE_AVX512)
  // #if defined(__AVX512F__) && defined(__AVX512BW__)  CURRENTLY I HAVE NO WAY TO TEST THIS AS I HAVE NO ACCESS TO AVX512 HARDWARE
    // AVX512BW is required for byte/word (16-bit) operations
    // #define USE_AVX512
  #if defined(__AVX2__)
    #define USE_AVX2
  #elif defined(__SSE4_2__)
    #define USE_SSE
  #else
    #define USE_SCALAR
  #endif
#endif

typedef int16_t SW_SCORE_INT_16;

#ifdef USE_SCALAR

typedef int16_t SIMD_INT_16;

static inline SIMD_INT_16 load_int_unaligned_16(const void* p) {
    SIMD_INT_16 v;
    std::memcpy(&v, p, sizeof v); 
    return v;
}
static inline void store_int_unaligned_16(void* p, SIMD_INT_16 v) {
    std::memcpy(p, &v, sizeof v);
}

#define SET1_INT_16(x) ((int16_t)(x))
#define LOAD_INT_16(p) (*(p))
#define STORE_INT_16(p, v) (*(p) = (v))
#define LOADU_INT_16(p) load_int_unaligned_16(p)
#define STOREU_INT_16(p, v) store_int_unaligned_16(p, v)

// Standard Add (Explicit cast to prevent int promotion)
#define ADD_INT_16(a, b) ((SIMD_INT_16)((a) + (b)))

// Saturated Add (Clamps at 32767) - Explicit casts for std::min/max deduction
#define ADDS_INT_16(a, b) ((SIMD_INT_16)std::max((int)-32768, std::min((int)32767, (int)((a) + (b)))))

#define MAX_INT_16(a, b) std::max((SIMD_INT_16)(a), (SIMD_INT_16)(b))
#define CMP_GT_INT_16(a, b) ((a) > (b) ? 1 : 0)
#define CMP_GT_INT_16_EXCEPT_LAST(a, b) 0

// Count bytes macro
#define COUNT_EQUAL_BYTES_16(a, b) (\
    ((a & 0xFF00) == (b & 0xFF00)) + \
    ((a & 0x00FF) == (b & 0x00FF)) )

#define INT_PER_BLOCK_16 1
#define BYTES_PER_BLOCK_16 (sizeof(SIMD_INT_16))

#elif defined(USE_SSE)
#include <immintrin.h>

typedef __m128i SIMD_INT_16;

#define SET1_INT_16 _mm_set1_epi16
#define LOAD_INT_16 _mm_load_si128
#define STORE_INT_16 _mm_store_si128
#define LOADU_INT_16 _mm_loadu_si128
#define STOREU_INT_16 _mm_storeu_si128

#define ADD_INT_16 _mm_add_epi16
#define ADDS_INT_16 _mm_adds_epi16 // Saturated Add

#define MAX_INT_16 _mm_max_epi16
// Returns non-zero mask if any element is greater
#define CMP_GT_INT_16(a, b) _mm_movemask_epi8(_mm_cmpgt_epi16(a, b))
#define CMP_GT_INT_16_EXCEPT_LAST(a, b) (CMP_GT_INT_16(a, b) & 0x3FFF)

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u32(_mm_movemask_epi8(_mm_cmpeq_epi8(a, b)))

#define INT_PER_BLOCK_16 8
#define BYTES_PER_BLOCK_16 16

#elif defined(USE_AVX2)
#include <immintrin.h>

typedef __m256i SIMD_INT_16;

#define SET1_INT_16 _mm256_set1_epi16
#define LOAD_INT_16 _mm256_load_si256
#define STORE_INT_16 _mm256_store_si256
#define LOADU_INT_16 _mm256_loadu_si256
#define STOREU_INT_16 _mm256_storeu_si256

#define ADD_INT_16 _mm256_add_epi16
#define ADDS_INT_16 _mm256_adds_epi16 // Saturated Add

#define MAX_INT_16 _mm256_max_epi16

// Movemask epi8 returns 32 bits. If any 16-bit compare was true (0xFFFF), 
// the corresponding 2 bits in the mask will be set. result != 0 is sufficient.
#define CMP_GT_INT_16(a, b) _mm256_movemask_epi8(_mm256_cmpgt_epi16((a),(b)))
#define CMP_GT_INT_16_EXCEPT_LAST(a, b) (CMP_GT_INT_16(a, b) & 0x3FFFFFFF)

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u32(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)))

#define INT_PER_BLOCK_16 16
#define BYTES_PER_BLOCK_16 32

#elif defined(USE_AVX512)
#include <immintrin.h>

typedef __m512i SIMD_INT_16;

#define SET1_INT_16 _mm512_set1_epi16
#define LOAD_INT_16 _mm512_load_si512
#define STORE_INT_16 _mm512_store_si512
#define LOADU_INT_16 _mm512_loadu_si512
#define STOREU_INT_16 _mm512_storeu_si512

#define ADD_INT_16 _mm512_add_epi16
#define ADDS_INT_16 _mm512_adds_epi16 // Saturated Add

#define MAX_INT_16 _mm512_max_epi16
// AVX512 mask comparison returns a bitmask directly (__mmask32)
#define CMP_GT_INT_16 _mm512_cmpgt_epi16_mask
#define CMP_GT_INT_16_EXCEPT_LAST(a, b) (CMP_GT_INT_16(a, b) & 0x7FFFFFFF)

#define COUNT_EQUAL_BYTES_16(a, b) _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(a, b))

#define INT_PER_BLOCK_16 32
#define BYTES_PER_BLOCK_16 64

#else
#error "Either USE_SCALAR, USE_SSE, USE_AVX2, or USE_AVX512 must be defined"
#endif

#endif /* SIMD_MACROS_H_ */
