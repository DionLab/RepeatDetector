#ifndef _CONFIG_H
#define _CONFIG_H
/*
 ********************************************************************************
 * SPECIFIC ARCHITECTURE
 ********************************************************************************
 */
#cmakedefine USE_WINAPI
#ifdef USE_WINAPI
#define __USE_WINAPI__
#endif

#cmakedefine USE_AFFINITY
#if defined(USE_AFFINITY)
# if defined(__GNUC__) && !defined(USE_WINAPI)
#  define _GNU_SOURCE
# endif
#endif

/*
 ********************************************************************************
 * Math Library
 ********************************************************************************
 */
#ifdef __INTEL_COMPILER
# include <mathimf.h>
#elif __IBMC__
# include <math.h>
# include <mass.h>
#else
# include <math.h>
#endif

/*
 ********************************************************************************
 * ALLOCA header file
 ********************************************************************************
 */
#cmakedefine HAVE_ALLOCA_H
#if !defined(HAVE_ALLOCA_H) && defined(__GNUC__)
/* This is not defined in MinGW */
#define alloca __builtin_alloca
#endif

/*
 ********************************************************************************
 * MM_MALLOC header file
 ********************************************************************************
 */
#cmakedefine HAVE_MM_MALLOC_H 1
#ifdef HAVE_MM_MALLOC_H
# include <mm_malloc.h>
#else
#include <stdlib.h>
/* We can't depend on <stdlib.h> since the prototype of posix_memalign
   may not be visible.  */
#ifndef __cplusplus
extern int posix_memalign (void **, size_t, size_t);
#else
extern "C" int posix_memalign (void **, size_t, size_t) throw ();
#endif

static inline void * __ALWAYS_INLINE
_mm_malloc (size_t size, size_t alignment)
{
  void *ptr;
  if (alignment == 1)
    return malloc (size);
  if (alignment == 2 || (sizeof (void *) == 8 && alignment == 4))
    alignment = sizeof (void *);
  if (posix_memalign (&ptr, alignment, size) == 0)
    return ptr;
  else
    return NULL;
}

static inline void __ALWAYS_INLINE
_mm_free (void * ptr)
{
  free (ptr);
}
#endif

/*
 ********************************************************************************
 * Executable build using inline functions
 ********************************************************************************
 */
#ifndef BUILD_LIBRARY
#define __USE_INLINE_FUNCTIONS__
#endif

#define SUCCESS 0


#endif /* _CONFIG_H */
