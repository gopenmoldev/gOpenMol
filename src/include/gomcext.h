#ifndef INC_GOPENMOL_GOM_CC
#define INC_GOPENMOL_GOM_CC

/* Extensions to the C language. */

/* CHECK_FORMAT_PRINTF */
#ifdef __GNUC__
#  define CHECK_FORMAT_PRINTF          __attribute__ ((format (printf, 1, 2)))
#  define CHECK_FORMAT_PRINTF_1_EXTRA  __attribute__ ((format (printf, 2, 3)))
#  define CHECK_FORMAT_PRINTF_2_EXTRA  __attribute__ ((format (printf, 3, 4)))
#  define CHECK_FORMAT_VPRINTF         __attribute__ ((format (printf, 1, 0)))
#  define CHECK_FORMAT_VPRINTF_1_EXTRA __attribute__ ((format (printf, 2, 0)))
#  define CHECK_FORMAT_VPRINTF_2_EXTRA __attribute__ ((format (printf, 3, 0)))
#else
#  define CHECK_FORMAT_PRINTF
#  define CHECK_FORMAT_PRINTF_1_EXTRA
#  define CHECK_FORMAT_PRINTF_2_EXTRA
#  define CHECK_FORMAT_VPRINTF
#  define CHECK_FORMAT_VPRINTF_1_EXTRA
#  define CHECK_FORMAT_VPRINTF_2_EXTRA
#endif

/* CHECK_FORMAT_SCANF */
#ifdef __GNUC__
#  define CHECK_FORMAT_SCANF          __attribute__ ((format (scanf, 1, 2)))
#  define CHECK_FORMAT_SCANF_1_EXTRA  __attribute__ ((format (scanf, 2, 3)))
#  define CHECK_FORMAT_SCANF_2_EXTRA  __attribute__ ((format (scanf, 3, 4)))
#  define CHECK_FORMAT_VSCANF         __attribute__ ((format (scanf, 1, 0)))
#  define CHECK_FORMAT_VSCANF_1_EXTRA __attribute__ ((format (scanf, 2, 0)))
#  define CHECK_FORMAT_VSCANF_2_EXTRA __attribute__ ((format (scanf, 3, 0)))
#else
#  define CHECK_FORMAT_SCANF
#  define CHECK_FORMAT_SCANF_1_EXTRA
#  define CHECK_FORMAT_SCANF_2_EXTRA
#  define CHECK_FORMAT_VSCANF
#  define CHECK_FORMAT_VSCANF_1_EXTRA
#  define CHECK_FORMAT_VSCANF_2_EXTRA
#endif

#endif /* INC_GOPENMOL_GOM_CC */
