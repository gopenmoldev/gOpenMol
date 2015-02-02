#define DELEGATE_TO_VFMT(name,paramlist,params) \
int gomp_##name##f paramlist \
{ \
    va_list args; \
    int len; \
    va_start(args, format); \
    len = gomp_v##name##f params; \
    va_end(args); \
    return len; \
}

#undef NEED_GOMP_VXPRINTF
#if !defined(HAVE_C99_XPRINTF)
#  define NEED_GOMP_VXPRINTF
#elif !defined(HAVE_VASPRINTF) && !defined(HAVE_VA_COPY)
#  define NEED_GOMP_VXPRINTF
#endif

#if defined(NEED_GOMP_VXPRINTF) || !defined(HAVE_C99_SCANF)

/**
 * The maximum length of the printf field.
 * Field format is %[<flags>][<width>][.<precision>][<size>]<type>
 * */
#define FORMAT_FIELD_LEN   64

/**
 * The string prodused by
 *     sprintf(buff, "%.*f", MAX_PRINTF_WIDTH, MAX_PRINTF_F_VALUE)
 * must fit into buff[PRINTF_BUFF_LEN].
 * sprintf(buff, "%.*f", 128, 1e30) will produce something like
 * 1000000...00000.000000...000000
 *  |< 30 zeros >| |< 128 zeros >|
 * Note that the result may depend on the current locale (i.e. there
 * can be thousand separators and a desimal point may vary).
 */
#define MAX_PRINTF_WIDTH   128
#define MAX_PRINTF_F_VALUE 1e30
#define PRINTF_BUFF_LEN    256

#if defined(WIN32)
#  define LLONG              __int64
#  define ULLONG    unsigned __int64
#  define PRINTF_FMT_LL "I64"
#elif SIZEOF_LONG_LONG_INT
#  define LLONG               long long int
#  define ULLONG     unsigned long long int
#  define PRINTF_FMT_LL "ll"
#else
#  define LLONG              long int
#  define ULLONG    unsigned long int
#  define PRINTF_FMT_LL "l"
#endif

extern int gomp_parse_format_size(
    char *, int *, const char **, va_list *, size_t);

enum printf_size {
    pfsz_hh,
    pfsz_h,
    pfsz,
    pfsz_l,
    pfsz_ll,
    pfsz_j,
    pfsz_t,
    pfsz_z,
    pfsz_L
};

#define PFSZ_DELEGATE_AS(type) \
    ( ( sizeof(type) >= sizeof(LLONG   ) ) ? pfsz_ll : \
      ( sizeof(type) >= sizeof(long int) ) ? pfsz_l  : \
      pfsz )


typedef struct vxprintf_data_s vxprintf_data;
typedef int (*append_func)(const char *, size_t, vxprintf_data *);
extern int gomp_vxprintf(vxprintf_data *);

struct vxprintf_data_s {
    append_func  append_text;
    size_t       index;
    const char  *format;
    va_list      args;
    FILE        *file;
};

typedef struct vsnprintf_data_s {
    vxprintf_data data;
    char *dst;
    size_t count;
} vsnprintf_data;

#define INT_TYPE \
         'd': \
    case 'i'

#define UINT_TYPE \
         'o': \
    case 'u': \
    case 'x': \
    case 'X'

#define FLOAT_TYPE \
         'f': \
    case 'F': \
    case 'e': \
    case 'E': \
    case 'g': \
    case 'G': \
    case 'a': \
    case 'A'

#endif /* NEED_GOMP_VXPRINTF */
