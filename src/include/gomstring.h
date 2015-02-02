extern int   gomp_IsStringAFloat(const char *);
extern char *gomp_CopyString(char *, const char *, size_t);
extern int   gomp_String2Lower(char *);
extern int   gomp_Indexo(const char * ,const char *);
extern const char *gomp_StrTok(char *, const char *);
#define STRTOK gomp_StrTok
extern int   gomp_StringTrim(char *);
extern int   gomp_SplitString(char *, const char * , const char *[BUFF_LEN]);
extern int   gomp_Sign_char(char *, int);
extern int   gomp_CompareStrings(const char *, const char *, int);
