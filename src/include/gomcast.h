static inline void *ConstCast(const volatile void *);

static inline void *ConstCast(const volatile void *v)
{
    union ConstCast_t
    {
        const volatile void *cp;
        void                *p;
    } cc;
    cc.cp = v;
    return cc.p;
}

#define CONST_CAST(type,v) ((void)(1||(v)==(const type)0),(type)ConstCast(v))
