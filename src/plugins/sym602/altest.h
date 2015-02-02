#include "align.h"
#include <gopenmolext.h>

namespace Plugin {
namespace Symmetry {

class atomrec
{
public:
    int NAtom;
    int Str1;
    int *AtomList1;
    int Str2;
    int *AtomList2;
    int FixedAtom1;
    int FixedAtom2;
};

void AlignStructs(atomrec* ar);

} // namespace Symmetry
} // namespace Plugin
