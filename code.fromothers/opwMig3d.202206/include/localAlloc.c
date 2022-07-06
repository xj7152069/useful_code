
#include "fbCommon.h"
#include "cmpFileIndex.h"

char	**alloc2char(size_t n1, size_t n2)
{
        return	(char**)alloc2(n1, n2, sizeof(char));
}

void	free2char(char **p)
{
        free2((void**)p);
}

cmpFileIndex    **alloc2cmpFileIndex(size_t n1, size_t n2)
{
        return	(cmpFileIndex **)alloc2(n1, n2, sizeof(cmpFileIndex));
}

void	free2cmpFileIndex(cmpFileIndex **p)
{
        free2((void**)p);
}
