#ifndef LOCAL_ALLOC_H
#define LOCAL_ALLOC_H

// function prototpyes:
char	**alloc2char(size_t n1, size_t n2);
void	free2char(char **p);

cmpFileIndex    **alloc2cmpFileIndex(size_t n1, size_t n2);
void	free2cmpFileIndex(cmpFileIndex **p);

#endif
