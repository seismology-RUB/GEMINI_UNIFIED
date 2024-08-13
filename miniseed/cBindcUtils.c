/*--------------------------------------------------------------------
 *  Collection of C-code used for interfacing to Fortran 95
 *-------------------------------------------------------------------*/
#include <string.h>
#include <dirent.h>
/*--------------------------------------------------------------
 *  unpack the contents of a struct dirent
*/
int
unpackDirentFileSystem(struct dirent* entry,char* fname,unsigned char* ftype)
{
	int ls = strlen(entry->d_name);
	strncpy(fname,entry->d_name,ls);       /* fname provides 133 characters including \0 */
	*ftype = entry->d_type;
	return ls;
}
