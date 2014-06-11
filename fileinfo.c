#define _fileinfo

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//-----------------------------------------------------------------------------------
// Returns 1 if the file exists and has read permission, 0 otherwise
//-----------------------------------------------------------------------------------
int checkReadAccess(char *file)
{
    if( access(file,R_OK) == -1 ) return(0);

    return(1);
}
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
// Returns 1 on write permission, 0 otherwise
//-----------------------------------------------------------------------------------
int checkWriteAccess(char *file)
{
    if( access(file,W_OK) == -1 ) return(0);

    return(1);
}
//-----------------------------------------------------------------------------------

int checkFileExistence(char *file)
{
    if( access(file,F_OK) == -1) return(0);

	return(1);
}

int checkFileReadOK(char *file)
{
    if( access(file,R_OK) == -1 ) return(0);

	return(1);
}

int checkFileWriteOK(char *file)
{

    if( access(file,F_OK)==0 && access(file,W_OK) == -1 ) return(0);

	return(1);
}

int check_F_R_permission(char *file)
{

	if( access(file,F_OK) == -1 )
		return(0);

	if( access(file,R_OK) == -1 )
		return(0);

	return(1);
}

int getFileSize(const char *file)
{
	struct stat fileinfo;   // structure into  which  information is placed concerning the file

	if( stat(file, &fileinfo) == -1 )
	{
		printf("Error: stat() failure. There is something wrong with %s!\n",file);
		return(0);
	}

	return(fileinfo.st_size);
}

int isregular(const char *file)
{
	struct stat fileinfo;   // structure into  which  information is placed concerning the file

	if( stat(file, &fileinfo) == -1 )
	{
		printf("Error: stat() failure. There is something wrong with %s!\n",file);
		return(0);
	}

	if( S_ISREG(fileinfo.st_mode) ) 
		return(1); 
	else
		return(0);
}
