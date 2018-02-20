#include <stdio.h>
#include <stdlib.h>

void memory_allocation_error(const char *s)
{
	printf("\nMemory allocation error for variable %s, aborting ...\n\n",s);
	exit(EXIT_FAILURE);
}

void file_open_error(const char *s)
{
	printf("\nError: Cannot open \"%s\", aborting ...\n\n",s);
	exit(EXIT_FAILURE);
}

void errorMessage(const char *s)
{
	printf("\nError: %s\n\n",s);
	exit(EXIT_FAILURE);
}
