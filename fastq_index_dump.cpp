// fastq_index_dump.c
//
// THIS SOFTWARE IS DISTRIBUTED UNDER
// CC0 1.0 Universal
// Public Domain Dedication
// http://creativecommons.org/publicdomain/zero/1.0/
//
// Developed by:
// TOTAI MITSUYAMA, PhD.
// Biotechnology Research Institute for Drug Discovery
// National Institute of Advanced Industrial Science and Technology (AIST)
//
// Version; 0.1
// August 3. 2015.
//
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#ifdef DEBUG_BUILD
#define DEBUG(x) do { std::cerr << "Debug: " << __FILE__ << ": " << __LINE__ << ": " << x << std::endl; } while(0)
#else
#define DEBUG(x)
#endif /* DEBUG_BUILD */


int main (int argc, char **argv)
{
	if (argc == 1) {
		printf("Dump a given index file.\n");
		printf("Usage: fastq_index_dump <fastq_index_file>\n");
		return(0);
	}

	int arg_opt;

	struct stat file_stat;
	int ret=stat(argv[1], &file_stat);
	if (ret!=0) {
		fprintf(stderr, "Error: main: Unable to read an index file \"%s\".\n", argv[1]);
		return(1);
	}

	long int *index_table = NULL;
	index_table = (long int*)calloc(file_stat.st_size, sizeof(char));
	if (index_table==NULL) {
		fprintf(stderr, "Error: main: Memory allocation failed.\n");
		return(1);
	}

	long int index_raw_size = file_stat.st_size / (sizeof(long int));
	DEBUG("index_raw_size= " << index_raw_size);

	FILE *fp=NULL;
	fp = fopen(argv[1], "rb");
	fread(index_table, sizeof(char), file_stat.st_size, fp);
	fclose(fp);

	long int index_size = 0;
	for ( ; index_size < index_raw_size && index_table[index_size] != 0L; ++index_size) {}

	DEBUG("index_size = " << index_size);

	for (int i=0; i<index_size; ++i) {
		printf("%08d %ld\n", i, index_table[i]);
	}

	return(0);
}

