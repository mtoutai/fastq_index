// fastq_index_size.c
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

#define BLOCK_SIZE 4096


void usage(FILE *fp)
{
	fprintf(fp, "Prints the size of a given index file.\n");
	fprintf(fp, "Usage: fastq_index_size <fastq_index_file>\n");
}


int main (int argc, char **argv)
{
	if (argc == 1) {
		usage(stdout);
		return(0);
	}

	struct stat file_stat;
	int ret=stat(argv[1], &file_stat);
	if (ret!=0) {
		fprintf(stderr, "Error: Unable to read the specified index file \"%s\".\n", argv[1]);
		return(1);
	}

	size_t index_table_size = BLOCK_SIZE/sizeof(long int);
	DEBUG("index_table_size=" << index_table_size);
	long int index_table[index_table_size];
	long int index_raw_size = file_stat.st_size;
	DEBUG("index_raw_size=" << index_raw_size);

	FILE *fp=NULL;
	fp = fopen(argv[1], "rb");
	fseek(fp, -BLOCK_SIZE, SEEK_END);
	while (index_table[0] == 0) {
		fseek(fp, -BLOCK_SIZE, SEEK_CUR);
		fread(index_table, sizeof(char), BLOCK_SIZE, fp);
	}
	DEBUG("ftell=" << ftell(fp));

	long int index_size;
	for (index_size=index_table_size-1; index_size >= 0 && index_table[index_size] == 0L; --index_size) {}
	index_size++;
	index_size += (ftell(fp) - BLOCK_SIZE)/sizeof(long int);

	printf("%d\n", index_size);

	fclose(fp);
	return(0);
}

