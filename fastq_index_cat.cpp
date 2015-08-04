// fastq_index_cat.c
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

#define LINE_BUFFER_SIZE	10240


void ctconv_f(char *buf, size_t size)
{
	size_t i;
	for (i=0; i<size && buf[i]!=0; ++i) {
		switch(buf[i]) {
		case 'C':
		case 'c':
			buf[i]='T';
			break;
		default:
			buf[i]=toupper(buf[i]);
		}
	}
}

void ctconv_r(char *buf, size_t size)
{
	size_t i;
	for (i=0; i<size && buf[i]!=0; ++i) {
		switch(buf[i]) {
		case 'G':
		case 'g':
			buf[i]='A';
			break;
		default:
			buf[i]=toupper(buf[i]);
		}
	}
}


void read_write(FILE *fp, char *line_buffer, long int &nline)
{
	if (fgets(line_buffer, sizeof(line_buffer), fp)!=NULL) {
		fputs(line_buffer, stdout);
		++nline;
	}
}


void read_write_ctconv_f(FILE *fp, char *line_buffer, long int &nline)
{
	if (fgets(line_buffer, sizeof(line_buffer), fp)!=NULL) {
			if (nline % 4 == 1) ctconv_f(line_buffer, LINE_BUFFER_SIZE);
		fputs(line_buffer, stdout);
		++nline;
	}
}


void read_write_ctconv_r(FILE *fp, char *line_buffer, long int &nline)
{
	if (fgets(line_buffer, sizeof(line_buffer), fp)!=NULL) {
			if (nline % 4 == 1) ctconv_r(line_buffer, LINE_BUFFER_SIZE);
		fputs(line_buffer, stdout);
		++nline;
	}
}


void usage(FILE *fp)
{
	fprintf(fp, "Print a speficied block from a given FASTQ file.\n");
	fprintf(fp, "   With -b or -r options supply conversion of nucleotides.\n");
	fprintf(fp, "Usage: fastq_index_cat [-b|-r] -n <int(0~)> [-i <fastq_index_file>] <fastq_file>\n\n");
	fprintf(fp, "   Options:\n");
	fprintf(fp, "      -b    Applies C->T conversion for output sequences. Cancelled with -r option.\n");
	fprintf(fp, "      -r    Applies G->A conversion for output sequences.\n");
	fprintf(fp, "      -n    Specifies the block to output (the first block is 0).\n");
	fprintf(fp, "      -i    Specifies the index file to read. Default index file name is \"<fastq_file>.index\".\n\n");
}


int main (int argc, char **argv)
{
	if (argc == 1) {
		usage(stdout);
		return(0);
	}

	int arg_opt;
	char index_file[2048] = { '\0' };
	char nth_str[16];
	unsigned int nth=0;
	void (*io_function)(FILE *, char*, long int&) = read_write;

	while ((arg_opt=getopt(argc, argv, "bri:n:"))!=-1) {
		switch(arg_opt) {
		case 'b':
			io_function = read_write_ctconv_f;
			break;
		case 'r':
			io_function = read_write_ctconv_r;
			break;
		case 'i':
			strncpy(index_file, optarg, sizeof(index_file)-1);
			break;
		case 'n':
			strncpy(nth_str, optarg, sizeof(nth_str)-1);
			break;
		default:
			fprintf(stderr, "Error: main: Unknown option speficied.\n");
			usage(stderr);
			return(1);
		}
	}

	if (sscanf(nth_str, "%d", &nth)!=1 || nth<0) {
		fprintf(stderr, "Error: main: The given argument \"%s\" for option -n <int(0~)> is invalid.\n", nth_str);
		return(1);
	}

	if (index_file[0]=='\0') {
		if (strlen(argv[optind]) + strlen(".index") >= sizeof(index_file)) {
			fprintf(stderr, "Error: main: Filename for the index file (length=%d) is long inter than the internal buffer (length=%d).\n", strlen(argv[optind])+strlen(".index"), sizeof(index_file));
			return(1);
		}
		sprintf(index_file, "%s.index", argv[optind]);
	}
	DEBUG("Index file is \"" << index_file << "\"");

	struct stat file_stat;
	int ret=stat(index_file, &file_stat);
	if (ret!=0) {
		fprintf(stderr, "Error: main: Unable to read the index file \"%s\".\n", index_file);
		return(1);
	}

	long int *index_table = NULL;
	index_table = (long int*)calloc(file_stat.st_size, sizeof(char));
	if (index_table==NULL) {
		fprintf(stderr, "Error: main: Memory allocation failed.\n");
		return(1);
	}

	long int index_raw_size = file_stat.st_size / (sizeof(long int));
	DEBUG("index_size = " << index_raw_size);

	FILE *fp=NULL;
	fp = fopen(index_file, "rb");
	fread(index_table, sizeof(char), file_stat.st_size, fp);
	fclose(fp);

    long int index_size = 0;
    for ( ; index_size < index_raw_size && index_table[index_size] != 0L; ++index_size) {}

    DEBUG("index_size = " << index_size);

	if (nth > index_size) {
		fprintf(stderr, "Error: main: The given argument for \"%d\" for option -n <int(0~)> is out of range for the given index (0~%ld).\n", nth, index_size);
		return(1);
	}

	if (nth > 0 && index_table[nth-1] == 0L) {
		fprintf(stderr ,"Error: main: The offset has invalid value (%ld)\n", index_table[nth-1]);
		return(1);
	}

	long int nline=0;
	char line_buffer[LINE_BUFFER_SIZE];
	fp = NULL;
	fp = fopen(argv[optind], "r");
	if (nth > 0) {
		fseek(fp, index_table[nth-1], SEEK_SET);
	}
	if (nth == index_size) {
		while (! feof(fp)) {
			io_function(fp, line_buffer, nline);
		}
	}
	else {
		DEBUG("ftell(fp)=" << ftell(fp) << " index_table[" << nth << "]=" << index_table[nth]);
		while (ftell(fp) < index_table[nth]) {
			io_function(fp, line_buffer, nline);
		}
	}
    fclose(fp);

	DEBUG(nline << " lines have been output.");

	return(0);
}

