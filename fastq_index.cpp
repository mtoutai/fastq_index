// fastq_index.c
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define max(x,y) (x>y)?x:y

#ifdef DEBUG_BUILD
#define DEBUG(x) do { std::cerr << "Debug: " << __FILE__ << ": " << __LINE__ << ": " << x << std::endl; } while(0)
#else
#define DEBUG(x)
#endif /* DEBUG_BUILD */

#define DEFAULT_STEP 1000000
#define BLOCK_SIZE 4096


int estimate_seq_len (FILE *fp, size_t *head_len, size_t *seq_len)
{
	rewind (fp);
	int c=0;
	for (c=fgetc (fp); !feof (fp) & c!='\n'; c=fgetc (fp), ++(*head_len));
	for (c=fgetc (fp); !feof (fp) & c!='\n'; c=fgetc (fp), ++(*seq_len));
	rewind (fp);
	++(*head_len);	// counting '\n'
	++(*seq_len);	// counting '\n'
	return (*seq_len);
}


void usage(FILE *fp)
{
	fprintf (fp, "Create the index file for a given FASTQ file.\n");
	fprintf (fp, "Usage: fastq_index [-i <fastq index file>] [-s <number of items per block> (default: 1000000)] <fastq_file>\n");
}

int main (int argc, char **argv)
{
	if (argc==1) {
		usage (stdout);
		return (0);
	}

	char index_file[4096] = { '\0' };
	char step_str[256] = { '\0' };
	int arg_opt;
	while ((arg_opt=getopt(argc, argv, "i:s:"))!=-1) {
		switch(arg_opt) {
		case 'i':
			strncpy(index_file, optarg, sizeof(index_file)-1);
			break;
		case 's':
			strncpy(step_str, optarg, sizeof(step_str)-1);
			break;
		default:
			fprintf(stderr, "Error: main: Unknown option speficied.\n");
			usage(stderr);
			return(1);
		}
	}

	FILE *fp=NULL;
	fp=fopen (argv[optind], "r");
	if (fp==NULL) {
		fprintf (stderr, "Error: main: \"%s\" is not found.\n", argv[optind]);
		return (1);
	}

	int step=DEFAULT_STEP;
	if (step_str[0] != '\0') {
		if (sscanf(step_str, "%d", &step) != 1) {
			fprintf (stderr, "Error: main: \"-s %s\" has invalid argument.\n", step_str);
			return(1);
		}
	}
	int step_lines=step * 4;

	struct stat file_stat;
	stat (argv[optind], &file_stat);

	size_t head_len=0;
	size_t seq_len=0;

	estimate_seq_len (fp, &head_len, &seq_len);
	DEBUG("header_length=" << head_len << " sequence_length=" << seq_len);

	float estimated_lines=file_stat.st_size / (head_len + seq_len * 2 + 2);
	DEBUG("estimated_lines=" << estimated_lines);

	if (estimated_lines < step_lines) {
		fprintf (stderr, "Error: main: \"-s %s\" has too larger value than estimated number of lines (%.1f).\n", step_str, estimated_lines);
		return (1);
	}

	int buffer_size=(max (head_len, seq_len))*2;
	char *buffer=NULL;
	buffer=(char *)calloc (buffer_size, sizeof (char));
	if (buffer==NULL) {
		fprintf (stderr, "Error: main: memory allocation error.\n");
		return (1);
	}
	DEBUG("buffer_size=" << buffer_size);
	long int *index=NULL;	// why "long int" because ftell returns a "long int" value.
	size_t index_lines=BLOCK_SIZE/sizeof (long int);
	index=(long int*)calloc (index_lines, sizeof (long int));
	if (index==NULL) {
		fprintf (stderr, "Error: main: memory allocation error.\n");
		return (1);
	}
	DEBUG("index_lines=" << index_lines);

	if (index_file[0] == '\0') {
		sprintf (index_file, "%s.index", argv[optind]);
	}
	DEBUG("index_file_name=" << index_file);
	FILE *outfp=NULL;
	outfp=fopen (index_file, "wb");
	size_t line_no=0;
	size_t ii=0;
	while (!feof (fp)) {
		fgets (buffer, buffer_size, fp);
		++line_no;
		if (buffer[buffer_size-1] != 0 & buffer[buffer_size-1] != '\n') {
			fprintf (stderr, "Error: main: the line no. %d is longer than the buffer size (%d).\n", line_no, buffer_size);
			return (1);
		}
		if (line_no % step_lines==0) {
			if (ii>=index_lines) {
				fprintf (stderr, "Error: main: ii=%d >= index_lines=%d\n", ii, index_lines);
				return (-1);
			}
			index[ii++] = ftell (fp);
			if (ii==index_lines) {
				ii=0;
				fwrite (index, sizeof (long int), index_lines, outfp);
				memset (index, index_lines, sizeof (long int));
			}
		}
	}
	fclose (fp);
	if (ii>0) {
		fwrite (index, sizeof (long int), index_lines, outfp);
	}
	fclose (outfp);

	return (0);
}

