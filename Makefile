CPP=c++
CFLAGS=-DDEBUG_BUILD -g

all: fastq_index fastq_index_cat fastq_index_dump fastq_index_size

fastq_index: fastq_index.cpp
	${CPP} ${CFLAGS} -o fastq_index fastq_index.cpp

fastq_index_cat: fastq_index_cat.cpp
	${CPP} ${CFLAGS} -o fastq_index_cat fastq_index_cat.cpp

fastq_index_dump: fastq_index_dump.cpp
	${CPP} ${CFLAGS} -o fastq_index_dump fastq_index_dump.cpp

fastq_index_size: fastq_index_size.cpp
	${CPP} ${CFLAGS} -o fastq_index_size fastq_index_size.cpp

clean:
	rm -f *.o fastq_index fastq_index_cat fastq_index_dump fastq_index_size

