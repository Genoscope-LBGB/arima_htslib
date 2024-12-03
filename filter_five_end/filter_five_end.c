#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

void set_unmapped_flag(bam1_t *record) {
    uint32_t flag = record->core.flag;
    flag |= BAM_FUNMAP;
    record->core.flag = flag;
}

int is_five_prime_match(bam1_t *record) {
    uint32_t *cigar = bam_get_cigar(record);
    uint32_t flag = record->core.flag;
    int n_cigar = record->core.n_cigar;

    // Determine if the alignment is forward or reverse
    int is_reverse = flag & BAM_FREVERSE;

    // Check the first or last CIGAR operation based on alignment direction
    if (!is_reverse) {
        // Forward strand: check if the first CIGAR operation is a match
        int op = bam_cigar_op(cigar[0]);
        return (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF);
    } else {
        // Reverse strand: check if the last CIGAR operation is a match
        int op = bam_cigar_op(cigar[n_cigar - 1]);
        return (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF);
    }
}

int is_three_prime_match(bam1_t *record) {
    uint32_t *cigar = bam_get_cigar(record);
    int n_cigar = record->core.n_cigar;
    int op = bam_cigar_op(cigar[n_cigar - 1]);

    // Check if the last CIGAR operation is a match
    return (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF);
}

void process_record_group(int counter, bam1_t *first_record, bam1_t *five_prime_record, samFile *out, bam_hdr_t *header) {
    if (counter == 1) {
        if (five_prime_record) {
            // If only one record in group and it has a five prime match, write it
            if (sam_write1(out, header, five_prime_record) < 0) {
                fprintf(stderr, "Error: unable to write record\n");
                exit(1);
            }
        } else if (first_record) {
            // Otherwise, write the first record with unmapped flag set
            set_unmapped_flag(first_record);
            if (sam_write1(out, header, first_record) < 0) {
                fprintf(stderr, "Error: unable to write record\n");
                exit(1);
            }
        }
    } else if (counter == 2 && five_prime_record) {
        // If two records and we have a five prime match, write it
        if (sam_write1(out, header, five_prime_record) < 0) {
            fprintf(stderr, "Error: unable to write record\n");
            exit(1);
        }
    } else {
        // Write the first record with unmapped flag set if no five prime match
        if (first_record) {
            set_unmapped_flag(first_record);
            if (sam_write1(out, header, first_record) < 0) {
                fprintf(stderr, "Error: unable to write record\n");
                exit(1);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    char *output_file = NULL;
    samFile *out = NULL;
    char prev_id[1024] = "";
    int counter = 0;

    bam1_t *first_record = NULL;
    bam1_t *five_prime_record = NULL;

    if (argc < 2 || argc > 4) {
        fprintf(stderr, "Usage: %s <in.bam> [-o <out.bam>]\n", argv[0]);
        return 1;
    }

    // Parse arguments
    int i;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (i == 1) {
            argv[1] = argv[i];
        }
    }

    // Open input BAM file
    samFile *in = sam_open(argv[1], "r");
    if (in == NULL) {
        fprintf(stderr, "Error: unable to open input BAM file %s\n", argv[1]);
        return 1;
    }

    // Read the header
    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL) {
        fprintf(stderr, "Error: unable to read header from %s\n", argv[1]);
        sam_close(in);
        return 1;
    }

    // Open output BAM file or use stdout
    if (output_file) {
        out = sam_open(output_file, "wb");
        if (out == NULL) {
            fprintf(stderr, "Error: unable to open output BAM file %s\n", output_file);
            bam_hdr_destroy(header);
            sam_close(in);
            return 1;
        }
    } else {
        out = sam_open("-", "w");
    }

    // Write the header to the output file
    if (sam_hdr_write(out, header) != 0) {
        fprintf(stderr, "Error: unable to write header to output\n");
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    // Initialize the record structure
    bam1_t *record = bam_init1();
    if (record == NULL) {
        fprintf(stderr, "Error: unable to initialize BAM record\n");
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);
        return 1;
    }

    while (sam_read1(in, header, record) >= 0) {
        char *current_id = bam_get_qname(record);

        if (strcmp(prev_id, current_id) != 0 && prev_id[0] != '\0') {
            // Process previous group of records
            process_record_group(counter, first_record, five_prime_record, out, header);

            // Clean up previous records
            if (first_record) {
                bam_destroy1(first_record);
                first_record = NULL;
            }
            if (five_prime_record) {
                bam_destroy1(five_prime_record);
                five_prime_record = NULL;
            }

            counter = 0;
        }

        // Store records for the current group
        if (counter == 0) {
            first_record = bam_dup1(record);
        }
        if (is_five_prime_match(record) && five_prime_record == NULL) {
            five_prime_record = bam_dup1(record);
        }

        strcpy(prev_id, current_id);
        counter++;
    }

    // Final record processing
    process_record_group(counter, first_record, five_prime_record, out, header);

    if (first_record) {
        bam_destroy1(first_record);
    }
    if (five_prime_record) {
        bam_destroy1(five_prime_record);
    }

    // Clean up
    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    return 0;
}
