#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <htslib/sam.h>

void print_usage() {
    fprintf(stderr, "Usage: ./two_read_bam_combiner <read 1 bam> <read 2 bam> <minimum map quality filter> [-o output.bam]\n");
}

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 6) {
        print_usage();
        return 1;
    }

    char *read1_bam = argv[1];
    char *read2_bam = argv[2];
    int min_map_quality = atoi(argv[3]);
    char *output_bam = NULL;

    if (argc == 6 && strcmp(argv[4], "-o") == 0) {
        output_bam = argv[5];
    }

    samFile *fp1 = sam_open(read1_bam, "r");
    samFile *fp2 = sam_open(read2_bam, "r");
    if (!fp1 || !fp2) {
        fprintf(stderr, "Failed to open BAM files.\n");
        return 1;
    }

    samFile *output = output_bam ? sam_open(output_bam, "w") : sam_open("-", "w");
    if (!output) {
        fprintf(stderr, "Failed to open output stream.\n");
        return 1;
    }

    bam_hdr_t *header1 = sam_hdr_read(fp1);
    bam_hdr_t *header2 = sam_hdr_read(fp2);
    if (!header1 || !header2) {
        fprintf(stderr, "Failed to read headers from BAM files.\n");
        return 1;
    }

    // Check if sequence headers (SQ lines) are identical
    int n_targets1 = header1->n_targets;
    int n_targets2 = header2->n_targets;
    if (n_targets1 != n_targets2) {
        fprintf(stderr, "Inconsistent number of sequences in BAM headers. BAM files must be aligned to the same reference.\n");
        return 1;
    }

    int i;
    for (i = 0; i < n_targets1; i++) {
        if (strcmp(header1->target_name[i], header2->target_name[i]) != 0) {
            fprintf(stderr, "Inconsistent sequence names in BAM headers. BAM files must be aligned to the same reference.\n");
            return 1;
        }
    }

    if (sam_hdr_write(output, header1) < 0) {
        fprintf(stderr, "Failed to write header to output stream.\n");
        return 1;
    }

    bam1_t *aln1 = bam_init1();
    bam1_t *aln2 = bam_init1();
    int counter = 0;
    int new_counter = 0;

    while (sam_read1(fp1, header1, aln1) >= 0 && sam_read1(fp2, header2, aln2) >= 0) {
        counter++;
        if (counter == (new_counter + 1000000)) {
            fprintf(stderr, "%d\n", counter);
            new_counter = counter;
        }

        // Check if read names are the same
        if (strcmp(bam_get_qname(aln1), bam_get_qname(aln2)) != 0) {
            fprintf(stderr, "The read ids of the two files do not match up at line number %d.\n", counter);
            return 1;
        }

        // Get mapping quality and flag
        int flag1 = aln1->core.flag;
        int flag2 = aln2->core.flag;
        int mapq1 = aln1->core.qual;
        int mapq2 = aln2->core.qual;

        // Skip reads that are unmapped or below minimum map quality
        if ((flag1 & BAM_FUNMAP) || (flag2 & BAM_FUNMAP) || (mapq1 < min_map_quality) || (mapq2 < min_map_quality)) {
            continue;
        }

        // Determine distance between reads if on the same chromosome
        int32_t dist1 = 0, dist2 = 0;
        if (aln1->core.tid == aln2->core.tid) {
            int32_t dist = abs(aln1->core.pos - aln2->core.pos);
            if (aln1->core.pos >= aln2->core.pos) {
                dist1 = -dist;
                dist2 = dist;
            } else {
                dist1 = dist;
                dist2 = -dist;
            }
        }

        // Update mate information in both alignments
        aln1->core.mtid = aln2->core.tid;
        aln1->core.mpos = aln2->core.pos;
        aln1->core.isize = dist1;

        aln2->core.mtid = aln1->core.tid;
        aln2->core.mpos = aln1->core.pos;
        aln2->core.isize = dist2;

        // Set proper pair flags, first/second in pair flags, mate strand flags, and clear supplementary alignment flag
        aln1->core.flag |= BAM_FPROPER_PAIR | BAM_FREAD1 | BAM_FPAIRED;
        aln2->core.flag |= BAM_FPROPER_PAIR | BAM_FREAD2 | BAM_FPAIRED;

        aln1->core.flag &= ~BAM_FSUPPLEMENTARY;
        aln2->core.flag &= ~BAM_FSUPPLEMENTARY;

        if (flag2 & BAM_FREVERSE) {
            aln1->core.flag |= BAM_FMREVERSE;
        }
        if (flag1 & BAM_FREVERSE) {
            aln2->core.flag |= BAM_FMREVERSE;
        }

        if (sam_write1(output, header1, aln1) < 0) {
            fprintf(stderr, "Failed to write alignment 1 to output stream.\n");
            return 1;
        }
        if (sam_write1(output, header2, aln2) < 0) {
            fprintf(stderr, "Failed to write alignment 2 to output stream.\n");
            return 1;
        }
    }

    bam_destroy1(aln1);
    bam_destroy1(aln2);
    bam_hdr_destroy(header1);
    bam_hdr_destroy(header2);
    sam_close(fp1);
    sam_close(fp2);
    sam_close(output);

    return 0;
}
