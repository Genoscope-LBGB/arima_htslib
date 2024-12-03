#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Please provide the <bamfile>\n");
        return 1;
    }
    const char *bam_file = argv[1];

    // Open the BAM file
    htsFile *fp = hts_open(bam_file, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open BAM file %s\n", bam_file);
        return 1;
    }

    // Read the header
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Failed to read header from BAM file %s\n", bam_file);
        hts_close(fp);
        return 1;
    }

    bam1_t *rec = bam_init1();

    // Declare counters
    long all = 0, intra = 0, inter = 0;
    long intra_1 = 0, intra_10 = 0, intra_15 = 0, intra_20 = 0;

    // Counters for valid pairs
    long total_valid_pairs = 0;
    long *valid_pairs_per_seq = calloc(hdr->n_targets, sizeof(long));
    if (!valid_pairs_per_seq) {
        fprintf(stderr, "Memory allocation failed\n");
        bam_destroy1(rec);
        bam_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    // Variables for sequence lengths in base pairs
    long total_assembly_length = 0;
    long *sequence_length_bp = calloc(hdr->n_targets, sizeof(long));
    if (!sequence_length_bp) {
        fprintf(stderr, "Memory allocation failed\n");
        free(valid_pairs_per_seq);
        bam_destroy1(rec);
        bam_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    // Loop variable declaration
    int i;

    // Calculate total assembly length and sequence lengths in base pairs
    for (i = 0; i < hdr->n_targets; i++) {
        total_assembly_length += hdr->target_len[i];
        sequence_length_bp[i] = hdr->target_len[i];
    }

    // Iterate over the records
    while (sam_read1(fp, hdr, rec) >= 0) {
        all++;

        int32_t tid = rec->core.tid;
        int32_t mtid = rec->core.mtid;
        int32_t isize = rec->core.isize;

        // Check if the read is properly paired and is the first in pair
        if ((rec->core.flag & BAM_FREAD1) && (rec->core.flag & BAM_FPROPER_PAIR)) {
            total_valid_pairs++;
            if (tid >= 0 && tid < hdr->n_targets) {
                valid_pairs_per_seq[tid]++;
            }
        }

        // Check if the mate is on the same chromosome
        if (mtid == tid && mtid != -1) {
            intra++;
            int32_t abs_isize;  // Declare variable at the beginning of the block
            abs_isize = isize >= 0 ? isize : -isize;  // Calculate absolute insert size
            if (abs_isize >= 1000) intra_1++;
            if (abs_isize >= 10000) intra_10++;
            if (abs_isize >= 15000) intra_15++;
            if (abs_isize >= 20000) intra_20++;
        } else {
            inter++;
        }
    }

    // Divide the counters by 2
    all      = all / 2;
    intra    = intra / 2;
    intra_1  = intra_1 / 2;
    intra_10 = intra_10 / 2;
    intra_15 = intra_15 / 2;
    intra_20 = intra_20 / 2;
    inter    = inter / 2;

    // Calculate valid pairs per Mb
    double valid_pairs_per_Mb_total = 0.0;
    if (total_assembly_length > 0) {
        valid_pairs_per_Mb_total = total_valid_pairs / ((double)total_assembly_length / 1e6);
    }
    double *valid_pairs_per_Mb_seq = calloc(hdr->n_targets, sizeof(double));
    if (!valid_pairs_per_Mb_seq) {
        fprintf(stderr, "Memory allocation failed\n");
        free(sequence_length_bp);
        free(valid_pairs_per_seq);
        bam_destroy1(rec);
        bam_hdr_destroy(hdr);
        hts_close(fp);
        return 1;
    }

    for (i = 0; i < hdr->n_targets; i++) {
        if (sequence_length_bp[i] > 0) {
            valid_pairs_per_Mb_seq[i] = valid_pairs_per_seq[i] / ((double)sequence_length_bp[i] / 1e6);
        } else {
            valid_pairs_per_Mb_seq[i] = 0.0;
        }
    }

    // Print the results
    printf("All\t%ld\n", all);
    printf("All intra\t%ld\n", intra);
    printf("All intra 1kb\t%ld\n", intra_1);
    printf("All intra 10kb\t%ld\n", intra_10);
    printf("All intra 15kb\t%ld\n", intra_15);
    printf("All intra 20kb\t%ld\n", intra_20);
    printf("All inter\t%ld\n", inter);

    // Print total number of valid pairs
    printf("Total valid pairs over the whole assembly: %ld\n", total_valid_pairs);
    printf("Valid pairs per Mb over the whole assembly: %.2f\n", valid_pairs_per_Mb_total);

    // Print header
    printf("Sequence\tValid Pairs\tLength (bp)\tValid Pairs per Mb\n");

    // Print number of valid pairs per sequence on a single line
    for (i = 0; i < hdr->n_targets; i++) {
        printf("%s\t%ld\t%ld\t%.2f\n",
               hdr->target_name[i],
               valid_pairs_per_seq[i],
               sequence_length_bp[i],
               valid_pairs_per_Mb_seq[i]);
    }

    // Clean up
    free(valid_pairs_per_Mb_seq);
    free(sequence_length_bp);
    free(valid_pairs_per_seq);
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    hts_close(fp);

    return 0;
}
