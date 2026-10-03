#include "quant.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double p_correct(int q) {
    return 1.0 - pow(10.0, -q / 10.0);
}

static int cmp_int(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

int qq_build_mapping(int preserve_less_than, const int *bins, size_t n_bins,
                     int round_down, uint8_t mapping[256],
                     char *err, size_t err_len) {
    if (n_bins == 0) {
        snprintf(err, err_len, "--static-quantized-quals requires at least one quality value");
        return -1;
    }
    if (preserve_less_than < 0 || preserve_less_than > QQ_MAX_PHRED) {
        snprintf(err, err_len, "--preserve-qscores-less-than must be between 0 and %d, got %d",
                 QQ_MAX_PHRED, preserve_less_than);
        return -1;
    }

    int *bounds = malloc((n_bins + 1) * sizeof *bounds);
    if (!bounds) {
        snprintf(err, err_len, "out of memory");
        return -1;
    }
    for (size_t i = 0; i < n_bins; i++) {
        if (bins[i] < 0 || bins[i] > QQ_MAX_PHRED) {
            snprintf(err, err_len, "--static-quantized-quals values must be between 0 and %d, got %d",
                     QQ_MAX_PHRED, bins[i]);
            free(bounds);
            return -1;
        }
        bounds[i + 1] = bins[i];
    }
    qsort(bounds + 1, n_bins, sizeof *bounds, cmp_int);
    if (preserve_less_than >= bounds[1]) {
        snprintf(err, err_len,
                 "--preserve-qscores-less-than (%d) must be strictly below the smallest "
                 "--static-quantized-quals value (%d)", preserve_less_than, bounds[1]);
        free(bounds);
        return -1;
    }
    bounds[0] = preserve_less_than;

    size_t n = 1;
    for (size_t i = 1; i <= n_bins; i++)
        if (bounds[i] != bounds[n - 1]) bounds[n++] = bounds[i];

    for (int q = 0; q < QQ_MISSING_QUAL; q++) {
        if (q < preserve_less_than) {
            mapping[q] = (uint8_t)q;
            continue;
        }
        int best = bounds[0];
        if (round_down) {
            for (size_t i = 0; i < n && bounds[i] <= q; i++) best = bounds[i];
        } else {
            double pq = p_correct(q), best_d = fabs(pq - p_correct(bounds[0]));
            for (size_t i = 1; i < n; i++) {
                double d = fabs(pq - p_correct(bounds[i]));
                if (d < best_d) {
                    best_d = d;
                    best = bounds[i];
                }
            }
        }
        mapping[q] = (uint8_t)best;
    }
    mapping[QQ_MISSING_QUAL] = QQ_MISSING_QUAL;

    free(bounds);
    return 0;
}

void qq_apply(const uint8_t mapping[256], uint8_t *qual, size_t len, uint64_t counts[256]) {
    if (len == 0 || qual[0] == QQ_MISSING_QUAL) return;
    if (counts) {
        for (size_t i = 0; i < len; i++) {
            counts[qual[i]]++;
            qual[i] = mapping[qual[i]];
        }
    } else {
        for (size_t i = 0; i < len; i++) qual[i] = mapping[qual[i]];
    }
}

int qq_describe_mapping(const uint8_t mapping[256], int preserve_less_than, char *buf, size_t buf_len) {
    size_t used = 0;
    buf[0] = '\0';
    int start = preserve_less_than;
    for (int q = preserve_less_than; q <= QQ_MAX_PHRED; q++) {
        if (q < QQ_MAX_PHRED && mapping[q + 1] == mapping[start]) continue;
        int n;
        if (q == QQ_MAX_PHRED)
            n = snprintf(buf + used, buf_len - used, "%s%d+->%d", used ? " " : "", start, mapping[start]);
        else if (q == start)
            n = snprintf(buf + used, buf_len - used, "%s%d->%d", used ? " " : "", start, mapping[start]);
        else
            n = snprintf(buf + used, buf_len - used, "%s%d-%d->%d", used ? " " : "", start, q, mapping[start]);
        if (n < 0 || (size_t)n >= buf_len - used) return -1;
        used += (size_t)n;
        start = q + 1;
    }
    return 0;
}
