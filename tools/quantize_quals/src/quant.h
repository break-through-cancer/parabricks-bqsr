#ifndef QUANTIZE_QUALS_QUANT_H
#define QUANTIZE_QUALS_QUANT_H

#include <stddef.h>
#include <stdint.h>

#define QQ_MAX_PHRED 93
#define QQ_MISSING_QUAL 0xff

int qq_build_mapping(int preserve_less_than, const int *bins, size_t n_bins,
                     int round_down, uint8_t mapping[256],
                     char *err, size_t err_len);

void qq_apply(const uint8_t mapping[256], uint8_t *qual, size_t len);

#endif
