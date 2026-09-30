#include <stdio.h>
#include <string.h>

#include "quant.h"

static int failures = 0;

#define CHECK(cond, ...) do { \
    if (!(cond)) { failures++; fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } } while (0)

static int expected_nearest(int q) {
    if (q < 6) return q;
    if (q <= 7) return 6;
    if (q <= 12) return 10;
    if (q <= 22) return 20;
    return 30;
}

static int expected_round_down(int q) {
    if (q < 6) return q;
    if (q <= 9) return 6;
    if (q <= 19) return 10;
    if (q <= 29) return 20;
    return 30;
}

static void test_nearest_probability_space(void) {
    const int bins[] = {10, 20, 30};
    uint8_t m[256];
    char err[256];
    CHECK(qq_build_mapping(6, bins, 3, 0, m, err, sizeof err) == 0, "build failed: %s", err);
    for (int q = 0; q <= 254; q++)
        CHECK(m[q] == expected_nearest(q), "nearest q=%d got %d want %d", q, m[q], expected_nearest(q));
}

static void test_round_down(void) {
    const int bins[] = {10, 20, 30};
    uint8_t m[256];
    char err[256];
    CHECK(qq_build_mapping(6, bins, 3, 1, m, err, sizeof err) == 0, "build failed: %s", err);
    for (int q = 0; q <= 254; q++)
        CHECK(m[q] == expected_round_down(q), "round-down q=%d got %d want %d", q, m[q], expected_round_down(q));
}

static void test_missing_qual_passthrough(void) {
    const int bins[] = {10, 20, 30};
    uint8_t m[256];
    char err[256];
    qq_build_mapping(6, bins, 3, 0, m, err, sizeof err);
    CHECK(m[QQ_MISSING_QUAL] == QQ_MISSING_QUAL, "0xff mapped to %d", m[QQ_MISSING_QUAL]);

    uint8_t missing[4] = {0xff, 0xff, 0xff, 0xff};
    qq_apply(m, missing, 4);
    for (int i = 0; i < 4; i++) CHECK(missing[i] == 0xff, "missing qual altered at %d", i);
}

static void test_unsorted_duplicate_bins(void) {
    const int sorted[] = {10, 20, 30};
    const int messy[] = {30, 10, 20, 20, 10};
    uint8_t a[256], b[256];
    char err[256];
    CHECK(qq_build_mapping(6, sorted, 3, 0, a, err, sizeof err) == 0, "%s", err);
    CHECK(qq_build_mapping(6, messy, 5, 0, b, err, sizeof err) == 0, "%s", err);
    CHECK(memcmp(a, b, sizeof a) == 0, "unsorted/duplicate bins changed the mapping");
}

static void test_apply(void) {
    const int bins[] = {10, 20, 30};
    uint8_t m[256];
    char err[256];
    qq_build_mapping(6, bins, 3, 0, m, err, sizeof err);
    uint8_t q[] = {2, 7, 8, 13, 23, 41};
    const uint8_t want[] = {2, 6, 10, 20, 30, 30};
    qq_apply(m, q, sizeof q);
    CHECK(memcmp(q, want, sizeof q) == 0, "qq_apply produced wrong values");
}

static void test_validation(void) {
    uint8_t m[256];
    char err[256];
    const int bins[] = {10, 20, 30};
    const int neg[] = {-1, 20};
    const int big[] = {10, 94};

    CHECK(qq_build_mapping(6, bins, 0, 0, m, err, sizeof err) == -1, "empty bins accepted");
    CHECK(strstr(err, "--static-quantized-quals") != NULL, "empty-bins message: %s", err);

    CHECK(qq_build_mapping(10, bins, 3, 0, m, err, sizeof err) == -1, "preserve == min bin accepted");
    CHECK(qq_build_mapping(11, bins, 3, 0, m, err, sizeof err) == -1, "preserve > min bin accepted");
    CHECK(strstr(err, "--preserve-qscores-less-than") != NULL, "preserve message: %s", err);

    CHECK(qq_build_mapping(-1, bins, 3, 0, m, err, sizeof err) == -1, "negative preserve accepted");
    CHECK(qq_build_mapping(0, neg, 2, 0, m, err, sizeof err) == -1, "negative bin accepted");
    CHECK(qq_build_mapping(6, big, 2, 0, m, err, sizeof err) == -1, "bin > 93 accepted");

    CHECK(qq_build_mapping(0, bins, 3, 0, m, err, sizeof err) == 0, "preserve 0 rejected: %s", err);
    CHECK(m[0] == 0 && m[2] == 0 && m[3] == 10, "preserve 0 mapping wrong: m[2]=%d m[3]=%d", m[2], m[3]);
}

int main(void) {
    test_nearest_probability_space();
    test_round_down();
    test_missing_qual_passthrough();
    test_unsorted_duplicate_bins();
    test_apply();
    test_validation();
    if (failures) {
        fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    printf("test_quant: all checks passed\n");
    return 0;
}
