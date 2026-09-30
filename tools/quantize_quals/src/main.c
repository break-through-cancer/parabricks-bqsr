#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "quant.h"

#define QQ_VERSION "0.1.0"
#define QQ_MAX_BINS 256

static const char *usage =
    "Usage: quantize_quals --in <sam|bam|cram> --out <bam|cram|sam> [--ref <fasta>]\n"
    "                      --static-quantized-quals <Q> [<Q> ...]\n"
    "                      [--preserve-qscores-less-than <Q>] [--round-down-quantized]\n"
    "                      [--threads <N>]\n"
    "\n"
    "Replicates GATK ApplyBQSR static quality-score quantization on an already\n"
    "recalibrated alignment file. Only per-base QUAL values are changed.\n"
    "\n"
    "  --in <path>                        Input SAM/BAM/CRAM.\n"
    "  --out <path>                       Output; format from extension (.bam, .cram, .sam).\n"
    "                                     BAM/CRAM output of coordinate-sorted input is\n"
    "                                     indexed (<out>.bai / <out>.crai).\n"
    "  --ref <fasta>                      Reference; required when input or output is CRAM.\n"
    "  --static-quantized-quals <Q...>    Static bins; values may follow one flag or the\n"
    "                                     flag may repeat. Required.\n"
    "  --preserve-qscores-less-than <Q>   Qualities below Q are unchanged (default: 6).\n"
    "  --round-down-quantized             Round down to the largest bin <= Q instead of\n"
    "                                     the nearest bin in probability space.\n"
    "  --threads <N>                      BGZF/CRAM worker threads (default: 1).\n"
    "  --help, --version\n";

typedef struct {
    const char *in, *out, *ref;
    int bins[QQ_MAX_BINS];
    size_t n_bins;
    int preserve;
    int round_down;
    int threads;
} opts_t;

static int die(const char *fmt, const char *arg) {
    fputs("quantize_quals: ", stderr);
    fprintf(stderr, fmt, arg);
    fputc('\n', stderr);
    return 1;
}

static int parse_int(const char *s, int *out) {
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || end == s || *end != '\0' || v < INT_MIN || v > INT_MAX) return -1;
    *out = (int)v;
    return 0;
}

static int ends_with(const char *s, const char *suffix) {
    size_t ls = strlen(s), lx = strlen(suffix);
    return ls >= lx && strcmp(s + ls - lx, suffix) == 0;
}

static int parse_args(int argc, char **argv, opts_t *o) {
    int saw_bins_flag = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        int needs_value = !strcmp(a, "--in") || !strcmp(a, "--out") || !strcmp(a, "--ref") ||
                          !strcmp(a, "--threads") || !strcmp(a, "--preserve-qscores-less-than");
        if (needs_value && (i + 1 >= argc || !strncmp(argv[i + 1], "--", 2)))
            return die("%s requires a value", a);

        if (!strcmp(a, "--help") || !strcmp(a, "-h")) {
            fputs(usage, stdout);
            exit(0);
        } else if (!strcmp(a, "--version")) {
            printf("quantize_quals %s (htslib %s)\n", QQ_VERSION, hts_version());
            exit(0);
        } else if (!strcmp(a, "--in")) {
            o->in = argv[++i];
        } else if (!strcmp(a, "--out")) {
            o->out = argv[++i];
        } else if (!strcmp(a, "--ref")) {
            o->ref = argv[++i];
        } else if (!strcmp(a, "--threads")) {
            if (parse_int(argv[++i], &o->threads) || o->threads < 1)
                return die("--threads must be a positive integer, got '%s'", argv[i]);
        } else if (!strcmp(a, "--preserve-qscores-less-than")) {
            if (parse_int(argv[++i], &o->preserve))
                return die("--preserve-qscores-less-than value '%s' is not an integer", argv[i]);
        } else if (!strcmp(a, "--round-down-quantized")) {
            o->round_down = 1;
        } else if (!strcmp(a, "--static-quantized-quals")) {
            saw_bins_flag = 1;
            size_t before = o->n_bins;
            while (i + 1 < argc && strncmp(argv[i + 1], "--", 2)) {
                if (o->n_bins == QQ_MAX_BINS) return die("too many --static-quantized-quals values%s", "");
                if (parse_int(argv[++i], &o->bins[o->n_bins]))
                    return die("--static-quantized-quals value '%s' is not an integer", argv[i]);
                o->n_bins++;
            }
            if (o->n_bins == before) return die("%s requires at least one quality value", a);
        } else {
            return die("unknown option '%s' (see --help)", a);
        }
    }
    if (!o->in) return die("missing required option %s", "--in");
    if (!o->out) return die("missing required option %s", "--out");
    if (!saw_bins_flag) return die("missing required option %s", "--static-quantized-quals");
    return 0;
}

static char *command_line(int argc, char **argv) {
    kstring_t ks = KS_INITIALIZE;
    for (int i = 0; i < argc; i++) {
        if (i) kputc(' ', &ks);
        kputs(argv[i], &ks);
    }
    return ks_release(&ks);
}

static int is_coordinate_sorted(sam_hdr_t *hdr) {
    kstring_t so = KS_INITIALIZE;
    int sorted = sam_hdr_find_tag_hd(hdr, "SO", &so) == 0 && !strcmp(ks_str(&so), "coordinate");
    ks_free(&so);
    return sorted;
}

int main(int argc, char **argv) {
    opts_t o = { .preserve = 6, .threads = 1 };
    if (argc == 1) {
        fputs(usage, stderr);
        return 1;
    }
    if (parse_args(argc, argv, &o)) return 1;

    uint8_t mapping[256];
    char err[512];
    if (qq_build_mapping(o.preserve, o.bins, o.n_bins, o.round_down, mapping, err, sizeof err))
        return die("%s", err);

    const char *mode;
    const char *idx_ext = NULL;
    int out_is_cram = 0;
    if (ends_with(o.out, ".bam")) {
        mode = "wb";
        idx_ext = ".bai";
    } else if (ends_with(o.out, ".cram")) {
        mode = "wc";
        idx_ext = ".crai";
        out_is_cram = 1;
    } else if (ends_with(o.out, ".sam")) {
        mode = "w";
    } else {
        return die("--out '%s' must have a .bam, .cram or .sam extension", o.out);
    }
    if (out_is_cram && !o.ref) return die("--ref is required for CRAM output '%s'", o.out);

    samFile *in = sam_open(o.in, "r");
    if (!in) return die("cannot open input '%s'", o.in);
    if (hts_get_format(in)->format == cram && !o.ref) {
        sam_close(in);
        return die("--ref is required for CRAM input '%s'", o.in);
    }

    int rc = 1;
    samFile *out = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = NULL;
    char *cl = NULL;
    char *idx_fn = NULL;
    hts_tpool *pool = NULL;

    if (o.threads > 1) {
        pool = hts_tpool_init(o.threads);
        if (!pool) {
            die("failed to start thread pool%s", "");
            goto done;
        }
        htsThreadPool tp = { pool, 0 };
        hts_set_opt(in, HTS_OPT_THREAD_POOL, &tp);
    }
    if (o.ref && hts_set_fai_filename(in, o.ref) < 0) {
        die("cannot load reference '%s'", o.ref);
        goto done;
    }

    hdr = sam_hdr_read(in);
    if (!hdr) {
        die("cannot read header from '%s'", o.in);
        goto done;
    }
    cl = command_line(argc, argv);
    if (sam_hdr_add_pg(hdr, "quantize_quals", "VN", QQ_VERSION, "CL", cl, NULL) < 0) {
        die("cannot add @PG line to header%s", "");
        goto done;
    }

    out = sam_open(o.out, mode);
    if (!out) {
        die("cannot open output '%s'", o.out);
        goto done;
    }
    if (pool) {
        htsThreadPool tp = { pool, 0 };
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &tp);
    }
    if (o.ref && hts_set_fai_filename(out, o.ref) < 0) {
        die("cannot load reference '%s'", o.ref);
        goto done;
    }
    if (sam_hdr_write(out, hdr) < 0) {
        die("cannot write header to '%s'", o.out);
        goto done;
    }

    if (idx_ext) {
        if (is_coordinate_sorted(hdr)) {
            idx_fn = malloc(strlen(o.out) + strlen(idx_ext) + 1);
            strcpy(idx_fn, o.out);
            strcat(idx_fn, idx_ext);
            if (sam_idx_init(out, hdr, 0, idx_fn) < 0) {
                die("cannot initialise index '%s'", idx_fn);
                goto done;
            }
        } else {
            fprintf(stderr, "quantize_quals: warning: input '%s' is not coordinate-sorted "
                            "(@HD SO:coordinate); output will not be indexed\n", o.in);
        }
    }

    b = bam_init1();
    int r;
    while ((r = sam_read1(in, hdr, b)) >= 0) {
        qq_apply(mapping, bam_get_qual(b), (size_t)b->core.l_qseq);
        if (sam_write1(out, hdr, b) < 0) {
            die("failed writing a record to '%s'", o.out);
            goto done;
        }
    }
    if (r < -1) {
        die("failed reading a record from '%s' (truncated or corrupt input?)", o.in);
        goto done;
    }
    if (idx_fn && sam_idx_save(out) < 0) {
        die("cannot write index '%s'", idx_fn);
        goto done;
    }
    rc = 0;

done:
    if (b) bam_destroy1(b);
    if (out && sam_close(out) < 0 && rc == 0) rc = die("failed closing output '%s'", o.out);
    if (in && sam_close(in) < 0 && rc == 0) rc = die("failed closing input '%s'", o.in);
    if (hdr) sam_hdr_destroy(hdr);
    if (pool) hts_tpool_destroy(pool);
    free(cl);
    free(idx_fn);
    return rc;
}
