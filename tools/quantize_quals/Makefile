CC      ?= cc
CFLAGS  ?= -O2 -g
CFLAGS  += -std=c99 -Wall -Wextra -Werror -Isrc

HTSLIB_CFLAGS ?= $(shell pkg-config --cflags htslib 2>/dev/null)
HTSLIB_LIBS   ?= $(shell pkg-config --libs htslib 2>/dev/null || echo -lhts)

BUILD  := build
PREFIX ?= /usr/local

.PHONY: all test clean install

all: $(BUILD)/quantize_quals

$(BUILD):
	mkdir -p $@

$(BUILD)/quant.o: src/quant.c src/quant.h | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/main.o: src/main.c src/quant.h | $(BUILD)
	$(CC) $(CFLAGS) $(HTSLIB_CFLAGS) -c $< -o $@

$(BUILD)/quantize_quals: $(BUILD)/main.o $(BUILD)/quant.o
	$(CC) $(CFLAGS) $^ -o $@ $(HTSLIB_LIBS) -lm

$(BUILD)/test_quant: tests/test_quant.c $(BUILD)/quant.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

test: $(BUILD)/test_quant $(BUILD)/quantize_quals
	$(BUILD)/test_quant
	QUANTIZE_QUALS=$(BUILD)/quantize_quals tests/test_cli.sh

install: $(BUILD)/quantize_quals
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $< $(DESTDIR)$(PREFIX)/bin/quantize_quals

clean:
	rm -rf $(BUILD)
