CC       ?= cc

DOTGENO_PREFIX ?= /usr/local

PREFIX ?= /usr/local
BINDIR  = $(PREFIX)/bin

TARGET = admix_io
SRCDIR = src
OUTDIR = bin

CPPFLAGS += -I$(SRCDIR) -I$(DOTGENO_PREFIX)/include -D_POSIX_C_SOURCE=200809L
CFLAGS   ?= -std=c17 -Wall -Wextra -Wpedantic -O2
LDFLAGS  += -L$(DOTGENO_PREFIX)/lib -Wl,-rpath,$(DOTGENO_PREFIX)/lib
LDLIBS   += -ldotgeno -lm

SRC  = $(SRCDIR)/admix_io.c
HDRS = $(wildcard $(SRCDIR)/*.h)
BIN  = $(OUTDIR)/$(TARGET)

.PHONY: all clean install uninstall

all: $(BIN)

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(BIN): $(SRC) $(HDRS) | $(OUTDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

install: all
	mkdir -p $(DESTDIR)$(BINDIR)
	install -m 755 $(BIN) $(DESTDIR)$(BINDIR)/

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(TARGET)

clean:
	rm -rf $(OUTDIR)
