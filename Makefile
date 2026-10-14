.PHONY: all clean build help
PREFIX ?= $(HOME)/local

all: build

clean:
	$(RM) -rf build

build: build/swatplus-modflow

build/swatplus-modflow:
	scripts/build.sh

install: build/swatplus-modflow
	install -d $(PREFIX)/bin
	install -m 755 $< $(PREFIX)/bin

help:
	@echo "Available targets:"
	@echo "  all      Build SWAT+MODFLOW (default)"
	@echo "  build    Build SWAT+MODFLOW"
	@echo "  clean    Remove the build directory"
	@echo "  install  Install to PREFIX/bin (default: ~/local/bin)"
	@echo "  help     Show this help message"
