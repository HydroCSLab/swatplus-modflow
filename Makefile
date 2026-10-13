.DEFAULT_GOAL := help

BUILD_DIR ?= build
MESON ?= meson

.PHONY: help build clean

help:
	@printf 'Targets:\n  build  Compile SWAT+MODFLOW incrementally\n  clean  Remove compiled files\n'

build:
	@if [ ! -f '$(BUILD_DIR)/build.ninja' ]; then $(MESON) setup '$(BUILD_DIR)'; fi
	$(MESON) compile -C '$(BUILD_DIR)'

clean:
	rm -rf -- '$(BUILD_DIR)'
