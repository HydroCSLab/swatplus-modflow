.PHONY: all clean build get-build

all: build

clean:
	$(RM) -rf build

build:
	scripts/build.sh

get-build:
	scripts/get-build.sh
