#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Build cellranger.
#

LIBPY=lib/python
VERSION=$(shell git describe --tags --always --dirty)

export ROOT_DIR=$(shell pwd)

.PHONY: test clean

#
# Targets for development builds.
#
all: tenkit-all bamtofastq websummary

tenkit-all:
	make -C tenkit all

websummary:
	make -C lib/python/websummary

lib/bin:
	mkdir -p lib/bin

bamtofastq: lib/bin
	cd lib/bamtofastq; CARGO_HOME=$(ROOT_DIR)/.cargo cargo build --release
	cp lib/bamtofastq/target/release/bamtofastq lib/bin

clean-bamtofastq:
	rm -Rf .cargo
	rm -Rf lib/bin/bamtofastq
	rm -Rf lib/bamtofastq/target

clean:
	if [[ -d "$(GOPATH)" ]] ; then rm -rf $(GOPATH)/bin ; fi
	if [[ -d "$(GOPATH)" ]] ; then rm -rf $(GOPATH)/pkg ; fi
	rm -rf lib/bin

clean-post-tenkit-submodule:
	rm -rf $(LIBPY)/striped_smith_waterman
	rm -rf $(LIBPY)/tenkit
	rm -rf $(LIBPY)/tenx_numerics
	rm -rf mro/stages/bcl_processor
	rm -rf mro/stages/preflight/bcl_processor

