


COMPONENTS=$(shell ls src/components/*.js)
dist=dist/tenx-websummary-script.min.js dist/tenx-websummary-styles.min.css

.PHONY: all
all:$(dist)

$(dist): $(COMPONENTS)
	cd src; npm install
	cd src; npm run build
