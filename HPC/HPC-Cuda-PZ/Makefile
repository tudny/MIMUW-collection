
BUILD_DIR = build
BUILD_NAME = gpugenv
INSTALL_DIR = .
ZIP_FILE = "$(BUILD_NAME).Aleksander.Tudruj.zip"
SOURCES = src CMakeLists.txt Makefile README.md report.pdf
REPORT_DIR = report
DATA_DIR = data

.PHONY: all run install clean pack report build

all: install

build:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=Release .. && make

run: build
	./$(BUILD_DIR)/$(BUILD_NAME)

install: build
	cp $(BUILD_DIR)/$(BUILD_NAME) $(INSTALL_DIR)

report:
	$(MAKE) -C $(REPORT_DIR) all

pack: report
	cp $(REPORT_DIR)/out/report.pdf .
	zip -r $(ZIP_FILE) $(SOURCES)

clean:
	$(MAKE) -C $(REPORT_DIR) clean
	rm -rf $(BUILD_DIR)
	rm -f $(INSTALL_DIR)/$(BUILD_NAME)
	rm -f *.zip
	rm -f report.pdf

pull-data:
	mkdir -p $(DATA_DIR) ; \
	cd $(DATA_DIR) ; \
	scp entropy:/home/krzadca/biomatch-data/dbNSFP4.5a_variant_100c_sorted.chr22 . ; \
	scp entropy:/home/krzadca/biomatch-data/homo_sapiens_GRCh38.vcf .
