SOURCES = report.tex macros.tex
OUT_DIR = out

.PHONY: all clean

all: $(OUT_DIR)/report.pdf

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(OUT_DIR)/report.pdf: report.tex macros.tex $(OUT_DIR)
	pdflatex -output-directory $(OUT_DIR) report.tex 2>&1 > /dev/null

clean:
	rm -rf $(OUT_DIR)
