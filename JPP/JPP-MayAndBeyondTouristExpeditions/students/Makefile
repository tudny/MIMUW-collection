SICSTUSHOME=/opt/sicstus
SICSTUSBIN= $(SICSTUSHOME)/bin
PL = $(SICSTUSBIN)/sicstus
SPLD = spld
SPLDFLAGS = --static --exechome=$(SICSTUSBIN)

BINFILES = at429630
BINFILE = wyprawy

all: $(BINFILES)
	cp $(BINFILES) $(BINFILE)

%: %.sav
	$(SPLD) $(SPLDFLAGS) $< -o $@ 

%.sav: %.pl
	echo "compile('$<'). save_program('$@')." | $(PL)

clean:
	rm -f $(ALL) *.sav