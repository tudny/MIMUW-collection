BNFC       = /home/students/inf/PUBLIC/MRJP/bin/bnfc
ifneq ($(shell which bnfc),)
	BNFC       = $(shell which bnfc)
endif

HAPPY      = happy
HAPPY_OPTS = --array --info --ghc --coerce
ALEX       = alex
ALEX_OPTS  = --ghc

.PHONY : all clean cleanall

all : Latte/Abs.hs Latte/Lex.hs Latte/Par.hs Latte/Print.hs

# Rules for building the parser.

Latte/Abs.hs Latte/Lex.x Latte/Par.y Latte/Print.hs Latte/Test.hs : Latte.cf
	$(BNFC) --haskell -d --functor Latte.cf

%.hs : %.y
	${HAPPY} ${HAPPY_OPTS} $<

%.hs : %.x
	${ALEX} ${ALEX_OPTS} $<

clean:
	-rm -rf Latte/

cleanall: clean
