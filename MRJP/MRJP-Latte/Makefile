
.PHONY: default clean cleanall pack

default: latc_llvm

debug: latc_llvm_debug

latc_llvm: build
	cabal install --install-method=copy --overwrite-policy=always --installdir=.

latc_llvm_debug: build
	# Those lines are for debugging purposes
	rm -f latc_llvm
	# ln -s dist-newstyle/build/x86_64-linux/ghc-8.10.7/MRJP-Latte-0.1.0.0/x/latc_llvm/build/latc_llvm/latc_llvm .
	cp dist-newstyle/build/x86_64-linux/ghc-8.10.7/MRJP-Latte-0.1.0.0/x/latc_llvm/build/latc_llvm/latc_llvm .

build:
	$(MAKE) -C lib -f lib.makefile
	$(MAKE) -C src -f latte.makefile
	cabal build -j12 latc_llvm

clean:
	$(MAKE) -C lib -f lib.makefile clean
	$(MAKE) -C src -f latte.makefile clean
	cabal clean

cleanall: clean
	$(MAKE) -C lib -f lib.makefile cleanall
	$(MAKE) -C src -f latte.makefile cleanall
	rm -f latc_llvm

pack:
	tar -czf at429630-Latte.tar.gz \
		lib/*.makefile \
		lib/*.ll \
		lib/*.c \
		src/*.hs \
		src/Backend/*.hs \
		src/Frontend/*.hs \
		src/*.makefile \
		src/Latte.cf \
		Makefile \
		MRJP-Latte.cabal \
		README.md

save: pack
	mkdir -p backups
	mv at429630-Latte.tar.gz backups/at429630-Latte-`date +%Y-%m-%d_%H-%M-%S`.tar.gz
