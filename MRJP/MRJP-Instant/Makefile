
all: jvm llvm

jvm: lib/Runtime.class
	$(MAKE) -C src Target/insc_jvm
	cp -u src/Target/insc_jvm insc_jvm  

llvm: lib/runtime.bc
	$(MAKE) -C src Target/insc_llvm
	cp -u src/Target/insc_llvm insc_llvm

lib/Runtime.class: lib/Runtime.java
	javac lib/Runtime.java

lib/runtime.bc: lib/runtime.ll
	llvm-as lib/runtime.ll -o lib/runtime.bc

# Rules for cleaning generated files (exclusing generated binaries).
clean:
	$(MAKE) -C src clean

# Rules for cleaning generated files (including generated binaries).
distclean:
	$(MAKE) -C src distclean
	rm -f examples/*.j examples/*.class examples/*.ll examples/*.bc
	rm -f insc_jvm insc_llvm lib/Runtime.class lib/runtime.bc

pack: distclean
	tar -czvf at429630-instant.tar.gz lib src Makefile README.md examples

pack-solution: distclean
	tar -czvf at429630-instant-solution.tar.gz lib src Makefile README.md
