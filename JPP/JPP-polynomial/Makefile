
test: polyTest polyBenkeTest
	./polyTest
	./polyBenkeTest
	doctest DensePoly.hs
	doctest SparsePoly.hs

polyTest: DensePoly.hs SparsePoly.hs THTestPoly.hs
	ghc -o polyTest THTestPoly.hs -odir 'obj' -hidir 'obj'

polyBenkeTest: original_THTestPoly.hs
	ghc -o polyBenkeTest original_THTestPoly.hs -odir 'obj' -hidir 'obj'

clean:
	rm -f *.o *.hi *.dyn_o *.dyn_hi polyTest

