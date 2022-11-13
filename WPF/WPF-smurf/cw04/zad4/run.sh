
# Sygnatury modułów (TUTAJ WPISUJESZ SWOJE SYGNATURY)
SIGNATURES=''

# Moduły (w tym moduł główny) (TUTAJ WPISUJESZ SWOJE MODULY)
MODULES='head_tail.ml'

# Nazwa programu wykonywalnego (TUTAJ WPISUJESZ NAZWE PILKU PO KOMPILACJI)
RUNNABLE='HT'

echo "====================================="
echo "Kompilowanie ocamlopt"
echo "SYGNATURY: \"$SIGNATURES\""
echo "MODUŁY: \"$MODULES\""
echo "NAZWA PILKU WYKONYWALNEGO: \"$RUNNABLE\""
echo "====================================="

# Funkcja kompilująca (TEGO RACZEJ NIE DOTYKAJ)
ocamlopt $SIGNATURES $MODULES -o $RUNNABLE


# Sprawdzenie czy kompilacja się powiodła
if test -f "$RUNNABLE"; then
    echo "Uruchamianie $RUNNABLE"
    
    # Uruchamianie programu
	./$RUNNABLE
	
else
	echo "Compilation failed!"
fi



