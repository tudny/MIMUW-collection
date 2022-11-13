# JNP-FuzzyTests

## Uruchomienie
Ścieżka `<path-to-fuzzy.h,fuzzy.cc>` powinna zawierać pliki `fuzzy.h` oraz `fuzzy.cc`
Argument test_name jest opcjonalny. W przypadku jego braku uruchamiane są wszystkie testy.
```bash
./test.sh <path-to-fuzzy.h,fuzzy.cc> [test_name]
```
Skrypt stworzy linki do plików fuzzy.h, fuzzy.cc.
Automatycznie skrypt kompiluje przy użyciu g++-10, chyba że zostanie wykryta nazwa maszyny `students`.
