# Instant

## Usage
Run 
```bash
make
``` 
to build jvm and llvm compilers. Both executables depend on `lib` contents.


Each compiler can be build on their own. Respectively for `jvm` and `llvm`:
```bash
make jvm
make llvm
```

## Dependencies
`Jvm` compiler depends on `Runtime.class` (`Runtime.java`).
- Source: own implementation.

`Llvm` compiler depends on `runtime.bc` (`runtime.ll`).
- Source: MRJP's Moodle.
