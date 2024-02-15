% MRJP - Assignment 1
% Marcin Benke
% 2021-09-30

Zadanie 1
=========

A program in the Instant language consists of a sequence of statements
separated by semicolons.

There are two kinds of statements:


* expression - prints its value on stdout,
* assignment of the form `variable = expression` - assigns value of
  the expression to he variable in the LHS; doe snot print anything.

Expressions are built from integer literals, variables and arithmetic
operators. Evaluation order within an expression is not predefined
(you can choose whatever order suits you best)

BNFC syntax:

~~~
Prog. Program ::= [Stmt] ;
SAss. Stmt ::= Ident "=" Exp;
SExp. Stmt ::= Exp ;
separator Stmt ";" ;

ExpAdd.            Exp1   ::= Exp2 "+"  Exp1 ;
ExpSub.            Exp2   ::= Exp2 "-"  Exp3 ;
ExpMul.            Exp3   ::= Exp3 "*"  Exp4 ;
ExpDiv.            Exp3   ::= Exp3 "/"  Exp4 ;
ExpLit.            Exp4   ::= Integer ;
ExpVar.            Exp4   ::= Ident ;
coercions Exp 4;
~~~

**Note:**

* addition binds **to the right**
* addition and multiplication are commutative but not associative

Your task is to write a compiler from Instant to JVM and LLVM.

In this assignment, the generated code should execute all the operations specified in the
input program. Hence it is not allowed to replace the expression `2+3`
by constant 5, omitting assignments to unused variables, etc.
Improving generated code will be a subject of later assignments.

The only allowed, indeed desirable improvement is choosing an
evaluation order so as to minimize the needed JVM stack size. In any
case needed stack size must be computed and declared. (`clever`
solutions like `.limit stack 1000` will not be appreciated). Similarly
you should compute and declare the number of needed locals.

Technical requirements
------------

1. Your solution should be submitted as a packed tar achive (.tar.gz
   or .tgz)
2. The **root** directory of this archive should contain at least
    * A text file README describing how to compile and run the
      program, used libraries and tools, project directory structure,
      possibly references to more detaild documentation.
    * A Makefile to build the project
    * src directory containing only source files of your solution (and
      possibly the Instant.cf supplied by us); auxiliary files should
      be placed in other directories.
3. All used libraries (apart from the standard library of the
programming language used) must be described in README
4. Your submission must compile on students by running `make`
5. After the build, the root of the project must contain executable files
`insc_jvm`  and `insc_llvm`
6. Executing `insc_jvm foo/bar/baz.ins` for a correct program
   `baz.ins` should create files `baz.j` (Jasmin) and `baz.class`
   (JVM) in the directory `foo/bar` (running Jasmin with `-d`
   may be helpful here). Needed runtime library methods should be placed in
   `Runtime.class` in the lib subdirectory.

Executing `insc_llvm foo/bar/baz.ins` for a correct program `baz.ins`
should create files `baz.ll` (text LLVM) and `baz.bc` (lli-exectutable
bitcode)  in the directory `foo/bar`.


Grading
---------

For this assignment you may get at most 6 points. Approximately

* LLVM 2p
* JVM 3p
* For JVM: optimizing order of evaluation for JVM 1p, redundant swap elimination, instruction choice - 1p

Late submissions
----------

Solutions submitted late will be penalised with 1p for every delay
week started. The final date after which no solutions will be accepted
is Dec 5.

Rules
------

The  assignment must be completed on your own. In particular:

* it is not allowed to look at other students' code, or make your code
  available in any way to other students
* all code that is not your own must be clearly marked with a source attribution

Example programs
----------

The archive [instant210930.tgz](instant210930.tgz) contains example
programs with their expected output, as well as `Instant.cf`
containing the BNFC grammar of Instant.

----
&copy; 2021 Marcin Benke
