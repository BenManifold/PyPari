PyPari
======

A Python Interface for the Pari/Gp Computer Algebra System

At present effective use of this tool requires some knowledge of both Ctypes as well as the use of PariLib in library mode.
I'm working on implementing gp functionality by linking up GP expressions with their appropriate pari-library function calls.
For now, assuming you have pari installed and have libpari.so in /usr/lib, simply import PyPari.py in the python shell.
Numeric Types, Strings, and Lists should all convert naturally to their pari analogs cia Ctypes. String arguements will be evaluated by the by the Pari syntax analyzer, and in this way we can naturally pass things like polynomials and matrices to parilib functions. Creating Pari Gen objects is handled by invocation of the PyPari.Gen class constructor, and passing string arguements to this constructor is the preferred way to directly create lists, matrices, polynomials, and so on. 
