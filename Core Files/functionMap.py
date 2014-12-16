import PyPari
import ctypes
PyPari.pari_init(400000,400000)

def bernvec(*argv):
  '''
  bernvec
  Class: basic
  Section: transcendental
  C-Name: bernvec
  Prototype: L
  Help: bernvec(x): Vector of rational Bernoulli numbers B_0, B_2,...up to
   B_(2x).
  Doc: creates a vector containing, as rational numbers,
   the \idx{Bernoulli numbers} $B_0$, $B_2$,\dots, $B_{2x}$.
   This routine is obsolete. Use \kbd{bernfrac} instead each time you need a
   Bernoulli number in exact form.
   
   \misctitle{Note} This routine is implemented using repeated independent
   calls to \kbd{bernfrac}, which is faster than the standard recursion in exact
   arithmetic. It is only kept for backward compatibility: it is not faster than
   individual calls to \kbd{bernfrac}, its output uses a lot of memory space,
   and coping with the index shift is awkward.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bernvec(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def listcreate(*argv):
  '''
  listcreate
  Class: basic
  Section: linear_algebra
  C-Name: listcreate
  Prototype: D0,L,
  Help: listcreate(): creates an empty list.
  Description: 
   (?gen):list        listcreate()
  Doc: creates an empty list. This routine used to have a mandatory argument,
   which is now ignored (for backward compatibility). In fact, this function
   has become redundant and obsolete; it will disappear in future versions of
   PARI: just use \kbd{List()}
   % \syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.listcreate(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def listkill(*argv):
  '''
  listkill
  Class: basic
  Section: linear_algebra
  C-Name: listkill
  Prototype: vG
  Help: listkill(L): obsolete, retained for backward compatibility.
  Doc: obsolete, retained for backward compatibility. Just use \kbd{L = List()}
   instead of \kbd{listkill(L)}. In most cases, you won't even need that, e.g.
   local variables are automatically cleared when a user function returns.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.listkill(*c_arg_tuple)

def factorcantor(*argv):
  '''
  factorcantor
  Class: basic
  Section: number_theoretical
  C-Name: factcantor
  Prototype: GG
  Help: factorcantor(x,p): factorization mod p of the polynomial x using
   Cantor-Zassenhaus.
  Doc: factors the polynomial $x$ modulo the
   prime $p$, using distinct degree plus
   \idx{Cantor-Zassenhaus}\sidx{Zassenhaus}. The coefficients of $x$ must be
   operation-compatible with $\Z/p\Z$. The result is a two-column matrix, the
   first column being the irreducible polynomials dividing $x$, and the second
   the exponents. If you want only the \emph{degrees} of the irreducible
   polynomials (for example for computing an $L$-function), use
   $\kbd{factormod}(x,p,1)$. Note that the \kbd{factormod} algorithm is
   usually faster than \kbd{factorcantor}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factcantor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factorff(*argv):
  '''
  factorff
  Class: basic
  Section: number_theoretical
  C-Name: factorff
  Prototype: GDGDG
  Help: factorff(x,{p},{a}): factorization of the polynomial x in the finite field
   F_p[X]/a(X)F_p[X].
  Doc: factors the polynomial $x$ in the field
   $\F_q$ defined by the irreducible polynomial $a$ over $\F_p$. The
   coefficients of $x$ must be operation-compatible with $\Z/p\Z$. The result
   is a two-column matrix: the first column contains the irreducible factors of
   $x$, and the second their exponents. If all the coefficients of $x$ are in
   $\F_p$, a much faster algorithm is applied, using the computation of
   isomorphisms between finite fields.
   
   Either $a$ or $p$ can omitted (in which case both are ignored) if x has
   \typ{FFELT} coefficients; the function then becomes identical to \kbd{factor}:
   \bprog
   ? factorff(x^2 + 1, 5, y^2+3)  \\ over F_5[y]/(y^2+3) ~ F_25
   %1 =
   [Mod(Mod(1, 5), Mod(1, 5)*y^2 + Mod(3, 5))*x
    + Mod(Mod(2, 5), Mod(1, 5)*y^2 + Mod(3, 5)) 1]
   
   [Mod(Mod(1, 5), Mod(1, 5)*y^2 + Mod(3, 5))*x
    + Mod(Mod(3, 5), Mod(1, 5)*y^2 + Mod(3, 5)) 1]
   ? t = ffgen(y^2 + Mod(3,5), 't); \\ a generator for F_25 as a t_FFELT
   ? factorff(x^2 + 1)   \\ not enough information to determine the base field
    ***   at top-level: factorff(x^2+1)
    ***                 ^---------------
    *** factorff: incorrect type in factorff.
   ? factorff(x^2 + t^0) \\ make sure a coeff. is a t_FFELT
   %3 =
   [x + 2 1]
   
   [x + 3 1]
   ? factorff(x^2 + t + 1)
   %11 =
   [x + (2*t + 1) 1]
   
   [x + (3*t + 4) 1]
   @eprog\noindent
   Notice that the second syntax is easier to use and much more readable.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factorff(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factormod(*argv):
  '''
  factormod
  Class: basic
  Section: number_theoretical
  C-Name: factormod0
  Prototype: GGD0,L,
  Help: factormod(x,p,{flag=0}): factors the polynomial x modulo the prime p, using Berlekamp. flag is optional, and can be 0: default or 1:
   only the degrees of the irreducible factors are given.
  Doc: factors the polynomial $x$ modulo the prime integer $p$, using
   \idx{Berlekamp}. The coefficients of $x$ must be operation-compatible with
   $\Z/p\Z$. The result is a two-column matrix, the first column being the
   irreducible polynomials dividing $x$, and the second the exponents. If $\fl$
   is non-zero, outputs only the \emph{degrees} of the irreducible polynomials
   (for example, for computing an $L$-function). A different algorithm for
   computing the mod $p$ factorization is \kbd{factorcantor} which is sometimes
   faster.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factormod0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polrootsff(*argv):
  '''
  polrootsff
  Class: basic
  Section: number_theoretical
  C-Name: polrootsff
  Prototype: GDGDG
  Help: polrootsff(x,{p},{a}): returns the roots of the polynomial x in the finite
   field F_p[X]/a(X)F_p[X]. a or p can be omitted if x has t_FFELT coefficients.
  Doc: returns the vector of distinct roots of the polynomial $x$ in the field
   $\F_q$ defined by the irreducible polynomial $a$ over $\F_p$. The
   coefficients of $x$ must be operation-compatible with $\Z/p\Z$.
   Either $a$ or $p$ can omitted (in which case both are ignored) if x has
   \typ{FFELT} coefficients:
   \bprog
   ? polrootsff(x^2 + 1, 5, y^2+3)  \\ over F_5[y]/(y^2+3) ~ F_25
   %1 = [Mod(Mod(3, 5), Mod(1, 5)*y^2 + Mod(3, 5)),
         Mod(Mod(2, 5), Mod(1, 5)*y^2 + Mod(3, 5))]
   ? t = ffgen(y^2 + Mod(3,5), 't); \\ a generator for F_25 as a t_FFELT
   ? polrootsff(x^2 + 1)   \\ not enough information to determine the base field
    ***   at top-level: polrootsff(x^2+1)
    ***                 ^-----------------
    *** polrootsff: incorrect type in factorff.
   ? polrootsff(x^2 + t^0) \\ make sure one coeff. is a t_FFELT
   %3 = [3, 2]
   ? polrootsff(x^2 + t + 1)
   %4 = [2*t + 1, 3*t + 4]
   @eprog\noindent
   Notice that the second syntax is easier to use and much more readable.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polrootsff(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polrootsmod(*argv):
  '''
  polrootsmod
  Class: basic
  Section: polynomials
  C-Name: rootmod0
  Prototype: GGD0,L,
  Help: polrootsmod(pol,p,{flag=0}): roots mod the prime p of the polynomial pol. flag is
   optional, and can be 0: default, or 1: use a naive search, useful for small p.
  Description: 
   (pol, int, ?0):vec           rootmod($1, $2)
   (pol, int, 1):vec            rootmod2($1, $2)
   (pol, int, #small):vec       $"Bad flag in polrootsmod"
   (pol, int, small):vec        rootmod0($1, $2, $3)
  Doc: row vector of roots modulo $p$ of the polynomial \var{pol}.
   Multiple roots are \emph{not} repeated.
   \bprog
   ? polrootsmod(x^2-1,2)
   %1 = [Mod(1, 2)]~
   @eprog\noindent
   If $p$ is very small, you may set $\fl=1$, which uses a naive search.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rootmod0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polhensellift(*argv):
  '''
  polhensellift
  Class: basic
  Section: polynomials
  C-Name: polhensellift
  Prototype: GGGL
  Help: polhensellift(A, B, p, e): lift the factorization B of A modulo p to a
   factorization modulo p^e using Hensel lift. The factors in B must be
   pairwise relatively prime modulo p.
  Doc: given a prime $p$, an integral polynomial $A$ whose leading coefficient
   is a $p$-unit, a vector $B$ of integral polynomials that are monic and
   pairwise relatively prime modulo $p$, and whose product is congruent to
   $A/\text{lc}(A)$ modulo $p$, lift the elements of $B$ to polynomials whose
   product is congruent to $A$ modulo $p^e$.
   
   More generally, if $T$ is an integral polynomial irreducible mod $p$, and
   $B$ is a factorization of $A$ over the finite field $\F_p[t]/(T)$, you can
   lift it to $\Z_p[t]/(T, p^e)$ by replacing the $p$ argument with $[p,T]$:
   \bprog
   ? { T = t^3 - 2; p = 7; A = x^2 + t + 1;
       B = [x + (3*t^2 + t + 1), x + (4*t^2 + 6*t + 6)];
       r = polhensellift(A, B, [p, T], 6) }
   %1 = [x + (20191*t^2 + 50604*t + 75783), x + (97458*t^2 + 67045*t + 41866)]
   ? liftall( r[1] * r[2] * Mod(Mod(1,p^6),T) )
   %2 = x^2 + (t + 1)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polhensellift(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polcyclofactors(*argv):
  '''
  polcyclofactors
  Class: basic
  Section: polynomials
  C-Name: polcyclofactors
  Prototype: G
  Help: polcyclofactors(f): returns a vector of polynomials, whose product is
   the product of distinct cyclotomic polynomials dividing f.
  Doc: returns a vector of polynomials, whose product is the product of
   distinct cyclotomic polynomials dividing $f$.
   \bprog
   ? f = x^10+5*x^8-x^7+8*x^6-4*x^5+8*x^4-3*x^3+7*x^2+3;
   ? v = polcyclofactors(f)
   %2 = [x^2 + 1, x^2 + x + 1, x^4 - x^3 + x^2 - x + 1]
   ? apply(poliscycloprod, v)
   %3 = [1, 1, 1]
   ? apply(poliscyclo, v)
   %4 = [4, 3, 10]
   @eprog\noindent In general, the poynomials are products of cyclotomic
   polynomials and not themselves irreducible:
   \bprog
   ? g = x^8+2*x^7+6*x^6+9*x^5+12*x^4+11*x^3+10*x^2+6*x+3;
   ? polcyclofactors(g)
   %2 = [x^6 + 2*x^5 + 3*x^4 + 3*x^3 + 3*x^2 + 2*x + 1]
   ? factor(%[1])
   %3 =
   [            x^2 + x + 1 1]
   
   [x^4 + x^3 + x^2 + x + 1 1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polcyclofactors(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poliscyclo(*argv):
  '''
  poliscyclo
  Class: basic
  Section: polynomials
  C-Name: poliscyclo
  Prototype: lG
  Help: poliscyclo(f): returns 0 if f is not a cyclotomic polynomial, and n
   > 0 if f = Phi_n, the n-th cyclotomic polynomial.
  Doc: returns 0 if $f$ is not a cyclotomic polynomial, and $n > 0$ if $f =
   \Phi_n$, the $n$-th cyclotomic polynomial.
   \bprog
   ? poliscyclo(x^4-x^2+1)
   %1 = 12
   ? polcyclo(12)
   %2 = x^4 - x^2 + 1
   ? poliscyclo(x^4-x^2-1)
   %3 = 0
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.poliscyclo(*c_arg_tuple)

def poliscycloprod(*argv):
  '''
  poliscycloprod
  Class: basic
  Section: polynomials
  C-Name: poliscycloprod
  Prototype: lG
  Help: poliscycloprod(f): returns 1 if f is a product of cyclotomic
   polynonials, and 0 otherwise.
  Doc: returns 1 if $f$ is a product of cyclotomic polynomial, and $0$
   otherwise.
   \bprog
   ? f = x^6+x^5-x^3+x+1;
   ? poliscycloprod(f)
   %2 = 1
   ? factor(f)
   %3 =
   [  x^2 + x + 1 1]
   
   [x^4 - x^2 + 1 1]
   ? [ poliscyclo(T) | T <- %[,1] ]
   %4 = [3, 12]
   ? polcyclo(3) * polcyclo(12)
   %5 = x^6 + x^5 - x^3 + x + 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.poliscycloprod(*c_arg_tuple)

def matisdiagonal(*argv):
  '''
  matisdiagonal
  Class: basic
  Section: linear_algebra
  C-Name: isdiagonal
  Prototype: iG
  Help: matisdiagonal(x): true(1) if x is a diagonal matrix, false(0)
   otherwise.
  Doc: returns true (1) if $x$ is a diagonal matrix, false (0) if not.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.isdiagonal(*c_arg_tuple)

def matid(*argv):
  '''
  matid
  Class: basic
  Section: linear_algebra
  C-Name: matid
  Prototype: L
  Help: matid(n): identity matrix of order n.
  Description: 
   (small):vec    matid($1)
  Doc: creates the $n\times n$ identity matrix.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matid(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matdet(*argv):
  '''
  matdet
  Class: basic
  Section: linear_algebra
  C-Name: det0
  Prototype: GD0,L,
  Help: matdet(x,{flag=0}): determinant of the matrix x using an appropriate
   algorithm depending on the coefficients. If (optional) flag is set to 1, use
   classical Gaussian elimination (usually worse than the default).
  Description: 
   (gen, ?0):gen           det($1)
   (gen, 1):gen            det2($1)
   (gen, #small):gen       $"incorrect flag in matdet"
   (gen, small):gen        det0($1, $2)
  Doc: determinant of the square matrix $x$.
   
   If $\fl=0$, uses an appropriate algorithm depending on the coefficients:
   
   \item integer entries: modular method due to Dixon, Pernet and Stein.
   
   \item real or $p$-adic entries: classical Gaussian elimination using maximal
   pivot.
   
   \item intmod entries: classical Gaussian elimination using first non-zero
   pivot.
   
   \item other cases: Gauss-Bareiss.
   
   If $\fl=1$, uses classical Gaussian elimination with appropriate pivoting
   strategy (maximal pivot for real or $p$-adic coefficients). This is usually
   worse than the default.
  Variant: Also available are \fun{GEN}{det}{GEN x} ($\fl=0$),
   \fun{GEN}{det2}{GEN x} ($\fl=1$) and \fun{GEN}{ZM_det}{GEN x} for integer
   entries.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.det0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matdetint(*argv):
  '''
  matdetint
  Class: basic
  Section: linear_algebra
  C-Name: detint
  Prototype: G
  Help: matdetint(B): some multiple of the determinant of the lattice
   generated by the columns of B (0 if not of maximal rank). Useful with
   mathnfmod.
  Doc: 
   Let $B$ be an $m\times n$ matrix with integer coefficients. The
   \emph{determinant} $D$ of the lattice generated by the columns of $B$ is
   the square root of $\det(B^T B)$ if $B$ has maximal rank $m$, and $0$
   otherwise.
   
   This function uses the Gauss-Bareiss algorithm to compute a positive
   \emph{multiple} of $D$. When $B$ is square, the function actually returns
   $D = |\det B|$.
   
   This function is useful in conjunction with \kbd{mathnfmod}, which needs to
   know such a multiple. If the rank is maximal and the matrix non-square,
   you can obtain $D$ exactly using
   \bprog
     matdet( mathnfmod(B, matdetint(B)) )
   @eprog\noindent
   Note that as soon as one of the dimensions gets large ($m$ or $n$ is larger
   than 20, say), it will often be much faster to use \kbd{mathnf(B, 1)} or
   \kbd{mathnf(B, 4)} directly.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.detint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matsolve(*argv):
  '''
  matsolve
  Class: basic
  Section: linear_algebra
  C-Name: gauss
  Prototype: GG
  Help: matsolve(M,B): solution of MX=B (M matrix, B column vector).
  Doc: $M$ being an invertible matrix and $B$ a column
   vector, finds the solution $X$ of $MX=B$, using Dixon $p$-adic lifting method
   if $M$ and $B$ are integral and Gaussian elimination otherwise. This
   has the same effect as, but is faster, than $M^{-1}*B$.
  Variant: For integral input, the function
   \fun{GEN}{ZM_gauss}{GEN M,GEN B} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gauss(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matimagecompl(*argv):
  '''
  matimagecompl
  Class: basic
  Section: linear_algebra
  C-Name: imagecompl
  Prototype: G
  Help: matimagecompl(x): vector of column indices not corresponding to the
   indices given by the function matimage.
  Description: 
   (gen):vecsmall                imagecompl($1)
  Doc: gives the vector of the column indices which
   are not extracted by the function \kbd{matimage}, as a permutation
   (\typ{VECSMALL}). Hence the number of
   components of \kbd{matimagecompl(x)} plus the number of columns of
   \kbd{matimage(x)} is equal to the number of columns of the matrix $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.imagecompl(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matindexrank(*argv):
  '''
  matindexrank
  Class: basic
  Section: linear_algebra
  C-Name: indexrank
  Prototype: G
  Help: matindexrank(x): gives two extraction vectors (rows and columns) for
   the matrix x such that the extracted matrix is square of maximal rank.
  Doc: $x$ being a matrix of rank $r$, returns a vector with two
   \typ{VECSMALL} components $y$ and $z$ of length $r$ giving a list of rows
   and columns respectively (starting from 1) such that the extracted matrix
   obtained from these two vectors using $\tet{vecextract}(x,y,z)$ is
   invertible.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.indexrank(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matinverseimage(*argv):
  '''
  matinverseimage
  Class: basic
  Section: linear_algebra
  C-Name: inverseimage
  Prototype: GG
  Help: matinverseimage(x,y): an element of the inverse image of the vector y
   by the matrix x if one exists, the empty vector otherwise.
  Doc: given a matrix $x$ and
   a column vector or matrix $y$, returns a preimage $z$ of $y$ by $x$ if one
   exists (i.e such that $x z = y$), an empty vector or matrix otherwise. The
   complete inverse image is $z + \text{Ker} x$, where a basis of the kernel of
   $x$ may be obtained by \kbd{matker}.
   \bprog
   ? M = [1,2;2,4];
   ? matinverseimage(M, [1,2]~)
   %2 = [1, 0]~
   ? matinverseimage(M, [3,4]~)
   %3 = []~    \\@com no solution
   ? matinverseimage(M, [1,3,6;2,6,12])
   %4 =
   [1 3 6]
   
   [0 0 0]
   ? matinverseimage(M, [1,2;3,4])
   %5 = [;]    \\@com no solution
   ? K = matker(M)
   %6 =
   [-2]
   
   [1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.inverseimage(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mateigen(*argv):
  '''
  mateigen
  Class: basic
  Section: linear_algebra
  C-Name: mateigen
  Prototype: GD0,L,p
  Help: mateigen(x,{flag=0}): complex eigenvectors of the matrix x given as
   columns of a matrix H. If flag=1, return [L,H], where L contains the
   eigenvalues and H the corresponding eigenvectors.
  Doc: returns the (complex) eigenvectors of $x$ as columns of a matrix.
   If $\fl=1$, return $[L,H]$, where $L$ contains the
   eigenvalues and $H$ the corresponding eigenvectors; multiple eigenvalues are
   repeated according to the eigenspace dimension (which may be less
   than the eigenvalue multiplicity in the characteristic polynomial).
   
   This function first computes the characteristic polynomial of $x$ and
   approximates its complex roots $(\lambda_i)$, then tries to compute the
   eigenspaces as kernels of the $x - \lambda_i$. This algorithm is
   ill-conditioned and is likely to miss kernel vectors if some roots of the
   characteristic polynomial are close, in particular if it has multiple roots.
   \bprog
   ? A = [13,2; 10,14]; mateigen(A)
   %1 =
   [-1/2 2/5]
   
   [   1   1]
   ? [L,H] = mateigen(A, 1);
   ? L
   %3 = [9, 18]
   ? H
   %4 =
   [-1/2 2/5]
   
   [   1   1]
   @eprog\noindent
   For symmetric matrices, use \tet{qfjacobi} instead; for Hermitian matrices,
   compute
   \bprog
    A = real(x);
    B = imag(x);
    y = matconcat([A, -B; B, A]);
   @eprog\noindent and apply \kbd{qfjacobi} to $y$.
  Variant: Also available is \fun{GEN}{eigen}{GEN x, long prec} ($\fl = 0$)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mateigen(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matimage(*argv):
  '''
  matimage
  Class: basic
  Section: linear_algebra
  C-Name: matimage0
  Prototype: GD0,L,
  Help: matimage(x,{flag=0}): basis of the image of the matrix x. flag is
   optional and can be set to 0 or 1, corresponding to two different algorithms.
  Description: 
   (gen, ?0):vec           image($1)
   (gen, 1):vec            image2($1)
   (gen, #small)           $"incorrect flag in matimage"
   (gen, small):vec        matimage0($1, $2)
  Doc: gives a basis for the image of the
   matrix $x$ as columns of a matrix. A priori the matrix can have entries of
   any type. If $\fl=0$, use standard Gauss pivot. If $\fl=1$, use
   \kbd{matsupplement} (much slower: keep the default flag!).
  Variant: Also available is \fun{GEN}{image}{GEN x} ($\fl=0$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matimage0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matker(*argv):
  '''
  matker
  Class: basic
  Section: linear_algebra
  C-Name: matker0
  Prototype: GD0,L,
  Help: matker(x,{flag=0}): basis of the kernel of the matrix x. flag is
   optional, and may be set to 0: default; non-zero: x is known to have
   integral entries.
  Description: 
   (gen, ?0):vec           ker($1)
   (gen, 1):vec            keri($1)
   (gen, #small)           $"incorrect flag in matker"
   (gen, small):vec        matker0($1, $2)
  Doc: gives a basis for the kernel of the matrix $x$ as columns of a matrix.
   The matrix can have entries of any type, provided they are compatible with
   the generic arithmetic operations ($+$, $\times$ and $/$).
   
   If $x$ is known to have integral entries, set $\fl=1$.
  Variant: Also available are \fun{GEN}{ker}{GEN x} ($\fl=0$),
   \fun{GEN}{keri}{GEN x} ($\fl=1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matker0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matsolvemod(*argv):
  '''
  matsolvemod
  Class: basic
  Section: linear_algebra
  C-Name: matsolvemod0
  Prototype: GGGD0,L,
  Help: matsolvemod(M,D,B,{flag=0}): one solution of system of congruences
   MX=B mod D (M matrix, B and D column vectors). If (optional) flag is
   non-null return all solutions.
  Doc: $M$ being any integral matrix,
   $D$ a column vector of non-negative integer moduli, and $B$ an integral
   column vector, gives a small integer solution to the system of congruences
   $\sum_i m_{i,j}x_j\equiv b_i\pmod{d_i}$ if one exists, otherwise returns
   zero. Shorthand notation: $B$ (resp.~$D$) can be given as a single integer,
   in which case all the $b_i$ (resp.~$d_i$) above are taken to be equal to $B$
   (resp.~$D$).
   \bprog
   ? M = [1,2;3,4];
   ? matsolvemod(M, [3,4]~, [1,2]~)
   %2 = [-2, 0]~
   ? matsolvemod(M, 3, 1) \\ M X = [1,1]~ over F_3
   %3 = [-1, 1]~
   ? matsolvemod(M, [3,0]~, [1,2]~) \\ x + 2y = 1 (mod 3), 3x + 4y = 2 (in Z)
   %4 = [6, -4]~
   @eprog
   If $\fl=1$, all solutions are returned in the form of a two-component row
   vector $[x,u]$, where $x$ is a small integer solution to the system of
   congruences and $u$ is a matrix whose columns give a basis of the homogeneous
   system (so that all solutions can be obtained by adding $x$ to any linear
   combination of columns of $u$). If no solution exists, returns zero.
  Variant: Also available are \fun{GEN}{gaussmodulo}{GEN M, GEN D, GEN B}
   ($\fl=0$) and \fun{GEN}{gaussmodulo2}{GEN M, GEN D, GEN B} ($\fl=1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matsolvemod0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matrank(*argv):
  '''
  matrank
  Class: basic
  Section: linear_algebra
  C-Name: rank
  Prototype: lG
  Help: matrank(x): rank of the matrix x.
  Doc: rank of the matrix $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.rank(*c_arg_tuple)

def matsupplement(*argv):
  '''
  matsupplement
  Class: basic
  Section: linear_algebra
  C-Name: suppl
  Prototype: G
  Help: matsupplement(x): supplement the columns of the matrix x to an
   invertible matrix.
  Doc: assuming that the columns of the matrix $x$
   are linearly independent (if they are not, an error message is issued), finds
   a square invertible matrix whose first columns are the columns of $x$,
   i.e.~supplement the columns of $x$ to a basis of the whole space.
   \bprog
   ? matsupplement([1;2])
   %1 =
   [1 0]
   
   [2 1]
   @eprog
   Raises an error if $x$ has 0 columns, since (due to a long standing design
   bug), the dimension of the ambient space (the number of rows) is unknown in
   this case:
   \bprog
   ? matsupplement(matrix(2,0))
     ***   at top-level: matsupplement(matrix
     ***                 ^--------------------
     *** matsupplement: sorry, suppl [empty matrix] is not yet implemented.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.suppl(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def charpoly(*argv):
  '''
  charpoly
  Class: basic
  Section: linear_algebra
  C-Name: charpoly0
  Prototype: GDnD5,L,
  Help: charpoly(A,{v='x},{flag=5}): det(v*Id-A)=characteristic polynomial of
   the matrix or polmod A. flag is optional and ignored unless A is a matrix;
   it may be set to 0 (Le Verrier), 1 (Lagrange interpolation),
   2 (Hessenberg form), 3 (Berkowitz), 4 (modular) if A is integral,
   or 5 (default, choose best method).
   Algorithms 0 (Le Verrier) and 1 (Lagrange) assume that n! is invertible,
   where n is the dimension of the matrix.
  Doc: 
   \idx{characteristic polynomial}
   of $A$ with respect to the variable $v$, i.e.~determinant of $v*I-A$ if $A$
   is a square matrix.
   \bprog
   ? charpoly([1,2;3,4]);
   %1 = x^2 - 5*x - 2
   ? charpoly([1,2;3,4],, 't)
   %2 = t^2 - 5*t - 2
   @eprog\noindent
   If $A$ is not a square matrix, the function returns the characteristic
   polynomial of the map ``multiplication by $A$'' if $A$ is a scalar:
   \bprog
   ? charpoly(Mod(x+2, x^3-2))
   %1 = x^3 - 6*x^2 + 12*x - 10
   ? charpoly(I)
   %2 = x^2 + 1
   ? charpoly(quadgen(5))
   %3 = x^2 - x - 1
   ? charpoly(ffgen(ffinit(2,4)))
   %4 = Mod(1, 2)*x^4 + Mod(1, 2)*x^3 + Mod(1, 2)*x^2 + Mod(1, 2)*x + Mod(1, 2)
   @eprog
   
   The value of $\fl$ is only significant for matrices, and we advise to stick
   to the default value. Let $n$ be the dimension of $A$.
   
   If $\fl=0$, same method (Le Verrier's) as for computing the adjoint matrix,
   i.e.~using the traces of the powers of $A$. Assumes that $n!$ is
   invertible; uses $O(n^4)$ scalar operations.
   
   If $\fl=1$, uses Lagrange interpolation which is usually the slowest method.
   Assumes that $n!$ is invertible; uses $O(n^4)$ scalar operations.
   
   If $\fl=2$, uses the Hessenberg form. Assumes that the base ring is a field.
   Uses $O(n^3)$ scalar operations, but suffers from coefficient explosion
   unless the base field is finite or $\R$.
   
   If $\fl=3$, uses Berkowitz's division free algorithm, valid over any
   ring (commutative, with unit). Uses $O(n^4)$ scalar operations.
   
   If $\fl=4$, $x$ must be integral. Uses a modular algorithm: Hessenberg form
   for various small primes, then Chinese remainders.
   
   If $\fl=5$ (default), uses the ``best'' method given $x$.
   This means we use Berkowitz unless the base ring is $\Z$ (use $\fl=4$)
   or a field where coefficient explosion does not occur,
   e.g.~a finite field or the reals (use $\fl=2$).
  Variant: Also available are
   \fun{GEN}{charpoly}{GEN x, long v} ($\fl=5$),
   \fun{GEN}{caract}{GEN A, long v} ($\fl=1$),
   \fun{GEN}{carhess}{GEN A, long v} ($\fl=2$),
   \fun{GEN}{carberkowitz}{GEN A, long v} ($\fl=3$) and
   \fun{GEN}{caradj}{GEN A, long v, GEN *pt}. In this
   last case, if \var{pt} is not \kbd{NULL}, \kbd{*pt} receives the address of
   the adjoint matrix of $A$ (see \tet{matadjoint}), so both can be obtained at
   once.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.charpoly0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def norm(*argv):
  '''
  norm
  Class: basic
  Section: conversions
  C-Name: gnorm
  Prototype: G
  Help: norm(x): norm of x.
  Doc: 
   algebraic norm of $x$, i.e.~the product of $x$ with
   its conjugate (no square roots are taken), or conjugates for polmods. For
   vectors and matrices, the norm is taken componentwise and hence is not the
   $L^2$-norm (see \kbd{norml2}). Note that the norm of an element of
   $\R$ is its square, so as to be compatible with the complex norm.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def normlp(*argv):
  '''
  normlp
  Class: basic
  Section: linear_algebra
  C-Name: gnormlp
  Prototype: GDGp
  Help: normlp(x,{p}): Lp-norm of x; sup norm if p is omitted.
  Doc: 
   $L^p$-norm of $x$; sup norm if $p$ is omitted. More precisely,
   if $x$ is a scalar, \kbd{normlp}$(x, p)$ is defined to be \kbd{abs}$(x)$.
   If $x$ is a polynomial, a (row or column) vector or a matrix:
   
   \item  if $p$ is omitted, \kbd{normlp($x$)} is defined recursively as
   $\max_i \kbd{normlp}(x_i))$, where $(x_i)$ run through the components of~$x$.
   In particular, this yields the usual sup norm if $x$ is a polynomial or
   vector with complex components.
   
   \item otherwise, \kbd{normlp($x$, $p$)} is defined recursively as $(\sum_i
   \kbd{normlp}^p(x_i,p))^{1/p}$. In particular, this yields the usual $(\sum
   |x_i|^p)^{1/p}$ if $x$ is a polynomial or vector with complex components.
   
   \bprog
   ? v = [1,-2,3]; normlp(v)      \\ vector
   %1 = 3
   ? M = [1,-2;-3,4]; normlp(M)   \\ matrix
   %2 = 4
   ? T = (1+I) + I*x^2; normlp(T)
   %3 = 1.4142135623730950488016887242096980786
   ? normlp([[1,2], [3,4], 5, 6])   \\ recursively defined
   %4 = 6
   
   ? normlp(v, 1)
   %5 = 6
   ? normlp(M, 1)
   %6 = 10
   ? normlp(T, 1)
   %7 = 2.4142135623730950488016887242096980786
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gnormlp(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def norml2(*argv):
  '''
  norml2
  Class: basic
  Section: linear_algebra
  C-Name: gnorml2
  Prototype: G
  Help: norml2(x): square of the L2-norm of x.
  Doc: square of the $L^2$-norm of $x$. More precisely,
   if $x$ is a scalar, $\kbd{norml2}(x)$ is defined to be the square
   of the complex modulus of $x$ (real \typ{QUAD}s are not supported).
   If $x$ is a polynomial, a (row or column) vector or a matrix, \kbd{norml2($x$)} is
   defined recursively as $\sum_i \kbd{norml2}(x_i)$, where $(x_i)$ run through
   the components of $x$. In particular, this yields the usual $\sum |x_i|^2$
   (resp.~$\sum |x_{i,j}|^2$) if $x$ is a polynomial or vector (resp.~matrix) with
   complex components.
   
   \bprog
   ? norml2( [ 1, 2, 3 ] )      \\ vector
   %1 = 14
   ? norml2( [ 1, 2; 3, 4] )   \\ matrix
   %2 = 30
   ? norml2( 2*I + x )
   %3 = 5
   ? norml2( [ [1,2], [3,4], 5, 6 ] )   \\ recursively defined
   %4 = 91
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gnorml2(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def trace(*argv):
  '''
  trace
  Class: basic
  Section: linear_algebra
  C-Name: gtrace
  Prototype: G
  Help: trace(x): trace of x.
  Doc: this applies to quite general $x$. If $x$ is not a
   matrix, it is equal to the sum of $x$ and its conjugate, except for polmods
   where it is the trace as an algebraic number.
   
   For $x$ a square matrix, it is the ordinary trace. If $x$ is a
   non-square matrix (but not a vector), an error occurs.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtrace(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathess(*argv):
  '''
  mathess
  Class: basic
  Section: linear_algebra
  C-Name: hess
  Prototype: G
  Help: mathess(x): Hessenberg form of x.
  Doc: returns a matrix similar to the square matrix $x$, which is in upper Hessenberg
   form (zero entries below the first subdiagonal).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hess(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matintersect(*argv):
  '''
  matintersect
  Class: basic
  Section: linear_algebra
  C-Name: intersect
  Prototype: GG
  Help: matintersect(x,y): intersection of the vector spaces whose bases are
   the columns of x and y.
  Doc: $x$ and $y$ being two matrices with the same
   number of rows each of whose columns are independent, finds a basis of the
   $\Q$-vector space equal to the intersection of the spaces spanned by the
   columns of $x$ and $y$ respectively. The faster function
   \tet{idealintersect} can be used to intersect fractional ideals (projective
   $\Z_K$ modules of rank $1$); the slower but much more general function
   \tet{nfhnf} can be used to intersect general $\Z_K$-modules.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.intersect(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfjacobi(*argv):
  '''
  qfjacobi
  Class: basic
  Section: linear_algebra
  C-Name: jacobi
  Prototype: Gp
  Help: qfjacobi(A): eigenvalues and orthogonal matrix of eigenvectors of the
   real symmetric matrix A.
  Doc: apply Jacobi's eigenvalue algorithm to the real symmetric matrix $A$.
   This returns $[L, V]$, where
   
   \item $L$ is the vector of (real) eigenvalues of $A$, sorted in increasing
   order,
   
   \item $V$ is the corresponding orthogonal matrix of eigenvectors of $A$.
   
   \bprog
   ? \p19
   ? A = [1,2;2,1]; mateigen(A)
   %1 =
   [-1 1]
   
   [ 1 1]
   ? [L, H] = qfjacobi(A);
   ? L
   %3 = [-1.000000000000000000, 3.000000000000000000]~
   ? H
   %4 =
   [ 0.7071067811865475245 0.7071067811865475244]
   
   [-0.7071067811865475244 0.7071067811865475245]
   ? norml2( (A-L[1])*H[,1] )       \\ approximate eigenvector
   %5 = 9.403954806578300064 E-38
   ? norml2(H*H~ - 1)
   %6 = 2.350988701644575016 E-38   \\ close to orthogonal
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.jacobi(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matadjoint(*argv):
  '''
  matadjoint
  Class: basic
  Section: linear_algebra
  C-Name: matadjoint0
  Prototype: GD0,L,
  Help: matadjoint(M,{flag=0}): adjoint matrix of M using Leverrier-Faddeev's
   algorithm. If flag is 1, compute the characteristic polynomial independently
   first.
  Doc: 
   \idx{adjoint matrix} of $M$, i.e.~a matrix $N$
   of cofactors of $M$, satisfying $M*N=\det(M)*\Id$. $M$ must be a
   (non-necessarily invertible) square matrix of dimension $n$.
   If $\fl$ is 0 or omitted, we try to use Leverrier-Faddeev's algorithm,
   which assumes that $n!$ invertible. If it fails or $\fl = 1$,
   compute $T = \kbd{charpoly}(M)$ independently first and return
   $(-1)^{n-1} (T(x)-T(0))/x$ evaluated at $M$.
   \bprog
   ? a = [1,2,3;3,4,5;6,7,8] * Mod(1,4);
   %2 =
   [Mod(1, 4) Mod(2, 4) Mod(3, 4)]
   
   [Mod(3, 4) Mod(0, 4) Mod(1, 4)]
   
   [Mod(2, 4) Mod(3, 4) Mod(0, 4)]
   @eprog\noindent
   Both algorithms use $O(n^4)$ operations in the base ring, and are usually
   slower than computing the characteristic polynomial or the inverse of $M$
   directly.
  Variant: Also available are
   \fun{GEN}{adj}{GEN x} (\fl=0) and
   \fun{GEN}{adjsafe}{GEN x} (\fl=1).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matadjoint0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matcompanion(*argv):
  '''
  matcompanion
  Class: basic
  Section: linear_algebra
  C-Name: matcompanion
  Prototype: G
  Help: matcompanion(x): companion matrix to polynomial x.
  Doc: 
   the left companion matrix to the non-zero polynomial $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matcompanion(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matrixqz(*argv):
  '''
  matrixqz
  Class: basic
  Section: linear_algebra
  C-Name: matrixqz0
  Prototype: GDG
  Help: matrixqz(A,{p=0}): if p>=0, transforms the rational or integral mxn (m>=n)
   matrix A into an integral matrix with gcd of maximal determinants coprime to
   p. If p=-1, finds a basis of the intersection with Z^n of the lattice spanned
   by the columns of A. If p=-2, finds a basis of the intersection with Z^n of
   the Q-vector space spanned by the columns of A.
  Doc: $A$ being an $m\times n$ matrix in $M_{m,n}(\Q)$, let
   $\text{Im}_\Q A$ (resp.~$\text{Im}_\Z A$) the $\Q$-vector space
   (resp.~the $\Z$-module) spanned by the columns of $A$. This function has
   varying behavior depending on the sign of $p$:
   
   If $p \geq 0$, $A$ is assumed to have maximal rank $n\leq m$. The function
   returns a matrix $B\in M_{m,n}(\Z)$, with $\text{Im}_\Q B = \text{Im}_\Q A$,
   such that the GCD of all its $n\times n$ minors is coprime to
   $p$; in particular, if $p = 0$ (default), this GCD is $1$.
   \bprog
   ? minors(x) = vector(#x[,1], i, matdet(x[^i,]));
   ? A = [3,1/7; 5,3/7; 7,5/7]; minors(A)
   %1 = [4/7, 8/7, 4/7]   \\ determinants of all 2x2 minors
   ? B = matrixqz(A)
   %2 =
   [3 1]
   
   [5 2]
   
   [7 3]
   ? minors(%)
   %3 = [1, 2, 1]   \\ B integral with coprime minors
   @eprog
   
   If $p=-1$, returns the HNF basis of the lattice $\Z^n \cap \text{Im}_\Z A$.
   
   If $p=-2$, returns the HNF basis of the lattice $\Z^n \cap \text{Im}_\Q A$.
   \bprog
   ? matrixqz(A,-1)
   %4 =
   [8 5]
   
   [4 3]
   
   [0 1]
   
   ? matrixqz(A,-2)
   %5 =
   [2 -1]
   
   [1 0]
   
   [0 1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matrixqz0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def minpoly(*argv):
  '''
  minpoly
  Class: basic
  Section: linear_algebra
  C-Name: minpoly
  Prototype: GDn
  Help: minpoly(A,{v='x}): minimal polynomial of the matrix or polmod A.
  Doc: \idx{minimal polynomial}
   of $A$ with respect to the variable $v$., i.e. the monic polynomial $P$
   of minimal degree (in the variable $v$) such that $P(A) = 0$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.minpoly(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfgaussred(*argv):
  '''
  qfgaussred
  Class: basic
  Section: linear_algebra
  C-Name: qfgaussred
  Prototype: G
  Help: qfgaussred(q): square reduction of the (symmetric) matrix q (returns a
   square matrix whose i-th diagonal term is the coefficient of the i-th square
   in which the coefficient of the i-th variable is 1).
  Doc: 
   \idx{decomposition into squares} of the
   quadratic form represented by the symmetric matrix $q$. The result is a
   matrix whose diagonal entries are the coefficients of the squares, and the
   off-diagonal entries on each line represent the bilinear forms. More
   precisely, if $(a_{ij})$ denotes the output, one has
   $$ q(x) = \sum_i a_{ii} (x_i + \sum_{j \neq i} a_{ij} x_j)^2 $$
   \bprog
   ? qfgaussred([0,1;1,0])
   %1 =
   [1/2 1]
   
   [-1 -1/2]
   @eprog\noindent This means that $2xy = (1/2)(x+y)^2 - (1/2)(x-y)^2$.
  Variant: \fun{GEN}{qfgaussred_positive}{GEN q} assumes that $q$ is
    positive definite and is a little faster; returns \kbd{NULL} if a vector
    with negative norm occurs (non positive matrix or too many rounding errors).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfgaussred(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfsign(*argv):
  '''
  qfsign
  Class: basic
  Section: linear_algebra
  C-Name: qfsign
  Prototype: G
  Help: qfsign(x): signature of the symmetric matrix x.
  Doc: 
   returns $[p,m]$ the signature of the quadratic form represented by the
   symmetric matrix $x$. Namely, $p$ (resp.~$m$) is the number of positive
   (resp.~negative) eigenvalues of $x$.The result is computed using Gaussian
   reduction.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfsign(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def apply(*argv):
  '''
  apply
  Class: basic
  Section: programming/specific
  C-Name: apply0
  Prototype: GG
  Help: apply(f, A): apply function f to each entry in A.
  Wrapper: (G)
  Description: 
    (closure,gen):gen    genapply(${1 cookie}, ${1 wrapper}, $2)
  Doc: Apply the \typ{CLOSURE} \kbd{f} to the entries of \kbd{A}. If \kbd{A}
   is a scalar, return \kbd{f(A)}. If \kbd{A} is a polynomial or power series,
   apply \kbd{f} on all coefficients. If \kbd{A} is a vector or list, return
   the elements $f(x)$ where $x$ runs through \kbd{A}. If \kbd{A} is a matrix,
   return the matrix whose entries are the $f(\kbd{A[i,j]})$.
   \bprog
   ? apply(x->x^2, [1,2,3,4])
   %1 = [1, 4, 9, 16]
   ? apply(x->x^2, [1,2;3,4])
   %2 =
   [1 4]
   
   [9 16]
   ? apply(x->x^2, 4*x^2 + 3*x+ 2)
   %3 = 16*x^2 + 9*x + 4
   @eprog\noindent Note that many functions already act componentwise on
   vectors or matrices, but they almost never act on lists; in this
   case, \kbd{apply} is a good solution:
   \bprog
   ? L = List([Mod(1,3), Mod(2,4)]);
   ? lift(L)
     ***   at top-level: lift(L)
     ***                 ^-------
     *** lift: incorrect type in lift.
   ? apply(lift, L);
   %2 = List([1, 2])
   @eprog
   \misctitle{Remark} For $v$ a \typ{VEC}, \typ{COL}, \typ{LIST} or \typ{MAT},
   the alternative set-notations
   \bprog
   [g(x) | x <- v, f(x)]
   [x | x <- v, f(x)]
   [g(x) | x <- v]
   @eprog\noindent
   are available as shortcuts for
   \bprog
   apply(g, select(f, Vec(v)))
   select(f, Vec(v))
   apply(g, Vec(v))
   @eprog\noindent respectively:
   \bprog
   ? L = List([Mod(1,3), Mod(2,4)]);
   ? [ lift(x) | x<-L ]
   %2 = [1, 2]
   @eprog
   
   \synt{genapply}{void *E, GEN (*fun)(void*,GEN), GEN a}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.apply0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matdiagonal(*argv):
  '''
  matdiagonal
  Class: basic
  Section: linear_algebra
  C-Name: diagonal
  Prototype: G
  Help: matdiagonal(x): creates the diagonal matrix whose diagonal entries are
   the entries of the vector x.
  Doc: $x$ being a vector, creates the diagonal matrix
   whose diagonal entries are those of $x$.
   \bprog
   ? matdiagonal([1,2,3]);
   %1 =
   [1 0 0]
   
   [0 2 0]
   
   [0 0 3]
   @eprog\noindent Block diagonal matrices are easily created using
   \tet{matconcat}:
   \bprog
   ? U=[1,2;3,4]; V=[1,2,3;4,5,6;7,8,9];
   ? matconcat(matdiagonal([U, V]))
   %1 =
   [1 2 0 0 0]
   
   [3 4 0 0 0]
   
   [0 0 1 2 3]
   
   [0 0 4 5 6]
   
   [0 0 7 8 9]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.diagonal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def vecextract(*argv):
  '''
  vecextract
  Class: basic
  Section: linear_algebra
  C-Name: extract0
  Prototype: GGDG
  Help: vecextract(x,y,{z}): extraction of the components of the matrix or
   vector x according to y and z. If z is omitted, y represents columns, otherwise
   y corresponds to rows and z to columns. y and z can be vectors (of indices),
   strings (indicating ranges as in "1..10") or masks (integers whose binary
   representation indicates the indices to extract, from left to right 1, 2, 4,
   8, etc.).
  Description: 
   (vec,gen,?gen):vec  extract0($1, $2, $3)
  Doc: extraction of components of the vector or matrix $x$ according to $y$.
   In case $x$ is a matrix, its components are the \emph{columns} of $x$. The
   parameter $y$ is a component specifier, which is either an integer, a string
   describing a range, or a vector.
   
   If $y$ is an integer, it is considered as a mask: the binary bits of $y$ are
   read from right to left, but correspond to taking the components from left to
   right. For example, if $y=13=(1101)_2$ then the components 1,3 and 4 are
   extracted.
   
   If $y$ is a vector (\typ{VEC}, \typ{COL} or \typ{VECSMALL}), which must have
   integer entries, these entries correspond to the component numbers to be
   extracted, in the order specified.
   
   If $y$ is a string, it can be
   
   \item a single (non-zero) index giving a component number (a negative
   index means we start counting from the end).
   
   \item a range of the form \kbd{"$a$..$b$"}, where $a$ and $b$ are
   indexes as above. Any of $a$ and $b$ can be omitted; in this case, we take
   as default values $a = 1$ and $b = -1$, i.e.~ the first and last components
   respectively. We then extract all components in the interval $[a,b]$, in
   reverse order if $b < a$.
   
   In addition, if the first character in the string is \kbd{\pow}, the
   complement of the given set of indices is taken.
   
   If $z$ is not omitted, $x$ must be a matrix. $y$ is then the \emph{row}
   specifier, and $z$ the \emph{column} specifier, where the component specifier
   is as explained above.
   
   \bprog
   ? v = [a, b, c, d, e];
   ? vecextract(v, 5)         \\@com mask
   %1 = [a, c]
   ? vecextract(v, [4, 2, 1]) \\@com component list
   %2 = [d, b, a]
   ? vecextract(v, "2..4")    \\@com interval
   %3 = [b, c, d]
   ? vecextract(v, "-1..-3")  \\@com interval + reverse order
   %4 = [e, d, c]
   ? vecextract(v, "^2")      \\@com complement
   %5 = [a, c, d, e]
   ? vecextract(matid(3), "2..", "..")
   %6 =
   [0 1 0]
   
   [0 0 1]
   @eprog
   The range notations \kbd{v[i..j]} and \kbd{v[\pow i]} (for \typ{VEC} or
   \typ{COL}) and \kbd{M[i..j, k..l]} and friends (for \typ{MAT}) implement a
   subset of the above, in a simpler and \emph{faster} way, hence should be
   preferred in most common situations. The following features are not
   implemented in the range notation:
   
   \item reverse order,
   
   \item omitting either $a$ or $b$ in \kbd{$a$..$b$}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.extract0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Mat(*argv):
  '''
  Mat
  Class: basic
  Section: conversions
  C-Name: gtomat
  Prototype: DG
  Help: Mat({x=[]}): transforms any GEN x into a matrix. Empty matrix if x is
   omitted.
  Doc: 
   transforms the object $x$ into a matrix.
   If $x$ is already a matrix, a copy of $x$ is created.
   If $x$ is a row (resp. column) vector, this creates a 1-row (resp.
   1-column) matrix, \emph{unless} all elements are column (resp.~row) vectors
   of the same length, in which case the vectors are concatenated sideways
   and the associated big matrix is returned.
   If $x$ is a binary quadratic form, creates the associated $2\times 2$
   matrix. Otherwise, this creates a $1\times 1$ matrix containing $x$.
   
   \bprog
   ? Mat(x + 1)
   %1 =
   [x + 1]
   ? Vec( matid(3) )
   %2 = [[1, 0, 0]~, [0, 1, 0]~, [0, 0, 1]~]
   ? Mat(%)
   %3 =
   [1 0 0]
   
   [0 1 0]
   
   [0 0 1]
   ? Col( [1,2; 3,4] )
   %4 = [[1, 2], [3, 4]]~
   ? Mat(%)
   %5 =
   [1 2]
   
   [3 4]
   ? Mat(Qfb(1,2,3))
   %6 =
   [1 1]
   
   [1 3]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtomat(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mattranspose(*argv):
  '''
  mattranspose
  Class: basic
  Section: linear_algebra
  C-Name: gtrans
  Prototype: G
  Help: mattranspose(x): x~ = transpose of x.
  Doc: transpose of $x$ (also $x\til$).
   This has an effect only on vectors and matrices.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtrans(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matmuldiagonal(*argv):
  '''
  matmuldiagonal
  Class: basic
  Section: linear_algebra
  C-Name: matmuldiagonal
  Prototype: GG
  Help: matmuldiagonal(x,d): product of matrix x by diagonal matrix whose
   diagonal coefficients are those of the vector d, equivalent but faster than
   x*matdiagonal(d).
  Doc: product of the matrix $x$ by the diagonal
   matrix whose diagonal entries are those of the vector $d$. Equivalent to,
   but much faster than $x*\kbd{matdiagonal}(d)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matmuldiagonal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matmultodiagonal(*argv):
  '''
  matmultodiagonal
  Class: basic
  Section: linear_algebra
  C-Name: matmultodiagonal
  Prototype: GG
  Help: matmultodiagonal(x,y): product of matrices x and y, knowing that the
   result will be a diagonal matrix. Much faster than general multiplication in
   that case.
  Doc: product of the matrices $x$ and $y$ assuming that the result is a
   diagonal matrix. Much faster than $x*y$ in that case. The result is
   undefined if $x*y$ is not diagonal.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matmultodiagonal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def parapply(*argv):
  '''
  parapply
  Class: basic
  Section: programming/parallel
  C-Name: parapply
  Prototype: GG
  Help: parapply(f, x): parallel evaluation of f on the elements of x.
  Doc: parallel evaluation of \kbd{f} on the elements of \kbd{x}.
   The function \kbd{f} must not access global variables or variables
   declared with local(), and must be free of side effects.
   \bprog
   parapply(factor,[2^256 + 1, 2^193 - 1])
   @eprog
   factors $2^{256} + 1$ and $2^{193} - 1$ in parallel.
   \bprog
   {
     my(E = ellinit([1,3]), V = vector(12,i,randomprime(2^200)));
     parapply(p->ellcard(E,p), V)
   }
   @eprog
   computes the order of $E(\F_p)$ for $12$ random primes of $200$ bits.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.parapply(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def parselect(*argv):
  '''
  parselect
  Class: basic
  Section: programming/parallel
  C-Name: parselect
  Prototype: GGD0,L,
  Help: parselect(f, A, {flag = 0}): (parallel select) selects elements of A
   according to the selection function f which is tested in parallel. If flag
   is 1, return the indices of those elements (indirect selection)
  Doc: selects elements of $A$ according to the selection function $f$, done in
   parallel.  If \fl is $1$, return the indices of those elements (indirect
   selection) The function \kbd{f} must not access global variables or
   variables declared with local(), and must be free of side effects.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.parselect(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def select(*argv):
  '''
  select
  Class: basic
  Section: programming/specific
  C-Name: select0
  Prototype: GGD0,L,
  Help: select(f, A, {flag = 0}): selects elements of A according to the selection
   function f. If flag is 1, return the indices of those elements (indirect
   selection)
  Wrapper: (bG)
  Description: 
    (gen,gen):gen    genselect(${1 cookie}, ${1 wrapper}, $2)
    (gen,gen,0):gen  genselect(${1 cookie}, ${1 wrapper}, $2)
    (gen,gen,1):gen  genindexselect(${1 cookie}, ${1 wrapper}, $2)
  Doc: We first describe the default behaviour, when $\fl$ is 0 or omitted.
   Given a vector or list \kbd{A} and a \typ{CLOSURE} \kbd{f}, \kbd{select}
   returns the elements $x$ of \kbd{A} such that $f(x)$ is non-zero. In other
   words, \kbd{f} is seen as a selection function returning a boolean value.
   \bprog
   ? select(x->isprime(x), vector(50,i,i^2+1))
   %1 = [2, 5, 17, 37, 101, 197, 257, 401, 577, 677, 1297, 1601]
   ? select(x->(x<100), %)
   %2 = [2, 5, 17, 37]
   @eprog\noindent returns the primes of the form $i^2+1$ for some $i\leq 50$,
   then the elements less than 100 in the preceding result. The \kbd{select}
   function also applies to a matrix \kbd{A}, seen as a vector of columns, i.e. it
   selects columns instead of entries, and returns the matrix whose columns are
   the selected ones.
   
   \misctitle{Remark} For $v$ a \typ{VEC}, \typ{COL}, \typ{LIST} or \typ{MAT},
   the alternative set-notations
   \bprog
   [g(x) | x <- v, f(x)]
   [x | x <- v, f(x)]
   [g(x) | x <- v]
   @eprog\noindent
   are available as shortcuts for
   \bprog
   apply(g, select(f, Vec(v)))
   select(f, Vec(v))
   apply(g, Vec(v))
   @eprog\noindent respectively:
   \bprog
   ? [ x | x <- vector(50,i,i^2+1), isprime(x) ]
   %1 = [2, 5, 17, 37, 101, 197, 257, 401, 577, 677, 1297, 1601]
   @eprog
   
   \noindent If $\fl = 1$, this function returns instead the \emph{indices} of
   the selected elements, and not the elements themselves (indirect selection):
   \bprog
   ? V = vector(50,i,i^2+1);
   ? select(x->isprime(x), V, 1)
   %2 = Vecsmall([1, 2, 4, 6, 10, 14, 16, 20, 24, 26, 36, 40])
   ? vecextract(V, %)
   %3 = [2, 5, 17, 37, 101, 197, 257, 401, 577, 677, 1297, 1601]
   @eprog\noindent
   The following function lists the elements in $(\Z/N\Z)^*$:
   \bprog
   ? invertibles(N) = select(x->gcd(x,N) == 1, [1..N])
   @eprog
   
   \noindent Finally
   \bprog
   ? select(x->x, M)
   @eprog\noindent selects the non-0 entries in \kbd{M}. If the latter is a
   \typ{MAT}, we extract the matrix of non-0 columns. Note that \emph{removing}
   entries instead of selecting them just involves replacing the selection
   function \kbd{f} with its negation:
   \bprog
   ? select(x->!isprime(x), vector(50,i,i^2+1))
   @eprog
   
   \synt{genselect}{void *E, long (*fun)(void*,GEN), GEN a}. Also available
   is \fun{GEN}{genindexselect}{void *E, long (*fun)(void*, GEN), GEN a},
   corresponding to $\fl = 1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.select0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def vecsum(*argv):
  '''
  vecsum
  Class: basic
  Section: linear_algebra
  C-Name: vecsum
  Prototype: G
  Help: vecsum(v): return the sum of the component of the vector v
  Doc: return the sum of the component of the vector $v$
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.vecsum(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def addhelp(*argv):
  '''
  addhelp
  Class: basic
  Section: programming/specific
  C-Name: addhelp
  Prototype: vrs
  Help: addhelp(sym,str): add/change help message for the symbol sym.
  Doc: changes the help message for the symbol \kbd{sym}. The string \var{str}
   is expanded on the spot and stored as the online help for \kbd{sym}. It is
   recommended to document global variables and user functions in this way,
   although \kbd{gp} will not protest if you don't.
   
   You can attach a help text to an alias, but it will never be
   shown: aliases are expanded by the \kbd{?} help operator and we get the help
   of the symbol the alias points to. Nothing prevents you from modifying the
   help of built-in PARI functions. But if you do, we would like to hear why you
   needed it!
   
   Without \tet{addhelp}, the standard help for user functions consists of its
   name and definition.
   \bprog
   gp> f(x) = x^2;
   gp> ?f
   f =
     (x)->x^2
   
   @eprog\noindent Once addhelp is applied to $f$, the function code is no
   longer included. It can still be consulted by typing the function name:
   \bprog
   gp> addhelp(f, "Square")
   gp> ?f
   Square
   
   gp> f
   %2 = (x)->x^2
   @eprog
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.addhelp(*c_arg_tuple)

def alias(*argv):
  '''
  alias
  Class: basic
  Section: programming/specific
  C-Name: alias0
  Prototype: vrr
  Help: alias(newsym,sym): defines the symbol newsym as an alias for the symbol
   sym.
  Doc: defines the symbol \var{newsym} as an alias for the the symbol \var{sym}:
   \bprog
   ? alias("det", "matdet");
   ? det([1,2;3,4])
   %1 = -2
   @eprog\noindent
   You are not restricted to ordinary functions, as in the above example:
   to alias (from/to) member functions, prefix them with `\kbd{\_.}';
   to alias operators, use their internal name, obtained by writing
   \kbd{\_} in lieu of the operators argument: for instance, \kbd{\_!} and
   \kbd{!\_} are the internal names of the factorial and the
   logical negation, respectively.
   \bprog
   ? alias("mod", "_.mod");
   ? alias("add", "_+_");
   ? alias("_.sin", "sin");
   ? mod(Mod(x,x^4+1))
   %2 = x^4 + 1
   ? add(4,6)
   %3 = 10
   ? Pi.sin
   %4 = 0.E-37
   @eprog
   Alias expansion is performed directly by the internal GP compiler.
   Note that since alias is performed at compilation-time, it does not
   require any run-time processing, however it only affects GP code
   compiled \emph{after} the alias command is evaluated. A slower but more
   flexible alternative is to use variables. Compare
   \bprog
   ? fun = sin;
   ? g(a,b) = intnum(t=a,b,fun(t));
   ? g(0, Pi)
   %3 = 2.0000000000000000000000000000000000000
   ? fun = cos;
   ? g(0, Pi)
   %5 = 1.8830410776607851098 E-39
   @eprog\noindent
   with
   \bprog
   ? alias(fun, sin);
   ? g(a,b) = intnum(t=a,b,fun(t));
   ? g(0,Pi)
   %2 = 2.0000000000000000000000000000000000000
   ? alias(fun, cos);  \\ Oops. Does not affect *previous* definition!
   ? g(0,Pi)
   %3 = 2.0000000000000000000000000000000000000
   ? g(a,b) = intnum(t=a,b,fun(t)); \\ Redefine, taking new alias into account
   ? g(0,Pi)
   %5 = 1.8830410776607851098 E-39
   @eprog
   
   A sample alias file \kbd{misc/gpalias} is provided with
   the standard distribution.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.alias0(*c_arg_tuple)

def kill(*argv):
  '''
  kill
  Class: basic
  Section: programming/specific
  C-Name: kill0
  Prototype: vr
  Help: kill(sym): restores the symbol sym to its ``undefined'' status and kill
   associated help messages.
  Doc: restores the symbol \kbd{sym} to its ``undefined'' status, and deletes any
   help messages associated to \kbd{sym} using \kbd{addhelp}. Variable names
   remain known to the interpreter and keep their former priority: you cannot
   make a variable ``less important" by killing it!
   \bprog
   ? z = y = 1; y
   %1 = 1
   ? kill(y)
   ? y            \\ restored to ``undefined'' status
   %2 = y
   ? variable()
   %3 = [x, y, z] \\ but the variable name y is still known, with y > z !
   @eprog\noindent
   For the same reason, killing a user function (which is an ordinary
   variable holding a \typ{CLOSURE}) does not remove its name from the list of
   variable names.
   
   If the symbol is associated to a variable --- user functions being an
   important special case ---, one may use the \idx{quote} operator
   \kbd{a = 'a} to reset variables to their starting values. However, this
   will not delete a help message associated to \kbd{a}, and is also slightly
   slower than \kbd{kill(a)}.
   \bprog
   ? x = 1; addhelp(x, "foo"); x
   %1 = 1
   ? x = 'x; x   \\ same as 'kill', except we don't delete help.
   %2 = x
   ? ?x
   foo
   @eprog\noindent
   On the other hand, \kbd{kill} is the only way to remove aliases and installed
   functions.
   \bprog
   ? alias(fun, sin);
   ? kill(fun);
   
   ? install(addii, GG);
   ? kill(addii);
   @eprog
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.kill0(*c_arg_tuple)

def type(*argv):
  '''
  type
  Class: basic
  Section: programming/specific
  C-Name: type0
  Prototype: G
  Help: type(x): return the type of the GEN x.
  Description: 
   (gen):typ              typ($1)
  Doc: this is useful only under \kbd{gp}. Returns the internal type name of
   the PARI object $x$ as a  string. Check out existing type names with the
   metacommand \b{t}. For example \kbd{type(1)} will return "\typ{INT}".
  Variant: The macro \kbd{typ} is usually simpler to use since it returns a
   \kbd{long} that can easily be matched with the symbols \typ{*}. The name
   \kbd{type} was avoided since it is a reserved identifier for some compilers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.type0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Qfb(*argv):
  '''
  Qfb
  Class: basic
  Section: conversions
  C-Name: Qfb0
  Prototype: GGGDGp
  Help: Qfb(a,b,c,{D=0.}): binary quadratic form a*x^2+b*x*y+c*y^2. D is
   optional (0.0 by default) and initializes Shanks's distance if b^2-4*a*c>0.
  Doc: creates the binary quadratic form\sidx{binary quadratic form}
   $ax^2+bxy+cy^2$. If $b^2-4ac>0$, initialize \idx{Shanks}' distance
   function to $D$. Negative definite forms are not implemented,
   use their positive definite counterpart instead.
  Variant: Also available are
   \fun{GEN}{qfi}{GEN a, GEN b, GEN c} (assumes $b^2-4ac<0$) and
   \fun{GEN}{qfr}{GEN a, GEN b, GEN c, GEN D} (assumes $b^2-4ac>0$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_params.append(argv[4])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Qfb0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbnucomp(*argv):
  '''
  qfbnucomp
  Class: basic
  Section: number_theoretical
  C-Name: nucomp
  Prototype: GGG
  Help: qfbnucomp(x,y,L): composite of primitive positive definite quadratic
   forms x and y using nucomp and nudupl, where L=[|D/4|^(1/4)] is precomputed.
  Doc: \idx{composition} of the primitive positive
   definite binary quadratic forms $x$ and $y$ (type \typ{QFI}) using the NUCOMP
   and NUDUPL algorithms of \idx{Shanks}, \`a la Atkin. $L$ is any positive
   constant, but for optimal speed, one should take $L=|D|^{1/4}$, where $D$ is
   the common discriminant of $x$ and $y$. When $x$ and $y$ do not have the same
   discriminant, the result is undefined.
   
   The current implementation is straightforward and in general \emph{slower}
   than the generic routine (since the latter takes advantage of asymptotically
   fast operations and careful optimizations).
  Variant: Also available is \fun{GEN}{nudupl}{GEN x, GEN L} when $x=y$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nucomp(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbnupow(*argv):
  '''
  qfbnupow
  Class: basic
  Section: number_theoretical
  C-Name: nupow
  Prototype: GG
  Help: qfbnupow(x,n): n-th power of primitive positive definite quadratic
   form x using nucomp and nudupl.
  Doc: $n$-th power of the primitive positive definite
   binary quadratic form $x$ using \idx{Shanks}'s NUCOMP and NUDUPL algorithms
   (see \kbd{qfbnucomp}, in particular the final warning).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nupow(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbprimeform(*argv):
  '''
  qfbprimeform
  Class: basic
  Section: number_theoretical
  C-Name: primeform
  Prototype: GGp
  Help: qfbprimeform(x,p): returns the prime form of discriminant x, whose
   first coefficient is p.
  Doc: prime binary quadratic form of discriminant
   $x$ whose first coefficient is $p$, where $|p|$ is a prime number.
   By abuse of notation,
   $p = \pm 1$ is also valid and returns the unit form. Returns an
   error if $x$ is not a quadratic residue mod $p$, or if $x < 0$ and $p < 0$.
   (Negative definite \typ{QFI} are not implemented.) In the case where $x>0$,
   the ``distance'' component of the form is set equal to zero according to the
   current precision.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.primeform(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbcompraw(*argv):
  '''
  qfbcompraw
  Class: basic
  Section: number_theoretical
  C-Name: qfbcompraw
  Prototype: GG
  Help: qfbcompraw(x,y): Gaussian composition without reduction of the binary
   quadratic forms x and y.
  Doc: \idx{composition} of the binary quadratic forms $x$ and $y$, without
   \idx{reduction} of the result. This is useful e.g.~to compute a generating
   element of an ideal. The result is undefined if $x$ and $y$ do not have the
   same discriminant.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbcompraw(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbpowraw(*argv):
  '''
  qfbpowraw
  Class: basic
  Section: number_theoretical
  C-Name: qfbpowraw
  Prototype: GL
  Help: qfbpowraw(x,n): n-th power without reduction of the binary quadratic
   form x.
  Doc: $n$-th power of the binary quadratic form
   $x$, computed without doing any \idx{reduction} (i.e.~using \kbd{qfbcompraw}).
   Here $n$ must be non-negative and $n<2^{31}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbpowraw(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbred(*argv):
  '''
  qfbred
  Class: basic
  Section: number_theoretical
  C-Name: qfbred0
  Prototype: GD0,L,DGDGDG
  Help: qfbred(x,{flag=0},{d},{isd},{sd}): reduction of the binary
   quadratic form x. All other args. are optional. The arguments d, isd and
   sd, if
   present, supply the values of the discriminant, floor(sqrt(d)) and sqrt(d)
   respectively. If d<0, its value is not used and all references to Shanks's
   distance hereafter are meaningless. flag can be any of 0: default, uses
   Shanks's distance function d; 1: use d, do a single reduction step; 2: do
   not use d; 3: do not use d, single reduction step.
  Doc: reduces the binary quadratic form $x$ (updating Shanks's distance function
   if $x$ is indefinite). The binary digits of $\fl$ are toggles meaning
   
   \quad 1: perform a single \idx{reduction} step
   
   \quad 2: don't update \idx{Shanks}'s distance
   
   The arguments $d$, \var{isd}, \var{sd}, if present, supply the values of the
   discriminant, $\floor{\sqrt{d}}$, and $\sqrt{d}$ respectively
   (no checking is done of these facts). If $d<0$ these values are useless,
   and all references to Shanks's distance are irrelevant.
  Variant: Also available are
   
   \fun{GEN}{redimag}{GEN x} (for definite $x$),
   
   \noindent and for indefinite forms:
   
   \fun{GEN}{redreal}{GEN x}
   
   \fun{GEN}{rhoreal}{GEN x} (= \kbd{qfbred(x,1)}),
   
   \fun{GEN}{redrealnod}{GEN x, GEN isd} (= \kbd{qfbred(x,2,,isd)}),
   
   \fun{GEN}{rhorealnod}{GEN x, GEN isd} (= \kbd{qfbred(x,3,,isd)}).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_params.append(argv[4].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbred0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbsolve(*argv):
  '''
  qfbsolve
  Class: basic
  Section: number_theoretical
  C-Name: qfbsolve
  Prototype: GG
  Help: qfbsolve(Q,p): Return [x,y] so that Q(x,y)=p where Q is a binary
   quadratic form and p a prime number, or 0 if there is no solution.
  Doc: Solve the equation $Q(x,y)=p$ over the integers,
   where $Q$ is a binary quadratic form and $p$ a prime number.
   
   Return $[x,y]$ as a two-components vector, or zero if there is no solution.
   Note that this function returns only one solution and not all the solutions.
   
   Let $D = \disc Q$. The algorithm used runs in probabilistic polynomial time
   in $p$ (through the computation of a square root of $D$ modulo $p$); it is
   polynomial time in $D$ if $Q$ is imaginary, but exponential time if $Q$ is
   real (through the computation of a full cycle of reduced forms). In the
   latter case, note that \tet{bnfisprincipal} provides a solution in heuristic
   subexponential time in $D$ assuming the GRH.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbsolve(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadgen(*argv):
  '''
  quadgen
  Class: basic
  Section: number_theoretical
  C-Name: quadgen
  Prototype: G
  Help: quadgen(D): standard generator of quadratic order of discriminant D.
  Doc: creates the quadratic
   number\sidx{omega} $\omega=(a+\sqrt{D})/2$ where $a=0$ if $D\equiv0\mod4$,
   $a=1$ if $D\equiv1\mod4$, so that $(1,\omega)$ is an integral basis for the
   quadratic order of discriminant $D$. $D$ must be an integer congruent to 0 or
   1 modulo 4, which is not a square.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadgen(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadpoly(*argv):
  '''
  quadpoly
  Class: basic
  Section: number_theoretical
  C-Name: quadpoly0
  Prototype: GDn
  Help: quadpoly(D,{v='x}): quadratic polynomial corresponding to the
   discriminant D, in variable v.
  Doc: creates the ``canonical'' quadratic
   polynomial (in the variable $v$) corresponding to the discriminant $D$,
   i.e.~the minimal polynomial of $\kbd{quadgen}(D)$. $D$ must be an integer
   congruent to 0 or 1 modulo 4, which is not a square.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadpoly0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bestappr(*argv):
  '''
  bestappr
  Class: basic
  Section: number_theoretical
  C-Name: bestappr
  Prototype: GDG
  Help: bestappr(x, {B}): returns a rational approximation to x, whose
   denominator is limited by B, if present. This function applies to reals,
   intmods, p-adics, and rationals of course. Otherwise it applies recursively
   to all components.
  Doc: using variants of the extended Euclidean algorithm, returns a rational
   approximation $a/b$ to $x$, whose denominator is limited
   by $B$, if present. If $B$ is omitted, return the best approximation
   affordable given the input accuracy; if you are looking for true rational
   numbers, presumably approximated to sufficient accuracy, you should first
   try that option. Otherwise, $B$ must be a positive real scalar (impose
   $0 < b \leq B$).
   
   \item If $x$ is a \typ{REAL} or a \typ{FRAC}, this function uses continued
   fractions.
   \bprog
   ? bestappr(Pi, 100)
   %1 = 22/7
   ? bestappr(0.1428571428571428571428571429)
   %2 = 1/7
   ? bestappr([Pi, sqrt(2) + 'x], 10^3)
   %3 = [355/113, x + 1393/985]
   @eprog
   By definition, $a/b$ is the best rational approximation to $x$ if
   $|b x - a| < |v x - u|$ for all integers $(u,v)$ with $0 < v \leq B$.
   (Which implies that $n/d$ is a convergent of the continued fraction of $x$.)
   
   \item If $x$ is a \typ{INTMOD} modulo $N$ or a \typ{PADIC} of precision $N =
   p^k$, this function performs rational modular reconstruction modulo $N$. The
   routine then returns the unique rational number $a/b$ in coprime integers
   $|a| < N/2B$ and $b\leq B$ which is congruent to $x$ modulo $N$. Omitting
   $B$ amounts to choosing it of the order of $\sqrt{N/2}$. If rational
   reconstruction is not possible (no suitable $a/b$ exists), returns $[]$.
   \bprog
   ? bestappr(Mod(18526731858, 11^10))
   %1 = 1/7
   ? bestappr(Mod(18526731858, 11^20))
   %2 = []
   ? bestappr(3 + 5 + 3*5^2 + 5^3 + 3*5^4 + 5^5 + 3*5^6 + O(5^7))
   %2 = -1/3
   @eprog\noindent In most concrete uses, $B$ is a prime power and we performed
   Hensel lifting to obtain $x$.
   
   The function applies recursively to components of complex objects
   (polynomials, vectors, \dots). If rational reconstruction fails for even a
   single entry, return $[]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bestappr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bestapprPade(*argv):
  '''
  bestapprPade
  Class: basic
  Section: number_theoretical
  C-Name: bestapprPade
  Prototype: GD-1,L,
  Help: bestappr(x, {B}): returns a rational function approximation to x.
   This function applies to series, polmods, and rational functions of course.
   Otherwise it applies recursively to all components.
  Doc: using variants of the extended Euclidean algorithm, returns a rational
   function approximation $a/b$ to $x$, whose denominator is limited
   by $B$, if present. If $B$ is omitted, return the best approximation
   affordable given the input accuracy; if you are looking for true rational
   functions, presumably approximated to sufficient accuracy, you should first
   try that option. Otherwise, $B$ must be a non-negative real (impose
   $0 \leq \text{degree}(b) \leq B$).
   
   \item If $x$ is a \typ{RFRAC} or \typ{SER}, this function uses continued
   fractions.
   \bprog
   ? bestapprPade((1-x^11)/(1-x)+O(x^11))
   %1 = 1/(-x + 1)
   ? bestapprPade([1/(1+x+O(x^10)), (x^3-2)/(x^3+1)], 1)
   %2 =  [1/(x + 1), -2]
   @eprog
   
   \item If $x$ is a \typ{POLMOD} modulo $N$ or a \typ{SER} of precision $N =
   t^k$, this function performs rational modular reconstruction modulo $N$. The
   routine then returns the unique rational function $a/b$ in coprime
   polynomials, with $\text{degree}(b)\leq B$ which is congruent to $x$ modulo
   $N$. Omitting $B$ amounts to choosing it of the order of $N/2$. If rational
   reconstruction is not possible (no suitable $a/b$ exists), returns $[]$.
   \bprog
   ? bestapprPade(Mod(1+x+x^2+x^3+x^4, x^4-2))
   %1 = (2*x - 1)/(x - 1)
   ? % * Mod(1,x^4-2)
   %2 = Mod(x^3 + x^2 + x + 3, x^4 - 2)
   ? bestapprPade(Mod(1+x+x^2+x^3+x^5, x^9))
   %2 = []
   ? bestapprPade(Mod(1+x+x^2+x^3+x^5, x^10))
   %3 = (2*x^4 + x^3 - x - 1)/(-x^5 + x^3 + x^2 - 1)
   @eprog\noindent
   The function applies recursively to components of complex objects
   (polynomials, vectors, \dots). If rational reconstruction fails for even a
   single entry, return $[]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bestapprPade(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def chinese(*argv):
  '''
  chinese
  Class: basic
  Section: number_theoretical
  C-Name: chinese
  Prototype: GDG
  Help: chinese(x,{y}): x,y being both intmods (or polmods) computes z in the
   same residue classes as x and y.
  Description: 
   (gen):gen      chinese1($1)
   (gen, gen):gen chinese($1, $2)
  Doc: if $x$ and $y$ are both intmods or both polmods, creates (with the same
   type) a $z$ in the same residue class as $x$ and in the same residue class as
   $y$, if it is possible.
   \bprog
   ? chinese(Mod(1,2), Mod(2,3))
   %1 = Mod(5, 6)
   ? chinese(Mod(x,x^2-1), Mod(x+1,x^2+1))
   %2 = Mod(-1/2*x^2 + x + 1/2, x^4 - 1)
   @eprog\noindent
   This function also allows vector and matrix arguments, in which case the
   operation is recursively applied to each component of the vector or matrix.
   \bprog
   ? chinese([Mod(1,2),Mod(1,3)], [Mod(1,5),Mod(2,7)])
   %3 = [Mod(1, 10), Mod(16, 21)]
   @eprog\noindent
   For polynomial arguments in the same variable, the function is applied to each
   coefficient; if the polynomials have different degrees, the high degree terms
   are copied verbatim in the result, as if the missing high degree terms in the
   polynomial of lowest degree had been \kbd{Mod(0,1)}. Since the latter
   behavior is usually \emph{not} the desired one, we propose to convert the
   polynomials to vectors of the same length first:
   \bprog
    ? P = x+1; Q = x^2+2*x+1;
    ? chinese(P*Mod(1,2), Q*Mod(1,3))
    %4 = Mod(1, 3)*x^2 + Mod(5, 6)*x + Mod(3, 6)
    ? chinese(Vec(P,3)*Mod(1,2), Vec(Q,3)*Mod(1,3))
    %5 = [Mod(1, 6), Mod(5, 6), Mod(4, 6)]
    ? Pol(%)
    %6 = Mod(1, 6)*x^2 + Mod(5, 6)*x + Mod(4, 6)
   @eprog
   
   If $y$ is omitted, and $x$ is a vector, \kbd{chinese} is applied recursively
   to the components of $x$, yielding a residue belonging to the same class as all
   components of $x$.
   
   Finally $\kbd{chinese}(x,x) = x$ regardless of the type of $x$; this allows
   vector arguments to contain other data, so long as they are identical in both
   vectors.
  Variant: \fun{GEN}{chinese1}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.chinese(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def contfrac(*argv):
  '''
  contfrac
  Class: basic
  Section: number_theoretical
  C-Name: contfrac0
  Prototype: GDGD0,L,
  Help: contfrac(x,{b},{nmax}): continued fraction expansion of x (x
   rational,real or rational function). b and nmax are both optional, where b
   is the vector of numerators of the continued fraction, and nmax is a bound
   for the number of terms in the continued fraction expansion.
  Doc: returns the row vector whose components are the partial quotients of the
   \idx{continued fraction} expansion of $x$. In other words, a result
   $[a_0,\dots,a_n]$ means that $x \approx a_0+1/(a_1+\dots+1/a_n)$. The
   output is normalized so that $a_n \neq 1$ (unless we also have $n = 0$).
   
   The number of partial quotients $n+1$ is limited by \kbd{nmax}. If
   \kbd{nmax} is omitted, the expansion stops at the last significant partial
   quotient.
   \bprog
   ? \p19
     realprecision = 19 significant digits
   ? contfrac(Pi)
   %1 = [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2]
   ? contfrac(Pi,, 3)  \\ n = 2
   %2 = [3, 7, 15]
   @eprog\noindent
   $x$ can also be a rational function or a power series.
   
   If a vector $b$ is supplied, the numerators are equal to the coefficients
   of $b$, instead of all equal to $1$ as above; more precisely, $x \approx
   (1/b_0)(a_0+b_1/(a_1+\dots+b_n/a_n))$; for a numerical continued fraction
   ($x$ real), the $a_i$ are integers, as large as possible; if $x$ is a
   rational function, they are polynomials with $\deg a_i = \deg b_i + 1$.
   The length of the result is then equal to the length of $b$, unless the next
   partial quotient cannot be reliably computed, in which case the expansion
   stops. This happens when a partial remainder is equal to zero (or too small
   compared to the available significant digits for $x$ a \typ{REAL}).
   
   A direct implementation of the numerical continued fraction
   \kbd{contfrac(x,b)} described above would be
   \bprog
   \\ "greedy" generalized continued fraction
   cf(x, b) =
   { my( a= vector(#b), t );
   
     x *= b[1];
     for (i = 1, #b,
       a[i] = floor(x);
       t = x - a[i]; if (!t || i == #b, break);
       x = b[i+1] / t;
     ); a;
   }
   @eprog\noindent There is some degree of freedom when choosing the $a_i$; the
   program above can easily be modified to derive variants of the standard
   algorithm. In the same vein, although no builtin
   function implements the related \idx{Engel expansion} (a special kind of
   \idx{Egyptian fraction} decomposition: $x = 1/a_1 + 1/(a_1a_2) + \dots$ ),
   it can be obtained as follows:
   \bprog
   \\ n terms of the Engel expansion of x
   engel(x, n = 10) =
   { my( u = x, a = vector(n) );
     for (k = 1, n,
       a[k] = ceil(1/u);
       u = u*a[k] - 1;
       if (!u, break);
     ); a
   }
   @eprog
   
   \misctitle{Obsolete hack} (don't use this): If $b$ is an integer, \var{nmax}
   is ignored and the command is understood as \kbd{contfrac($x,, b$)}.
  Variant: Also available are \fun{GEN}{gboundcf}{GEN x, long nmax},
   \fun{GEN}{gcf}{GEN x} and \fun{GEN}{gcf2}{GEN b, GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.contfrac0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def contfracpnqn(*argv):
  '''
  contfracpnqn
  Class: basic
  Section: number_theoretical
  C-Name: contfracpnqn
  Prototype: GD-1,L,
  Help: contfracpnqn(x, {n=-1}): [p_n,p_{n-1}; q_n,q_{n-1}] corresponding to the
   continued fraction x. If n >= 0 is present, returns all convergents up
   to p_n/q_n.
  Doc: when $x$ is a vector or a one-row matrix, $x$
   is considered as the list of partial quotients $[a_0,a_1,\dots,a_n]$ of a
   rational number, and the result is the 2 by 2 matrix
   $[p_n,p_{n-1};q_n,q_{n-1}]$ in the standard notation of continued fractions,
   so $p_n/q_n=a_0+1/(a_1+\dots+1/a_n)$. If $x$ is a matrix with two rows
   $[b_0,b_1,\dots,b_n]$ and $[a_0,a_1,\dots,a_n]$, this is then considered as a
   generalized continued fraction and we have similarly
   $p_n/q_n=(1/b_0)(a_0+b_1/(a_1+\dots+b_n/a_n))$. Note that in this case one
   usually has $b_0=1$.
   
   If $n \geq 0$ is present, returns all convergents from $p_0/q_0$ up to
   $p_n/q_n$. (All convergents if $x$ is too small to compute the $n+1$
   requested convergents.)
   \bprog
   ? a=contfrac(Pi,20)
   %1 = [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2]
   ? contfracpnqn(a,3)
   %2 =
   [3 22 333 355]
   
   [1  7 106 113]
   
   ? contfracpnqn(a,7)
   %3 =
   [3 22 333 355 103993 104348 208341 312689]
   
   [1  7 106 113  33102  33215  66317  99532]
   @eprog
  Variant: also available is \fun{GEN}{pnqn}{GEN x} for $n = -1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.contfracpnqn(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def fibonacci(*argv):
  '''
  fibonacci
  Class: basic
  Section: number_theoretical
  C-Name: fibo
  Prototype: L
  Help: fibonacci(x): fibonacci number of index x (x C-integer).
  Doc: $x^{\text{th}}$ Fibonacci number.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.fibo(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbhclassno(*argv):
  '''
  qfbhclassno
  Class: basic
  Section: number_theoretical
  C-Name: hclassno
  Prototype: G
  Help: qfbhclassno(x): Hurwitz-Kronecker class number of x>0.
  Doc: \idx{Hurwitz class number} of $x$, where
   $x$ is non-negative and congruent to 0 or 3 modulo 4. For $x > 5\cdot
   10^5$, we assume the GRH, and use \kbd{quadclassunit} with default
   parameters.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hclassno(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def hilbert(*argv):
  '''
  hilbert
  Class: basic
  Section: number_theoretical
  C-Name: hilbert
  Prototype: lGGDG
  Help: hilbert(x,y,{p}): Hilbert symbol at p of x,y.
  Doc: \idx{Hilbert symbol} of $x$ and $y$ modulo the prime $p$, $p=0$ meaning
   the place at infinity (the result is undefined if $p\neq 0$ is not prime).
   
   It is possible to omit $p$, in which case we take $p = 0$ if both $x$
   and $y$ are rational, or one of them is a real number. And take $p = q$
   if one of $x$, $y$ is a \typ{INTMOD} modulo $q$ or a $q$-adic. (Incompatible
   types will raise an error.)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.hilbert(*c_arg_tuple)

def isfundamental(*argv):
  '''
  isfundamental
  Class: basic
  Section: number_theoretical
  C-Name: isfundamental
  Prototype: lG
  Help: isfundamental(x): true(1) if x is a fundamental discriminant
   (including 1), false(0) if not.
  Description: 
   (int):bool       Z_isfundamental($1)
   (gen):bool       isfundamental($1)
  Doc: true (1) if $x$ is equal to 1 or to the discriminant of a quadratic
   field, false (0) otherwise.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.isfundamental(*c_arg_tuple)

def ispolygonal(*argv):
  '''
  ispolygonal
  Class: basic
  Section: number_theoretical
  C-Name: ispolygonal
  Prototype: lGGD&
  Help: ispolygonal(x,s,{&N}): true(1) if x is an s-gonal number, false(0) if
   not (s > 2). If N is given set it to n if x is the n-th s-gonal number.
  Doc: true (1) if the integer $x$ is an s-gonal number, false (0) if not.
   The parameter $s > 2$ must be a \typ{INT}. If $N$ is given, set it to $n$
   if $x$ is the $n$-th $s$-gonal number.
   \bprog
   ? ispolygonal(36, 3, &N)
   %1 = 1
   ? N
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.ispolygonal(*c_arg_tuple)

def ispower(*argv):
  '''
  ispower
  Class: basic
  Section: number_theoretical
  C-Name: ispower
  Prototype: lGDGD&
  Help: ispower(x,{k},{&n}): if k > 0 is given, return true (1) if x is a k-th
   power, false (0) if not. If k is omitted, return the maximal k >= 2 such
   that x = n^k is a perfect power, or 0 if no such k exist.
   If n is present, and the function returns a non-zero result, set n to the
   k-th root of x.
  Description: 
   (int):small       Z_isanypower($1, NULL)
   (int, &int):small Z_isanypower($1, &$2)
  Doc: if $k$ is given, returns true (1) if $x$ is a $k$-th power, false
   (0) if not.
   
   If $k$ is omitted, only integers and fractions are allowed for $x$ and the
   function returns the maximal $k \geq 2$ such that $x = n^k$ is a perfect
   power, or 0 if no such $k$ exist; in particular \kbd{ispower(-1)},
   \kbd{ispower(0)}, and \kbd{ispower(1)} all return $0$.
   
   If a third argument $\&n$ is given and $x$ is indeed a $k$-th power, sets
   $n$ to a $k$-th root of $x$.
   
   \noindent For a \typ{FFELT} \kbd{x}, instead of omitting \kbd{k} (which is
   not allowed for this type), it may be natural to set
   \bprog
   k = (x.p ^ poldegree(x.pol) - 1) / fforder(x)
   @eprog
  Variant: Also available is
   \fun{long}{gisanypower}{GEN x, GEN *pty} ($k$ omitted).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.ispower(*c_arg_tuple)

def isprimepower(*argv):
  '''
  isprimepower
  Class: basic
  Section: number_theoretical
  C-Name: isprimepower
  Prototype: lGD&
  Help: isprimepower(x,{&n}): if x = p^k is a prime power (p prime, k > 0),
   return k, else return 0. If n is present, and the function returns a non-zero
   result, set n to p, the k-th root of x.
  Doc: if $x = p^k$ is a prime power ($p$ prime, $k > 0$), return $k$, else
   return 0. If a second argument $\&n$ is given and $x$ is indeed
   the $k$-th power of a prime $p$, sets $n$ to $p$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.isprimepower(*c_arg_tuple)

def issquare(*argv):
  '''
  issquare
  Class: basic
  Section: number_theoretical
  C-Name: issquareall
  Prototype: lGD&
  Help: issquare(x,{&n}): true(1) if x is a square, false(0) if not. If n is
   given puts the exact square root there if it was computed.
  Description: 
   (int):bool        Z_issquare($1)
   (gen):bool        issquare($1)
   (int, &int):bool  Z_issquarerem($1, &$2)
   (gen, &gen):bool  issquareall($1, &$2)
  Doc: true (1) if $x$ is a square, false (0)
   if not. What ``being a square'' means depends on the type of $x$: all
   \typ{COMPLEX} are squares, as well as all non-negative \typ{REAL}; for
   exact types such as \typ{INT}, \typ{FRAC} and \typ{INTMOD}, squares are
   numbers of the form $s^2$ with $s$ in $\Z$, $\Q$ and $\Z/N\Z$ respectively.
   \bprog
   ? issquare(3)          \\ as an integer
   %1 = 0
   ? issquare(3.)         \\ as a real number
   %2 = 1
   ? issquare(Mod(7, 8))  \\ in Z/8Z
   %3 = 0
   ? issquare( 5 + O(13^4) )  \\ in Q_13
   %4 = 0
   @eprog
   If $n$ is given, a square root of $x$ is put into $n$.
   \bprog
   ? issquare(4, &n)
   %1 = 1
   ? n
   %2 = 2
   @eprog
   For polynomials, either we detect that the characteristic is 2 (and check
   directly odd and even-power monomials) or we assume that $2$ is invertible
   and check whether squaring the truncated power series for the square root
   yields the original input.
  Variant: Also available is \fun{long}{issquare}{GEN x}. Deprecated
   GP-specific functions \fun{GEN}{gissquare}{GEN x} and
   \fun{GEN}{gissquareall}{GEN x, GEN *pt} return \kbd{gen\_0} and \kbd{gen\_1}
   instead of a boolean value.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.issquareall(*c_arg_tuple)

def kronecker(*argv):
  '''
  kronecker
  Class: basic
  Section: number_theoretical
  C-Name: kronecker
  Prototype: lGG
  Help: kronecker(x,y): kronecker symbol (x/y).
  Description: 
   (small, small):small  kross($1, $2)
   (int, small):small    krois($1, $2)
   (small, int):small    krosi($1, $2)
   (gen, gen):small      kronecker($1, $2)
  Doc: 
   \idx{Kronecker symbol} $(x|y)$, where $x$ and $y$ must be of type integer. By
   definition, this is the extension of \idx{Legendre symbol} to $\Z \times \Z$
   by total multiplicativity in both arguments with the following special rules
   for $y = 0, -1$ or $2$:
   
   \item $(x|0) = 1$ if $|x| = 1$ and $0$ otherwise.
   
   \item $(x|-1) = 1$ if $x \geq 0$ and $-1$ otherwise.
   
   \item $(x|2) = 0$ if $x$ is even and $1$ if $x = 1,-1 \mod 8$ and $-1$
   if $x=3,-3 \mod 8$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.kronecker(*c_arg_tuple)

def logint(*argv):
  '''
  logint
  Class: basic
  Section: number_theoretical
  C-Name: logint0
  Prototype: lGGD&
  Help: logint(x,b,&z): return the largest integer e so that b^e <= x, where the
   parameters b > 1 and x > 0 are both integers. If the parameter z is present,
   set it to b^e.
  Description: 
   (gen,2):small        expi($1)
   (gen,gen,&int):small logint0($1, $2, &$3)
  Doc: Return the largest integer $e$ so that $b^e \leq x$, where the
   parameters $b > 1$ and $x > 0$ are both integers. If the parameter $z$ is
   present, set it to $b^e$.
   \bprog
   ? logint(1000, 2)
   %1 = 9
   ? 2^9
   %2 = 512
   ? logint(1000, 2, &z)
   %3 = 9
   ? z
   %4 = 512
   @eprog\noindent The number of digits used to write $b$ in base $x$ is
   \kbd{1 + logint(x,b)}:
   \bprog
   ? #digits(1000!, 10)
   %5 = 2568
   ? logint(1000!, 10)
   %6 = 2567
   @eprog\noindent This function may conveniently replace
   \bprog
     floor( log(x) / log(b) )
   @eprog\noindent which may not give the correct answer since PARI
   does not guarantee exact rounding.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.logint0(*c_arg_tuple)

def qfbclassno(*argv):
  '''
  qfbclassno
  Class: basic
  Section: number_theoretical
  C-Name: qfbclassno0
  Prototype: GD0,L,
  Help: qfbclassno(D,{flag=0}): class number of discriminant D using Shanks's
   method by default. If (optional) flag is set to 1, use Euler products.
  Doc: ordinary class number of the quadratic
   order of discriminant $D$. In the present version \vers, a $O(D^{1/2})$
   algorithm is used for $D > 0$ (using Euler product and the functional
   equation) so $D$ should not be too large, say $D < 10^8$, for the time to be
   reasonable. On the other hand, for $D < 0$ one can reasonably compute
   \kbd{qfbclassno($D$)} for $|D|<10^{25}$, since the routine uses
   \idx{Shanks}'s method which is in $O(|D|^{1/4})$. For larger values of $|D|$,
   see \kbd{quadclassunit}.
   
   If $\fl=1$, compute the class number using \idx{Euler product}s and the
   functional equation. However, it is in $O(|D|^{1/2})$.
   
   \misctitle{Important warning} For $D < 0$, this function may give incorrect
   results when the class group has many cyclic factors,
   because implementing \idx{Shanks}'s method in full generality slows it down
   immensely. It is therefore strongly recommended to double-check results using
   either the version with $\fl = 1$ or the function \kbd{quadclassunit}.
   
   \misctitle{Warning} Contrary to what its name implies, this routine does not
   compute the number of classes of binary primitive forms of discriminant $D$,
   which is equal to the \emph{narrow} class number. The two notions are the same
   when $D < 0$ or the fundamental unit $\varepsilon$ has negative norm; when $D
   > 0$ and $N\varepsilon > 0$, the number of classes of forms is twice the
   ordinary class number. This is a problem which we cannot fix for backward
   compatibility reasons. Use the following routine if you are only interested
   in the number of classes of forms:
   \bprog
   QFBclassno(D) =
   qfbclassno(D) * if (D < 0 || norm(quadunit(D)) < 0, 1, 2)
   @eprog\noindent
   Here are a few examples:
   \bprog
   ? qfbclassno(400000028)
   time = 3,140 ms.
   %1 = 1
   ? quadclassunit(400000028).no
   time = 20 ms. \\@com{ much faster}
   %2 = 1
   ? qfbclassno(-400000028)
   time = 0 ms.
   %3 = 7253 \\@com{ correct, and fast enough}
   ? quadclassunit(-400000028).no
   time = 0 ms.
   %4 = 7253
   @eprog\noindent
   See also \kbd{qfbhclassno}.
  Variant: The following functions are also available:
   
   \fun{GEN}{classno}{GEN D} ($\fl = 0$)
   
   \fun{GEN}{classno2}{GEN D} ($\fl = 1$).
   
   \noindent Finally
   
   \fun{GEN}{hclassno}{GEN D} computes the class number of an imaginary
   quadratic field by counting reduced forms, an $O(|D|)$ algorithm.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbclassno0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quaddisc(*argv):
  '''
  quaddisc
  Class: basic
  Section: number_theoretical
  C-Name: quaddisc
  Prototype: G
  Help: quaddisc(x): discriminant of the quadratic field Q(sqrt(x)).
  Doc: discriminant of the quadratic field $\Q(\sqrt{x})$, where $x\in\Q$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quaddisc(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadregulator(*argv):
  '''
  quadregulator
  Class: basic
  Section: number_theoretical
  C-Name: quadregulator
  Prototype: Gp
  Help: quadregulator(x): regulator of the real quadratic field of
   discriminant x.
  Doc: regulator of the quadratic field of positive discriminant $x$. Returns
   an error if $x$ is not a discriminant (fundamental or not) or if $x$ is a
   square. See also \kbd{quadclassunit} if $x$ is large.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadregulator(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadunit(*argv):
  '''
  quadunit
  Class: basic
  Section: number_theoretical
  C-Name: quadunit
  Prototype: G
  Help: quadunit(D): fundamental unit of the quadratic field of discriminant D
   where D must be positive.
  Doc: fundamental unit\sidx{fundamental units} of the
   real quadratic field $\Q(\sqrt D)$ where  $D$ is the positive discriminant
   of the field. If $D$ is not a fundamental discriminant, this probably gives
   the fundamental unit of the corresponding order. $D$ must be an integer
   congruent to 0 or 1 modulo 4, which is not a square; the result is a
   quadratic number (see \secref{se:quadgen}).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadunit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sqrtint(*argv):
  '''
  sqrtint
  Class: basic
  Section: number_theoretical
  C-Name: sqrtint
  Prototype: G
  Help: sqrtint(x): integer square root of x, where x is a non-negative integer.
  Description: 
   (gen):int sqrtint($1)
  Doc: returns the integer square root of $x$, i.e. the largest integer $y$
   such that $y^2 \leq x$, where $x$ a non-negative integer.
   \bprog
   ? N = 120938191237; sqrtint(N)
   %1 = 347761
   ? sqrt(N)
   %2 = 347761.68741970412747602130964414095216
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sqrtint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def znlog(*argv):
  '''
  znlog
  Class: basic
  Section: number_theoretical
  C-Name: znlog
  Prototype: GGDG
  Help: znlog(x,g,{o}): return the discrete logarithm of x in
   (Z/nZ)* in base g. If present, o represents the multiplicative
   order of g. Return [] if no solution exist.
  Doc: discrete logarithm of $x$ in $(\Z/N\Z)^*$ in base $g$.
   The result is $[]$ when $x$ is not a power of $g$.
   If present, $o$ represents the multiplicative order of $g$, see
   \secref{se:DLfun}; the preferred format for this parameter is
   \kbd{[ord, factor(ord)]}, where \kbd{ord} is the order of $g$.
   This provides a definite speedup when the discrete log problem is simple:
   \bprog
   ? p = nextprime(10^4); g = znprimroot(p); o = [p-1, factor(p-1)];
   ? for(i=1,10^4, znlog(i, g, o))
   time = 205 ms.
   ? for(i=1,10^4, znlog(i, g))
   time = 244 ms. \\ a little slower
   @eprog
   
   The result is undefined if $g$ is not invertible mod $N$ or if the supplied
   order is incorrect.
   
   This function uses
   
   \item a combination of generic discrete log algorithms (see below).
   
   \item in $(\Z/N\Z)^*$ when $N$ is prime: a linear sieve index calculus
   method, suitable for $N < 10^{50}$, say, is used for large prime divisors of
   the order.
   
   The generic discrete log algorithms are:
   
   \item Pohlig-Hellman algorithm, to reduce to groups of prime order $q$,
   where $q | p-1$ and $p$ is an odd prime divisor of $N$,
   
   \item Shanks baby-step/giant-step ($q < 2^{32}$ is small),
   
   \item Pollard rho method ($q > 2^{32}$).
   
   The latter two algorithms require $O(\sqrt{q})$ operations in the group on
   average, hence will not be able to treat cases where $q > 10^{30}$, say.
   In addition, Pollard rho is not able to handle the case where there are no
   solutions: it will enter an infinite loop.
   \bprog
   ? g = znprimroot(101)
   %1 = Mod(2,101)
   ? znlog(5, g)
   %2 = 24
   ? g^24
   %3 = Mod(5, 101)
   
   ? G = znprimroot(2 * 101^10)
   %4 = Mod(110462212541120451003, 220924425082240902002)
   ? znlog(5, G)
   %5 = 76210072736547066624
   ? G^% == 5
   %6 = 1
   ? N = 2^4*3^2*5^3*7^4*11; g = Mod(13, N); znlog(g^110, g)
   %7 = 110
   ? znlog(6, Mod(2,3))  \\ no solution
   %8 = []
   @eprog\noindent For convenience, $g$ is also allowed to be a $p$-adic number:
   \bprog
   ? g = 3+O(5^10); znlog(2, g)
   %1 = 1015243
   ? g^%
   %2 = 2 + O(5^10)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.znlog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def znorder(*argv):
  '''
  znorder
  Class: basic
  Section: number_theoretical
  C-Name: znorder
  Prototype: GDG
  Help: znorder(x,{o}): order of the integermod x in (Z/nZ)*.
   Optional o represents a multiple of the order of the element.
  Description: 
   (gen):int             order($1)
   (gen,):int            order($1)
   (gen,int):int         znorder($1, $2)
  Doc: $x$ must be an integer mod $n$, and the
   result is the order of $x$ in the multiplicative group $(\Z/n\Z)^*$. Returns
   an error if $x$ is not invertible.
   The parameter o, if present, represents a non-zero
   multiple of the order of $x$, see \secref{se:DLfun}; the preferred format for
   this parameter is \kbd{[ord, factor(ord)]}, where \kbd{ord = eulerphi(n)}
   is the cardinality of the group.
  Variant: Also available is \fun{GEN}{order}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.znorder(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def znprimroot(*argv):
  '''
  znprimroot
  Class: basic
  Section: number_theoretical
  C-Name: znprimroot
  Prototype: G
  Help: znprimroot(n): returns a primitive root of n when it exists.
  Doc: returns a primitive root (generator) of $(\Z/n\Z)^*$, whenever this
   latter group is cyclic ($n = 4$ or $n = 2p^k$ or $n = p^k$, where $p$ is an
   odd prime and $k \geq 0$). If the group is not cyclic, the result is
   undefined. If $n$ is a prime power, then the smallest positive primitive
   root is returned. This may not be true for $n = 2p^k$, $p$ odd.
   
   Note that this function requires factoring $p-1$ for $p$ as above,
   in order to determine the exact order of elements in
   $(\Z/n\Z)^*$: this is likely to be costly if $p$ is large.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.znprimroot(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def znstar(*argv):
  '''
  znstar
  Class: basic
  Section: number_theoretical
  C-Name: znstar
  Prototype: G
  Help: znstar(n): 3-component vector v, giving the structure of (Z/nZ)^*.
   v[1] is the order (i.e. eulerphi(n)), v[2] is a vector of cyclic components,
   and v[3] is a vector giving the corresponding generators.
  Doc: gives the structure of the multiplicative group
   $(\Z/n\Z)^*$ as a 3-component row vector $v$, where $v[1]=\phi(n)$ is the
   order of that group, $v[2]$ is a $k$-component row-vector $d$ of integers
   $d[i]$ such that $d[i]>1$ and $d[i]\mid d[i-1]$ for $i \ge 2$ and
   $(\Z/n\Z)^* \simeq \prod_{i=1}^k(\Z/d[i]\Z)$, and $v[3]$ is a $k$-component row
   vector giving generators of the image of the cyclic groups $\Z/d[i]\Z$.
   \bprog
   ? G = znstar(40)
   %1 = [16, [4, 2, 2], [Mod(17, 40), Mod(21, 40), Mod(11, 40)]]
   ? G.no   \\ eulerphi(40)
   %2 = 16
   ? G.cyc  \\ cycle structure
   %3 = [4, 2, 2]
   ? G.gen  \\ generators for the cyclic components
   %4 = [Mod(17, 40), Mod(21, 40), Mod(11, 40)]
   ? apply(znorder, G.gen)
   %5 = [4, 2, 2]
   @eprog\noindent According to the above definitions, \kbd{znstar(0)} is
   \kbd{[2, [2], [-1]]}, corresponding to $\Z^*$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.znstar(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def core(*argv):
  '''
  core
  Class: basic
  Section: number_theoretical
  C-Name: core0
  Prototype: GD0,L,
  Help: core(n,{flag=0}): unique squarefree integer d
   dividing n such that n/d is a square. If (optional) flag is non-null, output
   the two-component row vector [d,f], where d is the unique squarefree integer
   dividing n such that n/d=f^2 is a square.
  Doc: if $n$ is an integer written as
   $n=df^2$ with $d$ squarefree, returns $d$. If $\fl$ is non-zero,
   returns the two-element row vector $[d,f]$. By convention, we write $0 = 0
   \times 1^2$, so \kbd{core(0, 1)} returns $[0,1]$.
  Variant: Also available are \fun{GEN}{core}{GEN n} ($\fl = 0$) and
   \fun{GEN}{core2}{GEN n} ($\fl = 1$)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.core0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def coredisc(*argv):
  '''
  coredisc
  Class: basic
  Section: number_theoretical
  C-Name: coredisc0
  Prototype: GD0,L,
  Help: coredisc(n,{flag=0}): discriminant of the quadratic field Q(sqrt(n)).
   If (optional) flag is non-null, output a two-component row vector [d,f],
   where d is the discriminant of the quadratic field Q(sqrt(n)) and n=df^2. f
   may be a half integer.
  Doc: a \emph{fundamental discriminant} is an integer of the form $t\equiv 1
   \mod 4$ or $4t \equiv 8,12 \mod 16$, with $t$ squarefree (i.e.~$1$ or the
   discriminant of a quadratic number field). Given a non-zero integer
   $n$, this routine returns the (unique) fundamental discriminant $d$
   such that $n=df^2$, $f$ a positive rational number. If $\fl$ is non-zero,
   returns the two-element row vector $[d,f]$. If $n$ is congruent to
   0 or 1 modulo 4, $f$ is an integer, and a half-integer otherwise.
   
   By convention, \kbd{coredisc(0, 1))} returns $[0,1]$.
   
   Note that \tet{quaddisc}$(n)$ returns the same value as \kbd{coredisc}$(n)$,
   and also works with rational inputs $n\in\Q^*$.
  Variant: Also available are \fun{GEN}{coredisc}{GEN n} ($\fl = 0$) and
   \fun{GEN}{coredisc2}{GEN n} ($\fl = 1$)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.coredisc0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def digits(*argv):
  '''
  digits
  Class: basic
  Section: conversions
  C-Name: digits
  Prototype: GDG
  Help: digits(x,{b=10}): gives the vector formed by the digits of x in base b (x and b
   integers).
  Doc: 
   outputs the vector of the digits of $|x|$ in base $b$, where $x$ and $b$ are integers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.digits(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def divisors(*argv):
  '''
  divisors
  Class: basic
  Section: number_theoretical
  C-Name: divisors
  Prototype: G
  Help: divisors(x): gives a vector formed by the divisors of x in increasing
   order.
  Description: 
   (gen):vec        divisors($1)
  Doc: creates a row vector whose components are the
   divisors of $x$. The factorization of $x$ (as output by \tet{factor}) can
   be used instead.
   
   By definition, these divisors are the products of the irreducible
   factors of $n$, as produced by \kbd{factor(n)}, raised to appropriate
   powers (no negative exponent may occur in the factorization). If $n$ is
   an integer, they are the positive divisors, in increasing order.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.divisors(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sumdigits(*argv):
  '''
  sumdigits
  Class: basic
  Section: number_theoretical
  C-Name: sumdigits
  Prototype: G
  Help: sumdigits(n): sum of (decimal) digits in the integer n.
  Doc: sum of (decimal) digits in the integer $n$.
   \bprog
   ? sumdigits(123456789)
   %1 = 45
   @eprog\noindent Other bases that 10 are not supported. Note that the sum of
   bits in $n$ is returned by \tet{hammingweight}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sumdigits(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def zetak(*argv):
  '''
  zetak
  Class: basic
  Section: number_fields
  C-Name: gzetakall
  Prototype: GGD0,L,p
  Help: zetak(nfz,x,{flag=0}): Dedekind zeta function of the number field nfz
   at x, where nfz is the vector computed by zetakinit (NOT by nfinit); flag is
   optional, and can be 0: default, compute zetak, or non-zero: compute the
   lambdak function, i.e. with the gamma factors.
  Doc: \var{znf} being a number
   field initialized by \kbd{zetakinit} (\emph{not} by \kbd{nfinit}),
   computes the value of the \idx{Dedekind} zeta function of the number
   field at the complex number $x$. If $\fl=1$ computes Dedekind $\Lambda$
   function instead (i.e.~the product of the Dedekind zeta function by its gamma
   and exponential factors).
   
   \misctitle{CAVEAT} This implementation is not satisfactory and must be
   rewritten. In particular
   
   \item The accuracy of the result depends in an essential way on the
   accuracy of both the \kbd{zetakinit} program and the current accuracy.
   Be wary in particular that $x$ of large imaginary part or, on the
   contrary, very close to an ordinary integer will suffer from precision
   loss, yielding fewer significant digits than expected. Computing with 28
   digits of relative accuracy, we have
   \bprog
   ? zeta(3)
   %1 = 1.202056903159594285399738161
   ? zeta(3-1e-20)
   %2 = 1.202056903159594285401719424
   ? zetak(zetakinit(x), 3-1e-20)
   %3 = 1.2020569031595952919  \\ 5 digits are wrong
   ? zetak(zetakinit(x), 3-1e-28)
   %4 = -25.33411749           \\ junk
   @eprog
   
   \item As the precision increases, results become unexpectedly
   completely wrong:
   \bprog
   ? \p100
   ? zetak(zetakinit(x^2-5), -1) - 1/30
   %1 = 7.26691813 E-108    \\ perfect
   ? \p150
   ? zetak(zetakinit(x^2-5), -1) - 1/30
   %2 = -2.486113578 E-156  \\ perfect
   ? \p200
   ? zetak(zetakinit(x^2-5), -1) - 1/30
   %3 = 4.47... E-75        \\ more than half of the digits are wrong
   ? \p250
   ? zetak(zetakinit(x^2-5), -1) - 1/30
   %4 = 1.6 E43             \\ junk
   @eprog
  Variant: See also \fun{GEN}{glambdak}{GEN znf, GEN x, long prec} or
   \fun{GEN}{gzetak}{GEN znf, GEN x, long prec}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gzetakall(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def zetakinit(*argv):
  '''
  zetakinit
  Class: basic
  Section: number_fields
  C-Name: initzeta
  Prototype: Gp
  Help: zetakinit(bnf): compute number field information necessary to use zetak.
   bnf may also be an irreducible polynomial.
  Doc: computes a number of initialization data
   concerning the number field associated to \kbd{bnf} so as to be able
   to compute the \idx{Dedekind} zeta and lambda functions, respectively
   $\kbd{zetak}(x)$ and $\kbd{zetak}(x,1)$, at the current real precision. If
   you do not need the \kbd{bnfinit} data somewhere else, you may call it
   with an irreducible polynomial instead of a \var{bnf}: it will call
   \kbd{bnfinit} itself.
   
   The result is a 9-component vector $v$ whose components are very technical
   and cannot really be used except through the \kbd{zetak} function.
   
   This function is very inefficient and should be rewritten. It needs to
   computes millions of coefficients of the corresponding Dirichlet series if
   the precision is big. Unless the discriminant is small it will not be able
   to handle more than 9 digits of relative precision. For instance,
   \kbd{zetakinit(x\pow 8 - 2)} needs 440MB of memory at default precision.
   
   This function will fail with the message
   \bprog
    *** bnrL1: overflow in zeta_get_N0 [need too many primes].
   @eprog\noindent if the approximate functional equation requires us to sum
   too many terms (if the discriminant of the number field is too large).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.initzeta(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def dirzetak(*argv):
  '''
  dirzetak
  Class: basic
  Section: number_fields
  C-Name: dirzetak
  Prototype: GG
  Help: dirzetak(nf,b): Dirichlet series of the Dedekind zeta function of the
   number field nf up to the bound b-1.
  Doc: gives as a vector the first $b$
   coefficients of the \idx{Dedekind} zeta function of the number field $\var{nf}$
   considered as a \idx{Dirichlet series}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.dirzetak(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfgaloisapply(*argv):
  '''
  nfgaloisapply
  Class: basic
  Section: number_fields
  C-Name: galoisapply
  Prototype: GGG
  Help: nfgaloisapply(nf,aut,x): Apply the Galois automorphism aut to the object
   x (element or ideal) in the number field nf.
  Doc: let $\var{nf}$ be a
   number field as output by \kbd{nfinit}, and let \var{aut} be a \idx{Galois}
   automorphism of $\var{nf}$ expressed by its image on the field generator
   (such automorphisms can be found using \kbd{nfgaloisconj}). The function
   computes the action of the automorphism \var{aut} on the object $x$ in the
   number field; $x$ can be a number field element, or an ideal (possibly
   extended). Because of possible confusion with elements and ideals, other
   vector or matrix arguments are forbidden.
    \bprog
    ? nf = nfinit(x^2+1);
    ? L = nfgaloisconj(nf)
    %2 = [-x, x]~
    ? aut = L[1]; /* the non-trivial automorphism */
    ? nfgaloisapply(nf, aut, x)
    %4 = Mod(-x, x^2 + 1)
    ? P = idealprimedec(nf,5); /* prime ideals above 5 */
    ? nfgaloisapply(nf, aut, P[2]) == P[1]
    %7 = 0 \\ !!!!
    ? idealval(nf, nfgaloisapply(nf, aut, P[2]), P[1])
    %8 = 1
   @eprog\noindent The surprising failure of the equality test (\kbd{\%7}) is
   due to the fact that although the corresponding prime ideals are equal, their
   representations are not. (A prime ideal is specificed by a uniformizer, and
   there is no guarantee that applying automorphisms yields the same elements
   as a direct \kbd{idealprimedec} call.)
   
   The automorphism can also be given as a column vector, representing the
   image of \kbd{Mod(x, nf.pol)} as an algebraic number. This last
   representation is more efficient and should be preferred if a given
   automorphism must be used in many such calls.
   \bprog
    ? nf = nfinit(x^3 - 37*x^2 + 74*x - 37);
    ? l = nfgaloisconj(nf); aut = l[2] \\ @com automorphisms in basistoalg form
    %2 = -31/11*x^2 + 1109/11*x - 925/11
    ? L = matalgtobasis(nf, l); AUT = L[2] \\ @com same in algtobasis form
    %3 = [16, -6, 5]~
    ? v = [1, 2, 3]~; nfgaloisapply(nf, aut, v) == nfgaloisapply(nf, AUT, v)
    %4 = 1 \\ @com same result...
    ? for (i=1,10^5, nfgaloisapply(nf, aut, v))
    time = 1,451 ms.
    ? for (i=1,10^5, nfgaloisapply(nf, AUT, v))
    time = 1,045 ms.  \\ @com but the latter is faster
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisapply(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealfrobenius(*argv):
  '''
  idealfrobenius
  Class: basic
  Section: number_fields
  C-Name: idealfrobenius
  Prototype: GGG
  Help: idealfrobenius(nf,gal,pr): Returns the Frobenius element (pr|nf/Q)
   associated with the unramified prime ideal pr in prid format, in the Galois
   group gal of the number field nf.
  Doc: Let $K$ be the number field defined by $nf$ and assume $K/\Q$ be a
   Galois extension with Galois group given \kbd{gal=galoisinit(nf)},
   and that $pr$ is the prime ideal $\goth{P}$ in prid format, and that
   $\goth{P}$ is unramified.
   This function returns a permutation of \kbd{gal.group} which defines the
   automorphism $\sigma=\left(\goth{P}\over K/\Q \right)$, i.e the Frobenius
   element associated to $\goth{P}$. If $p$ is the unique prime number
   in $\goth{P}$, then $\sigma(x)\equiv x^p\mod\P$ for all $x\in\Z_K$.
   \bprog
   ? nf = nfinit(polcyclo(31));
   ? gal = galoisinit(nf);
   ? pr = idealprimedec(nf,101)[1];
   ? g = idealfrobenius(nf,gal,pr);
   ? galoispermtopol(gal,g)
   %5 = x^8
   @eprog\noindent This is correct since $101\equiv 8\mod{31}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealfrobenius(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealramgroups(*argv):
  '''
  idealramgroups
  Class: basic
  Section: number_fields
  C-Name: idealramgroups
  Prototype: GGG
  Help: idealramgroups(nf,gal,pr): let pr be a prime ideal in prid format, and
   gal the Galois group of the number field nf, return a vector g such that g[1]
   is the decomposition group of pr, g[2] is the inertia group, g[i] is the
   (i-2)th ramification group of pr, all trivial subgroups being omitted.
  Doc: Let $K$ be the number field defined by \var{nf} and assume that $K/\Q$ is
   Galois with Galois group $G$ given by \kbd{gal=galoisinit(nf)}.
   Let \var{pr} be the prime ideal $\goth{P}$ in prid format.
   This function returns a vector $g$ of subgroups of \kbd{gal}
   as follow:
   
   \item \kbd{g[1]} is the decomposition group of $\goth{P}$,
   
   \item \kbd{g[2]} is $G_0(\goth{P})$, the inertia group of $\goth{P}$,
   
   and for $i\geq 2$,
   
   \item \kbd{g[i]} is $G_{i-2}(\goth{P})$, the $i-2$-th \idx{ramification
   group} of $\goth{P}$.
   
   \noindent The length of $g$ is the number of non-trivial groups in the
   sequence, thus is $0$ if $e=1$ and $f=1$, and $1$ if $f>1$ and $e=1$.
   The following function computes the cardinality of a subgroup of $G$,
   as given by the components of $g$:
   \bprog
   card(H) =my(o=H[2]); prod(i=1,#o,o[i]);
   @eprog
   \bprog
   ? nf=nfinit(x^6+3); gal=galoisinit(nf); pr=idealprimedec(nf,3)[1];
   ? g = idealramgroups(nf, gal, pr);
   ? apply(card,g)
   %4 = [6, 6, 3, 3, 3] \\ cardinalities of the G_i
   @eprog
   
   \bprog
   ? nf=nfinit(x^6+108); gal=galoisinit(nf); pr=idealprimedec(nf,2)[1];
   ? iso=idealramgroups(nf,gal,pr)[2]
   %4 = [[Vecsmall([2, 3, 1, 5, 6, 4])], Vecsmall([3])]
   ? nfdisc(galoisfixedfield(gal,iso,1))
   %5 = -3
   @eprog\noindent The field fixed by the inertia group of $2$ is not ramified at
   $2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealramgroups(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfcertify(*argv):
  '''
  nfcertify
  Class: basic
  Section: number_fields
  C-Name: nfcertify
  Prototype: G
  Help: nfcertify(nf): returns a vector of composite integers used to certify
   nf.zk and nf.disc unconditionally (both are correct when the output
   is the empty vector).
  Doc: $\var{nf}$ being as output by
   \kbd{nfinit}, checks whether the integer basis is known unconditionally.
   This is in particular useful when the argument to \kbd{nfinit} was of the
   form $[T, \kbd{listP}]$, specifying a finite list of primes when
   $p$-maximality had to be proven.
   
   The function returns a vector of composite integers. If this vector is
   empty, then \kbd{nf.zk} and \kbd{nf.disc} are correct. Otherwise, the
   result is dubious. In order to obtain a certified result, one must
   completely factor each of the given integers, then \kbd{addprime} each of
   them, then check whether \kbd{nfdisc(nf.pol)} is equal to \kbd{nf.disc}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfcertify(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfinit(*argv):
  '''
  nfinit
  Class: basic
  Section: number_fields
  C-Name: nfinit0
  Prototype: GD0,L,p
  Help: nfinit(pol,{flag=0}): pol being a nonconstant irreducible polynomial,
   gives the vector: [pol,[r1,r2],discf,index,[M,MC,T2,T,different] (see
   manual),r1+r2 first roots, integral basis, matrix of power basis in terms of
   integral basis, multiplication table of basis]. flag is optional and can be
   set to 0: default; 1: do not compute different; 2: first use polred to find
   a simpler polynomial; 3: outputs a two-element vector [nf,Mod(a,P)], where
   nf is as in 2 and Mod(a,P) is a polmod equal to Mod(x,pol) and P=nf.pol.
  Description: 
   (gen, ?0):nf:prec       nfinit0($1, 0, prec)
   (gen, 1):nf:prec        nfinit0($1, 1, prec)
   (gen, 2):nf:prec        nfinit0($1, 2, prec)
   (gen, 3):gen:prec       nfinit0($1, 3, prec)
   (gen, 4):nf:prec        nfinit0($1, 4, prec)
   (gen, 5):gen:prec       nfinit0($1, 5, prec)
   (gen, #small):void      $"incorrect flag in nfinit"
   (gen, small):gen:prec   nfinit0($1, $2, prec)
  Doc: \var{pol} being a non-constant,
   preferably monic, irreducible polynomial in $\Z[X]$, initializes a
   \emph{number field} structure (\kbd{nf}) associated to the field $K$ defined
   by \var{pol}. As such, it's a technical object passed as the first argument
   to most \kbd{nf}\var{xxx} functions, but it contains some information which
   may be directly useful. Access to this information via \emph{member
   functions} is preferred since the specific data organization specified below
   may change in the future. Currently, \kbd{nf} is a row vector with 9
   components:
   
   $\var{nf}[1]$ contains the polynomial \var{pol} (\kbd{\var{nf}.pol}).
   
   $\var{nf}[2]$ contains $[r1,r2]$ (\kbd{\var{nf}.sign}, \kbd{\var{nf}.r1},
   \kbd{\var{nf}.r2}), the number of real and complex places of $K$.
   
   $\var{nf}[3]$ contains the discriminant $d(K)$ (\kbd{\var{nf}.disc}) of $K$.
   
   $\var{nf}[4]$ contains the index of $\var{nf}[1]$ (\kbd{\var{nf}.index}),
   i.e.~$[\Z_K : \Z[\theta]]$, where $\theta$ is any root of $\var{nf}[1]$.
   
   $\var{nf}[5]$ is a vector containing 7 matrices $M$, $G$, \var{roundG}, $T$,
   $MD$, $TI$, $MDI$ useful for certain computations in the number field $K$.
   
   \quad\item $M$ is the $(r1+r2)\times n$ matrix whose columns represent
   the numerical values of the conjugates of the elements of the integral
   basis.
   
   \quad\item $G$ is an $n\times n$ matrix such that $T2 = {}^t G G$,
   where $T2$ is the quadratic form $T_2(x) = \sum |\sigma(x)|^2$, $\sigma$
   running over the embeddings of $K$ into $\C$.
   
   \quad\item \var{roundG} is a rescaled copy of $G$, rounded to nearest
   integers.
   
   \quad\item $T$ is the $n\times n$ matrix whose coefficients are
   $\text{Tr}(\omega_i\omega_j)$ where the $\omega_i$ are the elements of the
   integral basis. Note also that $\det(T)$ is equal to the discriminant of the
   field $K$. Also, when understood as an ideal, the matrix $T^{-1}$
   generates the codifferent ideal.
   
   \quad\item The columns of $MD$ (\kbd{\var{nf}.diff}) express a $\Z$-basis
   of the different of $K$ on the integral basis.
   
   \quad\item $TI$ is equal to the primitive part of $T^{-1}$, which has integral
   coefficients.
   
   \quad\item Finally, $MDI$ is a two-element representation (for faster
   ideal product) of $d(K)$ times the codifferent ideal
   (\kbd{\var{nf}.disc$*$\var{nf}.codiff}, which is an integral ideal). $MDI$
   is only used in \tet{idealinv}.
   
   $\var{nf}[6]$ is the vector containing the $r1+r2$ roots
   (\kbd{\var{nf}.roots}) of $\var{nf}[1]$ corresponding to the $r1+r2$
   embeddings of the number field into $\C$ (the first $r1$ components are real,
   the next $r2$ have positive imaginary part).
   
   $\var{nf}[7]$ is an integral basis for $\Z_K$ (\kbd{\var{nf}.zk}) expressed
   on the powers of~$\theta$. Its first element is guaranteed to be $1$. This
   basis is LLL-reduced with respect to $T_2$ (strictly speaking, it is a
   permutation of such a basis, due to the condition that the first element be
   $1$).
   
   $\var{nf}[8]$ is the $n\times n$ integral matrix expressing the power
   basis in terms of the integral basis, and finally
   
   $\var{nf}[9]$ is the $n\times n^2$ matrix giving the multiplication table
   of the integral basis.
   
   If a non monic polynomial is input, \kbd{nfinit} will transform it into a
   monic one, then reduce it (see $\fl=3$). It is allowed, though not very
   useful given the existence of \tet{nfnewprec}, to input a \kbd{nf} or a
   \kbd{bnf} instead of a polynomial.
   
   \bprog
   ? nf = nfinit(x^3 - 12); \\ initialize number field Q[X] / (X^3 - 12)
   ? nf.pol   \\ defining polynomial
   %2 = x^3 - 12
   ? nf.disc  \\ field discriminant
   %3 = -972
   ? nf.index \\ index of power basis order in maximal order
   %4 = 2
   ? nf.zk    \\ integer basis, lifted to Q[X]
   %5 = [1, x, 1/2*x^2]
   ? nf.sign  \\ signature
   %6 = [1, 1]
   ? factor(abs(nf.disc ))  \\ determines ramified primes
   %7 =
   [2 2]
   
   [3 5]
   ? idealfactor(nf, 2)
   %8 =
   [[2, [0, 0, -1]~, 3, 1, [0, 1, 0]~] 3]  \\ @com $\goth{p}_2^3$
   @eprog
   
   \misctitle{Huge discriminants, helping nfdisc}
   
   In case \var{pol} has a huge discriminant which is difficult to factor,
   it is hard to compute from scratch the maximal order. The special input
   format $[\var{pol}, B]$ is also accepted where \var{pol} is a polynomial as
   above and $B$ has one of the following forms
   
   \item an integer basis, as would be computed by \tet{nfbasis}: a vector of
   polynomials with first element $1$. This is useful if the maximal order is
   known in advance.
   
   \item an argument \kbd{listP} which specifies a list of primes (see
   \tet{nfbasis}). Instead of the maximal order, \kbd{nfinit} then computes an
   order which is maximal at these particular primes as well as the primes
   contained in the private prime table (see \tet{addprimes}). The result is
   unconditionnaly correct when the discriminant \kbd{nf.disc} factors
   completely over this set of primes. The function \tet{nfcertify} automates
   this:
   \bprog
   ? pol = polcompositum(x^5 - 101, polcyclo(7))[1];
   ? nf = nfinit( [pol, 10^3] );
   ? nfcertify(nf)
   %3 = []
   @eprog\noindent A priori, \kbd{nf.zk} defines an order which is only known
   to be maximal at all primes $\leq 10^3$ (no prime $\leq 10^3$ divides
   \kbd{nf.index}). The certification step proves the correctness of the
   computation.
   \medskip
   
   If $\fl=2$: \var{pol} is changed into another polynomial $P$ defining the same
   number field, which is as simple as can easily be found using the
   \tet{polredbest} algorithm, and all the subsequent computations are done
   using this new polynomial. In particular, the first component of the result
   is the modified polynomial.
   
   If $\fl=3$, apply \kbd{polredbest} as in case 2, but outputs
   $[\var{nf},\kbd{Mod}(a,P)]$, where $\var{nf}$ is as before and
   $\kbd{Mod}(a,P)=\kbd{Mod}(x,\var{pol})$ gives the change of
   variables. This is implicit when \var{pol} is not monic: first a linear change
   of variables is performed, to get a monic polynomial, then \kbd{polredbest}.
  Variant: Also available are
   \fun{GEN}{nfinit}{GEN x, long prec} ($\fl = 0$),
   \fun{GEN}{nfinitred}{GEN x, long prec} ($\fl = 2$),
   \fun{GEN}{nfinitred2}{GEN x, long prec} ($\fl = 3$).
   Instead of the above hardcoded numerical flags in \kbd{nfinit0}, one should
   rather use
   
   \fun{GEN}{nfinitall}{GEN x, long flag, long prec}, where \fl\ is an
   or-ed combination of
   
   \item \tet{nf_RED}: find a simpler defining polynomial,
   
   \item \tet{nf_ORIG}: if \tet{nf_RED} set, also return the change of variable,
   
   \item \tet{nf_ROUND2}: \emph{Deprecated}. Slow down the routine by using an
   obsolete normalization algorithm (do not use this one!),
   
   \item \tet{nf_PARTIALFACT}: \emph{Deprecated}. Lazy factorization of the
   polynomial discriminant. Result is conditional unless \kbd{nfcertify}
   can certify it.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfinit0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfisincl(*argv):
  '''
  nfisincl
  Class: basic
  Section: number_fields
  C-Name: nfisincl
  Prototype: GG
  Help: nfisincl(x,y): tests whether the number field x is isomorphic to a
   subfield of y (where x and y are either polynomials or number fields as
   output by nfinit). Return 0 if not, and otherwise all the isomorphisms. If y
   is a number field, a faster algorithm is used.
  Doc: tests whether the number field $K$ defined
   by the polynomial $x$ is conjugate to a subfield of the field $L$ defined
   by $y$ (where $x$ and $y$ must be in $\Q[X]$). If they are not, the output
   is the number 0. If they are, the output is a vector of polynomials, each
   polynomial $a$ representing an embedding of $K$ into $L$, i.e.~being such
   that $y\mid x\circ a$.
   
   If $y$ is a number field (\var{nf}), a much faster algorithm is used
   (factoring $x$ over $y$ using \tet{nffactor}). Before version 2.0.14, this
   wasn't guaranteed to return all the embeddings, hence was triggered by a
   special flag. This is no more the case.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfisincl(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfisisom(*argv):
  '''
  nfisisom
  Class: basic
  Section: number_fields
  C-Name: nfisisom
  Prototype: GG
  Help: nfisisom(x,y): as nfisincl but tests whether x is isomorphic to y.
  Doc: as \tet{nfisincl}, but tests for isomorphism. If either $x$ or $y$ is a
   number field, a much faster algorithm will be used.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfisisom(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfnewprec(*argv):
  '''
  nfnewprec
  Class: basic
  Section: number_fields
  C-Name: nfnewprec
  Prototype: Gp
  Help: nfnewprec(nf): transform the number field data nf into new data using
   the current (usually larger) precision.
  Doc: transforms the number field $\var{nf}$
   into the corresponding data using current (usually larger) precision. This
   function works as expected if $\var{nf}$ is in fact a $\var{bnf}$ (update
   $\var{bnf}$ to current precision) but may be quite slow (many generators of
   principal ideals have to be computed).
  Variant: See also \fun{GEN}{bnfnewprec}{GEN bnf, long prec}
   and \fun{GEN}{bnrnewprec}{GEN bnr, long prec}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfnewprec(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polredord(*argv):
  '''
  polredord
  Class: basic
  Section: number_fields
  C-Name: polredord
  Prototype: G
  Help: polredord(x): reduction of the polynomial x, staying in the same order.
  Doc: finds polynomials with reasonably small
   coefficients and of the same degree as that of $x$ defining suborders of the
   order defined by $x$. One of the polynomials always defines $\Q$ (hence
   is equal to $(x-1)^n$, where $n$ is the degree), and another always defines
   the same order as $x$ if $x$ is irreducible. Useless function: try
   \kbd{polredbest}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polredord(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polgalois(*argv):
  '''
  polgalois
  Class: basic
  Section: number_fields
  C-Name: polgalois
  Prototype: Gp
  Help: polgalois(T): Galois group of the polynomial T (see manual for group
   coding). Return [n, s, k, name] where n is the group order, s the signature,
   k the index and name is the GAP4 name of the transitive group.
  Doc: \idx{Galois} group of the non-constant
   polynomial $T\in\Q[X]$. In the present version \vers, $T$ must be irreducible
   and the degree $d$ of $T$ must be less than or equal to 7. If the
   \tet{galdata} package has been installed, degrees 8, 9, 10 and 11 are also
   implemented. By definition, if $K = \Q[x]/(T)$, this computes the action of
   the Galois group of the Galois closure of $K$ on the $d$ distinct roots of
   $T$, up to conjugacy (corresponding to different root orderings).
   
   The output is a 4-component vector $[n,s,k,name]$ with the
   following meaning: $n$ is the cardinality of the group, $s$ is its signature
   ($s=1$ if the group is a subgroup of the alternating group $A_d$, $s=-1$
   otherwise) and name is a character string containing name of the transitive
   group according to the GAP 4 transitive groups library by Alexander Hulpke.
   
   $k$ is more arbitrary and the choice made up to version~2.2.3 of PARI is rather
   unfortunate: for $d > 7$, $k$ is the numbering of the group among all
   transitive subgroups of $S_d$, as given in ``The transitive groups of degree up
   to eleven'', G.~Butler and J.~McKay, \emph{Communications in Algebra}, vol.~11,
   1983,
   pp.~863--911 (group $k$ is denoted $T_k$ there). And for $d \leq 7$, it was ad
   hoc, so as to ensure that a given triple would denote a unique group.
   Specifically, for polynomials of degree $d\leq 7$, the groups are coded as
   follows, using standard notations
   \smallskip
   In degree 1: $S_1=[1,1,1]$.
   \smallskip
   In degree 2: $S_2=[2,-1,1]$.
   \smallskip
   In degree 3: $A_3=C_3=[3,1,1]$, $S_3=[6,-1,1]$.
   \smallskip
   In degree 4: $C_4=[4,-1,1]$, $V_4=[4,1,1]$, $D_4=[8,-1,1]$, $A_4=[12,1,1]$,
   $S_4=[24,-1,1]$.
   \smallskip
   In degree 5: $C_5=[5,1,1]$, $D_5=[10,1,1]$, $M_{20}=[20,-1,1]$,
   $A_5=[60,1,1]$, $S_5=[120,-1,1]$.
   \smallskip
   In degree 6: $C_6=[6,-1,1]$, $S_3=[6,-1,2]$, $D_6=[12,-1,1]$, $A_4=[12,1,1]$,
   $G_{18}=[18,-1,1]$, $S_4^-=[24,-1,1]$, $A_4\times C_2=[24,-1,2]$,
   $S_4^+=[24,1,1]$, $G_{36}^-=[36,-1,1]$, $G_{36}^+=[36,1,1]$,
   $S_4\times C_2=[48,-1,1]$, $A_5=PSL_2(5)=[60,1,1]$, $G_{72}=[72,-1,1]$,
   $S_5=PGL_2(5)=[120,-1,1]$, $A_6=[360,1,1]$, $S_6=[720,-1,1]$.
   \smallskip
   In degree 7: $C_7=[7,1,1]$, $D_7=[14,-1,1]$, $M_{21}=[21,1,1]$,
   $M_{42}=[42,-1,1]$, $PSL_2(7)=PSL_3(2)=[168,1,1]$, $A_7=[2520,1,1]$,
   $S_7=[5040,-1,1]$.
   \smallskip
   This is deprecated and obsolete, but for reasons of backward compatibility,
   we cannot change this behavior yet. So you can use the default
   \tet{new_galois_format} to switch to a consistent naming scheme, namely $k$ is
   always the standard numbering of the group among all transitive subgroups of
   $S_n$. If this default is in effect, the above groups will be coded as:
   \smallskip
   In degree 1: $S_1=[1,1,1]$.
   \smallskip
   In degree 2: $S_2=[2,-1,1]$.
   \smallskip
   In degree 3: $A_3=C_3=[3,1,1]$, $S_3=[6,-1,2]$.
   \smallskip
   In degree 4: $C_4=[4,-1,1]$, $V_4=[4,1,2]$, $D_4=[8,-1,3]$, $A_4=[12,1,4]$,
   $S_4=[24,-1,5]$.
   \smallskip
   In degree 5: $C_5=[5,1,1]$, $D_5=[10,1,2]$, $M_{20}=[20,-1,3]$,
   $A_5=[60,1,4]$, $S_5=[120,-1,5]$.
   \smallskip
   In degree 6: $C_6=[6,-1,1]$, $S_3=[6,-1,2]$, $D_6=[12,-1,3]$, $A_4=[12,1,4]$,
   $G_{18}=[18,-1,5]$, $A_4\times C_2=[24,-1,6]$, $S_4^+=[24,1,7]$,
   $S_4^-=[24,-1,8]$, $G_{36}^-=[36,-1,9]$, $G_{36}^+=[36,1,10]$,
   $S_4\times C_2=[48,-1,11]$, $A_5=PSL_2(5)=[60,1,12]$, $G_{72}=[72,-1,13]$,
   $S_5=PGL_2(5)=[120,-1,14]$, $A_6=[360,1,15]$, $S_6=[720,-1,16]$.
   \smallskip
   In degree 7: $C_7=[7,1,1]$, $D_7=[14,-1,2]$, $M_{21}=[21,1,3]$,
   $M_{42}=[42,-1,4]$, $PSL_2(7)=PSL_3(2)=[168,1,5]$, $A_7=[2520,1,6]$,
   $S_7=[5040,-1,7]$.
   \smallskip
   
   \misctitle{Warning} The method used is that of resolvent polynomials and is
   sensitive to the current precision. The precision is updated internally but,
   in very rare cases, a wrong result may be returned if the initial precision
   was not sufficient.
  Variant: To enable the new format in library mode,
   set the global variable \tet{new_galois_format} to $1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polgalois(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polred(*argv):
  '''
  polred
  Class: basic
  Section: number_fields
  C-Name: polred0
  Prototype: GD0,L,DG
  Help: polred(T,{flag=0}): Deprecated, use polredbest. Reduction of the
   polynomial T (gives minimal polynomials only). The following binary digits of
   (optional) flag are significant 1: partial reduction, 2: gives also elements.
  Doc: This function is \emph{deprecated}, use \tet{polredbest} instead.
   Finds polynomials with reasonably small coefficients defining subfields of
   the number field defined by $T$. One of the polynomials always defines $\Q$
   (hence is equal to $x-1$), and another always defines the same number field
   as $T$ if $T$ is irreducible.
   
   All $T$ accepted by \tet{nfinit} are also allowed here;
   in particular, the format \kbd{[T, listP]} is recommended, e.g. with
   $\kbd{listP} = 10^5$ or a vector containing all ramified primes. Otherwise,
   the maximal order of $\Q[x]/(T)$ must be computed.
   
   The following binary digits of $\fl$ are significant:
   
   1: Possibly use a suborder of the maximal order. The
   primes dividing the index of the order chosen are larger than
   \tet{primelimit} or divide integers stored in the \tet{addprimes} table.
   This flag is \emph{deprecated}, the \kbd{[T, listP]} format is more
   flexible.
   
   2: gives also elements. The result is a two-column matrix, the first column
   giving primitive elements defining these subfields, the second giving the
   corresponding minimal polynomials.
   \bprog
   ? M = polred(x^4 + 8, 2)
   %1 =
   [1 x - 1]
   
   [1/2*x^2 x^2 + 2]
   
   [1/4*x^3 x^4 + 2]
   
   [x x^4 + 8]
   ? minpoly(Mod(M[2,1], x^4+8))
   %2 = x^2 + 2
   @eprog
   
   \synt{polred}{GEN T} ($\fl = 0$). Also available is
   \fun{GEN}{polred2}{GEN T} ($\fl = 2$). The function \kbd{polred0} is
   deprecated, provided for backward compatibility.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polred0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polredabs(*argv):
  '''
  polredabs
  Class: basic
  Section: number_fields
  C-Name: polredabs0
  Prototype: GD0,L,
  Help: polredabs(T,{flag=0}): a smallest generating polynomial of the number
   field for the T2 norm on the roots, with smallest index for the minimal T2
   norm. flag is optional, whose binary digit mean 1: give the element whose
   characteristic polynomial is the given polynomial. 4: give all polynomials
   of minimal T2 norm (give only one of P(x) and P(-x)).
  Doc: returns a canonical defining polynomial $P$ for the number field
   $\Q[X]/(T)$ defined by $T$, such that the sum of the squares of the modulus
   of the roots (i.e.~the $T_2$-norm) is minimal. Different $T$ defining
   isomorphic number fields will yield the same $P$. All $T$ accepted by
   \tet{nfinit} are also allowed here, e.g. non-monic polynomials, or pairs
   \kbd{[T, listP]} specifying that a non-maximal order may be used.
   
   \misctitle{Warning 1} Using a \typ{POL} $T$ requires fully factoring the
   discriminant of $T$, which may be very hard. The format \kbd{[T, listP]}
   computes only a suborder of the maximal order and replaces this part of the
   algorithm by a polynomial time computation. In that case the polynomial $P$
   is a priori no longer canonical, and it may happen that it does not have
   minimal $T_2$ norm. The routine attempts to certify the result independently
   of this order computation (as per \tet{nfcertify}: we try to prove that the
   order is maximal); if it fails, the routine returns $0$ instead of $P$.
   In order to force an output in that case as well, you may either use
   \tet{polredbest}, or \kbd{polredabs(,16)}, or
   \bprog
     polredabs([T, nfbasis([T, listP])])
   @eprog\noindent (In all three cases, the result is no longer canonical.)
   
   \misctitle{Warning 2} Apart from the factorization of the discriminant of
   $T$, this routine runs in polynomial time for a \emph{fixed} degree.
   But the complexity is exponential in the degree: this routine
   may be exceedingly slow when the number field has many subfields, hence a
   lot of elements of small $T_2$-norm. If you do not need a canonical
   polynomial, the function \tet{polredbest} is in general much faster (it runs
   in polynomial time), and tends to return polynomials with smaller
   discriminants.
   
   The binary digits of $\fl$ mean
   
   1: outputs a two-component row vector $[P,a]$, where $P$ is the default
   output and \kbd{Mod(a, P)} is a root of the original $T$.
   
   4: gives \emph{all} polynomials of minimal $T_2$ norm; of the two polynomials
   $P(x)$ and $\pm P(-x)$, only one is given.
   
   16: Possibly use a suborder of the maximal order, \emph{without} attempting to
   certify the result as in Warning 1: we always return a polynomial and never
   $0$. The result is a priori not canonical.
   
   \bprog
   ? T = x^16 - 136*x^14 + 6476*x^12 - 141912*x^10 + 1513334*x^8 \
         - 7453176*x^6 + 13950764*x^4 - 5596840*x^2 + 46225
   ? T1 = polredabs(T); T2 = polredbest(T);
   ? [ norml2(polroots(T1)), norml2(polroots(T2)) ]
   %3 = [88.0000000, 120.000000]
   ? [ sizedigit(poldisc(T1)), sizedigit(poldisc(T2)) ]
   %4 = [75, 67]
   @eprog
  Variant: Instead of the above hardcoded numerical flags, one should use an
   or-ed combination of
   
   \item \tet{nf_PARTIALFACT}: possibly use a suborder of the maximal order,
   \emph{without} attempting to certify the result.
   
   \item \tet{nf_ORIG}: return $[P, a]$, where \kbd{Mod(a, P)} is a root of $T$.
   
   \item \tet{nf_RAW}: return $[P, b]$, where \kbd{Mod(b, T)} is a root of $P$.
   The algebraic integer $b$ is the raw result produced by the small vectors
   enumeration in the maximal order; $P$ was computed as the characteristic
   polynomial of \kbd{Mod(b, T)}. \kbd{Mod(a, P)} as in \tet{nf_ORIG}
   is obtained with \tet{modreverse}.
   
   \item \tet{nf_ADDZK}: if $r$ is the result produced with some of the above
   flags (of the form $P$ or $[P,c]$), return \kbd{[r,zk]}, where \kbd{zk} is a
   $\Z$-basis for the maximal order of $\Q[X]/(P)$.
   
   \item \tet{nf_ALL}: return a vector of results of the above form, for all
   polynomials of minimal $T_2$-norm.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polredabs0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polredbest(*argv):
  '''
  polredbest
  Class: basic
  Section: number_fields
  C-Name: polredbest
  Prototype: GD0,L,
  Help: polredbest(T,{flag=0}): reduction of the polynomial T (gives minimal
   polynomials only). If flag=1, gives also elements.
  Doc: finds a polynomial with reasonably
   small coefficients defining the same number field as $T$.
   All $T$ accepted by \tet{nfinit} are also allowed here (e.g. non-monic
   polynomials, \kbd{nf}, \kbd{bnf}, \kbd{[T,Z\_K\_basis]}). Contrary to
   \tet{polredabs}, this routine runs in polynomial time, but it offers no
   guarantee as to the minimality of its result.
   
   This routine computes an LLL-reduced basis for the ring of integers of
   $\Q[X]/(T)$, then examines small linear combinations of the basis vectors,
   computing their characteristic polynomials. It returns the \emph{separable}
   $P$ polynomial of smallest discriminant (the one with lexicographically
   smallest \kbd{abs(Vec(P))} in case of ties). This is a good candidate
   for subsequent number field computations, since it guarantees that
   the denominators of algebraic integers, when expressed in the power basis,
   are reasonably small. With no claim of minimality, though.
   
   It can happen that iterating this functions yields better and better
   polynomials, until it stabilizes:
   \bprog
   ? \p5
   ? P = X^12+8*X^8-50*X^6+16*X^4-3069*X^2+625;
   ? poldisc(P)*1.
   %2 = 1.2622 E55
   ? P = polredbest(P);
   ? poldisc(P)*1.
   %4 = 2.9012 E51
   ? P = polredbest(P);
   ? poldisc(P)*1.
   %6 = 8.8704 E44
   @eprog\noindent In this example, the initial polynomial $P$ is the one
   returned by \tet{polredabs}, and the last one is stable.
   
   If $\fl = 1$: outputs a two-component row vector $[P,a]$,  where $P$ is the
   default output and \kbd{Mod(a, P)} is a root of the original $T$.
   \bprog
   ? [P,a] = polredbest(x^4 + 8, 1)
   %1 = [x^4 + 2, Mod(x^3, x^4 + 2)]
   ? charpoly(a)
   %2 = x^4 + 8
   @eprog\noindent In particular, the map $\Q[x]/(T) \to \Q[x]/(P)$,
   $x\mapsto \kbd{Mod(a,P)}$ defines an isomorphism of number fields, which can
   be computed as
   \bprog
     subst(lift(Q), 'x, a)
   @eprog\noindent if $Q$ is a \typ{POLMOD} modulo $T$; \kbd{b = modreverse(a)}
   returns a \typ{POLMOD} giving the inverse of the above map (which should be
   useless since $\Q[x]/(P)$ is a priori a better representation for the number
   field and its elements).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polredbest(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfpolredabs(*argv):
  '''
  rnfpolredabs
  Class: basic
  Section: number_fields
  C-Name: rnfpolredabs
  Prototype: GGD0,L,
  Help: rnfpolredabs(nf,pol,{flag=0}): given a pol with coefficients in nf,
   finds a relative simpler polynomial defining the same field. Binary digits
   of flag mean: 1: return also the element whose characteristic polynomial is
   the given polynomial, 2: return an absolute polynomial, 16: partial
   reduction.
  Doc: THIS FUNCTION IS OBSOLETE: use \tet{rnfpolredbest} instead.
   Relative version of \kbd{polredabs}. Given a monic polynomial \var{pol}
   with coefficients in $\var{nf}$, finds a simpler relative polynomial defining
   the same field. The binary digits of $\fl$ mean
   
   The binary digits of $\fl$ correspond to $1$: add information to convert
   elements to the new representation, $2$: absolute polynomial, instead of
   relative, $16$: possibly use a suborder of the maximal order. More precisely:
   
   0: default, return $P$
   
   1: returns $[P,a]$ where $P$ is the default output and $a$,
   a \typ{POLMOD} modulo $P$, is a root of \var{pol}.
   
   2: returns \var{Pabs}, an absolute, instead of a relative, polynomial.
   Same as but faster than
   \bprog
     rnfequation(nf, rnfpolredabs(nf,pol))
   @eprog
   
   3: returns $[\var{Pabs},a,b]$, where \var{Pabs} is an absolute polynomial
   as above, $a$, $b$ are \typ{POLMOD} modulo \var{Pabs}, roots of \kbd{nf.pol}
   and \var{pol} respectively.
   
   16: possibly use a suborder of the maximal order. This is slower than the
   default when the relative discriminant is smooth, and much faster otherwise.
   See \secref{se:polredabs}.
   
   \misctitle{Warning} In the present implementation, \kbd{rnfpolredabs}
   produces smaller polynomials than \kbd{rnfpolred} and is usually
   faster, but its complexity is still exponential in the absolute degree.
   The function \tet{rnfpolredbest} runs in polynomial time, and  tends  to
   return polynomials with smaller discriminants.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfpolredabs(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfpolredbest(*argv):
  '''
  rnfpolredbest
  Class: basic
  Section: number_fields
  C-Name: rnfpolredbest
  Prototype: GGD0,L,
  Help: rnfpolredbest(nf,pol,{flag=0}): given a pol with coefficients in nf,
   finds a relative polynomial P defining the same field, hopefully simpler
   than pol; flag
   can be 0: default, 1: return [P,a], where a is a root of pol
   2: return an absolute polynomial Pabs, 3:
   return [Pabs, a,b], where a is a root of nf.pol and b is a root of pol.
  Doc: relative version of \kbd{polredbest}. Given a monic polynomial \var{pol}
   with coefficients in $\var{nf}$, finds a simpler relative polynomial $P$
   defining the same field. As opposed to \tet{rnfpolredabs} this function does
   not return a \emph{smallest} (canonical) polynomial with respect to some
   measure, but it does run in polynomial time.
   
   The binary digits of $\fl$ correspond to $1$: add information to convert
   elements to the new representation, $2$: absolute polynomial, instead of
   relative. More precisely:
   
   0: default, return $P$
   
   1: returns $[P,a]$ where $P$ is the default output and $a$,
   a \typ{POLMOD} modulo $P$, is a root of \var{pol}.
   
   2: returns \var{Pabs}, an absolute, instead of a relative, polynomial.
   Same as but faster than
   \bprog
     rnfequation(nf, rnfpolredbest(nf,pol))
   @eprog
   
   3: returns $[\var{Pabs},a,b]$, where \var{Pabs} is an absolute polynomial
   as above, $a$, $b$ are \typ{POLMOD} modulo \var{Pabs}, roots of \kbd{nf.pol}
   and \var{pol} respectively.
   
   \bprog
   ? K = nfinit(y^3-2); pol = x^2 +x*y + y^2;
   ? [P, a] = rnfpolredbest(K,pol,1);
   ? P
   %3 = x^2 - x + Mod(y - 1, y^3 - 2)
   ? a
   %4 = Mod(Mod(2*y^2+3*y+4,y^3-2)*x + Mod(-y^2-2*y-2,y^3-2),
            x^2 - x + Mod(y-1,y^3-2))
   ? subst(K.pol,y,a)
   %5 = 0
   ? [Pabs, a, b] = rnfpolredbest(K,pol,3);
   ? Pabs
   %7 = x^6 - 3*x^5 + 5*x^3 - 3*x + 1
   ? a
   %8 = Mod(-x^2+x+1, x^6-3*x^5+5*x^3-3*x+1)
   ? b
   %9 = Mod(2*x^5-5*x^4-3*x^3+10*x^2+5*x-5, x^6-3*x^5+5*x^3-3*x+1)
   ? subst(K.pol,y,a)
   %10 = 0
   ? substvec(pol,[x,y],[a,b])
   %11 = 0
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfpolredbest(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poltschirnhaus(*argv):
  '''
  poltschirnhaus
  Class: basic
  Section: number_fields
  C-Name: tschirnhaus
  Prototype: G
  Help: poltschirnhaus(x): random Tschirnhausen transformation of the
   polynomial x.
  Doc: applies a random Tschirnhausen
   transformation to the polynomial $x$, which is assumed to be non-constant
   and separable, so as to obtain a new equation for the \'etale algebra
   defined by $x$. This is for instance useful when computing resolvents,
   hence is used by the \kbd{polgalois} function.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.tschirnhaus(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfmodprinit(*argv):
  '''
  nfmodprinit
  Class: basic
  Section: number_fields
  C-Name: nfmodprinit
  Prototype: GG
  Help: nfmodprinit(nf,pr): transform the 5 element row vector pr representing
   a prime ideal into modpr format necessary for all operations mod pr in the
   number field nf (see manual for details about the format).
  Doc: transforms the prime ideal \var{pr} into \tet{modpr} format necessary
   for all operations modulo \var{pr} in the number field \var{nf}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfmodprinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltreducemodpr(*argv):
  '''
  nfeltreducemodpr
  Class: basic
  Section: number_fields
  C-Name: nfreducemodpr
  Prototype: GGG
  Help: nfeltreducemodpr(nf,x,pr): element x modulo pr in nf, where pr is in
   modpr format (see nfmodprinit).
  Doc: given an element $x$ of the number field $\var{nf}$ and a prime ideal
   \var{pr} in \kbd{modpr} format compute a canonical representative for the
   class of $x$ modulo \var{pr}.
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nf\_to\_Fq}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfreducemodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polcompositum(*argv):
  '''
  polcompositum
  Class: basic
  Section: number_fields
  C-Name: polcompositum0
  Prototype: GGD0,L,
  Help: polcompositum(P,Q,{flag=0}): vector of all possible compositums
   of the number fields defined by the polynomials P and Q. If (optional)
   flag is set (i.e non-null), output for each compositum, not only the
   compositum polynomial pol, but a vector [R,a,b,k] where a (resp. b) is a root
   of P (resp. Q) expressed as a polynomial modulo R,
   and a small integer k such that al2+k*al1 is the chosen root of R.
  Doc: \sidx{compositum} $P$ and $Q$
   being squarefree polynomials in $\Z[X]$ in the same variable, outputs
   the simple factors of the \'etale $\Q$-algebra $A = \Q(X, Y) / (P(X), Q(Y))$.
   The factors are given by a list of polynomials $R$ in $\Z[X]$, associated to
   the number field $\Q(X)/ (R)$, and sorted by increasing degree (with respect
   to lexicographic ordering for factors of equal degrees). Returns an error if
   one of the polynomials is not squarefree.
   
   Note that it is more efficient to reduce to the case where $P$ and $Q$ are
   irreducible first. The routine will not perform this for you, since it may be
   expensive, and the inputs are irreducible in most applications anyway. In
   this case, there will be a single factor $R$ if and only if the number
   fields defined by $P$ and $Q$ are disjoint.
   
   Assuming $P$ is irreducible (of smaller degree than $Q$ for efficiency), it
   is in general much faster to proceed as follows
   \bprog
   nf = nfinit(P); L = nffactor(nf, Q)[,1];
   vector(#L, i, rnfequation(nf, L[i]))
   @eprog\noindent
   to obtain the same result. If you are only interested in the degrees of the
   simple factors, the \kbd{rnfequation} instruction can be replaced by a
   trivial \kbd{poldegree(P) * poldegree(L[i])}.
   
   If $\fl=1$, outputs a vector of 4-component vectors $[R,a,b,k]$, where $R$
   ranges through the list of all possible compositums as above, and $a$
   (resp. $b$) expresses the root of $P$ (resp. $Q$) as an element of
   $\Q(X)/(R)$. Finally, $k$ is a small integer such that $b + ka = X$ modulo
   $R$.
   
   A compositum is often defined by a complicated polynomial, which it is
   advisable to reduce before further work. Here is an example involving
   the field $\Q(\zeta_5, 5^{1/5})$:
   \bprog
   ? L = polcompositum(x^5 - 5, polcyclo(5), 1); \\@com list of $[R,a,b,k]$
   ? [R, a] = L[1];  \\@com pick the single factor, extract $R,a$ (ignore $b,k$)
   ? R               \\@com defines the compositum
   %3 = x^20 + 5*x^19 + 15*x^18 + 35*x^17 + 70*x^16 + 141*x^15 + 260*x^14\
   + 355*x^13 + 95*x^12 - 1460*x^11 - 3279*x^10 - 3660*x^9 - 2005*x^8    \
   + 705*x^7 + 9210*x^6 + 13506*x^5 + 7145*x^4 - 2740*x^3 + 1040*x^2     \
   - 320*x + 256
   ? a^5 - 5         \\@com a fifth root of $5$
   %4 = 0
   ? [T, X] = polredbest(R, 1);
   ? T     \\@com simpler defining polynomial for $\Q[x]/(R)$
   %6 = x^20 + 25*x^10 + 5
   ? X     \\ @com root of $R$ in $\Q[y]/(T(y))$
   %7 = Mod(-1/11*x^15 - 1/11*x^14 + 1/22*x^10 - 47/22*x^5 - 29/11*x^4 + 7/22,\
   x^20 + 25*x^10 + 5)
   ? a = subst(a.pol, 'x, X)  \\@com \kbd{a} in the new coordinates
   %8 = Mod(1/11*x^14 + 29/11*x^4, x^20 + 25*x^10 + 5)
   ? a^5 - 5
   %9 = 0
   @eprog
  Variant: Also available are
   \fun{GEN}{compositum}{GEN P, GEN Q} ($\fl = 0$) and
   \fun{GEN}{compositum2}{GEN P, GEN Q} ($\fl = 1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polcompositum0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealprimedec(*argv):
  '''
  idealprimedec
  Class: basic
  Section: number_fields
  C-Name: idealprimedec
  Prototype: GG
  Help: idealprimedec(nf,p): prime ideal decomposition of the prime number p
   in the number field nf as a vector of 5 component vectors [p,a,e,f,b]
   representing the prime ideals pZ_K+a. Z_K, e,f as usual, a as vector of
   components on the integral basis, b Lenstra's constant.
  Doc: computes the prime ideal
   decomposition of the (positive) prime number $p$ in the number field $K$
   represented by \var{nf}. If a non-prime $p$ is given the result is undefined.
   
   The result is a vector of \tev{prid} structures, each representing one of the
   prime ideals above $p$ in the number field $\var{nf}$. The representation
   $\kbd{pr}=[p,a,e,f,\var{mb}]$ of a prime ideal means the following: $a$ and
   is an algebraic integer in the maximal order $\Z_K$ and the prime ideal is
   equal to $\goth{p} = p\Z_K + a\Z_K$;
   $e$ is the ramification index; $f$ is the residual index;
   finally, \var{mb} is the multiplication table associated to the algebraic
   integer $b$ is such that $\goth{p}^{-1}=\Z_K+ b/ p\Z_K$, which is used
   internally to compute valuations. In other words if $p$ is inert,
   then \var{mb} is the integer $1$, and otherwise it's a square \typ{MAT}
   whose $j$-th column is $b \cdot \kbd{nf.zk[j]}$.
   
   The algebraic number $a$ is guaranteed to have a
   valuation equal to 1 at the prime ideal (this is automatic if $e>1$).
   
   The components of \kbd{pr} should be accessed by member functions: \kbd{pr.p},
   \kbd{pr.e}, \kbd{pr.f}, and \kbd{pr.gen} (returns the vector $[p,a]$):
   \bprog
   ? K = nfinit(x^3-2);
   ? L = idealprimedec(K, 5);
   ? #L       \\ 2 primes above 5 in Q(2^(1/3))
   %3 = 2
   ? p1 = L[1]; p2 = L[2];
   ? [p1.e, p1.f]    \\ the first is unramified of degree 1
   %4 = [1, 1]
   ? [p2.e, p2.f]    \\ the second is unramified of degree 2
   %5 = [1, 2]
   ? p1.gen
   %6 = [5, [2, 1, 0]~]
   ? nfbasistoalg(K, %[2])  \\ a uniformizer for p1
   %7 = Mod(x + 2, x^3 - 2)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealprimedec(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfbasis(*argv):
  '''
  rnfbasis
  Class: basic
  Section: number_fields
  C-Name: rnfbasis
  Prototype: GG
  Help: rnfbasis(bnf,M): given a projective Z_K-module M as output by
   rnfpseudobasis or rnfsteinitz, gives either a basis of M if it is free, or an
   n+1-element generating set.
  Doc: let $K$ the field represented by
   \var{bnf}, as output by \kbd{bnfinit}. $M$ is a projective $\Z_K$-module
   of rank $n$ ($M\otimes K$ is an $n$-dimensional $K$-vector space), given by a
   pseudo-basis of size $n$. The routine returns either a true $\Z_K$-basis of
   $M$ (of size $n$) if it exists, or an $n+1$-element generating set of $M$ if
   not.
   
   It is allowed to use an irreducible polynomial $P$ in $K[X]$ instead of $M$,
   in which case, $M$ is defined as the ring of integers of $K[X]/(P)$, viewed
   as a $\Z_K$-module.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfbasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfdedekind(*argv):
  '''
  rnfdedekind
  Class: basic
  Section: number_fields
  C-Name: rnfdedekind
  Prototype: GGDGD0,L,
  Help: rnfdedekind(nf,pol,{pr},{flag=0}): relative Dedekind criterion over the
   number field K, represented by nf, applied to the order O_K[X]/(P),
   modulo the prime ideal pr (at all primes if pr omitted, in which case
   flag is automatically set to 1).
   P is assumed to be monic, irreducible, in O_K[X].
   Returns [max,basis,v], where basis is a pseudo-basis of the
   enlarged order, max is 1 iff this order is pr-maximal, and v is the
   valuation at pr of the order discriminant. If flag is set, just return 1 if
   the order is maximal, and 0 if not.
  Doc: given a number field $K$ coded by $\var{nf}$ and a monic
   polynomial $P\in \Z_K[X]$, irreducible over $K$ and thus defining a relative
   extension $L$ of $K$, applies \idx{Dedekind}'s criterion to the order
   $\Z_K[X]/(P)$, at the prime ideal \var{pr}. It is possible to set \var{pr}
   to a vector of prime ideals (test maximality at all primes in the vector),
   or to omit altogether, in which case maximality at \emph{all} primes is tested;
   in this situation \fl\ is automatically set to $1$.
   
   The default historic behavior (\fl\ is 0 or omitted and \var{pr} is a
   single prime ideal) is not so useful since
   \kbd{rnfpseudobasis} gives more information and is generally not that
   much slower. It returns a 3-component vector $[\var{max}, \var{basis}, v]$:
   
   \item \var{basis} is a pseudo-basis of an enlarged order $O$ produced by
   Dedekind's criterion, containing the original order $\Z_K[X]/(P)$
   with index a power of \var{pr}. Possibly equal to the original order.
   
   \item \var{max} is a flag equal to 1 if the enlarged order $O$
   could be proven to be \var{pr}-maximal and to 0 otherwise; it may still be
   maximal in the latter case if \var{pr} is ramified in $L$,
   
   \item $v$ is the valuation at \var{pr} of the order discriminant.
   
   If \fl\ is non-zero, on the other hand, we just return $1$ if the order
   $\Z_K[X]/(P)$ is \var{pr}-maximal (resp.~maximal at all relevant primes, as
   described above), and $0$ if not. This is much faster than the default,
   since the enlarged order is not computed.
   \bprog
   ? nf = nfinit(y^2-3); P = x^3 - 2*y;
   ? pr3 = idealprimedec(nf,3)[1];
   ? rnfdedekind(nf, P, pr3)
   %2 = [1, [[1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 1, 1]], 8]
   ? rnfdedekind(nf, P, pr3, 1)
   %3 = 1
   @eprog\noindent In this example, \kbd{pr3} is the ramified ideal above $3$,
   and the order generated by the cube roots of $y$ is already
   \kbd{pr3}-maximal. The order-discriminant has valuation $8$. On the other
   hand, the order is not maximal at the prime above 2:
   \bprog
   ? pr2 = idealprimedec(nf,2)[1];
   ? rnfdedekind(nf, P, pr2, 1)
   %5 = 0
   ? rnfdedekind(nf, P, pr2)
   %6 = [0, [[2, 0, 0; 0, 1, 0; 0, 0, 1], [[1, 0; 0, 1], [1, 0; 0, 1],
        [1, 1/2; 0, 1/2]]], 2]
   @eprog
   The enlarged order is not proven to be \kbd{pr2}-maximal yet. In fact, it
   is; it is in fact the maximal order:
   \bprog
   ? B = rnfpseudobasis(nf, P)
   %7 = [[1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 1, [1, 1/2; 0, 1/2]],
        [162, 0; 0, 162], -1]
   ? idealval(nf,B[3], pr2)
   %4 = 2
   @eprog\noindent
   It is possible to use this routine with non-monic
   $P = \sum_{i\leq n} a_i X^i \in \Z_K[X]$ if $\fl = 1$;
   in this case, we test maximality of Dedekind's order generated by
   $$1, a_n \alpha, a_n\alpha^2 + a_{n-1}\alpha, \dots,
   a_n\alpha^{n-1} + a_{n-1}\alpha^{n-2} + \cdots + a_1\alpha.$$
   The routine will fail if $P$ is $0$ on the projective line over the residue
   field $\Z_K/\kbd{pr}$ (FIXME).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfdedekind(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfdet(*argv):
  '''
  rnfdet
  Class: basic
  Section: number_fields
  C-Name: rnfdet
  Prototype: GG
  Help: rnfdet(nf,M): given a pseudo-matrix M, compute its determinant.
  Doc: given a pseudo-matrix $M$ over the maximal
   order of $\var{nf}$, computes its determinant.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfdet(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfdisc(*argv):
  '''
  rnfdisc
  Class: basic
  Section: number_fields
  C-Name: rnfdiscf
  Prototype: GG
  Help: rnfdisc(nf,pol): given a pol with coefficients in nf, gives a
   2-component vector [D,d], where D is the relative ideal discriminant, and d
   is the relative discriminant in nf^*/nf*^2.
  Doc: given a number field $\var{nf}$ as
   output by \kbd{nfinit} and a polynomial \var{pol} with coefficients in
   $\var{nf}$ defining a relative extension $L$ of $\var{nf}$, computes the
   relative discriminant of $L$. This is a two-element row vector $[D,d]$, where
   $D$ is the relative ideal discriminant and $d$ is the relative discriminant
   considered as an element of $\var{nf}^*/{\var{nf}^*}^2$. The main variable of
   $\var{nf}$ \emph{must} be of lower priority than that of \var{pol}, see
   \secref{se:priority}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfdiscf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfequation(*argv):
  '''
  rnfequation
  Class: basic
  Section: number_fields
  C-Name: rnfequation0
  Prototype: GGD0,L,
  Help: rnfequation(nf,pol,{flag=0}): given a pol with coefficients in nf,
   gives an absolute equation z of the number field defined by pol. flag is
   optional, and can be 0: default, or non-zero, gives [z,al,k], where
   z defines the absolute equation L/Q as in the default behavior,
   al expresses as an element of L a root of the polynomial
   defining the base field nf, and k is a small integer such that
   t = b + k al is a root of z, for b a root of pol.
  Doc: given a number field
   $\var{nf}$ as output by \kbd{nfinit} (or simply a polynomial) and a
   polynomial \var{pol} with coefficients in $\var{nf}$ defining a relative
   extension $L$ of $\var{nf}$, computes an absolute equation of $L$ over
   $\Q$.
   
   The main variable of $\var{nf}$ \emph{must} be of lower priority than that
   of \var{pol} (see \secref{se:priority}). Note that for efficiency, this does
   not check whether the relative equation is irreducible over $\var{nf}$, but
   only if it is squarefree. If it is reducible but squarefree, the result will
   be the absolute equation of the \'etale algebra defined by \var{pol}. If
   \var{pol} is not squarefree, raise an \kbd{e\_DOMAIN} exception.
   \bprog
   ? rnfequation(y^2+1, x^2 - y)
   %1 = x^4 + 1
   ? T = y^3-2; rnfequation(nfinit(T), (x^3-2)/(x-Mod(y,T)))
   %2 = x^6 + 108  \\ Galois closure of Q(2^(1/3))
   @eprog
   
   If $\fl$ is non-zero, outputs a 3-component row vector $[z,a,k]$, where
   
   \item $z$ is the absolute equation of $L$ over $\Q$, as in the default
   behavior,
   
   \item $a$ expresses as a \typ{POLMOD} modulo $z$ a root $\alpha$ of the
   polynomial defining the base field $\var{nf}$,
   
   \item $k$ is a small integer such that $\theta = \beta+k\alpha$
   is a root of $z$, where $\beta$ is a root of $\var{pol}$.
   \bprog
   ? T = y^3-2; pol = x^2 +x*y + y^2;
   ? [z,a,k] = rnfequation(T, pol, 1);
   ? z
   %4 = x^6 + 108
   ? subst(T, y, a)
   %5 = 0
   ? alpha= Mod(y, T);
   ? beta = Mod(x*Mod(1,T), pol);
   ? subst(z, x, beta + k*alpha)
   %8 = 0
   @eprog
  Variant: Also available are
   \fun{GEN}{rnfequation}{GEN nf, GEN pol} ($\fl = 0$) and
   \fun{GEN}{rnfequation2}{GEN nf, GEN pol} ($\fl = 1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfequation0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfhnfbasis(*argv):
  '''
  rnfhnfbasis
  Class: basic
  Section: number_fields
  C-Name: rnfhnfbasis
  Prototype: GG
  Help: rnfhnfbasis(bnf,x): given an order x as output by rnfpseudobasis,
   gives either a true HNF basis of the order if it exists, zero otherwise.
  Doc: given $\var{bnf}$ as output by
   \kbd{bnfinit}, and either a polynomial $x$ with coefficients in $\var{bnf}$
   defining a relative extension $L$ of $\var{bnf}$, or a pseudo-basis $x$ of
   such an extension, gives either a true $\var{bnf}$-basis of $L$ in upper
   triangular Hermite normal form, if it exists, and returns $0$ otherwise.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfhnfbasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfisfree(*argv):
  '''
  rnfisfree
  Class: basic
  Section: number_fields
  C-Name: rnfisfree
  Prototype: lGG
  Help: rnfisfree(bnf,x): given an order x as output by rnfpseudobasis or
   rnfsteinitz, outputs true (1) or false (0) according to whether the order is
   free or not.
  Doc: given $\var{bnf}$ as output by
   \kbd{bnfinit}, and either a polynomial $x$ with coefficients in $\var{bnf}$
   defining a relative extension $L$ of $\var{bnf}$, or a pseudo-basis $x$ of
   such an extension, returns true (1) if $L/\var{bnf}$ is free, false (0) if
   not.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.rnfisfree(*c_arg_tuple)

def rnflllgram(*argv):
  '''
  rnflllgram
  Class: basic
  Section: number_fields
  C-Name: rnflllgram
  Prototype: GGGp
  Help: rnflllgram(nf,pol,order): given a pol with coefficients in nf and an
   order as output by rnfpseudobasis or similar, gives [[neworder],U], where
   neworder is a reduced order and U is the unimodular transformation matrix.
  Doc: given a polynomial
   \var{pol} with coefficients in \var{nf} defining a relative extension $L$ and
   a suborder \var{order} of $L$ (of maximal rank), as output by
   \kbd{rnfpseudobasis}$(\var{nf},\var{pol})$ or similar, gives
   $[[\var{neworder}],U]$, where \var{neworder} is a reduced order and $U$ is
   the unimodular transformation matrix.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnflllgram(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfpolred(*argv):
  '''
  rnfpolred
  Class: basic
  Section: number_fields
  C-Name: rnfpolred
  Prototype: GGp
  Help: rnfpolred(nf,pol): given a pol with coefficients in nf, finds a list
   of relative polynomials defining some subfields, hopefully simpler.
  Doc: THIS FUNCTION IS OBSOLETE: use \tet{rnfpolredbest} instead.
   Relative version of \kbd{polred}. Given a monic polynomial \var{pol} with
   coefficients in $\var{nf}$, finds a list of relative polynomials defining some
   subfields, hopefully simpler and containing the original field. In the present
   version \vers, this is slower and less efficient than \kbd{rnfpolredbest}.
   
   \misctitle{Remark} this function is based on an incomplete reduction
   theory of lattices over number fields, implemented by \kbd{rnflllgram}, which
   deserves to be improved.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfpolred(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfpseudobasis(*argv):
  '''
  rnfpseudobasis
  Class: basic
  Section: number_fields
  C-Name: rnfpseudobasis
  Prototype: GG
  Help: rnfpseudobasis(nf,pol): given a pol with coefficients in nf, gives a
   4-component vector [A,I,D,d] where [A,I] is a pseudo basis of the maximal
   order in HNF on the power basis, D is the relative ideal discriminant, and d
   is the relative discriminant in nf^*/nf*^2.
  Doc: given a number field
   $\var{nf}$ as output by \kbd{nfinit} and a polynomial \var{pol} with
   coefficients in $\var{nf}$ defining a relative extension $L$ of $\var{nf}$,
   computes a pseudo-basis $(A,I)$ for the maximal order $\Z_L$ viewed as a
   $\Z_K$-module, and the relative discriminant of $L$. This is output as a
   four-element row vector $[A,I,D,d]$, where $D$ is the relative ideal
   discriminant and $d$ is the relative discriminant considered as an element of
   $\var{nf}^*/{\var{nf}^*}^2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfpseudobasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfsteinitz(*argv):
  '''
  rnfsteinitz
  Class: basic
  Section: number_fields
  C-Name: rnfsteinitz
  Prototype: GG
  Help: rnfsteinitz(nf,x): given an order x as output by rnfpseudobasis,
   gives [A,I,D,d] where (A,I) is a pseudo basis where all the ideals except
   perhaps the last are trivial.
  Doc: given a number field $\var{nf}$ as
   output by \kbd{nfinit} and either a polynomial $x$ with coefficients in
   $\var{nf}$ defining a relative extension $L$ of $\var{nf}$, or a pseudo-basis
   $x$ of such an extension as output for example by \kbd{rnfpseudobasis},
   computes another pseudo-basis $(A,I)$ (not in HNF in general) such that all
   the ideals of $I$ except perhaps the last one are equal to the ring of
   integers of $\var{nf}$, and outputs the four-component row vector $[A,I,D,d]$
   as in \kbd{rnfpseudobasis}. The name of this function comes from the fact
   that the ideal class of the last ideal of $I$, which is well defined, is the
   \idx{Steinitz class} of the $\Z_K$-module $\Z_L$ (its image in $SK_0(\Z_K)$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfsteinitz(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfalgtobasis(*argv):
  '''
  nfalgtobasis
  Class: basic
  Section: number_fields
  C-Name: algtobasis
  Prototype: GG
  Help: nfalgtobasis(nf,x): transforms the algebraic number x into a column
   vector on the integral basis nf.zk.
  Doc: Given an algebraic number $x$ in the number field $\var{nf}$,
   transforms it to a column vector on the integral basis \kbd{\var{nf}.zk}.
   \bprog
   ? nf = nfinit(y^2 + 4);
   ? nf.zk
   %2 = [1, 1/2*y]
   ? nfalgtobasis(nf, [1,1]~)
   %3 = [1, 1]~
   ? nfalgtobasis(nf, y)
   %4 = [0, 2]~
   ? nfalgtobasis(nf, Mod(y, y^2+4))
   %4 = [0, 2]~
   @eprog
   This is the inverse function of \kbd{nfbasistoalg}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.algtobasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfbasistoalg(*argv):
  '''
  nfbasistoalg
  Class: basic
  Section: number_fields
  C-Name: basistoalg
  Prototype: GG
  Help: nfbasistoalg(nf,x): transforms the column vector x on the integral
   basis into an algebraic number.
  Doc: Given an algebraic number $x$ in the number field \kbd{nf}, transforms it
   into \typ{POLMOD} form.
   \bprog
   ? nf = nfinit(y^2 + 4);
   ? nf.zk
   %2 = [1, 1/2*y]
   ? nfbasistoalg(nf, [1,1]~)
   %3 = Mod(1/2*y + 1, y^2 + 4)
   ? nfbasistoalg(nf, y)
   %4 = Mod(y, y^2 + 4)
   ? nfbasistoalg(nf, Mod(y, y^2+4))
   %4 = Mod(y, y^2 + 4)
   @eprog
   This is the inverse function of \kbd{nfalgtobasis}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.basistoalg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ideallist(*argv):
  '''
  ideallist
  Class: basic
  Section: number_fields
  C-Name: ideallist0
  Prototype: GLD4,L,
  Help: ideallist(nf,bound,{flag=4}): vector of vectors L of all idealstar of
   all ideals of norm<=bound. If (optional) flag is present, its binary digits
   are toggles meaning 1: give generators; 2: add units; 4: give only the
   ideals and not the bid.
  Doc: computes the list
   of all ideals of norm less or equal to \var{bound} in the number field
   \var{nf}. The result is a row vector with exactly \var{bound} components.
   Each component is itself a row vector containing the information about
   ideals of a given norm, in no specific order, depending on the value of
   $\fl$:
   
   The possible values of $\fl$ are:
   
   \quad 0: give the \var{bid} associated to the ideals, without generators.
   
   \quad 1: as 0, but include the generators in the \var{bid}.
   
   \quad 2: in this case, \var{nf} must be a \var{bnf} with units. Each
   component is of the form $[\var{bid},U]$, where \var{bid} is as case 0
   and $U$ is a vector of discrete logarithms of the units. More precisely, it
   gives the \kbd{ideallog}s with respect to \var{bid} of \kbd{bnf.tufu}.
   This structure is technical, and only meant to be used in conjunction with
   \tet{bnrclassnolist} or \tet{bnrdisclist}.
   
   \quad 3: as 2, but include the generators in the \var{bid}.
   
   \quad 4: give only the HNF of the ideal.
   
   \bprog
   ? nf = nfinit(x^2+1);
   ? L = ideallist(nf, 100);
   ? L[1]
   %3 = [[1, 0; 0, 1]]  \\@com A single ideal of norm 1
   ? #L[65]
   %4 = 4               \\@com There are 4 ideals of norm 4 in $\Z[i]$
   @eprog
   If one wants more information, one could do instead:
   \bprog
   ? nf = nfinit(x^2+1);
   ? L = ideallist(nf, 100, 0);
   ? l = L[25]; vector(#l, i, l[i].clgp)
   %3 = [[20, [20]], [16, [4, 4]], [20, [20]]]
   ? l[1].mod
   %4 = [[25, 18; 0, 1], []]
   ? l[2].mod
   %5 = [[5, 0; 0, 5], []]
   ? l[3].mod
   %6 = [[25, 7; 0, 1], []]
   @eprog\noindent where we ask for the structures of the $(\Z[i]/I)^*$ for all
   three ideals of norm $25$. In fact, for all moduli with finite part of norm
   $25$ and trivial Archimedean part, as the last 3 commands show. See
   \tet{ideallistarch} to treat general moduli.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ideallist0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ideallistarch(*argv):
  '''
  ideallistarch
  Class: basic
  Section: number_fields
  C-Name: ideallistarch
  Prototype: GGG
  Help: ideallistarch(nf,list,arch): list is a vector of vectors of of bid's as
   output by ideallist. Return a vector of vectors with the same number of
   components as the original list. The leaves give information about
   moduli whose finite part is as in original list, in the same order, and
   Archimedean part is now arch. The information contained is of the same kind
   as was present in the input.
  Doc: 
   \var{list} is a vector of vectors of bid's, as output by \tet{ideallist} with
   flag $0$ to $3$. Return a vector of vectors with the same number of
   components as the original \var{list}. The leaves give information about
   moduli whose finite part is as in original list, in the same order, and
   Archimedean part is now \var{arch} (it was originally trivial). The
   information contained is of the same kind as was present in the input; see
   \tet{ideallist}, in particular the meaning of \fl.
   
   \bprog
   ? bnf = bnfinit(x^2-2);
   ? bnf.sign
   %2 = [2, 0]                         \\@com two places at infinity
   ? L = ideallist(bnf, 100, 0);
   ? l = L[98]; vector(#l, i, l[i].clgp)
   %4 = [[42, [42]], [36, [6, 6]], [42, [42]]]
   ? La = ideallistarch(bnf, L, [1,1]); \\@com add them to the modulus
   ? l = La[98]; vector(#l, i, l[i].clgp)
   %6 = [[168, [42, 2, 2]], [144, [6, 6, 2, 2]], [168, [42, 2, 2]]]
   @eprog
   Of course, the results above are obvious: adding $t$ places at infinity will
   add $t$ copies of $\Z/2\Z$ to the ray class group. The following application
   is more typical:
   \bprog
   ? L = ideallist(bnf, 100, 2);        \\@com units are required now
   ? La = ideallistarch(bnf, L, [1,1]);
   ? H = bnrclassnolist(bnf, La);
   ? H[98];
   %6 = [2, 12, 2]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ideallistarch(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealprincipalunits(*argv):
  '''
  idealprincipalunits
  Class: basic
  Section: number_fields
  C-Name: idealprincipalunits
  Prototype: GGL
  Help: idealprincipalunits(nf,pr,k): returns the structure [no, cyc, gen]
   of the multiplicative group (1 + pr) / (1 + pr^k)^*.
  Doc: given a prime ideal in \tet{idealprimedec} format,
   returns the multiplicative group $(1 + \var{pr}) / (1 + \var{pr}^k)$ as an
   abelian group. This function is much faster than \tet{idealstar} when the
   norm of \var{pr} is large, since it avoids (useless) work in the
   multiplicative group of the residue field.
   \bprog
   ? K = nfinit(y^2+1);
   ? P = idealprimedec(K,2)[1];
   ? G = idealprincipalunits(K, P, 20);
   ? G.cyc
   [512, 256, 4]   \\ Z/512 x Z/256 x Z/4
   ? G.gen
   %5 = [[-1, -2]~, 1021, [0, -1]~] \\ minimal generators of given order
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealprincipalunits(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealstar(*argv):
  '''
  idealstar
  Class: basic
  Section: number_fields
  C-Name: idealstar0
  Prototype: GGD1,L,
  Help: idealstar(nf,I,{flag=1}): gives the structure of (Z_K/I)^*. flag is
   optional, and can be 0: simply gives the structure as a 3-component vector v
   such that v[1] is the order (i.e. eulerphi(I)), v[2] is a vector of cyclic
   components, and v[3] is a vector giving the corresponding generators. If
   flag=1 (default), gives idealstarinit, i.e. a 6-component vector
   [I,v,fa,f2,U,V] where v is as above without the generators, fa is the prime
   ideal factorisation of I and f2, U and V are technical but essential to work
   in (Z_K/I)^*. Finally if flag=2, same as with flag=1 except that the
   generators are also given.
  Doc: outputs a \var{bid} structure,
   necessary for computing in the finite abelian group $G = (\Z_K/I)^*$. Here,
   \var{nf} is a number field and $I$ is a \var{modulus}: either an ideal in any
   form, or a row vector whose first component is an ideal and whose second
   component is a row vector of $r_1$ 0 or 1. Ideals can also be given
   by a factorization into prime ideals, as produced by \tet{idealfactor}.
   
   This \var{bid} is used in \tet{ideallog} to compute discrete logarithms. It
   also contains useful information which can be conveniently retrieved as
   \kbd{\var{bid}.mod} (the modulus),
   \kbd{\var{bid}.clgp} ($G$ as a finite abelian group),
   \kbd{\var{bid}.no} (the cardinality of $G$),
   \kbd{\var{bid}.cyc} (elementary divisors) and
   \kbd{\var{bid}.gen} (generators).
   
   If $\fl=1$ (default), the result is a \var{bid} structure without
   generators.
   
   If $\fl=2$, as $\fl=1$, but including generators, which wastes some time.
   
   If $\fl=0$, only outputs $(\Z_K/I)^*$ as an abelian group,
   i.e as a 3-component vector $[h,d,g]$: $h$ is the order, $d$ is the vector of
   SNF\sidx{Smith normal form} cyclic components and $g$ the corresponding
   generators.
  Variant: Instead the above hardcoded numerical flags, one should rather use
   \fun{GEN}{Idealstar}{GEN nf, GEN ideal, long flag}, where \kbd{flag} is
   an or-ed combination of \tet{nf_GEN} (include generators) and \tet{nf_INIT}
   (return a full \kbd{bid}, not a group), possibly $0$. This offers
   one more combination: gen, but no init.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealstar0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matalgtobasis(*argv):
  '''
  matalgtobasis
  Class: basic
  Section: number_fields
  C-Name: matalgtobasis
  Prototype: GG
  Help: matalgtobasis(nf,x): nfalgtobasis applied to every element of the
   vector or matrix x.
  Doc: $\var{nf}$ being a number field in \kbd{nfinit} format, and $x$ a
   (row or column) vector or matrix, apply \tet{nfalgtobasis} to each entry
   of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matalgtobasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matbasistoalg(*argv):
  '''
  matbasistoalg
  Class: basic
  Section: number_fields
  C-Name: matbasistoalg
  Prototype: GG
  Help: matbasistoalg(nf,x): nfbasistoalg applied to every element of the
   matrix or vector x.
  Doc: $\var{nf}$ being a number field in \kbd{nfinit} format, and $x$ a
   (row or column) vector or matrix, apply \tet{nfbasistoalg} to each entry
   of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matbasistoalg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltadd(*argv):
  '''
  nfeltadd
  Class: basic
  Section: number_fields
  C-Name: nfadd
  Prototype: GGG
  Help: nfadd(nf,x,y): element x+y in nf.
  Doc: 
   given two elements $x$ and $y$ in
   \var{nf}, computes their sum $x+y$ in the number field $\var{nf}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfadd(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltdiv(*argv):
  '''
  nfeltdiv
  Class: basic
  Section: number_fields
  C-Name: nfdiv
  Prototype: GGG
  Help: nfdiv(nf,x,y): element x/y in nf.
  Doc: given two elements $x$ and $y$ in
   \var{nf}, computes their quotient $x/y$ in the number field $\var{nf}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfdiv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltdiveuc(*argv):
  '''
  nfeltdiveuc
  Class: basic
  Section: number_fields
  C-Name: nfdiveuc
  Prototype: GGG
  Help: nfdiveuc(nf,x,y): gives algebraic integer q such that x-by is small.
  Doc: given two elements $x$ and $y$ in
   \var{nf}, computes an algebraic integer $q$ in the number field $\var{nf}$
   such that the components of $x-qy$ are reasonably small. In fact, this is
   functionally identical to \kbd{round(nfdiv(\var{nf},x,y))}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfdiveuc(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltdivrem(*argv):
  '''
  nfeltdivrem
  Class: basic
  Section: number_fields
  C-Name: nfdivrem
  Prototype: GGG
  Help: nfeltdivrem(nf,x,y): gives [q,r] such that r=x-by is small.
  Doc: given two elements $x$ and $y$ in
   \var{nf}, gives a two-element row vector $[q,r]$ such that $x=qy+r$, $q$ is
   an algebraic integer in $\var{nf}$, and the components of $r$ are
   reasonably small.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfdivrem(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltmod(*argv):
  '''
  nfeltmod
  Class: basic
  Section: number_fields
  C-Name: nfmod
  Prototype: GGG
  Help: nfeltmod(nf,x,y): gives r such that r=x-by is small with q algebraic
   integer.
  Doc: 
   given two elements $x$ and $y$ in
   \var{nf}, computes an element $r$ of $\var{nf}$ of the form $r=x-qy$ with
   $q$ and algebraic integer, and such that $r$ is small. This is functionally
   identical to
   $$\kbd{x - nfmul(\var{nf},round(nfdiv(\var{nf},x,y)),y)}.$$
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfmod(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltmul(*argv):
  '''
  nfeltmul
  Class: basic
  Section: number_fields
  C-Name: nfmul
  Prototype: GGG
  Help: nfmul(nf,x,y): element x.y in nf.
  Doc: 
   given two elements $x$ and $y$ in
   \var{nf}, computes their product $x*y$ in the number field $\var{nf}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfmul(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltnorm(*argv):
  '''
  nfeltnorm
  Class: basic
  Section: number_fields
  C-Name: nfnorm
  Prototype: GG
  Help: nfeltnorm(nf,x): norm of x.
  Doc: returns the absolute norm of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltpow(*argv):
  '''
  nfeltpow
  Class: basic
  Section: number_fields
  C-Name: nfpow
  Prototype: GGG
  Help: nfeltpow(nf,x,k): element x^k in nf.
  Doc: given an element $x$ in \var{nf}, and a positive or negative integer $k$,
   computes $x^k$ in the number field $\var{nf}$.
  Variant: \fun{GEN}{nfinv}{GEN nf, GEN x} correspond to $k = -1$, and
   \fun{GEN}{nfsqr}{GEN nf,GEN x} to $k = 2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfpow(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfelttrace(*argv):
  '''
  nfelttrace
  Class: basic
  Section: number_fields
  C-Name: nftrace
  Prototype: GG
  Help: nfelttrace(nf,x): trace of x.
  Doc: returns the absolute trace of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nftrace(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltval(*argv):
  '''
  nfeltval
  Class: basic
  Section: number_fields
  C-Name: nfval
  Prototype: lGGG
  Help: nfeltval(nf,x,pr): valuation of element x at the prime pr as output by
   idealprimedec.
  Doc: given an element $x$ in
   \var{nf} and a prime ideal \var{pr} in the format output by
   \kbd{idealprimedec}, computes their the valuation at \var{pr} of the
   element $x$. The same result could be obtained using
   \kbd{idealval(\var{nf},x,\var{pr})} (since $x$ would then be converted to a
   principal ideal), but it would be less efficient.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.nfval(*c_arg_tuple)

def rnfalgtobasis(*argv):
  '''
  rnfalgtobasis
  Class: basic
  Section: number_fields
  C-Name: rnfalgtobasis
  Prototype: GG
  Help: rnfalgtobasis(rnf,x): relative version of nfalgtobasis, where rnf is a
   relative numberfield.
  Doc: expresses $x$ on the relative
   integral basis. Here, $\var{rnf}$ is a relative number field extension $L/K$
   as output by \kbd{rnfinit}, and $x$ an element of $L$ in absolute form, i.e.
   expressed as a polynomial or polmod with polmod coefficients, \emph{not} on
   the relative integral basis.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfalgtobasis(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfbasistoalg(*argv):
  '''
  rnfbasistoalg
  Class: basic
  Section: number_fields
  C-Name: rnfbasistoalg
  Prototype: GG
  Help: rnfbasistoalg(rnf,x): relative version of nfbasistoalg, where rnf is a
   relative numberfield.
  Doc: computes the representation of $x$
   as a polmod with polmods coefficients. Here, $\var{rnf}$ is a relative number
   field extension $L/K$ as output by \kbd{rnfinit}, and $x$ an element of
   $L$ expressed on the relative integral basis.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfbasistoalg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfeltnorm(*argv):
  '''
  rnfeltnorm
  Class: basic
  Section: number_fields
  C-Name: rnfeltnorm
  Prototype: GG
  Help: rnfeltnorm(rnf,x): returns the relative norm N_{L/K}(x), as an element
   of K
  Doc: $\var{rnf}$ being a relative number field extension $L/K$ as output by
   \kbd{rnfinit} and $x$ being an element of $L$, returns the relative norm
   $N_{L/K}(x)$ as an element of $K$.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? rnfeltnorm(L, Mod(x, L.pol))
   %2 = Mod(x, x^2 + Mod(-y, y^2 + 1))
   ? rnfeltnorm(L, 2)
   %3 = 4
   ? rnfeltnorm(L, Mod(x, x^2-y))
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfeltnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfelttrace(*argv):
  '''
  rnfelttrace
  Class: basic
  Section: number_fields
  C-Name: rnfelttrace
  Prototype: GG
  Help: rnfelttrace(rnf,x): returns the relative trace N_{L/K}(x), as an element
   of K
  Doc: $\var{rnf}$ being a relative number field extension $L/K$ as output by
   \kbd{rnfinit} and $x$ being an element of $L$, returns the relative trace
   $N_{L/K}(x)$ as an element of $K$.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? rnfelttrace(L, Mod(x, L.pol))
   %2 = 0
   ? rnfelttrace(L, 2)
   %3 = 4
   ? rnfelttrace(L, Mod(x, x^2-y))
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfelttrace(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ideallog(*argv):
  '''
  ideallog
  Class: basic
  Section: number_fields
  C-Name: ideallog
  Prototype: GGG
  Help: ideallog(nf,x,bid): if bid is a big ideal, as given by
   idealstar(nf,I,1) or idealstar(nf,I,2), gives the vector of exponents on the
   generators bid[2][3] (even if these generators have not been computed).
  Doc: $\var{nf}$ is a number field,
   \var{bid} is as output by \kbd{idealstar(nf, D, \dots)} and $x$ a
   non-necessarily integral element of \var{nf} which must have valuation
   equal to 0 at all prime ideals in the support of $\kbd{D}$. This function
   computes the discrete logarithm of $x$ on the generators given in
   \kbd{\var{bid}.gen}. In other words, if $g_i$ are these generators, of orders
   $d_i$ respectively, the result is a column vector of integers $(x_i)$ such
   that $0\le x_i<d_i$ and
   $$x \equiv \prod_i g_i^{x_i} \pmod{\ ^*D}\enspace.$$
   Note that when the support of \kbd{D} contains places at infinity, this
   congruence implies also sign conditions on the associated real embeddings.
   See \tet{znlog} for the limitations of the underlying discrete log algorithms.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ideallog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealtwoelt(*argv):
  '''
  idealtwoelt
  Class: basic
  Section: number_fields
  C-Name: idealtwoelt0
  Prototype: GGDG
  Help: idealtwoelt(nf,x,{a}): two-element representation of an ideal x in the
   number field nf. If (optional) a is non-zero, first element will be equal to a.
  Doc: computes a two-element
   representation of the ideal $x$ in the number field $\var{nf}$, combining a
   random search and an approximation theorem; $x$ is an ideal
   in any form (possibly an extended ideal, whose principal part is ignored)
   
   \item When called as \kbd{idealtwoelt(nf,x)}, the result is a row vector
   $[a,\alpha]$ with two components such that $x=a\Z_K+\alpha\Z_K$ and $a$ is
   chosen to be the positive generator of $x\cap\Z$, unless $x$ was given as a
   principal ideal (in which case we may choose $a = 0$). The algorithm
   uses a fast lazy factorization of $x\cap \Z$ and runs in randomized
   polynomial time.
   
   \item When called as \kbd{idealtwoelt(nf,x,a)} with an explicit non-zero $a$
   supplied as third argument, the function assumes that $a \in x$ and returns
   $\alpha\in x$ such that $x = a\Z_K + \alpha\Z_K$. Note that we must factor
   $a$ in this case, and the algorithm is generally much slower than the
   default variant.
  Variant: Also available are
   \fun{GEN}{idealtwoelt}{GEN nf, GEN x} and
   \fun{GEN}{idealtwoelt2}{GEN nf, GEN x, GEN a}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealtwoelt0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealadd(*argv):
  '''
  idealadd
  Class: basic
  Section: number_fields
  C-Name: idealadd
  Prototype: GGG
  Help: idealadd(nf,x,y): sum of two ideals x and y in the number field
   defined by nf.
  Doc: sum of the two ideals $x$ and $y$ in the number field $\var{nf}$. The
   result is given in HNF.
   \bprog
    ? K = nfinit(x^2 + 1);
    ? a = idealadd(K, 2, x + 1)  \\ ideal generated by 2 and 1+I
    %2 =
    [2 1]
   
    [0 1]
    ? pr = idealprimedec(K, 5)[1];  \\ a prime ideal above 5
    ? idealadd(K, a, pr)     \\ coprime, as expected
    %4 =
    [1 0]
   
    [0 1]
   @eprog\noindent
   This function cannot be used to add arbitrary $\Z$-modules, since it assumes
   that its arguments are ideals:
   \bprog
     ? b = Mat([1,0]~);
     ? idealadd(K, b, b)     \\ only square t_MATs represent ideals
     *** idealadd: non-square t_MAT in idealtyp.
     ? c = [2, 0; 2, 0]; idealadd(K, c, c)   \\ non-sense
     %6 =
     [2 0]
   
     [0 2]
     ? d = [1, 0; 0, 2]; idealadd(K, d, d)   \\ non-sense
     %7 =
     [1 0]
   
     [0 1]
   
   @eprog\noindent In the last two examples, we get wrong results since the
   matrices $c$ and $d$ do not correspond to an ideal: the $\Z$-span of their
   columns (as usual interpreted as coordinates with respect to the integer basis
   \kbd{K.zk}) is not an $O_K$-module. To add arbitrary $\Z$-modules generated
   by the columns of matrices $A$ and $B$, use \kbd{mathnf(concat(A,B))}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealadd(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealaddtoone(*argv):
  '''
  idealaddtoone
  Class: basic
  Section: number_fields
  C-Name: idealaddtoone0
  Prototype: GGDG
  Help: idealaddtoone(nf,x,{y}): if y is omitted, when the sum of the ideals
   in the number field K defined by nf and given in the vector x is equal to
   Z_K, gives a vector of elements of the corresponding ideals who sum to 1.
   Otherwise, x and y are ideals, and if they sum up to 1, find one element in
   each of them such that the sum is 1.
  Doc: $x$ and $y$ being two co-prime
   integral ideals (given in any form), this gives a two-component row vector
   $[a,b]$ such that $a\in x$, $b\in y$ and $a+b=1$.
   
   The alternative syntax $\kbd{idealaddtoone}(\var{nf},v)$, is supported, where
   $v$ is a $k$-component vector of ideals (given in any form) which sum to
   $\Z_K$. This outputs a $k$-component vector $e$ such that $e[i]\in x[i]$ for
   $1\le i\le k$ and $\sum_{1\le i\le k}e[i]=1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealaddtoone0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealappr(*argv):
  '''
  idealappr
  Class: basic
  Section: number_fields
  C-Name: idealappr0
  Prototype: GGD0,L,
  Help: idealappr(nf,x,{flag=0}): x being a fractional ideal, gives an element
   b such that v_p(b)=v_p(x) for all prime ideals p dividing x, and v_p(b)>=0
   for all other p. If (optional) flag is non-null x must be a prime ideal
   factorization with possibly zero exponents.
  Doc: if $x$ is a fractional ideal
   (given in any form), gives an element $\alpha$ in $\var{nf}$ such that for
   all prime ideals $\goth{p}$ such that the valuation of $x$ at $\goth{p}$ is
   non-zero, we have $v_{\goth{p}}(\alpha)=v_{\goth{p}}(x)$, and
   $v_{\goth{p}}(\alpha)\ge0$ for all other $\goth{p}$.
   
   If $\fl$ is non-zero, $x$ must be given as a prime ideal factorization, as
   output by \kbd{idealfactor}, but possibly with zero or negative exponents.
   This yields an element $\alpha$ such that for all prime ideals $\goth{p}$
   occurring in $x$, $v_{\goth{p}}(\alpha)$ is equal to the exponent of
   $\goth{p}$ in $x$, and for all other prime ideals,
   $v_{\goth{p}}(\alpha)\ge0$. This generalizes $\kbd{idealappr}(\var{nf},x,0)$
   since zero exponents are allowed. Note that the algorithm used is slightly
   different, so that
   \bprog
     idealappr(nf, idealfactor(nf,x))
   @eprog\noindent
   may not be the same as \kbd{idealappr(nf,x,1)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealappr0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealchinese(*argv):
  '''
  idealchinese
  Class: basic
  Section: number_fields
  C-Name: idealchinese
  Prototype: GGG
  Help: idealchinese(nf,x,y): x being a prime ideal factorization and y a
   vector of elements, gives an element b such that v_p(b-y_p)>=v_p(x) for all
   prime ideals p dividing x, and v_p(b)>=0 for all other p.
  Doc: $x$ being a prime ideal factorization
   (i.e.~a 2 by 2 matrix whose first column contains prime ideals, and the second
   column integral exponents), $y$ a vector of elements in $\var{nf}$ indexed by
   the ideals in $x$, computes an element $b$ such that
   
   $v_{\goth{p}}(b - y_{\goth{p}}) \geq v_{\goth{p}}(x)$ for all prime ideals
   in $x$ and $v_{\goth{p}}(b)\geq 0$ for all other $\goth{p}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealchinese(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealcoprime(*argv):
  '''
  idealcoprime
  Class: basic
  Section: number_fields
  C-Name: idealcoprime
  Prototype: GGG
  Help: idealcoprime(nf,x,y): gives an element b in nf such that b. x is an
   integral ideal coprime to the integral ideal y.
  Doc: given two integral ideals $x$ and $y$
   in the number field $\var{nf}$, returns a $\beta$ in the field,
   such that $\beta\cdot x$ is an integral ideal coprime to $y$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealcoprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealdiv(*argv):
  '''
  idealdiv
  Class: basic
  Section: number_fields
  C-Name: idealdiv0
  Prototype: GGGD0,L,
  Help: idealdiv(nf,x,y,{flag=0}): quotient x/y of two ideals x and y in HNF
   in the number field nf. If (optional) flag is non-null, the quotient is
   supposed to be an integral ideal (slightly faster).
  Description: 
   (gen, gen, gen, ?0):gen        idealdiv($1, $2, $3)
   (gen, gen, gen, 1):gen         idealdivexact($1, $2, $3)
   (gen, gen, gen, #small):gen    $"invalid flag in idealdiv"
   (gen, gen, gen, small):gen     idealdiv0($1, $2, $3, $4)
  Doc: quotient $x\cdot y^{-1}$ of the two ideals $x$ and $y$ in the number
   field $\var{nf}$. The result is given in HNF.
   
   If $\fl$ is non-zero, the quotient $x \cdot y^{-1}$ is assumed to be an
   integral ideal. This can be much faster when the norm of the quotient is
   small even though the norms of $x$ and $y$ are large.
  Variant: Also available are \fun{GEN}{idealdiv}{GEN nf, GEN x, GEN y}
   ($\fl=0$) and \fun{GEN}{idealdivexact}{GEN nf, GEN x, GEN y} ($\fl=1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealdiv0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealfactor(*argv):
  '''
  idealfactor
  Class: basic
  Section: number_fields
  C-Name: idealfactor
  Prototype: GG
  Help: idealfactor(nf,x): factorization of the ideal x given in HNF into
   prime ideals in the number field nf.
  Doc: factors into prime ideal powers the
   ideal $x$ in the number field $\var{nf}$. The output format is similar to the
   \kbd{factor} function, and the prime ideals are represented in the form
   output by the \kbd{idealprimedec} function, i.e.~as 5-element vectors.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealfactor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealhnf(*argv):
  '''
  idealhnf
  Class: basic
  Section: number_fields
  C-Name: idealhnf0
  Prototype: GGDG
  Help: idealhnf(nf,u,{v}): hermite normal form of the ideal u in the number
   field nf if v is omitted. If called as idealhnf(nf,u,v), the ideal
   is given as uZ_K + vZ_K in the number field K defined by nf.
  Doc: gives the \idx{Hermite normal form} of the ideal $u\Z_K+v\Z_K$, where $u$
   and $v$ are elements of the number field $K$ defined by \kbd{nf}.
   \bprog
   ? nf = nfinit(y^3 - 2);
   ? idealhnf(nf, 2, y+1)
   %2 =
   [1 0 0]
   
   [0 1 0]
   
   [0 0 1]
   ? idealhnf(nf, y/2, [0,0,1/3]~)
   %3 =
   [1/3 0 0]
   
   [0 1/6 0]
   
   [0 0 1/6]
   @eprog
   
   If $b$ is omitted, returns the HNF of the ideal defined by $u$: $u$ may be an
   algebraic number (defining a principal ideal), a maximal ideal (as given by
   \kbd{idealprimedec} or \kbd{idealfactor}), or a matrix whose columns give
   generators for the ideal. This last format is a little complicated, but
   useful to reduce general modules to the canonical form once in a while:
   
   \item if strictly less than $N = [K:\Q]$ generators are given, $u$
   is the $\Z_K$-module they generate,
   
   \item if $N$ or more are given, it is \emph{assumed} that they form a
   $\Z$-basis of the ideal, in particular that the matrix has maximal rank $N$.
   This acts as \kbd{mathnf} since the $\Z_K$-module structure is (taken for
   granted hence) not taken into account in this case.
   \bprog
   ? idealhnf(nf, idealprimedec(nf,2)[1])
   %4 =
   [2 0 0]
   
   [0 1 0]
   
   [0 0 1]
   ? idealhnf(nf, [1,2;2,3;3,4])
   %5 =
   [1 0 0]
   
   [0 1 0]
   
   [0 0 1]
   @eprog\noindent Finally, when $K$ is quadratic with discriminant $D_K$, we
   allow $u =$ \kbd{Qfb(a,b,c)}, provided $b^2 - 4ac = D_K$. As usual,
   this represents the ideal $a \Z + (1/2)(-b + \sqrt{D_K}) \Z$.
   \bprog
   ? K = nfinit(x^2 - 60); K.disc
   %1 = 60
   ? idealhnf(K, qfbprimeform(60,2))
   %2 =
   [2 1]
   
   [0 1]
   ? idealhnf(K, Qfb(1,2,3))
     ***   at top-level: idealhnf(K,Qfb(1,2,3
     ***                 ^--------------------
     *** idealhnf: Qfb(1, 2, 3) has discriminant != 60 in idealhnf.
   @eprog
  Variant: Also available is \fun{GEN}{idealhnf}{GEN nf, GEN a}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealhnf0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealintersect(*argv):
  '''
  idealintersect
  Class: basic
  Section: number_fields
  C-Name: idealintersect
  Prototype: GGG
  Help: idealintersect(nf,A,B): intersection of two ideals A and B in the
   number field defined by nf.
  Doc: intersection of the two ideals
   $A$ and $B$ in the number field $\var{nf}$. The result is given in HNF.
   \bprog
   ? nf = nfinit(x^2+1);
   ? idealintersect(nf, 2, x+1)
   %2 =
   [2 0]
   
   [0 2]
   @eprog
   
   This function does not apply to general $\Z$-modules, e.g.~orders, since its
   arguments are replaced by the ideals they generate. The following script
   intersects $\Z$-modules $A$ and $B$ given by matrices of compatible
   dimensions with integer coefficients:
   \bprog
   ZM_intersect(A,B) =
   { my(Ker = matkerint(concat(A,B)));
     mathnf( A * Ker[1..#A,] )
   }
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealintersect(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealinv(*argv):
  '''
  idealinv
  Class: basic
  Section: number_fields
  C-Name: idealinv
  Prototype: GG
  Help: idealinv(nf,x): inverse of the ideal x in the number field nf.
  Description: 
   (gen, gen):gen        idealinv($1, $2)
  Doc: inverse of the ideal $x$ in the
   number field $\var{nf}$, given in HNF. If $x$ is an extended
   ideal\sidx{ideal (extended)}, its principal part is suitably
   updated: i.e. inverting $[I,t]$, yields $[I^{-1}, 1/t]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealinv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealred(*argv):
  '''
  idealred
  Class: basic
  Section: number_fields
  C-Name: idealred0
  Prototype: GGDG
  Help: idealred(nf,I,{v=0}): LLL reduction of the ideal I in the number
   field nf along direction v, in HNF.
  Doc: \idx{LLL} reduction of
   the ideal $I$ in the number field \var{nf}, along the direction $v$.
   The $v$ parameter is best left omitted, but if it is present, it must
   be an $\kbd{nf.r1} + \kbd{nf.r2}$-component vector of \emph{non-negative}
   integers. (What counts is the relative magnitude of the entries: if all
   entries are equal, the effect is the same as if the vector had been omitted.)
   
   This function finds a ``small'' $a$ in $I$ (see the end for technical details).
   The result is the Hermite normal form of
   the ``reduced'' ideal $J = r I/a$, where $r$ is the unique rational number such
   that $J$ is integral and primitive. (This is usually not a reduced ideal in
   the sense of \idx{Buchmann}.)
   \bprog
   ? K = nfinit(y^2+1);
   ? P = idealprimedec(K,5)[1];
   ? idealred(K, P)
   %3 =
   [1 0]
   
   [0 1]
   @eprog\noindent More often than not, a \idx{principal ideal} yields the unit
   ideal as above. This is a quick and dirty way to check if ideals are principal,
   but it is not a necessary condition: a non-trivial result does not prove that
   the ideal is non-principal. For guaranteed results, see \kbd{bnfisprincipal},
   which requires the computation of a full \kbd{bnf} structure.
   
   If the input is an extended ideal $[I,s]$, the output is $[J,sa/r]$; this way,
   one can keep track of the principal ideal part:
   \bprog
   ? idealred(K, [P, 1])
   %5 = [[1, 0; 0, 1], [-2, 1]~]
   @eprog\noindent
   meaning that $P$ is generated by $[-2, 1]~$. The number field element in the
   extended part is an algebraic number in any form \emph{or} a factorization
   matrix (in terms of number field elements, not ideals!). In the latter case,
   elements stay in factored form, which is a convenient way to avoid
   coefficient explosion; see also \tet{idealpow}.
   
   \misctitle{Technical note} The routine computes an LLL-reduced
   basis for the lattice $I$ equipped with the quadratic form
   $$|| x ||_v^2 = \sum_{i=1}^{r_1+r_2} 2^{v_i}\varepsilon_i|\sigma_i(x)|^2,$$
   where as usual the $\sigma_i$ are the (real and) complex embeddings and
   $\varepsilon_i = 1$, resp.~$2$, for a real, resp.~complex place. The element
   $a$ is simply the first vector in the LLL basis. The only reason you may want
   to try to change some directions and set some $v_i\neq 0$ is to randomize
   the elements found for a fixed ideal, which is heuristically useful in index
   calculus algorithms like \tet{bnfinit} and \tet{bnfisprincipal}.
   
   \misctitle{Even more technical note} In fact, the above is a white lie.
   We do not use $||\cdot||_v$ exactly but a rescaled rounded variant which
   gets us faster and simpler LLLs. There's no harm since we are not using any
   theoretical property of $a$ after all, except that it belongs to $I$ and is
   ``expected to be small''.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealred0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealmul(*argv):
  '''
  idealmul
  Class: basic
  Section: number_fields
  C-Name: idealmul0
  Prototype: GGGD0,L,
  Help: idealmul(nf,x,y,{flag=0}): product of the two ideals x and y in the
   number field nf. If (optional) flag is non-nul, reduce the result.
  Description: 
   (gen, gen, gen, ?0):gen        idealmul($1, $2, $3)
   (gen, gen, gen, 1):gen         idealmulred($1, $2, $3)
   (gen, gen, gen, #small):gen    $"invalid flag in idealmul"
   (gen, gen, gen, small):gen     idealmul0($1, $2, $3, $4)
  Doc: ideal multiplication of the ideals $x$ and $y$ in the number field
   \var{nf}; the result is the ideal product in HNF. If either $x$ or $y$
   are extended ideals\sidx{ideal (extended)}, their principal part is suitably
   updated: i.e. multiplying $[I,t]$, $[J,u]$ yields $[IJ, tu]$; multiplying
   $I$ and $[J, u]$ yields $[IJ, u]$.
   \bprog
   ? nf = nfinit(x^2 + 1);
   ? idealmul(nf, 2, x+1)
   %2 =
   [4 2]
   
   [0 2]
   ? idealmul(nf, [2, x], x+1)        \\ extended ideal * ideal
   %4 = [[4, 2; 0, 2], x]
   ? idealmul(nf, [2, x], [x+1, x])   \\ two extended ideals
   %5 = [[4, 2; 0, 2], [-1, 0]~]
   @eprog\noindent
   If $\fl$ is non-zero, reduce the result using \kbd{idealred}.
  Variant: 
   \noindent See also
   \fun{GEN}{idealmul}{GEN nf, GEN x, GEN y} ($\fl=0$) and
   \fun{GEN}{idealmulred}{GEN nf, GEN x, GEN y} ($\fl\neq0$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealmul0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealnorm(*argv):
  '''
  idealnorm
  Class: basic
  Section: number_fields
  C-Name: idealnorm
  Prototype: GG
  Help: idealnorm(nf,x): norm of the ideal x in the number field nf.
  Doc: computes the norm of the ideal~$x$ in the number field~$\var{nf}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealnumden(*argv):
  '''
  idealnumden
  Class: basic
  Section: number_fields
  C-Name: idealnumden
  Prototype: GG
  Help: idealnumden(nf,x): returns [A,B], where A,B are coprime integer ideals
   such that x = A/B
  Doc: returns $[A,B]$, where $A,B$ are coprime integer ideals
   such that $x = A/B$, in the number field $\var{nf}$.
   \bprog
   ? nf = nfinit(x^2+1);
   ? idealnumden(nf, (x+1)/2)
   %2 = [[1, 0; 0, 1], [2, 1; 0, 1]]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealnumden(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealpow(*argv):
  '''
  idealpow
  Class: basic
  Section: number_fields
  C-Name: idealpow0
  Prototype: GGGD0,L,
  Help: idealpow(nf,x,k,{flag=0}): k-th power of the ideal x in HNF in the
   number field nf. If (optional) flag is non-null, reduce the result.
  Doc: computes the $k$-th power of
   the ideal $x$ in the number field $\var{nf}$; $k\in\Z$.
   If $x$ is an extended
   ideal\sidx{ideal (extended)}, its principal part is suitably
   updated: i.e. raising $[I,t]$ to the $k$-th power, yields $[I^k, t^k]$.
   
   If $\fl$ is non-zero, reduce the result using \kbd{idealred}, \emph{throughout
   the (binary) powering process}; in particular, this is \emph{not} the same as
   as $\kbd{idealpow}(\var{nf},x,k)$ followed by reduction.
  Variant: 
   \noindent See also
   \fun{GEN}{idealpow}{GEN nf, GEN x, GEN k} and
   \fun{GEN}{idealpows}{GEN nf, GEN x, long k} ($\fl = 0$).
   Corresponding to $\fl=1$ is \fun{GEN}{idealpowred}{GEN nf, GEN vp, GEN k}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealpow0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealval(*argv):
  '''
  idealval
  Class: basic
  Section: number_fields
  C-Name: idealval
  Prototype: lGGG
  Help: idealval(nf,x,pr): valuation at pr given in idealprimedec format of the
   ideal x in the number field nf.
  Doc: gives the valuation of the ideal $x$ at the prime ideal \var{pr} in the
   number field $\var{nf}$, where \var{pr} is in \kbd{idealprimedec} format.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.idealval(*c_arg_tuple)

def nfisideal(*argv):
  '''
  nfisideal
  Class: basic
  Section: number_fields
  C-Name: isideal
  Prototype: lGG
  Help: nfisideal(nf,x): true(1) if x is an ideal in the number field nf,
   false(0) if not.
  Doc: returns 1 if $x$ is an ideal in the number field $\var{nf}$, 0 otherwise.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.isideal(*c_arg_tuple)

def idealmin(*argv):
  '''
  idealmin
  Class: basic
  Section: number_fields
  C-Name: idealmin
  Prototype: GGDG
  Help: idealmin(nf,ix,{vdir}): pseudo-minimum of the ideal ix in the direction
   vdir in the number field nf.
  Doc: \emph{This function is useless and kept for backward compatibility only,
   use \kbd{idealred}}. Computes a pseudo-minimum of the ideal $x$ in the
   direction \var{vdir} in the number field \var{nf}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealmin(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfdetint(*argv):
  '''
  nfdetint
  Class: basic
  Section: number_fields
  C-Name: nfdetint
  Prototype: GG
  Help: nfdetint(nf,x): multiple of the ideal determinant of the pseudo
   generating set x.
  Doc: given a pseudo-matrix $x$, computes a
   non-zero ideal contained in (i.e.~multiple of) the determinant of $x$. This
   is particularly useful in conjunction with \kbd{nfhnfmod}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfdetint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltdivmodpr(*argv):
  '''
  nfeltdivmodpr
  Class: basic
  Section: number_fields
  C-Name: nfdivmodpr
  Prototype: GGGG
  Help: nfeltdivmodpr(nf,x,y,pr): element x/y modulo pr in nf, where pr is in
   modpr format (see nfmodprinit).
  Doc: given two elements $x$
   and $y$ in \var{nf} and \var{pr} a prime ideal in \kbd{modpr} format (see
   \tet{nfmodprinit}), computes their quotient $x / y$ modulo the prime ideal
   \var{pr}.
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nf\_to\_Fq}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfdivmodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfhnf(*argv):
  '''
  nfhnf
  Class: basic
  Section: number_fields
  C-Name: nfhnf
  Prototype: GG
  Help: nfhnf(nf,x): if x=[A,I], gives a pseudo-basis of the module sum A_jI_j
  Doc: given a pseudo-matrix $(A,I)$, finds a
   pseudo-basis in \idx{Hermite normal form} of the module it generates.
  Variant: Also available:
   
   \fun{GEN}{rnfsimplifybasis}{GEN bnf, GEN x} simplifies the pseudo-basis
   given by $x = (A,I)$. The ideals in the list $I$ are integral, primitive and
   either trivial (equal to the full ring of integer) or non-principal.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfhnf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfhnfmod(*argv):
  '''
  nfhnfmod
  Class: basic
  Section: number_fields
  C-Name: nfhnfmod
  Prototype: GGG
  Help: nfhnfmod(nf,x,detx): if x=[A,I], and detx is a multiple of the ideal
   determinant of x, gives a pseudo-basis of the module sum A_jI_j.
  Doc: given a pseudo-matrix $(A,I)$
   and an ideal \var{detx} which is contained in (read integral multiple of) the
   determinant of $(A,I)$, finds a pseudo-basis in \idx{Hermite normal form}
   of the module generated by $(A,I)$. This avoids coefficient explosion.
   \var{detx} can be computed using the function \kbd{nfdetint}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfhnfmod(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfkermodpr(*argv):
  '''
  nfkermodpr
  Class: basic
  Section: number_fields
  C-Name: nfkermodpr
  Prototype: GGG
  Help: nfkermodpr(nf,x,pr): kernel of the matrix x in Z_K/pr, where pr is in
   modpr format (see nfmodprinit).
  Doc: kernel of the matrix $a$ in $\Z_K/\var{pr}$, where \var{pr} is in
   \key{modpr} format (see \kbd{nfmodprinit}).
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nfM\_to\_FqM}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfkermodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltmulmodpr(*argv):
  '''
  nfeltmulmodpr
  Class: basic
  Section: number_fields
  C-Name: nfmulmodpr
  Prototype: GGGG
  Help: nfeltmulmodpr(nf,x,y,pr): element x.y modulo pr in nf, where pr is in
   modpr format (see nfmodprinit).
  Doc: given two elements $x$ and
   $y$ in \var{nf} and \var{pr} a prime ideal in \kbd{modpr} format (see
   \tet{nfmodprinit}), computes their product $x*y$ modulo the prime ideal
   \var{pr}.
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nf\_to\_Fq}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfmulmodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltpowmodpr(*argv):
  '''
  nfeltpowmodpr
  Class: basic
  Section: number_fields
  C-Name: nfpowmodpr
  Prototype: GGGG
  Help: nfeltpowmodpr(nf,x,k,pr): element x^k modulo pr in nf, where pr is in
   modpr format (see nfmodprinit).
  Doc: given an element $x$ in \var{nf}, an integer $k$ and a prime ideal
   \var{pr} in \kbd{modpr} format
   (see \tet{nfmodprinit}), computes $x^k$ modulo the prime ideal \var{pr}.
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nf\_to\_Fq}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfpowmodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfeltreduce(*argv):
  '''
  nfeltreduce
  Class: basic
  Section: number_fields
  C-Name: nfreduce
  Prototype: GGG
  Help: nfeltreduce(nf,a,id): gives r such that a-r is in the ideal id and r
   is small.
  Doc: given an ideal \var{id} in
   Hermite normal form and an element $a$ of the number field $\var{nf}$,
   finds an element $r$ in $\var{nf}$ such that $a-r$ belongs to the ideal
   and $r$ is small.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfreduce(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfsnf(*argv):
  '''
  nfsnf
  Class: basic
  Section: number_fields
  C-Name: nfsnf
  Prototype: GG
  Help: nfsnf(nf,x): if x=[A,I,J], outputs [c_1,...c_n] Smith normal form of x.
  Doc: given a $\Z_K$-module $x$ associated to the integral pseudo-matrix
   $(A,I,J)$, returns an ideal list $d_1,\dots,d_n$ which is the \idx{Smith
   normal form} of $x$. In other words, $x$ is isomorphic to
   $\Z_K/d_1\oplus\cdots\oplus\Z_K/d_n$ and $d_i$ divides $d_{i-1}$ for $i\ge2$.
   
   See \secref{se:ZKmodules} for the definition of integral pseudo-matrix;
   briefly, it is input as a 3-component row vector $[A,I,J]$ where
   $I = [b_1,\dots,b_n]$ and $J = [a_1,\dots,a_n]$ are two ideal lists,
   and $A$ is a square $n\times n$ matrix with columns $(A_1,\dots,A_n)$,
   seen as elements in $K^n$ (with canonical basis $(e_1,\dots,e_n)$).
   This data defines the $\Z_K$ module $x$ given by
   $$ (b_1e_1\oplus\cdots\oplus b_ne_n) / (a_1A_1\oplus\cdots\oplus a_nA_n)
   \enspace, $$
   The integrality condition is $a_{i,j} \in b_i a_j^{-1}$ for all $i,j$. If it
   is not satisfied, then the $d_i$ will not be integral. Note that every
   finitely generated torsion module is isomorphic to a module of this form and
   even with $b_i=Z_K$ for all $i$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfsnf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfsolvemodpr(*argv):
  '''
  nfsolvemodpr
  Class: basic
  Section: number_fields
  C-Name: nfsolvemodpr
  Prototype: GGGG
  Help: nfsolvemodpr(nf,a,b,P): solution of a*x=b in Z_K/P, where a is a
   matrix and b a column vector, and where P is in modpr format (see
   nfmodprinit).
  Doc: let $P$ be a prime ideal in \key{modpr} format (see \kbd{nfmodprinit}),
   let $a$ be a matrix, invertible over the residue field, and let $b$ be
   a column vector or matrix. This function returns a solution of $a\cdot x =
   b$; the coefficients of $x$ are lifted to \var{nf} elements.
   \bprog
   ? K = nfinit(y^2+1);
   ? P = idealprimedec(K, 3)[1];
   ? P = nfmodprinit(K, P);
   ? a = [y+1, y; y, 0]; b = [1, y]~
   ? nfsolvemodpr(K, a,b, P)
   %5 = [1, 2]~
   @eprog
  Variant: This function is normally useless in library mode. Project your
   inputs to the residue field using \kbd{nfM\_to\_FqM}, then work there.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfsolvemodpr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfeltabstorel(*argv):
  '''
  rnfeltabstorel
  Class: basic
  Section: number_fields
  C-Name: rnfeltabstorel
  Prototype: GG
  Help: rnfeltabstorel(rnf,x): transforms the element x from absolute to
   relative representation.
  Doc: $\var{rnf}$ being a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and $x$ being an
   element of $L$ expressed as a polynomial modulo the absolute equation
   \kbd{\var{rnf}.pol}, computes $x$ as an element of the relative extension
   $L/K$ as a polmod with polmod coefficients.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? L.pol
   %2 = x^4 + 1
   ? rnfeltabstorel(L, Mod(x, L.pol))
   %3 = Mod(x, x^2 + Mod(-y, y^2 + 1))
   ? rnfeltabstorel(L, Mod(2, L.pol))
   %4 = 2
   ? rnfeltabstorel(L, Mod(x, x^2-y))
    ***   at top-level: rnfeltabstorel(L,Mod
    ***                 ^--------------------
    *** rnfeltabstorel: inconsistent moduli in rnfeltabstorel: x^2-y != x^4+1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfeltabstorel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfeltdown(*argv):
  '''
  rnfeltdown
  Class: basic
  Section: number_fields
  C-Name: rnfeltdown
  Prototype: GG
  Help: rnfeltdown(rnf,x): expresses x on the base field if possible; returns
   an error otherwise.
  Doc: $\var{rnf}$ being a relative number
   field extension $L/K$ as output by \kbd{rnfinit} and $x$ being an element of
   $L$ expressed as a polynomial or polmod with polmod coefficients, computes
   $x$ as an element of $K$ as a polmod, assuming $x$ is in $K$ (otherwise a
   domain error occurs).
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? L.pol
   %2 = x^4 + 1
   ? rnfeltdown(L, Mod(x^2, L.pol))
   %3 = Mod(y, y^2 + 1)
   ? rnfeltdown(L, Mod(y, x^2-y))
   %4 = Mod(y, y^2 + 1)
   ? rnfeltdown(L, Mod(y,K.pol))
   %5 = Mod(y, y^2 + 1)
   ? rnfeltdown(L, Mod(x, L.pol))
    ***   at top-level: rnfeltdown(L,Mod(x,x
    ***                 ^--------------------
    *** rnfeltdown: domain error in rnfeltdown: element not in the base field
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfeltdown(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfeltreltoabs(*argv):
  '''
  rnfeltreltoabs
  Class: basic
  Section: number_fields
  C-Name: rnfeltreltoabs
  Prototype: GG
  Help: rnfeltreltoabs(rnf,x): transforms the element x from relative to
   absolute representation.
  Doc: $\var{rnf}$ being a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and $x$ being an
   element of $L$ expressed as a polynomial or polmod with polmod
   coefficients, computes $x$ as an element of the absolute extension $L/\Q$ as
   a polynomial modulo the absolute equation \kbd{\var{rnf}.pol}.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? L.pol
   %2 = x^4 + 1
   ? rnfeltreltoabs(L, Mod(x, L.pol))
   %3 = Mod(x, x^4 + 1)
   ? rnfeltreltoabs(L, Mod(y, x^2-y))
   %4 = Mod(x^2, x^4 + 1)
   ? rnfeltreltoabs(L, Mod(y,K.pol))
   %5 = Mod(x^2, x^4 + 1)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfeltreltoabs(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfeltup(*argv):
  '''
  rnfeltup
  Class: basic
  Section: number_fields
  C-Name: rnfeltup
  Prototype: GG
  Help: rnfeltup(rnf,x): expresses x (belonging to the base field) on the
   relative field.
  Doc: $\var{rnf}$ being a relative number field extension $L/K$ as output by
   \kbd{rnfinit} and $x$ being an element of $K$, computes $x$ as an element of
   the absolute extension $L/\Q$ as a polynomial modulo the absolute equation
   \kbd{\var{rnf}.pol}.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y);
   ? L.pol
   %2 = x^4 + 1
   ? rnfeltup(L, Mod(y, K.pol))
   %4 = Mod(x^2, x^4 + 1)
   ? rnfeltup(L, y)
   %5 = Mod(x^2, x^4 + 1)
   ? rnfeltup(L, [1,2]~) \\ in terms of K.zk
   %6 = Mod(2*x^2 + 1, x^4 + 1)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfeltup(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealabstorel(*argv):
  '''
  rnfidealabstorel
  Class: basic
  Section: number_fields
  C-Name: rnfidealabstorel
  Prototype: GG
  Help: rnfidealabstorel(rnf,x): transforms the ideal x from absolute to
   relative representation.
  Doc: let $\var{rnf}$ be a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and $x$ be an ideal of
   the absolute extension $L/\Q$ given by a $\Z$-basis of elements of $L$.
   Returns the relative pseudo-matrix in HNF giving the ideal $x$ considered as
   an ideal of the relative extension $L/K$, i.e.~as a $\Z_K$-module.
   
   The reason why the input does not use the customary HNF in terms of a fixed
   $\Z$-basis for $\Z_L$ is precisely that no such basis has been explicitly
   specified. On the other hand, if you already computed an (absolute) \var{nf}
   structure \kbd{Labs} associated to $L$, and $m$ is in HNF, defining
   an (absolute) ideal with respect to the $\Z$-basis \kbd{Labs.zk}, then
   \kbd{Labs.zk * m} is a suitable $\Z$-basis for the ideal, and
   \bprog
     rnfidealabstorel(rnf, Labs.zk * m)
   @eprog\noindent converts $m$ to a relative ideal.
   \bprog
   ? K = nfinit(y^2+1); L = rnfinit(K, x^2-y); Labs = nfinit(L.pol);
   ? m = idealhnf(Labs, 17, x^3+2);
   ? B = rnfidealabstorel(L, Labs.zk * m)
   %3 = [[1, 8; 0, 1], [[17, 4; 0, 1], 1]]  \\ pseudo-basis for m as Z_K-module
   ? A = rnfidealreltoabs(L, B)
   %4 = [17, x^2 + 4, x + 8, x^3 + 8*x^2]   \\ Z-basis for m in Q[x]/(L.pol)
   ? mathnf(matalgtobasis(Labs, A))
   %5 =
   [17 8 4 2]
   
   [ 0 1 0 0]
   
   [ 0 0 1 0]
   
   [ 0 0 0 1]
   ? % == m
   %6 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealabstorel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealdown(*argv):
  '''
  rnfidealdown
  Class: basic
  Section: number_fields
  C-Name: rnfidealdown
  Prototype: GG
  Help: rnfidealdown(rnf,x): finds the intersection of the ideal x with the
   base field.
  Doc: let $\var{rnf}$ be a relative number
   field extension $L/K$ as output by \kbd{rnfinit}, and $x$ an ideal of
   $L$, given either in relative form or by a $\Z$-basis of elements of $L$
   (see \secref{se:rnfidealabstorel}). This function returns the ideal of $K$
   below $x$, i.e.~the intersection of $x$ with $K$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealdown(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealhnf(*argv):
  '''
  rnfidealhnf
  Class: basic
  Section: number_fields
  C-Name: rnfidealhnf
  Prototype: GG
  Help: rnfidealhnf(rnf,x): relative version of idealhnf, where rnf is a
   relative numberfield.
  Doc: $\var{rnf}$ being a relative number
   field extension $L/K$ as output by \kbd{rnfinit} and $x$ being a relative
   ideal (which can be, as in the absolute case, of many different types,
   including of course elements), computes the HNF pseudo-matrix associated to
   $x$, viewed as a $\Z_K$-module.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealhnf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealmul(*argv):
  '''
  rnfidealmul
  Class: basic
  Section: number_fields
  C-Name: rnfidealmul
  Prototype: GGG
  Help: rnfidealmul(rnf,x,y): relative version of idealmul, where rnf is a
   relative numberfield.
  Doc: $\var{rnf}$ being a relative number
   field extension $L/K$ as output by \kbd{rnfinit} and $x$ and $y$ being ideals
   of the relative extension $L/K$ given by pseudo-matrices, outputs the ideal
   product, again as a relative ideal.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealmul(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealnormabs(*argv):
  '''
  rnfidealnormabs
  Class: basic
  Section: number_fields
  C-Name: rnfidealnormabs
  Prototype: GG
  Help: rnfidealnormabs(rnf,x): absolute norm of the ideal x.
  Doc: let $\var{rnf}$ be a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and let $x$ be a
   relative ideal (which can be, as in the absolute case, of many different
   types, including of course elements). This function computes the norm of the
   $x$ considered as an ideal of the absolute extension $L/\Q$. This is
   identical to
   \bprog
      idealnorm(rnf, rnfidealnormrel(rnf,x))
   @eprog\noindent but faster.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealnormabs(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealnormrel(*argv):
  '''
  rnfidealnormrel
  Class: basic
  Section: number_fields
  C-Name: rnfidealnormrel
  Prototype: GG
  Help: rnfidealnormrel(rnf,x): relative norm of the ideal x.
  Doc: let $\var{rnf}$ be a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and let $x$ be a
   relative ideal (which can be, as in the absolute case, of many different
   types, including of course elements). This function computes the relative
   norm of $x$ as an ideal of $K$ in HNF.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealnormrel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealreltoabs(*argv):
  '''
  rnfidealreltoabs
  Class: basic
  Section: number_fields
  C-Name: rnfidealreltoabs
  Prototype: GG
  Help: rnfidealreltoabs(rnf,x): transforms the ideal x from relative to
   absolute representation.
  Doc: let $\var{rnf}$ be a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and let $x$ be a
   relative ideal, given as a $\Z_K$-module by a pseudo matrix $[A,I]$.
   This function returns the ideal $x$ as an absolute ideal of $L/\Q$ in
   the form of a $\Z$-basis, given by a vector of polynomials (modulo
   \kbd{rnf.pol}).
   
   The reason why we do not return the customary HNF in terms of a fixed
   $\Z$-basis for $\Z_L$ is precisely that no such basis has been explicitly
   specified. On the other hand, if you already computed an (absolute) \var{nf}
   structure \kbd{Labs} associated to $L$, then
   \bprog
     xabs = rnfidealreltoabs(L, x);
     xLabs = mathnf(matalgtobasis(Labs, xabs));
   @eprog\noindent computes a traditional HNF \kbd{xLabs} for $x$ in terms of
   the fixed $\Z$-basis \kbd{Labs.zk}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealreltoabs(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealtwoelt(*argv):
  '''
  rnfidealtwoelt
  Class: basic
  Section: number_fields
  C-Name: rnfidealtwoelement
  Prototype: GG
  Help: rnfidealtwoelt(rnf,x): relative version of idealtwoelt, where rnf
   is a relative numberfield.
  Doc: $\var{rnf}$ being a relative
   number field extension $L/K$ as output by \kbd{rnfinit} and $x$ being an
   ideal of the relative extension $L/K$ given by a pseudo-matrix, gives a
   vector of two generators of $x$ over $\Z_L$ expressed as polmods with polmod
   coefficients.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealtwoelement(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfidealup(*argv):
  '''
  rnfidealup
  Class: basic
  Section: number_fields
  C-Name: rnfidealup
  Prototype: GG
  Help: rnfidealup(rnf,x): lifts the ideal x (of the base field) to the
   relative field.
  Doc: let $\var{rnf}$ be a relative number
   field extension $L/K$ as output by \kbd{rnfinit} and let $x$ be an ideal of
   $K$. This function returns the ideal $x\Z_L$ as an absolute ideal of $L/\Q$,
   in the form of a $\Z$-basis, given by a vector of polynomials (modulo
   \kbd{rnf.pol}).
   
   The reason why we do not return the customary HNF in terms of a fixed
   $\Z$-basis for $\Z_L$ is precisely that no such basis has been explicitly
   specified. On the other hand, if you already computed an (absolute) \var{nf}
   structure \kbd{Labs} associated to $L$, then
   \bprog
     xabs = rnfidealup(L, x);
     xLabs = mathnf(matalgtobasis(Labs, xabs));
   @eprog\noindent computes a traditional HNF \kbd{xLabs} for $x$ in terms of
   the fixed $\Z$-basis \kbd{Labs.zk}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfidealup(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfinit(*argv):
  '''
  rnfinit
  Class: basic
  Section: number_fields
  C-Name: rnfinit
  Prototype: GG
  Help: rnfinit(nf,pol): pol being an irreducible polynomial
   defined over the number field nf, initializes a vector of data necessary for
   working in relative number fields (rnf functions). See manual for technical
   details.
  Doc: $\var{nf}$ being a number field in \kbd{nfinit}
   format considered as base field, and \var{pol} a polynomial defining a relative
   extension over $\var{nf}$, this computes data to work in the
   relative extension. The main variable of \var{pol} must be of higher priority
   (see \secref{se:priority}) than that of $\var{nf}$, and the coefficients of
   \var{pol} must be in $\var{nf}$.
   
   The result is a row vector, whose components are technical. In the following
   description, we let $K$ be the base field defined by $\var{nf}$ and $L/K$
   the large field associated to the \var{rnf}. Furthermore, we let
   $m = [K:\Q]$ the degree of the base field, $n = [L:K]$ the relative degree,
   $r_1$ and $r_2$ the number of real and complex places of $K$. Acces to this
   information via \emph{member functions} is preferred since the specific
   data organization specified below will change in the future.
   
   $\var{rnf}[1]$(\kbd{rnf.pol}) contains the relative polynomial \var{pol}.
   
   $\var{rnf}[2]$ contains the integer basis $[A,d]$ of $K$, as
   (integral) elements of $L/\Q$. More precisely, $A$ is a vector of
   polynomial with integer coefficients, $d$ is a denominator, and the integer
   basis is given by $A/d$.
   
   $\var{rnf}[3]$ (\kbd{rnf.disc}) is a two-component row vector
   $[\goth{d}(L/K),s]$ where $\goth{d}(L/K)$ is the relative ideal discriminant
   of $L/K$ and $s$ is the discriminant of $L/K$ viewed as an element of
   $K^*/(K^*)^2$, in other words it is the output of \kbd{rnfdisc}.
   
   $\var{rnf}[4]$(\kbd{rnf.index}) is the ideal index $\goth{f}$, i.e.~such
   that $d(pol)\Z_K=\goth{f}^2\goth{d}(L/K)$.
   
   $\var{rnf}[5]$ is currently unused.
   
   $\var{rnf}[6]$ is currently unused.
   
   $\var{rnf}[7]$ (\kbd{rnf.zk}) is the pseudo-basis $(A,I)$ for the maximal
   order $\Z_L$ as a $\Z_K$-module: $A$ is the relative integral pseudo basis
   expressed as polynomials (in the variable of $pol$) with polmod coefficients
   in $\var{nf}$, and the second component $I$ is the ideal list of the
   pseudobasis in HNF.
   
   $\var{rnf}[8]$ is the inverse matrix of the integral basis matrix, with
   coefficients polmods in $\var{nf}$.
   
   $\var{rnf}[9]$ is currently unused.
   
   $\var{rnf}[10]$ (\kbd{rnf.nf}) is $\var{nf}$.
   
   $\var{rnf}[11]$ is the output of \kbd{rnfequation(K, pol, 1)}. Namely, a
   vector $[P, a, k]$ describing the \emph{absolute} extension
   $L/\Q$: $P$ is an absolute equation, more conveniently obtained
   as \kbd{rnf.polabs}; $a$ expresses the generator $\alpha = y \mod \kbd{K.pol}$
   of the number field $K$ as an element of $L$, i.e.~a polynomial modulo the
   absolute equation $P$;
   
   $k$ is a small integer such that, if $\beta$ is an abstract root of \var{pol}
   and $\alpha$ the generator of $K$ given above, then $P(\beta + k\alpha) = 0$.
   
   \misctitle{Caveat.} Be careful if $k\neq0$ when dealing simultaneously with
   absolute and relative quantities since $L = \Q(\beta + k\alpha) =
   K(\alpha)$, and the generator chosen for the absolute extension is not the
   same as for the relative one. If this happens, one can of course go on
   working, but we advise to change the relative polynomial so that its root
   becomes $\beta + k \alpha$. Typical GP instructions would be
   \bprog
     [P,a,k] = rnfequation(K, pol, 1);
     if (k, pol = subst(pol, x, x - k*Mod(y, K.pol)));
     L = rnfinit(K, pol);
   @eprog
   
   $\var{rnf}[12]$ is by default unused and set equal to 0. This field is used
   to store further information about the field as it becomes available (which
   is rarely needed, hence would be too expensive to compute during the initial
   \kbd{rnfinit} call).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def algdep(*argv):
  '''
  algdep
  Class: basic
  Section: linear_algebra
  C-Name: algdep0
  Prototype: GLD0,L,
  Help: algdep(x,k,{flag=0}): algebraic relations up to degree n of x, using
   lindep([1,x,...,x^(k-1)], flag).
  Doc: \sidx{algebraic dependence}
   $x$ being real/complex, or $p$-adic, finds a polynomial of degree at most
   $k$ with integer coefficients having $x$ as approximate root. Note that the
   polynomial which is obtained is not necessarily the ``correct'' one. In fact
   it is not even guaranteed to be irreducible. One can check the closeness
   either by a polynomial evaluation (use \tet{subst}), or by computing the
   roots of the polynomial given by \kbd{algdep} (use \tet{polroots}).
   
   Internally, \tet{lindep}$([1,x,\ldots,x^k], \fl)$ is used.
   A non-zero value of $\fl$ may improve on the default behavior
   if the input number is known to a \emph{huge} accuracy, and you suspect the
   last bits are incorrect  (this truncates the number, throwing away the least
   significant bits), but default values are usually sufficient:
   \bprog
   ? \p200
   ? algdep(2^(1/6)+3^(1/5), 30);      \\ wrong in 0.8s
   ? algdep(2^(1/6)+3^(1/5), 30, 100); \\ wrong in 0.4s
   ? algdep(2^(1/6)+3^(1/5), 30, 170); \\ right in 0.8s
   ? algdep(2^(1/6)+3^(1/5), 30, 200); \\ wrong in 1.0s
   ? \p250
   ? algdep(2^(1/6)+3^(1/5), 30);      \\ right in 1.0s
   ? algdep(2^(1/6)+3^(1/5), 30, 200); \\ right in 1.0s
   ? \p500
   ? algdep(2^(1/6)+3^(1/5), 30);      \\ right in 2.9s
   ? \p1000
   ? algdep(2^(1/6)+3^(1/5), 30);      \\ right in 10.6s
   @eprog\noindent
   The changes in \kbd{defaultprecision} only affect the quality of the
   initial approximation to $2^{1/6} + 3^{1/5}$, \kbd{algdep} itself uses
   exact operations (the size of its operands depend on the accuracy of the
   input of course: more accurate input means slower operations).
   
   Proceeding by increments of 5 digits of accuracy, \kbd{algdep} with default
   flag produces its first correct result at 205 digits, and from then on a
   steady stream of correct results.
   
   The above example is the test case studied in a 2000 paper by Borwein and
   Lisonek: Applications of integer relation algorithms, \emph{Discrete Math.},
   {\bf 217}, p.~65--82. The version of PARI tested there was 1.39, which
   succeeded reliably from precision 265 on, in about 200 as much time as the
   current version.
  Variant: Also available is \fun{GEN}{algdep}{GEN x, long k} ($\fl=0$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.algdep0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def forqfvec(*argv):
  '''
  forqfvec
  Class: basic
  Section: linear_algebra
  C-Name: forqfvec0
  Prototype: vVGDGI
  Help: forqfvec(v,q,b,expr): q being a square and symmetric matrix
   representing a positive definite quadratic form, evaluate expr for all
   vector v such that q(v)<=b.
  Doc: $q$ being a square and symmetric matrix representing a positive definite
   quadratic form, evaluate \kbd{expr} for all vector $v$ such that $q(v)\leq b$.
   The formal variable $v$ runs through all such vectors in turn.
   \bprog
   ? forqfvec(v, [3,2;2,3], 3, print(v))
   [0, 1]~
   [1, 0]~
   [-1, 1]~
   @eprog
  Variant: The following function is also available:
   \fun{void}{forqfvec}{void *E, long (*fun)(void *, GEN, double), GEN q, GEN b}:
   Evaluate \kbd{fun(E,v,m)} on all $v$ such that $q(v)<b$, where $v$ is a
   \typ{VECSMALL} and $m=q(v)$ is a C double. The function \kbd{fun} must
   return $0$, unless \kbd{forqfvec} should stop, in which case, it should
   return $1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.forqfvec0(*c_arg_tuple)

def lindep(*argv):
  '''
  lindep
  Class: basic
  Section: linear_algebra
  C-Name: lindep0
  Prototype: GD0,L,
  Help: lindep(v,{flag=0}): integral linear dependencies between components of v.
   flag is optional, and can be 0: default, guess a suitable
   accuracy, or positive: accuracy to use for the computation, in decimal
   digits.
  Doc: \sidx{linear dependence} finds a small non-trivial integral linear
   combination between components of $v$. If none can be found return an empty
   vector.
   
   If $v$ is a vector with real/complex entries we use a floating point
   (variable precision) LLL algorithm. If $\fl = 0$ the accuracy is chosen
   internally using a crude heuristic. If $\fl > 0$ the computation is done with
   an accuracy of $\fl$ decimal digits. To get meaningful results in the latter
   case, the parameter $\fl$ should be smaller than the number of correct
   decimal digits in the input.
   
   \bprog
   ? lindep([sqrt(2), sqrt(3), sqrt(2)+sqrt(3)])
   %1 = [-1, -1, 1]~
   @eprog
   
   If $v$ is $p$-adic, $\fl$ is ignored and the algorithm LLL-reduces a
   suitable (dual) lattice.
   \bprog
   ? lindep([1, 2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)])
   %2 = [1, -2]~
   @eprog
   
   If $v$ is a matrix, $\fl$ is ignored and the function returns a non trivial
   kernel vector (combination of the columns).
   \bprog
   ? lindep([1,2,3;4,5,6;7,8,9])
   %3 = [1, -2, 1]~
   @eprog
   
   If $v$ contains polynomials or power series over some base field, finds a
   linear relation with coefficients in the field.
   \bprog
   ? lindep([x*y, x^2 + y, x^2*y + x*y^2, 1])
   %4 = [y, y, -1, -y^2]~
   @eprog\noindent For better control, it is preferable to use \typ{POL} rather
   than \typ{SER} in the input, otherwise one gets a linear combination which is
   $t$-adically small, but not necessarily $0$. Indeed, power series are first
   converted to the minimal absolute accuracy occuring among the entries of $v$
   (which can cause some coefficients to be ignored), then truncated to
   polynomials:
   \bprog
   ? v = [t^2+O(t^4), 1+O(t^2)]; L=lindep(v)
   %1 = [1, 0]~
   ? v*L
   %2 = t^2+O(t^4)  \\ small but not 0
   @eprog
  Variant: Also available are \fun{GEN}{lindep}{GEN v} (real/complex entries,
   $\fl=0$), \fun{GEN}{lindep2}{GEN v, long flag} (real/complex entries)
   \fun{GEN}{padic_lindep}{GEN v} ($p$-adic entries) and
   \fun{GEN}{Xadic_lindep}{GEN v} (polynomial entries).
   Finally \fun{GEN}{deplin}{GEN v} returns a non-zero kernel vector for a
   \typ{MAT} input.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.lindep0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathouseholder(*argv):
  '''
  mathouseholder
  Class: basic
  Section: linear_algebra
  C-Name: mathouseholder
  Prototype: GG
  Help: mathouseholder(Q,v): applies a sequence Q of Householder transforms
   to the vector or matrix v.
  Doc: \sidx{Householder transform}applies a sequence $Q$ of Householder
   transforms, as returned by \kbd{matqr}$(M,1)$ to the vector or matrix $v$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mathouseholder(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matqr(*argv):
  '''
  matqr
  Class: basic
  Section: linear_algebra
  C-Name: matqr
  Prototype: GD0,L,p
  Help: matqr(M,{flag=0}): returns [Q,R], the QR-decomposition of the square
   invertible matrix M. If flag=1, Q is given as a sequence of Householder
   transforms (faster and stabler).
  Doc: returns $[Q,R]$, the \idx{QR-decomposition} of the square invertible
   matrix $M$ with real entries: $Q$ is orthogonal and $R$ upper triangular. If
   $\fl=1$, the orthogonal matrix is returned as a sequence of Householder
   transforms: applying such a sequence is stabler and faster than
   multiplication by the corresponding $Q$ matrix.\sidx{Householder transform}
   More precisely, if
   \bprog
     [Q,R] = matqr(M);
     [q,r] = matqr(M, 1);
   @eprog\noindent then $r = R$ and \kbd{mathouseholder}$(q, M)$ is $R$;
   furthermore
   \bprog
     mathouseholder(q, matid(#M)) == Q~
   @eprog\noindent the inverse of $Q$. This function raises an error if the
   precision is too low or $x$ is singular.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matqr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfperfection(*argv):
  '''
  qfperfection
  Class: basic
  Section: linear_algebra
  C-Name: perf
  Prototype: G
  Help: qfperfection(G): rank of matrix of xx~ for x minimal vectors of a gram
   matrix G.
  Doc: 
   $G$ being a square and symmetric matrix with
   integer entries representing a positive definite quadratic form, outputs the
   perfection rank of the form. That is, gives the rank of the family of the $s$
   symmetric matrices $v_iv_i^t$, where $s$ is half the number of minimal
   vectors and the $v_i$ ($1\le i\le s$) are the minimal vectors.
   
   Since this requires computing the minimal vectors, the computations can
   become very lengthy as the dimension of $x$ grows.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.perf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfrep(*argv):
  '''
  qfrep
  Class: basic
  Section: linear_algebra
  C-Name: qfrep0
  Prototype: GGD0,L,
  Help: qfrep(q,B,{flag=0}): vector of (half) the number of vectors of norms
   from 1 to B for the integral and definite quadratic form q. If flag is 1,
   count vectors of even norm from 1 to 2B.
  Doc: 
   $q$ being a square and symmetric matrix with integer entries representing a
   positive definite quadratic form, count the vectors representing successive
   integers.
   
   \item If $\fl = 0$, count all vectors. Outputs the vector whose $i$-th
   entry, $1 \leq i \leq B$ is half the number of vectors $v$ such that $q(v)=i$.
   
   \item If $\fl = 1$, count vectors of even norm. Outputs the vector
   whose $i$-th entry, $1 \leq i \leq B$ is half the number of vectors such
   that $q(v) = 2i$.
   
   \bprog
   ? q = [2, 1; 1, 3];
   ? qfrep(q, 5)
   %2 = Vecsmall([0, 1, 2, 0, 0]) \\ 1 vector of norm 2, 2 of norm 3, etc.
   ? qfrep(q, 5, 1)
   %3 = Vecsmall([1, 0, 0, 1, 0]) \\ 1 vector of norm 2, 0 of norm 4, etc.
   @eprog\noindent
   This routine uses a naive algorithm based on \tet{qfminim}, and
   will fail if any entry becomes larger than $2^{31}$ (or $2^{63}$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfrep0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfminim(*argv):
  '''
  qfminim
  Class: basic
  Section: linear_algebra
  C-Name: qfminim0
  Prototype: GDGDGD0,L,p
  Help: qfminim(x,{b},{m},{flag=0}): x being a square and symmetric
   matrix representing a positive definite quadratic form, this function
   deals with the vectors of x whose norm is less than or equal to b,
   enumerated using the Fincke-Pohst algorithm, storing at most m vectors (no
   limit if m is omitted). The function searches for
   the minimal non-zero vectors if b is omitted. The precise behavior
   depends on flag. 0: seeks at most 2m vectors (unless m omitted), returns
   [N,M,mat] where N is the number of vectors found, M the maximum norm among
   these, and mat lists half the vectors (the other half is given by -mat). 1:
   ignores m and returns the first vector whose norm is less than b. 2: as 0
   but uses a more robust, slower implementation, valid for non integral
   quadratic forms.
  Doc: $x$ being a square and symmetric matrix representing a positive definite
   quadratic form, this function deals with the vectors of $x$ whose norm is
   less than or equal to $b$, enumerated using the Fincke-Pohst algorithm,
   storing at most $m$ vectors (no limit if $m$ is omitted). The function
   searches for the minimal non-zero vectors if $b$ is omitted. The behavior is
   undefined if $x$ is not positive definite (a ``precision too low'' error is
   most likely, although more precise error messages are possible). The precise
   behavior depends on $\fl$.
   
   If $\fl=0$ (default), seeks at most $2m$ vectors. The result is a
   three-component vector, the first component being the number of vectors
   found, the second being the maximum norm found, and the last vector is a
   matrix whose columns are the vectors found, only one being given for each
   pair $\pm v$ (at most $m$ such pairs, unless $m$ was omitted). The vectors
   are returned in no particular order.
   
   If $\fl=1$, ignores $m$ and returns $[N,v]$, where $v$ is a non-zero vector
   of length $N \leq b$, or $[]$ if no non-zero vector has length $\leq b$.
   If no explicit $b$ is provided, return a vector of smallish norm
   (smallest vector in an LLL-reduced basis).
   
   In these two cases, $x$ must have \emph{integral} entries. The
   implementation uses low precision floating point computations for maximal
   speed, which gives incorrect result when $x$ has large entries. (The
   condition is checked in the code and the routine raises an error if
   large rounding errors occur.) A more robust, but much slower,
   implementation is chosen if the following flag is used:
   
   If $\fl=2$, $x$ can have non integral real entries. In this case, if $b$
   is omitted, the ``minimal'' vectors only have approximately the same norm.
   If $b$ is omitted, $m$ is an upper bound for the number of vectors that
   will be stored and returned, but all minimal vectors are nevertheless
   enumerated. If $m$ is omitted, all vectors found are stored and returned;
   note that this may be a huge vector!
   
   \bprog
   ? x = matid(2);
   ? qfminim(x)  \\@com 4 minimal vectors of norm 1: $\pm[0,1]$, $\pm[1,0]$
   %2 = [4, 1, [0, 1; 1, 0]]
   ? { x =
   [4, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 1, 0,-1, 0, 0, 0,-2;
    2, 4,-2,-2, 0,-2, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 1,-1,-1;
    0,-2, 4, 0,-2, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0;
    0,-2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 1,-1, 0, 1,-1, 1, 0;
    0, 0,-2, 0, 4, 0, 0, 0, 1,-1, 0, 0, 1, 0, 0, 0,-2, 0, 0,-1, 1, 1, 0, 0;
   -2, -2,0, 0, 0, 4,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,-1, 1, 1;
    0, 0, 0, 0, 0,-2, 4,-2, 0, 0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 1,-1, 0;
    0, 0, 0, 0, 0, 0,-2, 4, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1,-1,-1, 0, 1, 0;
    0, 0, 0, 0, 1,-1, 0, 0, 4, 0,-2, 0, 1, 1, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0,-1, 0, 0, 0, 0, 4, 0, 0, 1, 1,-1, 1, 0, 0, 0, 1, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 4,-2, 0,-1, 0, 0, 0,-1, 0,-1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 4,-1, 1, 0, 0,-1, 1, 0, 1, 1, 1,-1, 0;
    1, 0,-1, 1, 1, 0, 0,-1, 1, 1, 0,-1, 4, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1,-1;
   -1,-1, 1,-1, 0, 0, 1, 0, 1, 1,-1, 1, 0, 4, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1, 1, 0, 4, 0, 0, 0, 0, 1, 1, 0, 0;
    0, 0, 1, 0,-2, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 4, 1, 1, 1, 0, 0, 1, 1;
    1, 0, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1, 0, 0, 0, 1, 4, 0, 1, 1, 0, 1, 0;
    0, 0, 0,-1, 0, 1, 0,-1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 4, 0, 1, 1, 0, 1;
   -1, -1,1, 0,-1, 1, 0,-1, 0, 1,-1, 1, 0, 1, 0, 0, 1, 1, 0, 4, 0, 0, 1, 1;
    0, 0,-1, 1, 1, 0, 0,-1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 4, 1, 0, 1;
    0, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 4, 0, 1;
    0,-1, 0, 1, 0, 1,-1, 1, 0, 1, 0,-1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 4, 1;
   -2,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 4]; }
   ? qfminim(x,,0)  \\ the Leech lattice has 196560 minimal vectors of norm 4
   time = 648 ms.
   %4 = [196560, 4, [;]]
   ? qfminim(x,,0,2); \\ safe algorithm. Slower and unnecessary here.
   time = 18,161 ms.
   %5 = [196560, 4.000061035156250000, [;]]
   @eprog\noindent\sidx{Leech lattice}\sidx{minimal vector}
   In the last example, we store 0 vectors to limit memory use. All minimal
   vectors are nevertheless enumerated. Provided \kbd{parisize} is about 50MB,
   \kbd{qfminim(x)} succeeds in 2.5 seconds.
  Variant: Also available are
   \fun{GEN}{minim}{GEN x, GEN b = NULL, GEN m = NULL} ($\fl=0$),
   \fun{GEN}{minim2}{GEN x, GEN b = NULL, GEN m = NULL} ($\fl=1$).
   \fun{GEN}{minim_raw}{GEN x, GEN b = NULL, GEN m = NULL} (do not perform LLL
   reduction on x).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_params.append(argv[4])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfminim0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def seralgdep(*argv):
  '''
  seralgdep
  Class: basic
  Section: linear_algebra
  C-Name: seralgdep
  Prototype: GLL
  Help: seralgdep(s,p,r): find a linear relation between powers (1,s, ..., s^p)
   of the series s, with polynomial coefficients of degree <= r.
  Doc: \sidx{algebraic dependence} finds a linear relation between powers $(1,s,
   \dots, s^p)$ of the series $s$, with polynomial coefficients of degree
   $\leq r$. In case no relation is found, return $0$.
   \bprog
   ? s = 1 + 10*y - 46*y^2 + 460*y^3 - 5658*y^4 + 77740*y^5 + O(y^6);
   ? seralgdep(s, 2, 2)
   %2 = -x^2 + (8*y^2 + 20*y + 1)
   ? subst(%, x, s)
   %3 = O(y^6)
   ? seralgdep(s, 1, 3)
   %4 = (-77*y^2 - 20*y - 1)*x + (310*y^3 + 231*y^2 + 30*y + 1)
   ? seralgdep(s, 1, 2)
   %5 = 0
   @eprog\noindent The series main variable must not be $x$, so as to be able
   to express the result as a polynomial in $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.seralgdep(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def zncoppersmith(*argv):
  '''
  zncoppersmith
  Class: basic
  Section: number_theoretical
  C-Name: zncoppersmith
  Prototype: GGGDG
  Help: zncoppersmith(P, N, X, {B=N}): finds all integers x
   with |x| <= X such that  gcd(N, P(x)) >= B. X should be smaller than
   exp((log B)^2 / (deg(P) log N)).
  Doc: $N$ being an integer and $P\in \Z[X]$, finds all integers $x$ with
   $|x| \leq X$ such that
   $$\gcd(N, P(x)) \geq B,$$
   using \idx{Coppersmith}'s algorithm (a famous application of the \idx{LLL}
   algorithm). $X$ must be smaller than $\exp(\log^2 B / (\deg(P) \log N))$:
   for $B = N$, this means $X < N^{1/\deg(P)}$. Some $x$ larger than $X$ may
   be returned if you are very lucky. The smaller $B$ (or the larger $X$), the
   slower the routine will be. The strength of Coppersmith method is the
   ability to find roots modulo a general \emph{composite} $N$: if $N$ is a prime
   or a prime power, \tet{polrootsmod} or \tet{polrootspadic} will be much
   faster.
   
   We shall now present two simple applications. The first one is
   finding non-trivial factors of $N$, given some partial information on the
   factors; in that case $B$ must obviously be smaller than the largest
   non-trivial divisor of $N$.
   \bprog
   setrand(1); \\ to make the example reproducible
   p = nextprime(random(10^30));
   q = nextprime(random(10^30)); N = p*q;
   p0 = p % 10^20; \\ assume we know 1) p > 10^29, 2) the last 19 digits of p
   p1 = zncoppersmith(10^19*x + p0, N, 10^12, 10^29)
   
   \\ result in 10ms.
   %1 = [35023733690]
   ? gcd(p1[1] * 10^19 + p0, N) == p
   %2 = 1
   @eprog\noindent and we recovered $p$, faster than by trying all
   possibilities $ < 10^{12}$.
   
   The second application is an attack on RSA with low exponent, when the
   message $x$ is short and the padding $P$ is known to the attacker. We use
   the same RSA modulus $N$ as in the first example:
   \bprog
   setrand(1);
   P = random(N);    \\ known padding
   e = 3;            \\ small public encryption exponent
   X = floor(N^0.3); \\ N^(1/e - epsilon)
   x0 = random(X);   \\ unknown short message
   C = lift( (Mod(x0,N) + P)^e ); \\ known ciphertext, with padding P
   zncoppersmith((P + x)^3 - C, N, X)
   
   \\ result in 244ms.
   %3 = [265174753892462432]
   ? %[1] == x0
   %4 = 1
   @eprog\noindent
   We guessed an integer of the order of $10^{18}$, almost instantly.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.zncoppersmith(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def binomial(*argv):
  '''
  binomial
  Class: basic
  Section: number_theoretical
  C-Name: binomial
  Prototype: GL
  Help: binomial(x,y): binomial coefficient x*(x-1)...*(x-y+1)/y! defined for
   y in Z and any x.
  Doc: \idx{binomial coefficient} $\binom{x}{y}$.
   Here $y$ must be an integer, but $x$ can be any PARI object.
  Variant: The function
   \fun{GEN}{binomialuu}{ulong n, ulong k} is also available, and so is
   \fun{GEN}{vecbinome}{long n}, which returns a vector $v$
   with $n+1$ components such that $v[k+1] = \kbd{binomial}(n,k)$ for $k$ from
   $0$ up to $n$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.binomial(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def cmp(*argv):
  '''
  cmp
  Class: basic
  Section: operators
  C-Name: cmp_universal
  Prototype: iGG
  Help: cmp(x,y): compare two arbitrary objects x and y (1 if x>y, 0 if x=y, -1
   if x<y). The function is used to implement sets, and has no useful
   mathematical meaning.
  Doc: gives the result of a comparison between arbitrary objects $x$ and $y$
   (as $-1$, $0$ or $1$). The underlying order relation is transitive,
   the function returns $0$ if and only if $x~\kbd{===}~y$, and its
   restriction to integers coincides with the customary one. Besides that,
   it has no useful mathematical meaning.
   
   In case all components are equal up to the smallest length of the operands,
   the more complex is considered to be larger. More precisely, the longest is
   the largest; when lengths are equal, we have matrix $>$ vector $>$ scalar.
   For example:
   \bprog
   ? cmp(1, 2)
   %1 = -1
   ? cmp(2, 1)
   %2 = 1
   ? cmp(1, 1.0)   \\ note that 1 == 1.0, but (1===1.0) is false.
   %3 = -1
   ? cmp(x + Pi, [])
   %4 = -1
   @eprog\noindent This function is mostly useful to handle sorted lists or
   vectors of arbitrary objets. For instance, if $v$ is a vector, the
   construction \kbd{vecsort(v, cmp)} is equivalent to \kbd{Set(v)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.cmp_universal(*c_arg_tuple)

def serconvol(*argv):
  '''
  serconvol
  Class: basic
  Section: polynomials
  C-Name: convol
  Prototype: GG
  Help: serconvol(x,y): convolution (or Hadamard product) of two power series.
  Doc: convolution (or \idx{Hadamard product}) of the
   two power series $x$ and $y$; in other words if $x=\sum a_k*X^k$ and $y=\sum
   b_k*X^k$ then $\kbd{serconvol}(x,y)=\sum a_k*b_k*X^k$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.convol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polcyclo(*argv):
  '''
  polcyclo
  Class: basic
  Section: polynomials
  C-Name: polcyclo_eval
  Prototype: LDG
  Help: polcyclo(n,{a = 'x}): n-th cyclotomic polynomial evaluated at a.
  Description: 
    (small,?var):gen     polcyclo($1,$2)
    (small,gen):gen      polcyclo_eval($1,$2)
  Doc: $n$-th cyclotomic polynomial, evaluated at $a$ (\kbd{'x} by default). The
   integer $n$ must be positive.
   
   Algorithm used: reduce to the case where $n$ is squarefree; to compute the
   cyclotomic polynomial, use $\Phi_{np}(x)=\Phi_n(x^p)/\Phi(x)$; to compute
   it evaluated, use $\Phi_n(x) = \prod_{d\mid n} (x^d-1)^{\mu(n/d)}$. In the
   evaluated case, the algorithm assumes that $a^d - 1$ is either $0$ or
   invertible, for all $d\mid n$. If this is not the case (the base ring has
   zero divisors), use \kbd{subst(polcyclo(n),x,a)}.
  Variant: The variant \fun{GEN}{polcyclo}{long n, long v} returns the $n$-th
   cyclotomic polynomial in variable $v$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polcyclo_eval(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def dirdiv(*argv):
  '''
  dirdiv
  Class: basic
  Section: number_theoretical
  C-Name: dirdiv
  Prototype: GG
  Help: dirdiv(x,y): division of the Dirichlet series x by the Dirichlet
   series y.
  Doc: $x$ and $y$ being vectors of perhaps different
   lengths but with $y[1]\neq 0$ considered as \idx{Dirichlet series}, computes
   the quotient of $x$ by $y$, again as a vector.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.dirdiv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def dirmul(*argv):
  '''
  dirmul
  Class: basic
  Section: number_theoretical
  C-Name: dirmul
  Prototype: GG
  Help: dirmul(x,y): multiplication of the Dirichlet series x by the Dirichlet
   series y.
  Doc: $x$ and $y$ being vectors of perhaps different lengths representing
   the \idx{Dirichlet series} $\sum_n x_n n^{-s}$ and $\sum_n y_n n^{-s}$,
   computes the product of $x$ by $y$, again as a vector.
   \bprog
   ? dirmul(vector(10,n,1), vector(10,n,moebius(n)))
   %1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
   @eprog\noindent
   The product
   length is the minimum of $\kbd{\#}x\kbd{*}v(y)$ and $\kbd{\#}y\kbd{*}v(x)$,
   where $v(x)$ is the index of the first non-zero coefficient.
   \bprog
   ? dirmul([0,1], [0,1]);
   %2 = [0, 0, 0, 1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.dirmul(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def getstack(*argv):
  '''
  getstack
  Class: basic
  Section: programming/specific
  C-Name: getstack
  Prototype: l
  Help: getstack(): current value of stack pointer avma.
  Doc: returns the current value of $\kbd{top}-\kbd{avma}$, i.e.~the number of
   bytes used up to now on the stack. Useful mainly for debugging purposes.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.getstack(*c_arg_tuple)

def gettime(*argv):
  '''
  gettime
  Class: basic
  Section: programming/specific
  C-Name: gettime
  Prototype: l
  Help: gettime(): time (in milliseconds) since last call to gettime.
  Doc: returns the time (in milliseconds) elapsed since either the last call to
   \kbd{gettime}, or to the beginning of the containing GP instruction (if
   inside \kbd{gp}), whichever came last.
   
   For a reentrant version, see \tet{getabstime}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.gettime(*c_arg_tuple)

def getabstime(*argv):
  '''
  getabstime
  Class: basic
  Section: programming/specific
  C-Name: getabstime
  Prototype: l
  Help: getabstime(): time (in milliseconds) since startup.
  Doc: returns the time (in milliseconds) elapsed since \kbd{gp} startup. This
   provides a reentrant version of \kbd{gettime}:
   \bprog
   my (t = getabstime());
   ...
   print("Time: ", getabstime() - t);
   @eprog
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.getabstime(*c_arg_tuple)

def Set(*argv):
  '''
  Set
  Class: basic
  Section: conversions
  C-Name: gtoset
  Prototype: DG
  Help: Set({x=[]}): convert x into a set, i.e. a row vector with strictly
   increasing coefficients. Empty set if x is omitted.
  Description: 
   ():vec           cgetg(1,t_VEC)
   (gen):vec        gtoset($1)
  Doc: 
   converts $x$ into a set, i.e.~into a row vector, with strictly increasing
   entries with respect to the (somewhat arbitrary) universal comparison function
   \tet{cmp}. Standard container types \typ{VEC}, \typ{COL}, \typ{LIST} and
   \typ{VECSMALL} are converted to the set with corresponding elements. All
   others are converted to a set with one element.
   \bprog
   ? Set([1,2,4,2,1,3])
   %1 = [1, 2, 3, 4]
   ? Set(x)
   %2 = [x]
   ? Set(Vecsmall([1,3,2,1,3]))
   %3 = [1, 2, 3]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtoset(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def serlaplace(*argv):
  '''
  serlaplace
  Class: basic
  Section: polynomials
  C-Name: laplace
  Prototype: G
  Help: serlaplace(x): replaces the power series sum of a_n*x^n/n! by sum of
   a_n*x^n. For the reverse operation, use serconvol(x,exp(X)).
  Doc: $x$ must be a power series with non-negative
   exponents. If $x=\sum (a_k/k!)*X^k$ then the result is $\sum a_k*X^k$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.laplace(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathilbert(*argv):
  '''
  mathilbert
  Class: basic
  Section: linear_algebra
  C-Name: mathilbert
  Prototype: L
  Help: mathilbert(n): Hilbert matrix of order n.
  Doc: $x$ being a \kbd{long}, creates the
   \idx{Hilbert matrix}of order $x$, i.e.~the matrix whose coefficient
   ($i$,$j$) is $1/ (i+j-1)$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mathilbert(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matpascal(*argv):
  '''
  matpascal
  Class: basic
  Section: linear_algebra
  C-Name: matqpascal
  Prototype: LDG
  Help: matpascal(n,{q}): Pascal triangle of order n if q is omited. q-Pascal
   triangle otherwise.
  Doc: creates as a matrix the lower triangular
   \idx{Pascal triangle} of order $x+1$ (i.e.~with binomial coefficients
   up to $x$). If $q$ is given, compute the $q$-Pascal triangle (i.e.~using
   $q$-binomial coefficients).
  Variant: Also available is \fun{GEN}{matpascal}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matqpascal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def modreverse(*argv):
  '''
  modreverse
  Class: basic
  Section: number_fields
  C-Name: modreverse
  Prototype: G
  Help: modreverse(z): reverse polmod of the polmod z, if it exists.
  Doc: let $z = \kbd{Mod(A, T)}$ be a polmod, and $Q$ be its minimal
   polynomial, which must satisfy $\text{deg}(Q) = \text{deg}(T)$.
   Returns a ``reverse polmod'' \kbd{Mod(B, Q)}, which is a root of $T$.
   
   This is quite useful when one changes the generating element in algebraic
   extensions:
   \bprog
   ? u = Mod(x, x^3 - x -1); v = u^5;
   ? w = modreverse(v)
   %2 = Mod(x^2 - 4*x + 1, x^3 - 5*x^2 + 4*x - 1)
   @eprog\noindent
   which means that $x^3 - 5x^2 + 4x -1$ is another defining polynomial for the
   cubic field
   $$\Q(u) = \Q[x]/(x^3 - x - 1) = \Q[x]/(x^3 - 5x^2 + 4x - 1) = \Q(v),$$
   and that $u \to v^2 - 4v + 1$ gives an explicit isomorphism. From this, it is
   easy to convert elements between the $A(u)\in \Q(u)$ and $B(v)\in \Q(v)$
   representations:
   \bprog
   ? A = u^2 + 2*u + 3; subst(lift(A), 'x, w)
   %3 = Mod(x^2 - 3*x + 3, x^3 - 5*x^2 + 4*x - 1)
   ? B = v^2 + v + 1;   subst(lift(B), 'x, v)
   %4 = Mod(26*x^2 + 31*x + 26, x^3 - x - 1)
   @eprog
   If the minimal polynomial of $z$ has lower degree than expected, the routine
   fails
   \bprog
   ? u = Mod(-x^3 + 9*x, x^4 - 10*x^2 + 1)
   ? modreverse(u)
    *** modreverse: domain error in modreverse: deg(minpoly(z)) < 4
    ***   Break loop: type 'break' to go back to GP prompt
   break> Vec( dbg_err() ) \\ ask for more info
   ["e_DOMAIN", "modreverse", "deg(minpoly(z))", "<", 4,
     Mod(-x^3 + 9*x, x^4 - 10*x^2 + 1)]
   break> minpoly(u)
   x^2 - 8
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.modreverse(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def numtoperm(*argv):
  '''
  numtoperm
  Class: basic
  Section: conversions
  C-Name: numtoperm
  Prototype: LG
  Help: numtoperm(n,k): permutation number k (mod n!) of n letters (n
   C-integer).
  Doc: generates the $k$-th permutation (as a row vector of length $n$) of the
   numbers $1$ to $n$. The number $k$ is taken modulo $n!\,$, i.e.~inverse
   function of \tet{permtonum}. The numbering used is the standard lexicographic
   ordering, starting at $0$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.numtoperm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def permtonum(*argv):
  '''
  permtonum
  Class: basic
  Section: conversions
  C-Name: permtonum
  Prototype: G
  Help: permtonum(x): ordinal (between 1 and n!) of permutation x.
  Doc: given a permutation $x$ on $n$ elements, gives the number $k$ such that
   $x=\kbd{numtoperm(n,k)}$, i.e.~inverse function of \tet{numtoperm}.
   The numbering used is the standard lexicographic ordering, starting at $0$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.permtonum(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polhermite(*argv):
  '''
  polhermite
  Class: basic
  Section: polynomials
  C-Name: polhermite_eval
  Prototype: LDG
  Help: polhermite(n,{a='x}): Hermite polynomial H(n,v) of degree n, evaluated
   at a.
  Description: 
    (small,?var):gen    polhermite($1,$2)
    (small,gen):gen     polhermite_eval($1,$2)
  Doc: $n^{\text{th}}$ \idx{Hermite} polynomial $H_n$ evaluated at $a$
   (\kbd{'x} by default), i.e.
   $$ H_n(x) = (-1)^n\*e^{x^2} \dfrac{d^n}{dx^n}e^{-x^2}.$$
  Variant: The variant \fun{GEN}{polhermite}{long n, long v} returns the $n$-th
   Hermite polynomial in variable $v$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polhermite_eval(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def pollegendre(*argv):
  '''
  pollegendre
  Class: basic
  Section: polynomials
  C-Name: pollegendre_eval
  Prototype: LDG
  Help: pollegendre(n,{a='x}): legendre polynomial of degree n evaluated at a.
  Description: 
    (small,?var):gen    pollegendre($1,$2)
    (small,gen):gen     pollegendre_eval($1,$2)
  Doc: $n^{\text{th}}$ \idx{Legendre polynomial} evaluated at $a$ (\kbd{'x} by
   default).
  Variant: To obtain the $n$-th Legendre polynomial in variable $v$,
   use \fun{GEN}{pollegendre}{long n, long v}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.pollegendre_eval(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polinterpolate(*argv):
  '''
  polinterpolate
  Class: basic
  Section: polynomials
  C-Name: polint
  Prototype: GDGDGD&
  Help: polinterpolate(X,{Y},{x},{&e}): polynomial interpolation at x
   according to data vectors X, Y (ie return P such that P(X[i]) = Y[i] for
   all i). If Y is omitted, return P such that P(i) = X[i]. If present, e
   will contain an error estimate on the returned value.
  Doc: given the data vectors
   $X$ and $Y$ of the same length $n$ ($X$ containing the $x$-coordinates,
   and $Y$ the corresponding $y$-coordinates), this function finds the
   \idx{interpolating polynomial} passing through these points and evaluates it
   at~$x$. If $Y$ is omitted, return the polynomial interpolating the
   $(i,X[i])$. If present, $e$ will contain an error estimate on the returned
   value.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polchebyshev(*argv):
  '''
  polchebyshev
  Class: basic
  Section: polynomials
  C-Name: polchebyshev_eval
  Prototype: LD1,L,DG
  Help: polchebyshev(n,{flag=1},{a='x}): Chebychev polynomial of the first (flag
   = 1) or second (flag = 2) kind, of degree n, evaluated at a.
  Description: 
   (small,?1,?var):gen polchebyshev1($1,$3)
   (small,2,?var):gen  polchebyshev2($1,$3)
   (small,small,?var):gen polchebyshev($1,$2,$3)
  Doc: returns the $n^{\text{th}}$
   \idx{Chebyshev} polynomial of the first kind $T_n$ ($\fl=1$) or the second
   kind $U_n$ ($\fl=2$), evaluated at $a$ (\kbd{'x} by default). Both series of
   polynomials satisfy the 3-term relation
   $$ P_{n+1} = 2xP_n - P_{n-1}, $$
   and are determined by the initial conditions $U_0 = T_0 = 1$, $T_1 = x$,
   $U_1 = 2x$. In fact $T_n' = n U_{n-1}$ and, for all complex numbers $z$, we
   have $T_n(\cos z) = \cos (nz)$ and $U_{n-1}(\cos z) = \sin(nz)/\sin z$.
   If $n \geq 0$, then these polynomials have degree $n$.  For $n < 0$,
   $T_n$ is equal to $T_{-n}$ and $U_n$ is equal to $-U_{-2-n}$.
   In particular, $U_{-1} = 0$.
  Variant: Also available are
   \fun{GEN}{polchebyshev}{long n, long \fl, long v},
   \fun{GEN}{polchebyshev1}{long n, long v} and
   \fun{GEN}{polchebyshev2}{long n, long v} for $T_n$ and $U_n$ respectively.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polchebyshev_eval(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poltchebi(*argv):
  '''
  poltchebi
  Class: basic
  Section: polynomials
  C-Name: polchebyshev1
  Prototype: LDn
  Help: poltchebi(n,{v='x}): deprecated alias for polchebyshev
  Doc: deprecated alias for \kbd{polchebyshev}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polchebyshev1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polrecip(*argv):
  '''
  polrecip
  Class: basic
  Section: polynomials
  C-Name: polrecip
  Prototype: G
  Help: polrecip(pol): reciprocal polynomial of pol.
  Doc: reciprocal polynomial of \var{pol}, i.e.~the coefficients are in
   reverse order. \var{pol} must be a polynomial.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polrecip(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def setbinop(*argv):
  '''
  setbinop
  Class: basic
  Section: linear_algebra
  C-Name: setbinop
  Prototype: GGDG
  Help: setbinop(f,X,{Y}): the set {f(x,y), x in X, y in Y}. If Y is omitted,
   assume that X = Y and that f is symmetric.
  Doc: the set whose elements are the f(x,y), where x,y run through X,Y.
   respectively. If $Y$ is omitted, assume that $X = Y$ and that $f$ is symmetric:
   $f(x,y) = f(y,x)$ for all $x,y$ in $X$.
   \bprog
   ? X = [1,2,3]; Y = [2,3,4];
   ? setbinop((x,y)->x+y, X,Y) \\ set X + Y
   %2 = [3, 4, 5, 6, 7]
   ? setbinop((x,y)->x-y, X,Y) \\ set X - Y
   %3 = [-3, -2, -1, 0, 1]
   ? setbinop((x,y)->x+y, X)   \\ set 2X = X + X
   %2 = [2, 3, 4, 5, 6]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.setbinop(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def setintersect(*argv):
  '''
  setintersect
  Class: basic
  Section: linear_algebra
  C-Name: setintersect
  Prototype: GG
  Help: setintersect(x,y): intersection of the sets x and y.
  Description: 
   (vec, vec):vec        setintersect($1, $2)
  Doc: intersection of the two sets $x$ and $y$ (see \kbd{setisset}).
   If $x$ or $y$ is not a set, the result is undefined.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.setintersect(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def setisset(*argv):
  '''
  setisset
  Class: basic
  Section: linear_algebra
  C-Name: setisset
  Prototype: lG
  Help: setisset(x): true(1) if x is a set (row vector with strictly
   increasing entries), false(0) if not.
  Doc: 
   returns true (1) if $x$ is a set, false (0) if
   not. In PARI, a set is a row vector whose entries are strictly
   increasing with respect to a (somewhat arbitray) universal comparison
   function. To convert any object into a set (this is most useful for
   vectors, of course), use the function \kbd{Set}.
   \bprog
   ? a = [3, 1, 1, 2];
   ? setisset(a)
   %2 = 0
   ? Set(a)
   %3 = [1, 2, 3]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.setisset(*c_arg_tuple)

def setminus(*argv):
  '''
  setminus
  Class: basic
  Section: linear_algebra
  C-Name: setminus
  Prototype: GG
  Help: setminus(x,y): set of elements of x not belonging to y.
  Description: 
   (vec, vec):vec        setminus($1, $2)
  Doc: difference of the two sets $x$ and $y$ (see \kbd{setisset}),
   i.e.~set of elements of $x$ which do not belong to $y$.
   If $x$ or $y$ is not a set, the result is undefined.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.setminus(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def setsearch(*argv):
  '''
  setsearch
  Class: basic
  Section: linear_algebra
  C-Name: setsearch
  Prototype: lGGD0,L,
  Help: setsearch(S,x,{flag=0}): determines whether x belongs to the set (or
   sorted list) S.
   If flag is 0 or omitted, returns 0 if it does not, otherwise returns the index
   j such that x==S[j]. If flag is non-zero, return 0 if x belongs to S,
   otherwise the index j where it should be inserted.
  Doc: determines whether $x$ belongs to the set $S$ (see \kbd{setisset}).
   
   We first describe the default behaviour, when $\fl$ is zero or omitted. If $x$
   belongs to the set $S$, returns the index $j$ such that $S[j]=x$, otherwise
   returns 0.
   \bprog
   ? T = [7,2,3,5]; S = Set(T);
   ? setsearch(S, 2)
   %2 = 1
   ? setsearch(S, 4)      \\ not found
   %3 = 0
   ? setsearch(T, 7)      \\ search in a randomly sorted vector
   %4 = 0 \\ WRONG !
   @eprog\noindent
   If $S$ is not a set, we also allow sorted lists with
   respect to the \tet{cmp} sorting function, without repeated entries,
   as per \tet{listsort}$(L,1)$; otherwise the result is undefined.
   \bprog
   ? L = List([1,4,2,3,2]); setsearch(L, 4)
   %1 = 0 \\ WRONG !
   ? listsort(L, 1); L    \\ sort L first
   %2 = List([1, 2, 3, 4])
   ? setsearch(L, 4)
   %3 = 4                 \\ now correct
   @eprog\noindent
   If $\fl$ is non-zero, this function returns the index $j$ where $x$ should be
   inserted, and $0$ if it already belongs to $S$. This is meant to be used for
   dynamically growing (sorted) lists, in conjunction with \kbd{listinsert}.
   \bprog
   ? L = List([1,5,2,3,2]); listsort(L,1); L
   %1 = List([1,2,3,5])
   ? j = setsearch(L, 4, 1)  \\ 4 should have been inserted at index j
   %2 = 4
   ? listinsert(L, 4, j); L
   %3 = List([1, 2, 3, 4, 5])
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.setsearch(*c_arg_tuple)

def setunion(*argv):
  '''
  setunion
  Class: basic
  Section: linear_algebra
  C-Name: setunion
  Prototype: GG
  Help: setunion(x,y): union of the sets x and y.
  Description: 
   (vec, vec):vec        setunion($1, $2)
  Doc: union of the two sets $x$ and $y$ (see \kbd{setisset}).
   If $x$ or $y$ is not a set, the result is undefined.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.setunion(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def stirling(*argv):
  '''
  stirling
  Class: basic
  Section: number_theoretical
  C-Name: stirling
  Prototype: LLD1,L,
  Help: stirling(n,k,{flag=1}): If flag=1 (default) return the Stirling number
   of the first kind s(n,k), if flag=2, return the Stirling number of the second
   kind S(n,k).
  Doc: \idx{Stirling number} of the first kind $s(n,k)$ ($\fl=1$, default) or
   of the second kind $S(n,k)$ (\fl=2), where $n$, $k$ are non-negative
   integers. The former is $(-1)^{n-k}$ times the
   number of permutations of $n$ symbols with exactly $k$ cycles; the latter is
   the number of ways of partitioning a set of $n$ elements into $k$ non-empty
   subsets. Note that if all $s(n,k)$ are needed, it is much faster to compute
   $$\sum_k s(n,k) x^k = x(x-1)\dots(x-n+1).$$
   Similarly, if a large number of $S(n,k)$ are needed for the same $k$,
   one should use
   $$\sum_n S(n,k) x^n = \dfrac{x^k}{(1-x)\dots(1-kx)}.$$
   (Should be implemented using a divide and conquer product.) Here are
   simple variants for $n$ fixed:
   \bprog
   /* list of s(n,k), k = 1..n */
   vecstirling(n) = Vec( factorback(vector(n-1,i,1-i*'x)) )
   
   /* list of S(n,k), k = 1..n */
   vecstirling2(n) =
   { my(Q = x^(n-1), t);
     vector(n, i, t = divrem(Q, x-i); Q=t[1]; t[2]);
   }
   @eprog
  Variant: Also available are \fun{GEN}{stirling1}{ulong n, ulong k}
   ($\fl=1$) and \fun{GEN}{stirling2}{ulong n, ulong k} ($\fl=2$).
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.stirling(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def vecsearch(*argv):
  '''
  vecsearch
  Class: basic
  Section: linear_algebra
  C-Name: vecsearch
  Prototype: lGGDG
  Help: vecsearch(v,x,{cmpf}): determines whether x belongs to the sorted
   vector v. If the comparison function cmpf is explicitly given, assume
   that v was sorted according to vecsort(, cmpf).
  Doc: determines whether $x$ belongs to the sorted vector or list $v$: return
   the (positive) index where $x$ was found, or $0$ if it does not belong to
   $v$.
   
   If the comparison function cmpf is omitted, we assume that $v$ is sorted in
   increasing order, according to the standard comparison function $<$, thereby
   restricting the possible types for $x$ and the elements of $v$ (integers,
   fractions or reals).
   
   If \kbd{cmpf} is present, it is understood as a comparison function and we
   assume that $v$ is sorted according to it, see \tet{vecsort} for how to
   encode comparison functions.
   \bprog
   ? v = [1,3,4,5,7];
   ? vecsearch(v, 3)
   %2 = 2
   ? vecsearch(v, 6)
   %3 = 0 \\ not in the list
   ? vecsearch([7,6,5], 5) \\ unsorted vector: result undefined
   %4 = 0
   @eprog
   
   By abuse of notation, $x$ is also allowed to be a matrix, seen as a vector
   of its columns; again by abuse of notation, a \typ{VEC} is considered
   as part of the matrix, if its transpose is one of the matrix columns.
   \bprog
   ? v = vecsort([3,0,2; 1,0,2]) \\ sort matrix columns according to lex order
   %1 =
   [0 2 3]
   
   [0 2 1]
   ? vecsearch(v, [3,1]~)
   %2 = 3
   ? vecsearch(v, [3,1])  \\ can search for x or x~
   %3 = 3
   ? vecsearch(v, [1,2])
   %4 = 0 \\ not in the list
   @eprog\noindent
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.vecsearch(*c_arg_tuple)

def vecsort(*argv):
  '''
  vecsort
  Class: basic
  Section: linear_algebra
  C-Name: vecsort0
  Prototype: GDGD0,L,
  Help: vecsort(x,{cmpf},{flag=0}): sorts the vector of vectors (or matrix) x in
   ascending order, according to the comparison function cmpf, if not omitted.
   (If cmpf is an integer, sort according to the value of the k-th component
   of each entry.) Binary digits of flag (if present) mean: 1: indirect sorting,
   return the permutation instead of the permuted vector, 2: sort using
   lexicographic order, 4: use descending instead of ascending order, 8: remove
   duplicate entries.
  Description: 
   (vecsmall,?gen):vecsmall       vecsort0($1, $2, 0)
   (vecsmall,?gen,small):vecsmall vecsort0($1, $2, $3)
   (vec, , ?0):vec                sort($1)
   (vec, , 1):vecsmall            indexsort($1)
   (vec, , 2):vec                 lexsort($1)
   (vec, gen):vec                 vecsort0($1, $2, 0)
   (vec, ?gen, 1):vecsmall        vecsort0($1, $2, 1)
   (vec, ?gen, 3):vecsmall        vecsort0($1, $2, 3)
   (vec, ?gen, 5):vecsmall        vecsort0($1, $2, 5)
   (vec, ?gen, 7):vecsmall        vecsort0($1, $2, 7)
   (vec, ?gen, 9):vecsmall        vecsort0($1, $2, 9)
   (vec, ?gen, 11):vecsmall       vecsort0($1, $2, 11)
   (vec, ?gen, 13):vecsmall       vecsort0($1, $2, 13)
   (vec, ?gen, 15):vecsmall       vecsort0($1, $2, 15)
   (vec, ?gen, #small):vec        vecsort0($1, $2, $3)
   (vec, ?gen, small):gen         vecsort0($1, $2, $3)
  Doc: sorts the vector $x$ in ascending order, using a mergesort method.
   $x$ must be a list, vector or matrix (seen as a vector of its columns).
   Note that mergesort is stable, hence the initial ordering of ``equal''
   entries (with respect to the sorting criterion) is not changed.
   
   If \kbd{cmpf} is omitted, we use the standard comparison function
   \kbd{lex}, thereby restricting the possible types for the elements of $x$
   (integers, fractions or reals and vectors of those). If \kbd{cmpf} is
   present, it is understood as a comparison function and we sort according to
   it. The following possibilities exist:
   
   \item an integer $k$: sort according to the value of the $k$-th
   subcomponents of the components of~$x$.
   
   \item a vector: sort lexicographically according to the components listed in
   the vector. For example, if $\kbd{cmpf}=\kbd{[2,1,3]}$, sort with respect to
   the second component, and when these are equal, with respect to the first,
   and when these are equal, with respect to the third.
   
   \item a comparison function (\typ{CLOSURE}), with two arguments $x$ and $y$,
   and returning an integer which is $<0$, $>0$ or $=0$ if $x<y$, $x>y$ or
   $x=y$ respectively. The \tet{sign} function is very useful in this context:
   \bprog
   ? vecsort([3,0,2; 1,0,2]) \\ sort columns according to lex order
   %1 =
   [0 2 3]
   
   [0 2 1]
   ? vecsort(v, (x,y)->sign(y-x))            \\@com reverse sort
   ? vecsort(v, (x,y)->sign(abs(x)-abs(y)))  \\@com sort by increasing absolute value
   ? cmpf(x,y) = my(dx = poldisc(x), dy = poldisc(y)); sign(abs(dx) - abs(dy))
   ? vecsort([x^2+1, x^3-2, x^4+5*x+1], cmpf)
   @eprog\noindent
   The last example used the named \kbd{cmpf} instead of an anonymous function,
   and sorts polynomials with respect to the absolute value of their
   discriminant. A more efficient approach would use precomputations to ensure
   a given discriminant is computed only once:
   \bprog
   ? DISC = vector(#v, i, abs(poldisc(v[i])));
   ? perm = vecsort(vector(#v,i,i), (x,y)->sign(DISC[x]-DISC[y]))
   ? vecextract(v, perm)
   @eprog\noindent Similar ideas apply whenever we sort according to the values
   of a function which is expensive to compute.
   
   \noindent The binary digits of \fl\ mean:
   
   \item 1: indirect sorting of the vector $x$, i.e.~if $x$ is an
   $n$-component vector, returns a permutation of $[1,2,\dots,n]$ which
   applied to the components of $x$ sorts $x$ in increasing order.
   For example, \kbd{vecextract(x, vecsort(x,,1))} is equivalent to
   \kbd{vecsort(x)}.
   
   \item 4: use descending instead of ascending order.
   
   \item 8: remove ``duplicate'' entries with respect to the sorting function
   (keep the first occurring entry).  For example:
   \bprog
     ? vecsort([Pi,Mod(1,2),z], (x,y)->0, 8)   \\@com make everything compare equal
     %1 = [3.141592653589793238462643383]
     ? vecsort([[2,3],[0,1],[0,3]], 2, 8)
     %2 = [[0, 1], [2, 3]]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.vecsort0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def binary(*argv):
  '''
  binary
  Class: basic
  Section: conversions
  C-Name: binaire
  Prototype: G
  Help: binary(x): gives the vector formed by the binary digits of x (x
   integer).
  Doc: 
   outputs the vector of the binary digits of $|x|$.
   Here $x$ can be an integer, a real number (in which case the result has two
   components, one for the integer part, one for the fractional part) or a
   vector/matrix.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.binaire(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bitand(*argv):
  '''
  bitand
  Class: basic
  Section: conversions
  C-Name: gbitand
  Prototype: GG
  Help: bitand(x,y): bitwise "and" of two integers x and y. Negative numbers
   behave as if modulo big power of 2.
  Description: 
   (small, small):small:parens        $(1)&$(2)
   (gen, gen):int        gbitand($1, $2)
  Doc: 
   bitwise \tet{and}
   \sidx{bitwise and}of two integers $x$ and $y$, that is the integer
   $$\sum_i (x_i~\kbd{and}~y_i) 2^i$$
   
   Negative numbers behave $2$-adically, i.e.~the result is the $2$-adic limit
   of \kbd{bitand}$(x_n,y_n)$, where $x_n$ and $y_n$ are non-negative integers
   tending to $x$ and $y$ respectively. (The result is an ordinary integer,
   possibly negative.)
   
   \bprog
   ? bitand(5, 3)
   %1 = 1
   ? bitand(-5, 3)
   %2 = 3
   ? bitand(-5, -3)
   %3 = -7
   @eprog
  Variant: Also available is
   \fun{GEN}{ibitand}{GEN x, GEN y}, which returns the bitwise \emph{and}
   of $|x|$ and $|y|$, two integers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbitand(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bitneg(*argv):
  '''
  bitneg
  Class: basic
  Section: conversions
  C-Name: gbitneg
  Prototype: GD-1,L,
  Help: bitneg(x,{n=-1}): bitwise negation of an integers x truncated to n
   bits. n=-1 means represent infinite sequences of bit 1 as negative numbers.
   Negative numbers behave as if modulo big power of 2.
  Doc: 
   \idx{bitwise negation} of an integer $x$,
   truncated to $n$ bits, $n\geq 0$, that is the integer
   $$\sum_{i=0}^{n-1} \kbd{not}(x_i) 2^i.$$
   The special case $n=-1$ means no truncation: an infinite sequence of
   leading $1$ is then represented as a negative number.
   
   See \secref{se:bitand} for the behavior for negative arguments.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbitneg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bitnegimply(*argv):
  '''
  bitnegimply
  Class: basic
  Section: conversions
  C-Name: gbitnegimply
  Prototype: GG
  Help: bitnegimply(x,y): bitwise "negated imply" of two integers x and y,
   in other words, x BITAND BITNEG(y). Negative numbers behave as if modulo big
   power of 2.
  Description: 
   (small, small):small:parens        $(1)&~$(2)
   (gen, gen):int        gbitnegimply($1, $2)
  Doc: 
   bitwise negated imply of two integers $x$ and
   $y$ (or \kbd{not} $(x \Rightarrow y)$), that is the integer $$\sum
   (x_i~\kbd{and not}(y_i)) 2^i$$
   
   See \secref{se:bitand} for the behavior for negative arguments.
  Variant: Also available is
   \fun{GEN}{ibitnegimply}{GEN x, GEN y}, which returns the bitwise negated
   imply of $|x|$ and $|y|$, two integers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbitnegimply(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bitor(*argv):
  '''
  bitor
  Class: basic
  Section: conversions
  C-Name: gbitor
  Prototype: GG
  Help: bitor(x,y): bitwise "or" of two integers x and y. Negative numbers
   behave as if modulo big power of 2.
  Description: 
   (small, small):small:parens        $(1)|$(2)
   (gen, gen):int        gbitor($1, $2)
  Doc: 
   \sidx{bitwise inclusive or}bitwise (inclusive)
   \tet{or} of two integers $x$ and $y$, that is the integer $$\sum
   (x_i~\kbd{or}~y_i) 2^i$$
   
   See \secref{se:bitand} for the behavior for negative arguments.
  Variant: Also available is
   \fun{GEN}{ibitor}{GEN x, GEN y}, which returns the bitwise \emph{ir}
   of $|x|$ and $|y|$, two integers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbitor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bittest(*argv):
  '''
  bittest
  Class: basic
  Section: conversions
  C-Name: gbittest
  Prototype: GL
  Help: bittest(x,n): gives bit number n (coefficient of 2^n) of the integer x.
   Negative numbers behave as if modulo big power of 2.
  Description: 
   (small, small):bool:parens     ($(1)>>$(2))&1
   (int, small):bool              bittest($1, $2)
   (gen, small):gen               gbittest($1, $2)
  Doc: 
   outputs the $n^{\text{th}}$ bit of $x$ starting
   from the right (i.e.~the coefficient of $2^n$ in the binary expansion of $x$).
   The result is 0 or 1.
   \bprog
   ? bittest(7, 3)
   %1 = 1 \\ the 3rd bit is 1
   ? bittest(7, 4)
   %2 = 0 \\ the 4th bit is 0
   @eprog\noindent
   See \secref{se:bitand} for the behavior at negative arguments.
  Variant: For a \typ{INT} $x$, the variant \fun{long}{bittest}{GEN x, long n} is
   generally easier to use, and if furthermore $n\ge 0$ the low-level function
   \fun{ulong}{int_bit}{GEN x, long n} returns \kbd{bittest(abs(x),n)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbittest(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bitxor(*argv):
  '''
  bitxor
  Class: basic
  Section: conversions
  C-Name: gbitxor
  Prototype: GG
  Help: bitxor(x,y): bitwise "exclusive or" of two integers x and y.
   Negative numbers behave as if modulo big power of 2.
  Description: 
   (small, small):small:parens        $(1)^$(2)
   (gen, gen):int        gbitxor($1, $2)
  Doc: 
   bitwise (exclusive) \tet{or}
   \sidx{bitwise exclusive or}of two integers $x$ and $y$, that is the integer
   $$\sum (x_i~\kbd{xor}~y_i) 2^i$$
   
   See \secref{se:bitand} for the behavior for negative arguments.
  Variant: Also available is
   \fun{GEN}{ibitxor}{GEN x, GEN y}, which returns the bitwise \emph{xor}
   of $|x|$ and $|y|$, two integers.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gbitxor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def hammingweight(*argv):
  '''
  hammingweight
  Class: basic
  Section: conversions
  C-Name: hammingweight
  Prototype: lG
  Help: hammingweight(x): returns the Hamming weight of x.
  Doc: 
   If $x$ is a \typ{INT}, return the binary Hamming weight of $|x|$. Otherwise
   $x$ must be of type \typ{POL}, \typ{VEC}, \typ{COL}, \typ{VECSMALL}, or
   \typ{MAT} and the function returns the number of non-zero coefficients of
   $x$.
   \bprog
   ? hammingweight(15)
   %1 = 4
   ? hammingweight(x^100 + 2*x + 1)
   %2 = 3
   ? hammingweight([Mod(1,2), 2, Mod(0,3)])
   %3 = 2
   ? hammingweight(matid(100))
   %4 = 100
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.hammingweight(*c_arg_tuple)

def quadclassunit(*argv):
  '''
  quadclassunit
  Class: basic
  Section: number_theoretical
  C-Name: quadclassunit0
  Prototype: GD0,L,DGp
  Help: quadclassunit(D,{flag=0},{tech=[]}): compute the structure of the
   class group and the regulator of the quadratic field of discriminant D.
   See manual for the optional technical parameters.
  Doc: \idx{Buchmann-McCurley}'s sub-exponential algorithm for computing the
   class group of a quadratic order of discriminant $D$.
   
   This function should be used instead of \tet{qfbclassno} or \tet{quadregula}
   when $D<-10^{25}$, $D>10^{10}$, or when the \emph{structure} is wanted. It
   is a special case of \tet{bnfinit}, which is slower, but more robust.
   
   The result is a vector $v$ whose components should be accessed using member
   functions:
   
   \item \kbd{$v$.no}: the class number
   
   \item \kbd{$v$.cyc}: a vector giving the structure of the class group as a
   product of cyclic groups;
   
   \item \kbd{$v$.gen}: a vector giving generators of those cyclic groups (as
   binary quadratic forms).
   
   \item \kbd{$v$.reg}: the regulator, computed to an accuracy which is the
   maximum of an internal accuracy determined by the program and the current
   default (note that once the regulator is known to a small accuracy it is
   trivial to compute it to very high accuracy, see the tutorial).
   
   The $\fl$ is obsolete and should be left alone. In older versions,
   it supposedly computed the narrow class group when $D>0$, but this did not
   work at all; use the general function \tet{bnfnarrow}.
   
   Optional parameter \var{tech} is a row vector of the form $[c_1, c_2]$,
   where $c_1 \leq c_2$ are non-negative real numbers which control the execution
   time and the stack size, see \ref{se:GRHbnf}. The parameter is used as a
   threshold to balance the relation finding phase against the final linear
   algebra. Increasing the default $c_1$ means that relations are easier
   to find, but more relations are needed and the linear algebra will be
   harder. The default value for $c_1$ is $0$ and means that it is taken equal
   to $c_2$. The parameter $c_2$ is mostly obsolete and should not be changed,
   but we still document it for completeness: we compute a tentative class
   group by generators and relations using a factorbase of prime ideals
   $\leq c_1 (\log |D|)^2$, then prove that ideals of norm
   $\leq c_2 (\log |D|)^2$ do
   not generate a larger group. By default an optimal $c_2$ is chosen, so that
   the result is provably correct under the GRH --- a famous result of Bach
   states that $c_2 = 6$ is fine, but it is possible to improve on this
   algorithmically. You may provide a smaller $c_2$, it will be ignored
   (we use the provably correct
   one); you may provide a larger $c_2$ than the default value, which results
   in longer computing times for equally correct outputs (under GRH).
  Variant: If you really need to experiment with the \var{tech} parameter, it is
   usually more convenient to use
   \fun{GEN}{Buchquad}{GEN D, double c1, double c2, long prec}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadclassunit0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadhilbert(*argv):
  '''
  quadhilbert
  Class: basic
  Section: number_theoretical
  C-Name: quadhilbert
  Prototype: Gp
  Help: quadhilbert(D): relative equation for the Hilbert class field
   of the quadratic field of discriminant D (which can also be a bnf).
  Doc: relative equation defining the
   \idx{Hilbert class field} of the quadratic field of discriminant $D$.
   
   If $D < 0$, uses complex multiplication (\idx{Schertz}'s variant).
   
   If $D > 0$ \idx{Stark units} are used and (in rare cases) a
   vector of extensions may be returned whose compositum is the requested class
   field. See \kbd{bnrstark} for details.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadhilbert(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def quadray(*argv):
  '''
  quadray
  Class: basic
  Section: number_theoretical
  C-Name: quadray
  Prototype: GGp
  Help: quadray(D,f): relative equation for the ray class field of
   conductor f for the quadratic field of discriminant D (which can also be a
   bnf).
  Doc: relative equation for the ray
   class field of conductor $f$ for the quadratic field of discriminant $D$
   using analytic methods. A \kbd{bnf} for $x^2 - D$ is also accepted in place
   of $D$.
   
   For $D < 0$, uses the $\sigma$ function and Schertz's method.
   
   For $D>0$, uses Stark's conjecture, and a vector of relative equations may be
   returned. See \tet{bnrstark} for more details.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.quadray(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfcompress(*argv):
  '''
  bnfcompress
  Class: basic
  Section: number_fields
  C-Name: bnfcompress
  Prototype: G
  Help: bnfcompress(bnf): converts bnf to a much smaller sbnf, containing the
   same information. Use bnfinit(sbnf) to recover a true bnf.
  Doc: computes a compressed version of \var{bnf} (from \tet{bnfinit}), a
   ``small Buchmann's number field'' (or \var{sbnf} for short) which contains
   enough information to recover a full $\var{bnf}$ vector very rapidly, but
   which is much smaller and hence easy to store and print. Calling
   \kbd{bnfinit} on the result recovers a true \kbd{bnf}, in general different
   from the original. Note that an \tev{snbf} is useless for almost all
   purposes besides storage, and must be converted back to \tev{bnf} form
   before use; for instance, no \kbd{nf*}, \kbd{bnf*} or member function
   accepts them.
   
   An \var{sbnf} is a 12 component vector $v$, as follows. Let \kbd{bnf} be
   the result of a full \kbd{bnfinit}, complete with units. Then $v[1]$ is
   \kbd{bnf.pol}, $v[2]$ is the number of real embeddings \kbd{bnf.sign[1]},
   $v[3]$ is \kbd{bnf.disc}, $v[4]$ is \kbd{bnf.zk}, $v[5]$ is the list of roots
   \kbd{bnf.roots}, $v[7]$ is the matrix $\kbd{W} = \kbd{bnf[1]}$,
   $v[8]$ is the matrix $\kbd{matalpha}=\kbd{bnf[2]}$,
   $v[9]$ is the prime ideal factor base \kbd{bnf[5]} coded in a compact way,
   and ordered according to the permutation \kbd{bnf[6]}, $v[10]$ is the
   2-component vector giving the number of roots of unity and a generator,
   expressed on the integral basis, $v[11]$ is the list of fundamental units,
   expressed on the integral basis, $v[12]$ is a vector containing the algebraic
   numbers alpha corresponding to the columns of the matrix \kbd{matalpha},
   expressed on the integral basis.
   
   All the components are exact (integral or rational), except for the roots in
   $v[5]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfcompress(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfinit(*argv):
  '''
  bnfinit
  Class: basic
  Section: number_fields
  C-Name: bnfinit0
  Prototype: GD0,L,DGp
  Help: bnfinit(P,{flag=0},{tech=[]}): compute the necessary data for future
   use in ideal and unit group computations, including fundamental units if
   they are not too large. flag and tech are both optional. flag can be any of
   0: default, 1: insist on having fundamental units.
   See manual for details about tech.
  Description: 
   (gen):bnf:prec           Buchall($1, 0, prec)
   (gen, 0):bnf:prec        Buchall($1, 0, prec)
   (gen, 1):bnf:prec        Buchall($1, nf_FORCE, prec)
   (gen, ?small, ?gen):bnf:prec        bnfinit0($1, $2, $3, prec)
  Doc: initializes a
   \var{bnf} structure. Used in programs such as \kbd{bnfisprincipal},
   \kbd{bnfisunit} or \kbd{bnfnarrow}. By default, the results are conditional
   on the GRH, see \ref{se:GRHbnf}. The result is a
   10-component vector \var{bnf}.
   
   This implements \idx{Buchmann}'s sub-exponential algorithm for computing the
   class group, the regulator and a system of \idx{fundamental units} of the
   general algebraic number field $K$ defined by the irreducible polynomial $P$
   with integer coefficients.
   
   If the precision becomes insufficient, \kbd{gp} does not strive to compute
   the units by default ($\fl=0$).
   
   When $\fl=1$, we insist on finding the fundamental units exactly. Be
   warned that this can take a very long time when the coefficients of the
   fundamental units on the integral basis are very large. If the fundamental
   units are simply too large to be represented in this form, an error message
   is issued. They could be obtained using the so-called compact representation
   of algebraic numbers as a formal product of algebraic integers. The latter is
   implemented internally but not publicly accessible yet.
   
   $\var{tech}$ is a technical vector (empty by default, see \ref{se:GRHbnf}).
   Careful use of this parameter may speed up your computations,
   but it is mostly obsolete and you should leave it alone.
   
   \smallskip
   
   The components of a \var{bnf} or \var{sbnf} are technical and never used by
   the casual user. In fact: \emph{never access a component directly, always use
   a proper member function.} However, for the sake of completeness and internal
   documentation, their description is as follows. We use the notations
   explained in the book by H. Cohen, \emph{A Course in Computational Algebraic
   Number Theory}, Graduate Texts in Maths \key{138}, Springer-Verlag, 1993,
   Section 6.5, and subsection 6.5.5 in particular.
   
   $\var{bnf}[1]$ contains the matrix $W$, i.e.~the matrix in Hermite normal
   form giving relations for the class group on prime ideal generators
   $(\goth{p}_i)_{1\le i\le r}$.
   
   $\var{bnf}[2]$ contains the matrix $B$, i.e.~the matrix containing the
   expressions of the prime ideal factorbase in terms of the $\goth{p}_i$.
   It is an $r\times c$ matrix.
   
   $\var{bnf}[3]$ contains the complex logarithmic embeddings of the system of
   fundamental units which has been found. It is an $(r_1+r_2)\times(r_1+r_2-1)$
   matrix.
   
   $\var{bnf}[4]$ contains the matrix $M''_C$ of Archimedean components of the
   relations of the matrix $(W|B)$.
   
   $\var{bnf}[5]$ contains the prime factor base, i.e.~the list of prime
   ideals used in finding the relations.
   
   $\var{bnf}[6]$ used to contain a permutation of the prime factor base, but
   has been obsoleted. It contains a dummy $0$.
   
   $\var{bnf}[7]$ or \kbd{\var{bnf}.nf} is equal to the number field data
   $\var{nf}$ as would be given by \kbd{nfinit}.
   
   $\var{bnf}[8]$ is a vector containing the classgroup \kbd{\var{bnf}.clgp}
   as a finite abelian group, the regulator \kbd{\var{bnf}.reg}, a $1$ (used to
   contain an obsolete ``check number''), the number of roots of unity and a
   generator \kbd{\var{bnf}.tu}, the fundamental units \kbd{\var{bnf}.fu}.
   
   $\var{bnf}[9]$ is a 3-element row vector used in \tet{bnfisprincipal} only
   and obtained as follows. Let $D = U W V$ obtained by applying the
   \idx{Smith normal form} algorithm to the matrix $W$ (= $\var{bnf}[1]$) and
   let $U_r$ be the reduction of $U$ modulo $D$. The first elements of the
   factorbase are given (in terms of \kbd{bnf.gen}) by the columns of $U_r$,
   with Archimedean component $g_a$; let also $GD_a$ be the Archimedean
   components of the generators of the (principal) ideals defined by the
   \kbd{bnf.gen[i]\pow bnf.cyc[i]}. Then $\var{bnf}[9]=[U_r, g_a, GD_a]$.
   
   $\var{bnf}[10]$ is by default unused and set equal to 0. This field is used
   to store further information about the field as it becomes available, which
   is rarely needed, hence would be too expensive to compute during the initial
   \kbd{bnfinit} call. For instance, the generators of the principal ideals
   \kbd{bnf.gen[i]\pow bnf.cyc[i]} (during a call to \tet{bnrisprincipal}), or
   those corresponding to the relations in $W$ and $B$ (when the \kbd{bnf}
   internal precision needs to be increased).
  Variant: 
   Also available is \fun{GEN}{Buchall}{GEN P, long flag, long prec},
   corresponding to \kbd{tech = NULL}, where
   \kbd{flag} is either $0$ (default) or \tet{nf_FORCE} (insist on finding
   fundamental units). The function
   \fun{GEN}{Buchall_param}{GEN P, double c1, double c2, long nrpid, long flag, long prec} gives direct access to the technical parameters.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfinit0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfisprincipal(*argv):
  '''
  bnfisprincipal
  Class: basic
  Section: number_fields
  C-Name: bnfisprincipal0
  Prototype: GGD1,L,
  Help: bnfisprincipal(bnf,x,{flag=1}): bnf being output by bnfinit (with
   flag<=2), gives [v,alpha], where v is the vector of exponents on
   the class group generators and alpha is the generator of the resulting
   principal ideal. In particular x is principal if and only if v is the zero
   vector. flag is optional, whose binary digits mean 1: output [v,alpha] (only v
   if unset); 2: increase precision until alpha can be computed (do not insist
   if unset).
  Doc: $\var{bnf}$ being the \sidx{principal ideal}
   number field data output by \kbd{bnfinit}, and $x$ being an ideal, this
   function tests whether the ideal is principal or not. The result is more
   complete than a simple true/false answer and solves general discrete
   logarithm problem. Assume the class group is $\oplus (\Z/d_i\Z)g_i$
   (where the generators $g_i$ and their orders $d_i$ are respectively given by
   \kbd{bnf.gen} and \kbd{bnf.cyc}). The routine returns a row vector $[e,t]$,
   where $e$ is a vector of exponents $0 \leq e_i < d_i$, and $t$ is a number
   field element such that
   $$ x = (t) \prod_i  g_i^{e_i}.$$
   For \emph{given} $g_i$ (i.e. for a given \kbd{bnf}), the $e_i$ are unique,
   and $t$ is unique modulo units.
   
   In particular, $x$ is principal if and only if $e$ is the zero vector. Note
   that the empty vector, which is returned when the class number is $1$, is
   considered to be a zero vector (of dimension $0$).
   \bprog
   ? K = bnfinit(y^2+23);
   ? K.cyc
   %2 = [3]
   ? K.gen
   %3 = [[2, 0; 0, 1]]          \\ a prime ideal above 2
   ? P = idealprimedec(K,3)[1]; \\ a prime ideal above 3
   ? v = bnfisprincipal(K, P)
   %5 = [[2]~, [3/4, 1/4]~]
   ? idealmul(K, v[2], idealfactorback(K, K.gen, v[1]))
   %6 =
   [3 0]
   
   [0 1]
   ? % == idealhnf(K, P)
   %7 = 1
   @eprog
   
   \noindent The binary digits of \fl mean:
   
   \item $1$: If set, outputs $[e,t]$ as explained above, otherwise returns
   only $e$, which is much easier to compute. The following idiom only tests
   whether an ideal is principal:
   \bprog
     is_principal(bnf, x) = !bnfisprincipal(bnf,x,0);
   @eprog
   
   \item $2$: It may not be possible to recover $t$, given the initial accuracy
   to which \kbd{bnf} was computed. In that case, a warning is printed and $t$ is
   set equal to the empty vector \kbd{[]\til}. If this bit is set,
   increase the precision and recompute needed quantities until $t$ can be
   computed. Warning: setting this may induce \emph{very} lengthy computations.
  Variant: Instead of the above hardcoded numerical flags, one should
   rather use an or-ed combination of the symbolic flags \tet{nf_GEN} (include
   generators, possibly a place holder if too difficult) and \tet{nf_FORCE}
   (insist on finding the generators).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfisprincipal0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfisunit(*argv):
  '''
  bnfisunit
  Class: basic
  Section: number_fields
  C-Name: bnfisunit
  Prototype: GG
  Help: bnfisunit(bnf,x): bnf being output by bnfinit, gives
   the column vector of exponents of x on the fundamental units and the roots
   of unity if x is a unit, the empty vector otherwise.
  Doc: \var{bnf} being the number field data
   output by \kbd{bnfinit} and $x$ being an algebraic number (type integer,
   rational or polmod), this outputs the decomposition of $x$ on the fundamental
   units and the roots of unity if $x$ is a unit, the empty vector otherwise.
   More precisely, if $u_1$,\dots,$u_r$ are the fundamental units, and $\zeta$
   is the generator of the group of roots of unity (\kbd{bnf.tu}), the output is
   a vector $[x_1,\dots,x_r,x_{r+1}]$ such that $x=u_1^{x_1}\cdots
   u_r^{x_r}\cdot\zeta^{x_{r+1}}$. The $x_i$ are integers for $i\le r$ and is an
   integer modulo the order of $\zeta$ for $i=r+1$.
   
   Note that \var{bnf} need not contain the fundamental unit explicitly:
   \bprog
   ? setrand(1); bnf = bnfinit(x^2-x-100000);
   ? bnf.fu
     ***   at top-level: bnf.fu
     ***                     ^--
     *** _.fu: missing units in .fu.
   ? u = [119836165644250789990462835950022871665178127611316131167, \
          379554884019013781006303254896369154068336082609238336]~;
   ? bnfisunit(bnf, u)
   %3 = [-1, Mod(0, 2)]~
   @eprog\noindent The given $u$ is the inverse of the fundamental unit
   implicitly stored in \var{bnf}. In this case, the fundamental unit was not
   computed and stored in algebraic form since the default accuracy was too
   low. (Re-run the command at \bs g1 or higher to see such diagnostics.)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfisunit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfsignunit(*argv):
  '''
  bnfsignunit
  Class: basic
  Section: number_fields
  C-Name: signunits
  Prototype: G
  Help: bnfsignunit(bnf): matrix of signs of the real embeddings of the system
   of fundamental units found by bnfinit.
  Doc: $\var{bnf}$ being as output by
   \kbd{bnfinit}, this computes an $r_1\times(r_1+r_2-1)$ matrix having $\pm1$
   components, giving the signs of the real embeddings of the fundamental units.
   The following functions compute generators for the totally positive units:
   
   \bprog
   /* exponents of totally positive units generators on bnf.tufu */
   tpuexpo(bnf)=
   { my(S,d,K);
   
     S = bnfsignunit(bnf); d = matsize(S);
     S = matrix(d[1],d[2], i,j, if (S[i,j] < 0, 1,0));
     S = concat(vectorv(d[1],i,1), S);   \\ add sign(-1)
     K = lift(matker(S * Mod(1,2)));
     if (K, mathnfmodid(K, 2), 2*matid(d[1]))
   }
   
   /* totally positive units */
   tpu(bnf)=
   { my(vu = bnf.tufu, ex = tpuexpo(bnf));
   
     vector(#ex-1, i, factorback(vu, ex[,i+1]))  \\ ex[,1] is 1
   }
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.signunits(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrclassno(*argv):
  '''
  bnrclassno
  Class: basic
  Section: number_fields
  C-Name: bnrclassno0
  Prototype: GDGDG
  Help: bnrclassno(A,{B},{C}): relative degree of the class field defined by
   A,B,C. [A,{B},{C}] is of type [bnr], [bnr,subgroup], [bnf,modulus],
   or [bnf,modulus,subgroup].
   Faster than bnrinit if only the ray class number is wanted.
  Doc: 
    let $A$, $B$, $C$ define a class field $L$ over a ground field $K$
   (of type \kbd{[\var{bnr}]},
   \kbd{[\var{bnr}, \var{subgroup}]},
   or \kbd{[\var{bnf}, \var{modulus}]},
   or \kbd{[\var{bnf}, \var{modulus},\var{subgroup}]},
   \secref{se:CFT}); this function returns the relative degree $[L:K]$.
   
   In particular if $A$ is a \var{bnf} (with units), and $B$ a modulus,
   this function returns the corresponding ray class number modulo $B$.
   One can input the associated \var{bid} (with generators if the subgroup
   $C$ is non trivial) for $B$ instead of the module itself, saving some time.
   
   This function is faster than \kbd{bnrinit} and should be used if only the
   ray class number is desired. See \tet{bnrclassnolist} if you need ray class
   numbers for all moduli less than some bound.
  Variant: Also available is
   \fun{GEN}{bnrclassno}{GEN bnf,GEN f} to compute the ray class number
   modulo~$f$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrclassno0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrclassnolist(*argv):
  '''
  bnrclassnolist
  Class: basic
  Section: number_fields
  C-Name: bnrclassnolist
  Prototype: GG
  Help: bnrclassnolist(bnf,list): if list is as output by ideallist or
   similar, gives list of corresponding ray class numbers.
  Doc: $\var{bnf}$ being as
   output by \kbd{bnfinit}, and \var{list} being a list of moduli (with units) as
   output by \kbd{ideallist} or \kbd{ideallistarch}, outputs the list of the
   class numbers of the corresponding ray class groups. To compute a single
   class number, \tet{bnrclassno} is more efficient.
   
   \bprog
   ? bnf = bnfinit(x^2 - 2);
   ? L = ideallist(bnf, 100, 2);
   ? H = bnrclassnolist(bnf, L);
   ? H[98]
   %4 = [1, 3, 1]
   ? l = L[1][98]; ids = vector(#l, i, l[i].mod[1])
   %5 = [[98, 88; 0, 1], [14, 0; 0, 7], [98, 10; 0, 1]]
   @eprog
   The weird \kbd{l[i].mod[1]}, is the first component of \kbd{l[i].mod}, i.e.
   the finite part of the conductor. (This is cosmetic: since by construction
   the Archimedean part is trivial, I do not want to see it). This tells us that
   the ray class groups modulo the ideals of norm 98 (printed as \kbd{\%5}) have
   respectively order $1$, $3$ and $1$. Indeed, we may check directly:
   \bprog
   ? bnrclassno(bnf, ids[2])
   %6 = 3
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrclassnolist(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrconductor(*argv):
  '''
  bnrconductor
  Class: basic
  Section: number_fields
  C-Name: bnrconductor0
  Prototype: GDGDGD0,L,
  Help: bnrconductor(A,{B},{C},{flag=0}): conductor f of the subfield of
   the ray class field given by A,B,C. flag is optional and
   can be 0: default, 1: returns [f, Cl_f, H], H subgroup of the ray class
   group modulo f defining the extension, 2: returns [f, bnr(f), H].
  Doc: conductor $f$ of the subfield of a ray class field as defined by $[A,B,C]$
   (of type \kbd{[\var{bnr}]},
   \kbd{[\var{bnr}, \var{subgroup}]},
   \kbd{[\var{bnf}, \var{modulus}]} or
   \kbd{[\var{bnf}, \var{modulus}, \var{subgroup}]},
   \secref{se:CFT})
   
   If $\fl = 0$, returns $f$.
   
   If $\fl = 1$, returns $[f, Cl_f, H]$, where $Cl_f$ is the ray class group
   modulo $f$, as a finite abelian group; finally $H$ is the subgroup of $Cl_f$
   defining the extension.
   
   If $\fl = 2$, returns $[f, \var{bnr}(f), H]$, as above except $Cl_f$ is
   replaced by a \kbd{bnr} structure, as output by $\tet{bnrinit}(,f,1)$.
  Variant: 
   Also available is \fun{GEN}{bnrconductor}{GEN bnr, GEN H, long flag}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrconductor0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrconductorofchar(*argv):
  '''
  bnrconductorofchar
  Class: basic
  Section: number_fields
  C-Name: bnrconductorofchar
  Prototype: GG
  Help: bnrconductorofchar(bnr,chi): conductor of the character chi on the ray
   class group bnr.
  Doc: \var{bnr} being a big
   ray number field as output by \kbd{bnrinit}, and \var{chi} being a row vector
   representing a \idx{character} as expressed on the generators of the ray
   class group, gives the conductor of this character as a modulus.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrconductorofchar(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrdisc(*argv):
  '''
  bnrdisc
  Class: basic
  Section: number_fields
  C-Name: bnrdisc0
  Prototype: GDGDGD0,L,
  Help: bnrdisc(A,{B},{C},{flag=0}): absolute or relative [N,R1,discf] of
   the field defined by A,B,C. [A,{B},{C}] is of type [bnr],
   [bnr,subgroup], [bnf, modulus] or [bnf,modulus,subgroup], where bnf is as
   output by bnfinit, bnr by bnrinit, and
   subgroup is the HNF matrix of a subgroup of the corresponding ray class
   group (if omitted, the trivial subgroup). flag is optional whose binary
   digits mean 1: give relative data; 2: return 0 if modulus is not the
   conductor.
  Doc: $A$, $B$, $C$ defining a class field $L$ over a ground field $K$
   (of type \kbd{[\var{bnr}]},
   \kbd{[\var{bnr}, \var{subgroup}]},
   \kbd{[\var{bnf}, \var{modulus}]} or
   \kbd{[\var{bnf}, \var{modulus}, \var{subgroup}]},
   \secref{se:CFT}), outputs data $[N,r_1,D]$ giving the discriminant and
   signature of $L$, depending on the binary digits of \fl:
   
   \item 1: if this bit is unset, output absolute data related to $L/\Q$:
   $N$ is the absolute degree $[L:\Q]$, $r_1$ the number of real places of $L$,
   and $D$ the discriminant of $L/\Q$. Otherwise, output relative data for $L/K$:
   $N$ is the relative degree $[L:K]$, $r_1$ is the number of real places of $K$
   unramified in $L$ (so that the number of real places of $L$ is equal to $r_1$
   times $N$), and $D$ is the relative discriminant ideal of $L/K$.
   
   \item 2: if this bit is set and if the modulus is not the conductor of $L$,
   only return 0.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrdisc0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrdisclist(*argv):
  '''
  bnrdisclist
  Class: basic
  Section: number_fields
  C-Name: bnrdisclist0
  Prototype: GGDG
  Help: bnrdisclist(bnf,bound,{arch}): gives list of discriminants of
   ray class fields of all conductors up to norm bound, in a long vector
   The ramified Archimedean places are given by arch; all possible values are
   taken if arch is omitted. Supports the alternative syntax
   bnrdisclist(bnf,list), where list is as output by ideallist or ideallistarch
   (with units).
  Doc: $\var{bnf}$ being as output by \kbd{bnfinit} (with units), computes a
   list of discriminants of Abelian extensions of the number field by increasing
   modulus norm up to bound \var{bound}. The ramified Archimedean places are
   given by \var{arch}; all possible values are taken if \var{arch} is omitted.
   
   The alternative syntax $\kbd{bnrdisclist}(\var{bnf},\var{list})$ is
   supported, where \var{list} is as output by \kbd{ideallist} or
   \kbd{ideallistarch} (with units), in which case \var{arch} is disregarded.
   
   The output $v$ is a vector of vectors, where $v[i][j]$ is understood to be in
   fact $V[2^{15}(i-1)+j]$ of a unique big vector $V$. (This awkward scheme
   allows for larger vectors than could be otherwise represented.)
   
   $V[k]$ is itself a vector $W$, whose length is the number of ideals of norm
   $k$. We consider first the case where \var{arch} was specified. Each
   component of $W$ corresponds to an ideal $m$ of norm $k$, and
   gives invariants associated to the ray class field $L$ of $\var{bnf}$ of
   conductor $[m, \var{arch}]$. Namely, each contains a vector $[m,d,r,D]$ with
   the following meaning: $m$ is the prime ideal factorization of the modulus,
   $d = [L:\Q]$ is the absolute degree of $L$, $r$ is the number of real places
   of $L$, and $D$ is the factorization of its absolute discriminant. We set $d
   = r = D = 0$ if $m$ is not the finite part of a conductor.
   
   If \var{arch} was omitted, all $t = 2^{r_1}$ possible values are taken and a
   component of $W$ has the form $[m, [[d_1,r_1,D_1], \dots, [d_t,r_t,D_t]]]$,
   where $m$ is the finite part of the conductor as above, and
   $[d_i,r_i,D_i]$ are the invariants of the ray class field of conductor
   $[m,v_i]$, where $v_i$ is the $i$-th Archimedean component, ordered by
   inverse lexicographic order; so $v_1 = [0,\dots,0]$, $v_2 = [1,0\dots,0]$,
   etc. Again, we set $d_i = r_i = D_i = 0$ if $[m,v_i]$ is not a conductor.
   
   Finally, each prime ideal $pr = [p,\alpha,e,f,\beta]$ in the prime
   factorization $m$ is coded as the integer $p\cdot n^2+(f-1)\cdot n+(j-1)$,
   where $n$ is the degree of the base field and $j$ is such that
   
   \kbd{pr = idealprimedec(\var{nf},p)[j]}.
   
   \noindent $m$ can be decoded using \tet{bnfdecodemodule}.
   
   Note that to compute such data for a single field, either \tet{bnrclassno}
   or \tet{bnrdisc} is more efficient.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrdisclist0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrinit(*argv):
  '''
  bnrinit
  Class: basic
  Section: number_fields
  C-Name: bnrinit0
  Prototype: GGD0,L,
  Help: bnrinit(bnf,f,{flag=0}): given a bnf as output by
   bnfinit and a modulus f, initializes data
   linked to the ray class group structure corresponding to this module. flag
   is optional, and can be 0: default, 1: compute also the generators.
  Description: 
   (gen,gen,?small):bnr       bnrinit0($1, $2, $3)
  Doc: $\var{bnf}$ is as
   output by \kbd{bnfinit}, $f$ is a modulus, initializes data linked to
   the ray class group structure corresponding to this module, a so-called
   \var{bnr} structure. One can input the associated \var{bid} with generators
   for $f$ instead of the module itself, saving some time.
   (As in \tet{idealstar}, the finite part of the conductor may be given
   by a factorization into prime ideals, as produced by \tet{idealfactor}.)
   
   The following member functions are available
   on the result: \kbd{.bnf} is the underlying \var{bnf},
   \kbd{.mod} the modulus, \kbd{.bid} the \var{bid} structure associated to the
   modulus; finally, \kbd{.clgp}, \kbd{.no}, \kbd{.cyc}, \kbd{.gen} refer to the
   ray class group (as a finite abelian group), its cardinality, its elementary
   divisors, its generators (only computed if $\fl = 1$).
   
   The last group of functions are different from the members of the underlying
   \var{bnf}, which refer to the class group; use \kbd{\var{bnr}.bnf.\var{xxx}}
   to access these, e.g.~\kbd{\var{bnr}.bnf.cyc} to get the cyclic decomposition
   of the class group.
   
   They are also different from the members of the underlying \var{bid}, which
   refer to $(\Z_K/f)^*$; use \kbd{\var{bnr}.bid.\var{xxx}} to access these,
   e.g.~\kbd{\var{bnr}.bid.no} to get $\phi(f)$.
   
   If $\fl=0$ (default), the generators of the ray class group are not computed,
   which saves time. Hence \kbd{\var{bnr}.gen} would produce an error.
   
   If $\fl=1$, as the default, except that generators are computed.
  Variant: Instead the above  hardcoded  numerical flags,  one should rather use
   \fun{GEN}{Buchray}{GEN bnf, GEN module, long flag}
   where flag is an or-ed combination of \kbd{nf\_GEN} (include generators)
   and \kbd{nf\_INIT} (if omitted, return just the cardinal of the ray class group
   and its structure), possibly 0.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrinit0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrisconductor(*argv):
  '''
  bnrisconductor
  Class: basic
  Section: number_fields
  C-Name: bnrisconductor0
  Prototype: lGDGDG
  Help: bnrisconductor(A,{B},{C}): returns 1 if the modulus is the
   conductor of the subfield of the ray class field given by A,B,C (see
   bnrdisc), and 0 otherwise. Slightly faster than bnrconductor if this is the
   only desired result.
  Doc: $A$, $B$, $C$ represent
   an extension of the base field, given by class field theory
   (see~\secref{se:CFT}). Outputs 1 if this modulus is the conductor, and 0
   otherwise. This is slightly faster than \kbd{bnrconductor}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.bnrisconductor0(*c_arg_tuple)

def bnrisprincipal(*argv):
  '''
  bnrisprincipal
  Class: basic
  Section: number_fields
  C-Name: bnrisprincipal
  Prototype: GGD1,L,
  Help: bnrisprincipal(bnr,x,{flag=1}): bnr being output by bnrinit, gives
   [v,alpha], where v is the vector of exponents on the class group
   generators and alpha is the generator of the resulting principal ideal. In
   particular x is principal if and only if v is the zero vector. If (optional)
   flag is set to 0, output only v.
  Doc: \var{bnr} being the
   number field data which is output by \kbd{bnrinit}$(,,1)$ and $x$ being an
   ideal in any form, outputs the components of $x$ on the ray class group
   generators in a way similar to \kbd{bnfisprincipal}. That is a 2-component
   vector $v$ where $v[1]$ is the vector of components of $x$ on the ray class
   group generators, $v[2]$ gives on the integral basis an element $\alpha$ such
   that $x=\alpha\prod_ig_i^{x_i}$.
   
   If $\fl=0$, outputs only $v_1$. In that case, \var{bnr} need not contain the
   ray class group generators, i.e.~it may be created with \kbd{bnrinit}$(,,0)$
   If $x$ is not coprime to the modulus of \var{bnr} the result is undefined.
  Variant: Instead of hardcoded  numerical flags,  one should rather
   use
   \fun{GEN}{isprincipalray}{GEN bnr, GEN x} for $\kbd{flag} = 0$, and if you
   want generators:
   \bprog
     bnrisprincipal(bnr, x, nf_GEN)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrisprincipal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfnarrow(*argv):
  '''
  bnfnarrow
  Class: basic
  Section: number_fields
  C-Name: buchnarrow
  Prototype: G
  Help: bnfnarrow(bnf): given a big number field as output by bnfinit, gives
   as a 3-component vector the structure of the narrow class group.
  Doc: $\var{bnf}$ being as output by
   \kbd{bnfinit}, computes the narrow class group of $\var{bnf}$. The output is
   a 3-component row vector $v$ analogous to the corresponding class group
   component \kbd{\var{bnf}.clgp} (\kbd{\var{bnf}[8][1]}): the first component
   is the narrow class number \kbd{$v$.no}, the second component is a vector
   containing the SNF\sidx{Smith normal form} cyclic components \kbd{$v$.cyc} of
   the narrow class group, and the third is a vector giving the generators of
   the corresponding \kbd{$v$.gen} cyclic groups. Note that this function is a
   special case of \kbd{bnrinit}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.buchnarrow(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfcertify(*argv):
  '''
  bnfcertify
  Class: basic
  Section: number_fields
  C-Name: bnfcertify0
  Prototype: lGD0,L,
  Help: bnfcertify(bnf,{flag = 0}): certify the correctness (i.e. remove the GRH) of the bnf data output by bnfinit. If flag is present, only certify that the class group is a quotient of the one computed in bnf (much simpler in general).
  Doc: $\var{bnf}$ being as output by
   \kbd{bnfinit}, checks whether the result is correct, i.e.~whether it is
   possible to remove the assumption of the Generalized Riemann
   Hypothesis\sidx{GRH}. It is correct if and only if the answer is 1. If it is
   incorrect, the program may output some error message, or loop indefinitely.
   You can check its progress by increasing the debug level. The \var{bnf}
   structure must contain the fundamental units:
   \bprog
   ? K = bnfinit(x^3+2^2^3+1); bnfcertify(K)
     ***   at top-level: K=bnfinit(x^3+2^2^3+1);bnfcertify(K)
     ***                                        ^-------------
     *** bnfcertify: missing units in bnf.
   ? K = bnfinit(x^3+2^2^3+1, 1); \\ include units
   ? bnfcertify(K)
   %3 = 1
   @eprog
   
   If flag is present, only certify that the class group is a quotient of the
   one computed in bnf (much simpler in general); likewise, the computed units
   may form a subgroup of the full unit group. In this variant, the units are
   no longer needed:
   \bprog
   ? K = bnfinit(x^3+2^2^3+1); bnfcertify(K, 1)
   %4 = 1
   @eprog
  Variant: Also available is  \fun{GEN}{bnfcertify}{GEN bnf} ($\fl=0$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.bnfcertify0(*c_arg_tuple)

def bnfdecodemodule(*argv):
  '''
  bnfdecodemodule
  Class: basic
  Section: number_fields
  C-Name: decodemodule
  Prototype: GG
  Help: bnfdecodemodule(nf,m): given a coded module m as in bnrdisclist,
   gives the true module.
  Doc: if $m$ is a module as output in the
   first component of an extension given by \kbd{bnrdisclist}, outputs the
   true module.
   \bprog
   ? K = bnfinit(x^2+23); L = bnrdisclist(K, 10); s = L[1][2]
   %1 = [[Mat([8, 1]), [[0, 0, 0]]], [Mat([9, 1]), [[0, 0, 0]]]]
   ? bnfdecodemodule(K, s[1][1])
   %2 =
   [2 0]
   
   [0 1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.decodemodule(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfconductor(*argv):
  '''
  rnfconductor
  Class: basic
  Section: number_fields
  C-Name: rnfconductor
  Prototype: GG
  Help: rnfconductor(bnf,pol): conductor of the Abelian extension
   of bnf defined by pol. The result is [conductor,rayclassgroup,subgroup],
   where conductor is the conductor itself, rayclassgroup the structure of the
   corresponding full ray class group, and subgroup the HNF defining the norm
   group (Artin or Takagi group) on the given generators rayclassgroup[3].
  Doc: given $\var{bnf}$
   as output by \kbd{bnfinit}, and \var{pol} a relative polynomial defining an
   \idx{Abelian extension}, computes the class field theory conductor of this
   Abelian extension. The result is a 3-component vector
   $[\var{conductor},\var{rayclgp},\var{subgroup}]$, where \var{conductor} is
   the conductor of the extension given as a 2-component row vector
   $[f_0,f_\infty]$, \var{rayclgp} is the full ray class group corresponding to
   the conductor given as a 3-component vector [h,cyc,gen] as usual for a group,
   and \var{subgroup} is a matrix in HNF defining the subgroup of the ray class
   group on the given generators gen.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfconductor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfisabelian(*argv):
  '''
  rnfisabelian
  Class: basic
  Section: number_fields
  C-Name: rnfisabelian
  Prototype: lGG
  Help: rnfisabelian(nf,T): T being a relative polynomial with coefficients
   in nf, return 1 if it defines an abelian extension, and 0 otherwise.
  Doc: $T$ being a relative polynomial with coefficients
   in \var{nf}, return 1 if it defines an abelian extension, and 0 otherwise.
   \bprog
   ? K = nfinit(y^2 + 23);
   ? rnfisabelian(K, x^3 - 3*x - y)
   %2 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.rnfisabelian(*c_arg_tuple)

def rnfnormgroup(*argv):
  '''
  rnfnormgroup
  Class: basic
  Section: number_fields
  C-Name: rnfnormgroup
  Prototype: GG
  Help: rnfnormgroup(bnr,pol): norm group (or Artin or Takagi group)
   corresponding to the Abelian extension of bnr.bnf defined by pol, where
   the module corresponding to bnr is assumed to be a multiple of the
   conductor. The result is the HNF defining the norm group on the
   generators in bnr.gen.
  Doc: 
   \var{bnr} being a big ray
   class field as output by \kbd{bnrinit} and \var{pol} a relative polynomial
   defining an \idx{Abelian extension}, computes the norm group (alias Artin
   or Takagi group) corresponding to the Abelian extension of
   $\var{bnf}=$\kbd{bnr.bnf}
   defined by \var{pol}, where the module corresponding to \var{bnr} is assumed
   to be a multiple of the conductor (i.e.~\var{pol} defines a subextension of
   bnr). The result is the HNF defining the norm group on the given generators
   of \kbd{bnr.gen}. Note that neither the fact that \var{pol} defines an
   Abelian extension nor the fact that the module is a multiple of the conductor
   is checked. The result is undefined if the assumption is not correct.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfnormgroup(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def subgrouplist(*argv):
  '''
  subgrouplist
  Class: basic
  Section: number_fields
  C-Name: subgrouplist0
  Prototype: GDGD0,L,
  Help: subgrouplist(bnr,{bound},{flag=0}): bnr being as output by bnrinit or
   a list of cyclic components of a finite Abelian group G, outputs the list of
   subgroups of G (of index bounded by bound, if not omitted), given as HNF
   left divisors of the SNF matrix corresponding to G. If flag=0 (default) and
   bnr is as output by bnrinit, gives only the subgroups for which the modulus
   is the conductor.
  Doc: \var{bnr} being as output by \kbd{bnrinit} or a list of cyclic components
   of a finite Abelian group $G$, outputs the list of subgroups of $G$. Subgroups
   are given as HNF left divisors of the SNF matrix corresponding to $G$.
   
   If $\fl=0$ (default) and \var{bnr} is as output by \kbd{bnrinit}, gives
   only the subgroups whose modulus is the conductor. Otherwise, the modulus is
   not taken into account.
   
   If \var{bound} is present, and is a positive integer, restrict the output to
   subgroups of index less than \var{bound}. If \var{bound} is a vector
   containing a single positive integer $B$, then only subgroups of index
   exactly equal to $B$ are computed. For instance
   \bprog
   ? subgrouplist([6,2])
   %1 = [[6, 0; 0, 2], [2, 0; 0, 2], [6, 3; 0, 1], [2, 1; 0, 1], [3, 0; 0, 2],
   [1, 0; 0, 2], [6, 0; 0, 1], [2, 0; 0, 1], [3, 0; 0, 1], [1, 0; 0, 1]]
   ? subgrouplist([6,2],3)    \\@com index less than 3
   %2 = [[2, 1; 0, 1], [1, 0; 0, 2], [2, 0; 0, 1], [3, 0; 0, 1], [1, 0; 0, 1]]
   ? subgrouplist([6,2],[3])  \\@com index 3
   %3 = [[3, 0; 0, 1]]
   ? bnr = bnrinit(bnfinit(x), [120,[1]], 1);
   ? L = subgrouplist(bnr, [8]);
   @eprog\noindent
   In the last example, $L$ corresponds to the 24 subfields of
   $\Q(\zeta_{120})$, of degree $8$ and conductor $120\infty$ (by setting \fl,
   we see there are a total of $43$ subgroups of degree $8$).
   \bprog
   ? vector(#L, i, galoissubcyclo(bnr, L[i]))
   @eprog\noindent
   will produce their equations. (For a general base field, you would
   have to rely on \tet{bnrstark}, or \tet{rnfkummer}.)
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.subgrouplist0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfisnorm(*argv):
  '''
  bnfisnorm
  Class: basic
  Section: number_fields
  C-Name: bnfisnorm
  Prototype: GGD1,L,
  Help: bnfisnorm(bnf,x,{flag=1}): Tries to tell whether x (in Q) is the norm
   of some fractional y (in bnf). Returns a vector [a,b] where x=Norm(a)*b.
   Looks for a solution which is a S-unit, with S a certain list of primes (in
   bnf) containing (among others) all primes dividing x. If bnf is known to be
   Galois, set flag=0 (in this case, x is a norm iff b=1). If flag is non zero
   the program adds to S all the primes: dividing flag if flag<0, or less than
   flag if flag>0. The answer is guaranteed (i.e x norm iff b=1) under GRH, if
   S contains all primes less than 12.log(disc(Bnf))^2, where Bnf is the Galois
   closure of bnf.
  Doc: tries to tell whether the
   rational number $x$ is the norm of some element y in $\var{bnf}$. Returns a
   vector $[a,b]$ where $x=Norm(a)*b$. Looks for a solution which is an $S$-unit,
   with $S$ a certain set of prime ideals containing (among others) all primes
   dividing $x$. If $\var{bnf}$ is known to be \idx{Galois}, set $\fl=0$ (in
   this case, $x$ is a norm iff $b=1$). If $\fl$ is non zero the program adds to
   $S$ the following prime ideals, depending on the sign of $\fl$. If $\fl>0$,
   the ideals of norm less than $\fl$. And if $\fl<0$ the ideals dividing $\fl$.
   
   Assuming \idx{GRH}, the answer is guaranteed (i.e.~$x$ is a norm iff $b=1$),
   if $S$ contains all primes less than $12\log(\disc(\var{Bnf}))^2$, where
   $\var{Bnf}$ is the Galois closure of $\var{bnf}$.
   
   See also \tet{bnfisintnorm}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfisnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfisnorm(*argv):
  '''
  rnfisnorm
  Class: basic
  Section: number_fields
  C-Name: rnfisnorm
  Prototype: GGD0,L,
  Help: rnfisnorm(T,a,{flag=0}): T is as output by rnfisnorminit applied to
   L/K. Tries to tell whether a is a norm from L/K. Returns a vector [x,q]
   where a=Norm(x)*q. Looks for a solution which is a S-integer, with S a list
   of places in K containing the ramified primes, generators of the class group
   of ext, as well as those primes dividing a. If L/K is Galois, omit flag,
   otherwise it is used to add more places to S: all the places above the
   primes p <= flag (resp. p | flag) if flag > 0 (resp. flag < 0). The answer
   is guaranteed (i.e a is a norm iff q=1) if L/K is Galois or, under GRH, if S
   contains all primes less than 12.log(disc(M))^2, where M is the normal
   closure of L/K.
  Doc: similar to
   \kbd{bnfisnorm} but in the relative case. $T$ is as output by
   \tet{rnfisnorminit} applied to the extension $L/K$. This tries to decide
   whether the element $a$ in $K$ is the norm of some $x$ in the extension
   $L/K$.
   
   The output is a vector $[x,q]$, where $a = \Norm(x)*q$. The
   algorithm looks for a solution $x$ which is an $S$-integer, with $S$ a list
   of places of $K$ containing at least the ramified primes, the generators of
   the class group of $L$, as well as those primes dividing $a$. If $L/K$ is
   Galois, then this is enough; otherwise, $\fl$ is used to add more primes to
   $S$: all the places above the primes $p \leq \fl$ (resp.~$p|\fl$) if $\fl>0$
   (resp.~$\fl<0$).
   
   The answer is guaranteed (i.e.~$a$ is a norm iff $q = 1$) if the field is
   Galois, or, under \idx{GRH}, if $S$ contains all primes less than
   $12\log^2\left|\disc(M)\right|$, where $M$ is the normal
   closure of $L/K$.
   
   If \tet{rnfisnorminit} has determined (or was told) that $L/K$ is
   \idx{Galois}, and $\fl \neq 0$, a Warning is issued (so that you can set
   $\fl = 1$ to check whether $L/K$ is known to be Galois, according to $T$).
   Example:
   
   \bprog
   bnf = bnfinit(y^3 + y^2 - 2*y - 1);
   p = x^2 + Mod(y^2 + 2*y + 1, bnf.pol);
   T = rnfisnorminit(bnf, p);
   rnfisnorm(T, 17)
   @eprog\noindent
   checks whether $17$ is a norm in the Galois extension $\Q(\beta) /
   \Q(\alpha)$, where $\alpha^3 + \alpha^2 - 2\alpha - 1 = 0$ and $\beta^2 +
   \alpha^2 + 2\alpha + 1 = 0$ (it is).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfisnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfisnorminit(*argv):
  '''
  rnfisnorminit
  Class: basic
  Section: number_fields
  C-Name: rnfisnorminit
  Prototype: GGD2,L,
  Help: rnfisnorminit(pol,polrel,{flag=2}): let K be defined by a root of pol,
   L/K the extension defined by polrel. Compute technical data needed by
   rnfisnorm to solve norm equations Nx = a, for x in L, and a in K. If flag=0,
   do not care whether L/K is Galois or not; if flag = 1, assume L/K is Galois;
   if flag = 2, determine whether L/K is Galois.
  Doc: let $K$ be defined by a root of \var{pol}, and $L/K$ the extension defined
   by the polynomial \var{polrel}. As usual, \var{pol} can in fact be an \var{nf},
   or \var{bnf}, etc; if \var{pol} has degree $1$ (the base field is $\Q$),
   polrel is also allowed to be an \var{nf}, etc. Computes technical data needed
   by \tet{rnfisnorm} to solve norm equations $Nx = a$, for $x$ in $L$, and $a$
   in $K$.
   
   If $\fl = 0$, do not care whether $L/K$ is Galois or not.
   
   If $\fl = 1$, $L/K$ is assumed to be Galois (unchecked), which speeds up
   \tet{rnfisnorm}.
   
   If $\fl = 2$, let the routine determine whether $L/K$ is Galois.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfisnorminit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfissunit(*argv):
  '''
  bnfissunit
  Class: basic
  Section: number_fields
  C-Name: bnfissunit
  Prototype: GGG
  Help: bnfissunit(bnf,sfu,x): bnf being output by bnfinit (with flag<=2), sfu
   by bnfsunit, gives the column vector of exponents of x on the fundamental
   S-units and the roots of unity if x is a unit, the empty vector otherwise.
  Doc: $\var{bnf}$ being output by
   \kbd{bnfinit}, \var{sfu} by \kbd{bnfsunit}, gives the column vector of
   exponents of $x$ on the fundamental $S$-units and the roots of unity.
   If $x$ is not a unit, outputs an empty vector.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfissunit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfsunit(*argv):
  '''
  bnfsunit
  Class: basic
  Section: number_fields
  C-Name: bnfsunit
  Prototype: GGp
  Help: bnfsunit(bnf,S): compute the fundamental S-units of the number field
   bnf output by bnfinit, S being a list of prime ideals. res[1] contains the
   S-units, res[5] the S-classgroup. See manual for details.
  Doc: computes the fundamental $S$-units of the
   number field $\var{bnf}$ (output by \kbd{bnfinit}), where $S$ is a list of
   prime ideals (output by \kbd{idealprimedec}). The output is a vector $v$ with
   6 components.
   
   $v[1]$ gives a minimal system of (integral) generators of the $S$-unit group
   modulo the unit group.
   
   $v[2]$ contains technical data needed by \kbd{bnfissunit}.
   
   $v[3]$ is an empty vector (used to give the logarithmic embeddings of the
   generators in $v[1]$ in version 2.0.16).
   
   $v[4]$ is the $S$-regulator (this is the product of the regulator, the
   determinant of $v[2]$ and the natural logarithms of the norms of the ideals
   in $S$).
   
   $v[5]$ gives the $S$-class group structure, in the usual format
   (a row vector whose three components give in order the $S$-class number,
   the cyclic components and the generators).
   
   $v[6]$ is a copy of $S$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfsunit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfhilbert(*argv):
  '''
  nfhilbert
  Class: basic
  Section: number_fields
  C-Name: nfhilbert0
  Prototype: lGGGDG
  Help: nfhilbert(nf,a,b,{pr}): if pr is omitted, global Hilbert symbol (a,b) in
   nf, that is 1 if X^2-aY^2-bZ^2 has a non-trivial solution (X,Y,Z) in nf, -1
   otherwise. Otherwise compute the local symbol modulo the prime ideal pr.
  Doc: if \var{pr} is omitted,
   compute the global quadratic \idx{Hilbert symbol} $(a,b)$ in $\var{nf}$, that
   is $1$ if $x^2 - a y^2 - b z^2$ has a non trivial solution $(x,y,z)$ in
   $\var{nf}$, and $-1$ otherwise. Otherwise compute the local symbol modulo
   the prime ideal \var{pr}, as output by \kbd{idealprimedec}.
  Variant: 
   Also available is \fun{long}{nfhilbert}{GEN bnf,GEN a,GEN b} (global
   quadratic Hilbert symbol).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.nfhilbert0(*c_arg_tuple)

def concat(*argv):
  '''
  concat
  Class: basic
  Section: linear_algebra
  C-Name: concat
  Prototype: GDG
  Help: concat(x,{y}): concatenation of x and y, which can be scalars, vectors
   or matrices, or lists (in this last case, both x and y have to be lists). If
   y is omitted, x has to be a list or row vector and its elements are
   concatenated.
  Description: 
   (mp,mp):vec           concat($1, $2)
   (vec,mp):vec          concat($1, $2)
   (mp,vec):vec          concat($1, $2)
   (vec,vec):vec         concat($1, $2)
   (list,list):list      concat($1, $2)
   (genstr,gen):genstr   concat($1, $2)
   (gen,genstr):genstr   concat($1, $2)
   (gen,?gen):gen        concat($1, $2)
  Doc: concatenation of $x$ and $y$. If $x$ or $y$ is
   not a vector or matrix, it is considered as a one-dimensional vector. All
   types are allowed for $x$ and $y$, but the sizes must be compatible. Note
   that matrices are concatenated horizontally, i.e.~the number of rows stays
   the same. Using transpositions, one can concatenate them vertically,
   but it is often simpler to use \tet{matconcat}.
   \bprog
   ? x = matid(2); y = 2*matid(2);
   ? concat(x,y)
   %2 =
   [1 0 2 0]
   
   [0 1 0 2]
   ? concat(x~,y~)~
   %3 =
   [1 0]
   
   [0 1]
   
   [2 0]
   
   [0 2]
   ? matconcat([x;y])
   %4 =
   [1 0]
   
   [0 1]
   
   [2 0]
   
   [0 2]
   @eprog\noindent
   To concatenate vectors sideways (i.e.~to obtain a two-row or two-column
   matrix), use \tet{Mat} instead, or \tet{matconcat}:
   \bprog
   ? x = [1,2];
   ? y = [3,4];
   ? concat(x,y)
   %3 = [1, 2, 3, 4]
   
   ? Mat([x,y]~)
   %4 =
   [1 2]
   
   [3 4]
   ? matconcat([x;y])
   %5 =
   [1 2]
   
   [3 4]
   @eprog
   Concatenating a row vector to a matrix having the same number of columns will
   add the row to the matrix (top row if the vector is $x$, i.e.~comes first, and
   bottom row otherwise).
   
   The empty matrix \kbd{[;]} is considered to have a number of rows compatible
   with any operation, in particular concatenation. (Note that this is
   \emph{not} the case for empty vectors \kbd{[~]} or \kbd{[~]\til}.)
   
   If $y$ is omitted, $x$ has to be a row vector or a list, in which case its
   elements are concatenated, from left to right, using the above rules.
   \bprog
   ? concat([1,2], [3,4])
   %1 = [1, 2, 3, 4]
   ? a = [[1,2]~, [3,4]~]; concat(a)
   %2 =
   [1 3]
   
   [2 4]
   
   ? concat([1,2; 3,4], [5,6]~)
   %3 =
   [1 2 5]
   
   [3 4 6]
   ? concat([%, [7,8]~, [1,2,3,4]])
   %5 =
   [1 2 5 7]
   
   [3 4 6 8]
   
   [1 2 3 4]
   @eprog
  Variant: \fun{GEN}{concat1}{GEN x} is a shortcut for \kbd{concat(x,NULL)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.concat(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matconcat(*argv):
  '''
  matconcat
  Class: basic
  Section: linear_algebra
  C-Name: matconcat
  Prototype: G
  Help: matconcat(v): concatenate the entries of v and return the resulting matrix
  Doc: returns a \typ{MAT} built from the entries of $v$, which may
   be a \typ{VEC} (concatenate horizontally), a \typ{COL} (concatenate
   vertically), or a \typ{MAT} (concatenate vertically each column, and
   concatenate vertically the resulting matrices). The entries of $v$ are always
   considered as matrices: they can themselves be \typ{VEC} (seen as a row
   matrix), a \typ{COL} seen as a column matrix), a \typ{MAT}, or a scalar (seen
   as an $1 \times 1$ matrix).
   \bprog
   ? A=[1,2;3,4]; B=[5,6]~; C=[7,8]; D=9;
   ? matconcat([A, B]) \\ horizontal
   %1 =
   [1 2 5]
   
   [3 4 6]
   ? matconcat([A, C]~) \\ vertical
   %2 =
   [1 2]
   
   [3 4]
   
   [7 8]
   ? matconcat([A, B; C, D]) \\ block matrix
   %3 =
   [1 2 5]
   
   [3 4 6]
   
   [7 8 9]
   @eprog\noindent
   If the dimensions of the entries to concatenate do not match up, the above
   rules are extended as follows:
   
   \item each entry $v_{i,j}$ of $v$ has a natural length and height: $1 \times
   1$ for a scalar, $1 \times n$ for a \typ{VEC} of length $n$, $n \times 1$
   for a \typ{COL}, $m \times n$ for an $m\times n$ \typ{MAT}
   
   \item let $H_i$ be the maximum over $j$ of the lengths of the $v_{i,j}$,
   let $L_j$ be the maximum over $i$ of the heights of the $v_{i,j}$.
   The dimensions of the $(i,j)$-th block in the concatenated matrix are
   $H_i \times L_j$.
   
   \item a scalar $s = v_{i,j}$ is considered as $s$ times an identity matrix
   of the block dimension $\min (H_i,L_j)$
   
   \item blocks are extended by 0 columns on the right and 0 rows at the
   bottom, as needed.
   
   \bprog
   ? matconcat([1, [2,3]~, [4,5,6]~]) \\ horizontal
   %4 =
   [1 2 4]
   
   [0 3 5]
   
   [0 0 6]
   ? matconcat([1, [2,3], [4,5,6]]~) \\ vertical
   %5 =
   [1 0 0]
   
   [2 3 0]
   
   [4 5 6]
   ? matconcat([B, C; A, D]) \\ block matrix
   %6 =
   [5 0 7 8]
   
   [6 0 0 0]
   
   [1 2 9 0]
   
   [3 4 0 9]
   ? U=[1,2;3,4]; V=[1,2,3;4,5,6;7,8,9];
   ? matconcat(matdiagonal([U, V])) \\ block diagonal
   %7 =
   [1 2 0 0 0]
   
   [3 4 0 0 0]
   
   [0 0 1 2 3]
   
   [0 0 4 5 6]
   
   [0 0 7 8 9]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matconcat(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def default(*argv):
  '''
  default
  Class: basic
  Section: programming/specific
  C-Name: default0
  Prototype: DrDs
  Help: default({key},{val}): returns the current value of the
   default key. If val is present, set opt to val first. If no argument is
   given, print a list of all defaults as well as their values.
  Description: 
   ("realprecision"):small:prec              getrealprecision()
   ("realprecision",small):small:prec        setrealprecision($2, &prec)
   ("seriesprecision"):small                 precdl
   ("seriesprecision",small):small:parens    precdl = $2
   ("debug"):small                           DEBUGLEVEL
   ("debug",small):small:parens              DEBUGLEVEL = $2
   ("debugmem"):small                        DEBUGMEM
   ("debugmem",small):small:parens           DEBUGMEM = $2
   ("debugfiles"):small                      DEBUGFILES
   ("debugfiles",small):small:parens         DEBUGFILES = $2
   ("factor_add_primes"):small               factor_add_primes
   ("factor_add_primes",small):small         factor_add_primes = $2
   ("factor_proven"):small                   factor_proven
   ("factor_proven",small):small             factor_proven = $2
   ("new_galois_format"):small               new_galois_format
   ("new_galois_format",small):small         new_galois_format = $2
  Doc: returns the default corresponding to keyword \var{key}. If \var{val} is
   present, sets the default to \var{val} first (which is subject to string
   expansion first). Typing \kbd{default()} (or \b{d}) yields the complete
   default list as well as their current values. See \secref{se:defaults} for an
   introduction to GP defaults, \secref{se:gp_defaults} for a
   list of available defaults, and \secref{se:meta} for some shortcut
   alternatives. Note that the shortcuts are meant for interactive use and
   usually display more information than \kbd{default}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.default0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellanalyticrank(*argv):
  '''
  ellanalyticrank
  Class: basic
  Section: elliptic_curves
  C-Name: ellanalyticrank
  Prototype: GDGp
  Help: ellanalyticrank(e, {eps}): returns the order of vanishing at s=1
   of the L-function of the elliptic curve e and the value of the first
   non-zero derivative. To determine this order, it is assumed that any
   value less than eps is zero. If no value of eps is given, a value of
   half the current precision is used.
  Doc: returns the order of vanishing at $s=1$ of the $L$-function of the
   elliptic curve $e$ and the value of the first non-zero derivative. To
   determine this order, it is assumed that any value less than \kbd{eps} is
   zero. If no value of \kbd{eps} is given, a value of half the current
   precision is used.
   \bprog
   ? e = ellinit("11a1"); \\ rank 0
   ? ellanalyticrank(e)
   %2 = [0, 0.2538418608559106843377589233]
   ? e = ellinit("37a1"); \\ rank 1
   ? ellanalyticrank(e)
   %4 = [1, 0.3059997738340523018204836835]
   ? e = ellinit("389a1"); \\ rank 2
   ? ellanalyticrank(e)
   %6 = [2, 1.518633000576853540460385214]
   ? e = ellinit("5077a1"); \\ rank 3
   ? ellanalyticrank(e)
   %8 = [3, 10.39109940071580413875185035]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellanalyticrank(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellL1(*argv):
  '''
  ellL1
  Class: basic
  Section: elliptic_curves
  C-Name: ellL1
  Prototype: GLp
  Help: ellL1(e, r): returns the value at s=1 of the derivative of order r of the L-function of the elliptic curve e assuming that r is at most the order of vanishing of the function at s=1.
  Doc: returns the value at $s=1$ of the derivative of order $r$ of the
   $L$-function of the elliptic curve $e$ assuming that $r$ is at most the order
   of vanishing of the $L$-function at $s=1$. (The result is wrong if $r$ is
   strictly larger than the order of vanishing at 1.)
   \bprog
   ? e = ellinit("11a1"); \\ order of vanishing is 0
   ? ellL1(e, 0)
   %2 = 0.2538418608559106843377589233
   ? e = ellinit("389a1");  \\ order of vanishing is 2
   ? ellL1(e, 0)
   %4 = -5.384067311837218089235032414 E-29
   ? ellL1(e, 1)
   %5 = 0
   ? ellL1(e, 2)
   %6 = 1.518633000576853540460385214
   @eprog\noindent
   The main use of this function, after computing at \emph{low} accuracy the
   order of vanishing using \tet{ellanalyticrank}, is to compute the
   leading term at \emph{high} accuracy to check (or use) the Birch and
   Swinnerton-Dyer conjecture:
   \bprog
   ? \p18
     realprecision = 18 significant digits
   ? ellanalyticrank(ellinit([0, 0, 1, -7, 6]))
   time = 32 ms.
   %1 = [3, 10.3910994007158041]
   ? \p200
     realprecision = 202 significant digits (200 digits displayed)
   ? ellL1(e, 3)
   time = 23,113 ms.
   %3 = 10.3910994007158041387518505103609170697263563756570092797@com$[\dots]$
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellL1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellconvertname(*argv):
  '''
  ellconvertname
  Class: basic
  Section: elliptic_curves
  C-Name: ellconvertname
  Prototype: G
  Help: ellconvertname(name): convert an elliptic curve name (as found in
   the elldata database) from a string to a triplet [conductor, isogeny class,
   index]. It will also convert a triplet back to a curve name.
  Doc: 
   converts an elliptic curve name, as found in the \tet{elldata} database,
   from a string to a triplet $[\var{conductor}, \var{isogeny class},
   \var{index}]$. It will also convert a triplet back to a curve name.
   Examples:
   \bprog
   ? ellconvertname("123b1")
   %1 = [123, 1, 1]
   ? ellconvertname(%)
   %2 = "123b1"
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellconvertname(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellidentify(*argv):
  '''
  ellidentify
  Class: basic
  Section: elliptic_curves
  C-Name: ellidentify
  Prototype: G
  Help: ellidentify(E): look up the elliptic curve E in the elldata database and
   return [[N, M, ...], C] where N is the name of the curve in Cremona's
   database, M the minimal model and C the coordinates change (see
   ellchangecurve).
  Doc: look up the elliptic curve $E$, defined by an arbitrary model over $\Q$,
   in the \tet{elldata} database.
   Return \kbd{[[N, M, G], C]}  where $N$ is the curve name in Cremona's
   elliptic curve database, $M$ is the minimal model, $G$ is a $\Z$-basis of
   the free part of the \idx{Mordell-Weil group} $E(\Q)$ and $C$ is the
   change of coordinates change, suitable for \kbd{ellchangecurve}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellidentify(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellsearch(*argv):
  '''
  ellsearch
  Class: basic
  Section: elliptic_curves
  C-Name: ellsearch
  Prototype: G
  Help: ellsearch(N): returns all curves in the elldata database matching
   constraint N:  given name (N = "11a1" or [11,0,1]),
   given isogeny class (N = "11a" or [11,0]), or
   given conductor (N = 11, "11", or [11]).
  Doc: This function finds all curves in the \tet{elldata} database satisfying
   the constraint defined by the argument $N$:
   
   \item if $N$ is a character string, it selects a given curve, e.g.
   \kbd{"11a1"}, or curves in the given isogeny class, e.g. \kbd{"11a"}, or
   curves with given condutor, e.g. \kbd{"11"};
   
   \item if $N$ is a vector of integers, it encodes the same constraints
   as the character string above, according to the \tet{ellconvertname}
   correspondance, e.g. \kbd{[11,0,1]} for \kbd{"11a1"}, \kbd{[11,0]} for
   \kbd{"11a"} and \kbd{[11]} for \kbd{"11"};
   
   \item if $N$ is an integer, curves with conductor $N$ are selected.
   
   If $N$ is a full curve name, e.g. \kbd{"11a1"} or \kbd{[11,0,1]},
   the output format is $[N, [a_1,a_2,a_3,a_4,a_6], G]$ where
   $[a_1,a_2,a_3,a_4,a_6]$ are the coefficients of the Weierstrass equation of
   the curve and $G$ is a $\Z$-basis of the free part of the \idx{Mordell-Weil
   group} associated to the curve.
   \bprog
   ? ellsearch("11a3")
   %1 = ["11a3", [0, -1, 1, 0, 0], []]
   ? ellsearch([11,0,3])
   %2 = ["11a3", [0, -1, 1, 0, 0], []]
   @eprog\noindent
   
   If $N$ is not a full curve name, then the output is a vector of all matching
   curves in the above format:
   \bprog
   ? ellsearch("11a")
   %1 = [["11a1", [0, -1, 1, -10, -20], []],
         ["11a2", [0, -1, 1, -7820, -263580], []],
         ["11a3", [0, -1, 1, 0, 0], []]]
   ? ellsearch("11b")
   %2 = []
   @eprog
  Variant: Also available is \fun{GEN}{ellsearchcurve}{GEN N} that only
   accepts complete curve names (as \typ{STR}).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellsearch(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellak(*argv):
  '''
  ellak
  Class: basic
  Section: elliptic_curves
  C-Name: akell
  Prototype: GG
  Help: ellak(E,n): computes the n-th Fourier coefficient of the L-function of
   the elliptic curve E (assumed E is an integral model).
  Doc: 
   computes the coefficient $a_n$ of the $L$-function of the elliptic curve
   $E/\Q$, i.e.~coefficients of a newform of weight 2 by the modularity theorem
   (\idx{Taniyama-Shimura-Weil conjecture}). $E$ must be an \var{ell} structure
   over $\Q$ as output by \kbd{ellinit}. $E$ must be given by an integral model,
   not necessarily minimal, although a minimal model will make the function
   faster.
   \bprog
   ? E = ellinit([0,1]);
   ? ellak(E, 10)
   %2 = 0
   ? e = ellinit([5^4,5^6]); \\ not minimal at 5
   ? ellak(e, 5) \\ wasteful but works
   %3 = -3
   ? E = ellminimalmodel(e); \\ now minimal
   ? ellak(E, 5)
   %5 = -3
   @eprog\noindent If the model is not minimal at a number of bad primes, then
   the function will be slower on those $n$ divisible by the bad primes.
   The speed should be comparable for other $n$:
   \bprog
   ? for(i=1,10^6, ellak(E,5))
   time = 820 ms.
   ? for(i=1,10^6, ellak(e,5)) \\ 5 is bad, markedly slower
   time = 1,249 ms.
   
   ? for(i=1,10^5,ellak(E,5*i))
   time = 977 ms.
   ? for(i=1,10^5,ellak(e,5*i)) \\ still slower but not so much on average
   time = 1,008 ms.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.akell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellan(*argv):
  '''
  ellan
  Class: basic
  Section: elliptic_curves
  C-Name: anell
  Prototype: GL
  Help: ellan(E,n): computes the first n Fourier coefficients of the
   L-function of the elliptic curve E (n<2^24 on a 32-bit machine).
  Doc: computes the vector of the first $n$ Fourier coefficients $a_k$
   corresponding to the elliptic curve $E$. The curve must be given by an
   integral model, not necessarily minimal, although a minimal model will make
   the function faster.
  Variant: Also available is \fun{GEN}{anellsmall}{GEN e, long n}, which
   returns a \typ{VECSMALL} instead of a \typ{VEC}, saving on memory.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.anell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellbil(*argv):
  '''
  ellbil
  Class: basic
  Section: elliptic_curves
  C-Name: bilhell
  Prototype: GGGp
  Help: ellbil(E,z1,z2): canonical bilinear form for the points z1,z2 on the
   elliptic curve E (assumed to be minimal). Either z1 or z2 can also be a
   vector/matrix of points.
  Doc: 
   if $z1$ and $z2$ are points on the elliptic
   curve $E$, assumed to be integral given by a minimal model, this function
   computes the value of the canonical bilinear form on $z1$, $z2$:
   $$ ( h(E,z1\kbd{+}z2) - h(E,z1) - h(E,z2) ) / 2 $$
   where \kbd{+} denotes of course addition on $E$. In addition, $z1$ or $z2$
   (but not both) can be vectors or matrices.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bilhell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elladd(*argv):
  '''
  elladd
  Class: basic
  Section: elliptic_curves
  C-Name: elladd
  Prototype: GGG
  Help: elladd(E,z1,z2): sum of the points z1 and z2 on elliptic curve E.
  Doc: 
   sum of the points $z1$ and $z2$ on the
   elliptic curve corresponding to $E$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elladd(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellap(*argv):
  '''
  ellap
  Class: basic
  Section: elliptic_curves
  C-Name: ellap
  Prototype: GDG
  Help: ellap(E,{p}): computes the trace of Frobenius a_p for the elliptic
   curve E, defined over Q or a finite field.
  Doc: 
   Let $E$ be an \var{ell} structure as output by \kbd{ellinit}, defined over
   $\Q$ or a finite field $\F_q$. The argument $p$ is best left omitted if the
   curve is defined over a finite field, and must be a prime number otherwise.
   This function computes the trace of Frobenius $t$ for the elliptic curve $E$,
   defined by the equation $\#E(\F_q) = q+1 - t$.
   
   If the curve is defined over $\Q$, $p$ must be explicitly given and the
   function computes the trace of the reduction over $\F_p$.
   The trace of Frobenius is also the $a_p$ coefficient in the curve $L$-series
   $L(E,s) = \sum_n a_n n^{-s}$, whence the function name. The equation must be
   integral at $p$ but need not be minimal at $p$; of course, a minimal model
   will be more efficient.
   \bprog
   ? E = ellinit([0,1]);  \\ y^2 = x^3 + 0.x + 1, defined over Q
   ? ellap(E, 7) \\ 7 necessary here
   %2 = -4       \\ #E(F_7) = 7+1-(-4) = 12
   ? ellcard(E, 7)
   %3 = 12       \\ OK
   
   ? E = ellinit([0,1], 11);  \\ defined over F_11
   ? ellap(E)       \\ no need to repeat 11
   %4 = 0
   ? ellap(E, 11)   \\ ... but it also works
   %5 = 0
   ? ellgroup(E, 13) \\ ouch, inconsistent input!
      ***   at top-level: ellap(E,13)
      ***                 ^-----------
      *** ellap: inconsistent moduli in Rg_to_Fp:
        11
        13
   
   ? Fq = ffgen(ffinit(11,3), 'a); \\ defines F_q := F_{11^3}
   ? E = ellinit([a+1,a], Fq);  \\ y^2 = x^3 + (a+1)x + a, defined over F_q
   ? ellap(E)
   %8 = -3
   @eprog
   
   \misctitle{Algorithms used} If $E/\F_q$ has CM by a principal imaginary
   quadratic order we use a fast explicit formula (involving essentially Kronecker
   symbols and Cornacchia's algorithm), in $O(\log q)^2$.
   Otherwise, we use Shanks-Mestre's baby-step/giant-step method, which runs in
   time $q(p^{1/4})$ using $O(q^{1/4})$ storage, hence becomes unreasonable when
   $q$ has about 30~digits. If the \tet{seadata} package is installed, the
   \tet{SEA} algorithm becomes available, heuristically in $\tilde{O}(\log
   q)^4$, and primes of the order of 200~digits become feasible. In very small
   characteristic (2,3,5,7 or $13$), we use Harley's algorithm.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellap(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellcard(*argv):
  '''
  ellcard
  Class: basic
  Section: elliptic_curves
  C-Name: ellcard
  Prototype: GDG
  Help: ellcard(E,{p}): computes the order of the group E(Fp)
   for the elliptic curve E, defined over Q or a finite field.
  Doc: Let $E$ be an \var{ell} structure as output by \kbd{ellinit}, defined over
   $\Q$ or a finite field $\F_q$. The argument $p$ is best left omitted if the
   curve is defined over a finite field, and must be a prime number otherwise.
   This function computes the order of the group $E(\F_q)$ (as would be
   computed by \tet{ellgroup}).
   
   If the curve is defined over $\Q$, $p$ must be explicitly given and the
   function computes the cardinal of the reduction over $\F_p$; the
   equation need not be minimal at $p$, but a minimal model will be more
   efficient. The reduction is allowed to be singular, and we return the order
   of the group of non-singular points in this case.
  Variant: Also available is \fun{GEN}{ellcard}{GEN E, GEN p} where $p$ is not
   \kbd{NULL}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellcard(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellchangecurve(*argv):
  '''
  ellchangecurve
  Class: basic
  Section: elliptic_curves
  C-Name: ellchangecurve
  Prototype: GG
  Help: ellchangecurve(E,v): change data on elliptic curve according to
   v=[u,r,s,t].
  Description: 
   (gen, gen):ell        ellchangecurve($1, $2)
  Doc: 
   changes the data for the elliptic curve $E$
   by changing the coordinates using the vector \kbd{v=[u,r,s,t]}, i.e.~if $x'$
   and $y'$ are the new coordinates, then $x=u^2x'+r$, $y=u^3y'+su^2x'+t$.
   $E$ must be an \var{ell} structure as output by \kbd{ellinit}. The special
   case $v = 1$ is also used instead of $[1,0,0,0]$ to denote the
   trivial coordinate change.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellchangecurve(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellchangepoint(*argv):
  '''
  ellchangepoint
  Class: basic
  Section: elliptic_curves
  C-Name: ellchangepoint
  Prototype: GG
  Help: ellchangepoint(x,v): change data on point or vector of points x on an
   elliptic curve according to v=[u,r,s,t].
  Doc: 
   changes the coordinates of the point or
   vector of points $x$ using the vector \kbd{v=[u,r,s,t]}, i.e.~if $x'$ and
   $y'$ are the new coordinates, then $x=u^2x'+r$, $y=u^3y'+su^2x'+t$ (see also
   \kbd{ellchangecurve}).
   \bprog
   ? E0 = ellinit([1,1]); P0 = [0,1]; v = [1,2,3,4];
   ? E = ellchangecurve(E0, v);
   ? P = ellchangepoint(P0,v)
   %3 = [-2, 3]
   ? ellisoncurve(E, P)
   %4 = 1
   ? ellchangepointinv(P,v)
   %5 = [0, 1]
   @eprog
  Variant: The reciprocal function \fun{GEN}{ellchangepointinv}{GEN x, GEN ch}
   inverts the coordinate change.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellchangepoint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellchangepointinv(*argv):
  '''
  ellchangepointinv
  Class: basic
  Section: elliptic_curves
  C-Name: ellchangepointinv
  Prototype: GG
  Help: ellchangepointinv(x,v): change data on point or vector of points x on an
   elliptic curve according to v=[u,r,s,t], inverse of ellchangepoint.
  Doc: 
   changes the coordinates of the point or vector of points $x$ using
   the inverse of the isomorphism associated to \kbd{v=[u,r,s,t]},
   i.e.~if $x'$ and $y'$ are the old coordinates, then $x=u^2x'+r$,
   $y=u^3y'+su^2x'+t$ (inverse of \kbd{ellchangepoint}).
   \bprog
   ? E0 = ellinit([1,1]); P0 = [0,1]; v = [1,2,3,4];
   ? E = ellchangecurve(E0, v);
   ? P = ellchangepoint(P0,v)
   %3 = [-2, 3]
   ? ellisoncurve(E, P)
   %4 = 1
   ? ellchangepointinv(P,v)
   %5 = [0, 1]  \\ we get back P0
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellchangepointinv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elldivpol(*argv):
  '''
  elldivpol
  Class: basic
  Section: elliptic_curves
  C-Name: elldivpol
  Prototype: GLDn
  Help: elldivpol(E,n,{v='x}): n-division polynomial f_n for the curve E in the
   variable v.
  Doc: $n$-division polynomial $f_n$ for the curve $E$ in the
   variable $v$. In standard notation, for any affine point $P = (X,Y)$ on the
   curve, we have
   $$[n]P = (\phi_n(P)\psi_n(P) : \omega_n(P) : \psi_n(P)^3)$$
   for some polynomials $\phi_n,\omega_n,\psi_n$ in
   $\Z[a_1,a_2,a_3,a_4,a_6][X,Y]$. We have $f_n(X) = \psi_n(X)$ for $n$ odd, and
   $f_n(X) = \psi_n(X,Y) (2Y + a_1X+a_3)$ for $n$ even. We have
   $$ f_1  = 1,\quad f_2 = 4X^3 + b_2X^2 + 2b_4 X + b_6, \quad f_3 = 3 X^4 + b_2 X^3 + 3b_4 X^2 + 3 b_6 X + b8, $$
   $$ f_4 = f_2(2X^6 + b_2 X^5 + 5b_4 X^4 + 10 b_6 X^3 + 10 b_8 X^2 +
   (b_2b_8-b_4b_6)X + (b_8b_4 - b_6^2)), \dots $$
   For $n \geq 2$, the roots of $f_n$ are the $X$-coordinates of points in $E[n]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elldivpol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elleisnum(*argv):
  '''
  elleisnum
  Class: basic
  Section: elliptic_curves
  C-Name: elleisnum
  Prototype: GLD0,L,p
  Help: elleisnum(w,k,{flag=0}): k being an even positive integer, computes the
   numerical value of the Eisenstein series of weight k at the lattice
   w, as given by ellperiods. When flag is non-zero and k=4 or 6, this gives the
   elliptic invariants g2 or g3 with the correct normalization.
  Doc: $k$ being an even positive integer, computes the numerical value of the
   Eisenstein series of weight $k$ at the lattice $w$, as given by
   \tet{ellperiods}, namely
   $$
   (2i \pi/\omega_2)^k
   \Big(1 + 2/\zeta(1-k) \sum_{n\geq 0} n^{k-1}q^n / (1-q^n)\Big),
   $$
   where $q = \exp(2i\pi \tau)$ and $\tau:=\omega_1/\omega_2$ belongs to the
   complex upper half-plane. It is also possible to directly input $w =
   [\omega_1,\omega_2]$, or an elliptic curve $E$ as given by \kbd{ellinit}.
   \bprog
   ? w = ellperiods([1,I]);
   ? elleisnum(w, 4)
   %2 = 2268.8726415508062275167367584190557607
   ? elleisnum(w, 6)
   %3 = -3.977978632282564763 E-33
   ? E = ellinit([1, 0]);
   ? elleisnum(E, 4, 1)
   %5 = -47.999999999999999999999999999999999998
   @eprog
   
   When \fl\ is non-zero and $k=4$ or 6, returns the elliptic invariants $g_2$
   or $g_3$, such that
   $$y^2 = 4x^3 - g_2 x - g_3$$
   is a Weierstrass equation for $E$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elleisnum(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elleta(*argv):
  '''
  elleta
  Class: basic
  Section: elliptic_curves
  C-Name: elleta
  Prototype: Gp
  Help: elleta(w): w=[w1,w2], returns the vector [eta1,eta2] of quasi-periods
   associated to [w1,w2].
  Doc: returns the quasi-periods $[\eta_1,\eta_2]$
   associated to the lattice basis $\var{w} = [\omega_1, \omega_2]$.
   Alternatively, \var{w} can be an elliptic curve $E$ as output by
   \kbd{ellinit}, in which case, the quasi periods associated to the period
   lattice basis \kbd{$E$.omega} (namely, \kbd{$E$.eta}) are returned.
   \bprog
   ? elleta([1, I])
   %1 = [3.141592653589793238462643383, 9.424777960769379715387930149*I]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elleta(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellfromj(*argv):
  '''
  ellfromj
  Class: basic
  Section: elliptic_curves
  C-Name: ellfromj
  Prototype: G
  Help: ellfromj(j): returns the coefficients [a1,a2,a3,a4,a6] of a fixed
   elliptic curve with j-invariant j.
  Doc: returns the coefficients $[a_1,a_2,a_3,a_4,a_6]$ of a fixed elliptic curve
   with $j$-invariant $j$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellfromj(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellgenerators(*argv):
  '''
  ellgenerators
  Class: basic
  Section: elliptic_curves
  C-Name: ellgenerators
  Prototype: G
  Help: ellgenerators(E): If E is an elliptic curve over the rationals,
   return the generators of the Mordell-Weil group associated to the curve.
   This relies on the curve being referenced in the elldata database.
   If E is an elliptic curve over a finite field Fq as output by ellinit(),
   return a minimal set of generators for the group E(Fq).
  Doc: 
   If $E$ is an elliptic curve over the rationals, return a $\Z$-basis of the
   free part of the \idx{Mordell-Weil group} associated to $E$.  This relies on
   the \tet{elldata} database being installed and referencing the curve, and so
   is only available for curves over $\Z$ of small conductors.
   If $E$ is an elliptic curve over a finite field $\F_q$ as output by
   \tet{ellinit}, return a minimal set of generators for the group $E(\F_q)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellgenerators(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellglobalred(*argv):
  '''
  ellglobalred
  Class: basic
  Section: elliptic_curves
  C-Name: ellglobalred
  Prototype: G
  Help: ellglobalred(E): E being an elliptic curve, returns [N,[u,r,s,t],c,
   faN,L], where N is the conductor of E, [u,r,s,t] leads to the standard model
   for E, c is the product of the local Tamagawa numbers c_p, faN is factor(N)
   and L[i] is elllocalred(E, faN[i,1]).
  Description: 
   (gen):gen        ellglobalred($1)
  Doc: 
   calculates the arithmetic conductor, the global
   minimal model of $E$ and the global \idx{Tamagawa number} $c$.
   $E$ must be an \var{ell} structure as output by \kbd{ellinit}, defined over
   $\Q$. The result is a vector $[N,v,c,F,L]$, where
   
   \item $N$ is the arithmetic conductor of the curve,
   
   \item $v$ gives the coordinate change for $E$ over $\Q$ to the minimal
   integral model (see \tet{ellminimalmodel}),
   
   \item $c$ is the product of the local Tamagawa numbers $c_p$, a quantity
   which enters in the \idx{Birch and Swinnerton-Dyer conjecture},\sidx{minimal model}
   
   \item $F$ is the factorization of $N$ over $\Z$.
   
   \item $L$ is a vector, whose $i$-th entry contains the local data
   at the $i$-th prime divisor of $N$, i.e. \kbd{L[i] = elllocalred(E,F[i,1])},
   where the local coordinate change has been deleted, and replaced by a $0$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellglobalred(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellgroup(*argv):
  '''
  ellgroup
  Class: basic
  Section: elliptic_curves
  C-Name: ellgroup0
  Prototype: GDGD0,L,
  Help: ellgroup(E,{p},{flag}): computes the structure of the group E(Fp)
   If flag is 1, return also generators.
  Doc: Let $E$ be an \var{ell} structure as output by \kbd{ellinit}, defined over
   $\Q$ or a finite field $\F_q$. The argument $p$ is best left omitted if the
   curve is defined over a finite field, and must be a prime number otherwise.
   This function computes the structure of the group $E(\F_q) \sim \Z/d_1\Z
   \times \Z/d_2\Z$, with $d_2\mid d_1$.
   
   If the curve is defined over $\Q$, $p$ must be explicitly given and the
   function computes the structure of the reduction over $\F_p$; the
   equation need not be minimal at $p$, but a minimal model will be more
   efficient. The reduction is allowed to be singular, and we return the
   structure of the (cyclic) group of non-singular points in this case.
   
   If the flag is $0$ (default), return $[d_1]$ or $[d_1, d_2]$, if $d_2>1$.
   If the flag is $1$, return a triple $[h,\var{cyc},\var{gen}]$, where
   $h$ is the curve cardinality, \var{cyc} gives the group structure as a
   product of cyclic groups (as per $\fl = 0$). More precisely, if $d_2 > 1$,
   the output is $[d_1d_2, [d_1,d_2],[P,Q]]$ where $P$ is
   of order $d_1$ and $[P,Q]$ generates the curve.
   \misctitle{Caution} It is not guaranteed that $Q$ has order $d_2$, which in
   the worst case requires an expensive discrete log computation. Only that
   \kbd{ellweilpairing(E, P, Q, d1)} has order $d_2$.
   \bprog
   ? E = ellinit([0,1]);  \\ y^2 = x^3 + 0.x + 1, defined over Q
   ? ellgroup(E, 7)
   %2 = [6, 2] \\ Z/6 x Z/2, non-cyclic
   ? E = ellinit([0,1] * Mod(1,11));  \\ defined over F_11
   ? ellgroup(E)   \\ no need to repeat 11
   %4 = [12]
   ? ellgroup(E, 11)   \\ ... but it also works
   %5 = [12]
   ? ellgroup(E, 13) \\ ouch, inconsistent input!
      ***   at top-level: ellgroup(E,13)
      ***                 ^--------------
      *** ellgroup: inconsistent moduli in Rg_to_Fp:
        11
        13
   ? ellgroup(E, 7, 1)
   %6 = [12, [6, 2], [[Mod(2, 7), Mod(4, 7)], [Mod(4, 7), Mod(4, 7)]]]
   @eprog\noindent
   If $E$ is defined over $\Q$, we allow singular reduction and in this case we
   return the structure of the group of non-singular points, satisfying
   $\#E_{ns}(\F_p) = p - a_p$.
   \bprog
   ? E = ellinit([0,5]);
   ? ellgroup(E, 5, 1)
   %2 = [5, [5], [[Mod(4, 5), Mod(2, 5)]]]
   ? ellap(E, 5)
   %3 = 0 \\ additive reduction at 5
   ? E = ellinit([0,-1,0,35,0]);
   ? ellgroup(E, 5, 1)
   %5 = [4, [4], [[Mod(2, 5), Mod(2, 5)]]]
   ? ellap(E, 5)
   %6 = 1 \\ split multiplicative reduction at 5
   ? ellgroup(E, 7, 1)
   %7 = [8, [8], [[Mod(3, 7), Mod(5, 7)]]]
   ? ellap(E, 7)
   %8 = -1 \\ non-split multiplicative reduction at 7
   @eprog
  Variant: Also available is \fun{GEN}{ellgroup}{GEN E, GEN p}, corresponding
   to \fl = 0.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellgroup0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellheight(*argv):
  '''
  ellheight
  Class: basic
  Section: elliptic_curves
  C-Name: ellheight0
  Prototype: GGD2,L,p
  Help: ellheight(E,x,{flag=2}): canonical height of point x on elliptic curve
   E. flag is optional and selects the algorithm
   used to compute the Archimedean local height. Its meaning is 0: use
   theta-functions, 1: use Tate's method, 2: use Mestre's AGM.
  Doc: global N\'eron-Tate height of the point $z$ on the elliptic curve
   $E$ (defined over $\Q$), using the normalization in Cremona's
   \emph{Algorithms for modular elliptic curves}. $E$
   must be an \kbd{ell} as output by \kbd{ellinit}; it needs not be given by a
   minimal model although the computation will be faster if it is. \fl\ selects
   the algorithm used to compute the Archimedean local height. If $\fl=0$,
   we use sigma and theta-functions and Silverman's trick (Computing
   heights on elliptic curves, \emph{Math.~Comp.} {\bf 51}; note that
   Silverman's height is twice ours). If
   $\fl=1$, use Tate's $4^n$ algorithm. If $\fl=2$, use Mestre's AGM algorithm.
   The latter converges quadratically and is much faster than the other two.
  Variant: Also available is \fun{GEN}{ghell}{GEN E, GEN x, long prec}
   ($\fl=2$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellheight0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellheegner(*argv):
  '''
  ellheegner
  Class: basic
  Section: elliptic_curves
  C-Name: ellheegner
  Prototype: G
  Help: ellheegner(E): return a rational non-torsion point on the elliptic curve E
   assumed to be of rank 1
  Doc: Let $E$ be an elliptic curve over the rationals, assumed to be of
   (analytic) rank $1$. This returns a non-torsion rational point on the curve,
   whose canonical height is equal to the product of the elliptic regulator by the
   analytic Sha.
   
   This uses the Heegner point method, described in Cohen GTM 239; the complexity
   is proportional to the product of the square root of the conductor and the
   height of the point (thus, it is preferable to apply it to strong Weil curves).
   \bprog
   ? E = ellinit([-157^2,0]);
   ? u = ellheegner(E); print(u[1], "\n", u[2])
   69648970982596494254458225/166136231668185267540804
   538962435089604615078004307258785218335/67716816556077455999228495435742408
   ? ellheegner(ellinit([0,1]))         \\ E has rank 0 !
    ***   at top-level: ellheegner(E=ellinit
    ***                 ^--------------------
    *** ellheegner: The curve has even analytic rank.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellheegner(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellinit(*argv):
  '''
  ellinit
  Class: basic
  Section: elliptic_curves
  C-Name: ellinit
  Prototype: GDGp
  Help: ellinit(x,{D=1}): let x be a vector [a1,a2,a3,a4,a6], or [a4,a6] if
   a1=a2=a3=0, defining the curve Y^2 + a1.XY + a3.Y = X^3 + a2.X^2 + a4.X +
   a6; x can also be a string, in which case the curve with matching name is
   retrieved from the elldata database, if available. This function initializes
   an elliptic curve over the domain D (inferred from coefficients if omitted).
  Description: 
   (gen, gen, small):ell:prec  ellinit($1, $2, prec)
  Doc: 
   initialize an \tet{ell} structure, associated to the elliptic curve $E$.
   $E$ is either
   
   \item a $5$-component vector $[a_1,a_2,a_3,a_4,a_6]$ defining the elliptic
   curve with Weierstrass equation
   $$ Y^2 + a_1 XY + a_3 Y = X^3 + a_2 X^2 + a_4 X + a_6, $$
   
   \item a $2$-component vector $[a_4,a_6]$ defining the elliptic
   curve with short Weierstrass equation
   $$ Y^2 = X^3 + a_4 X + a_6, $$
   
   \item a character string in Cremona's notation, e.g. \kbd{"11a1"}, in which
   case the curve is retrieved from the \tet{elldata} database if available.
   
   The optional argument $D$ describes the domain over which the curve is
   defined:
   
   \item the \typ{INT} $1$ (default): the field of rational numbers $\Q$.
   
   \item a \typ{INT} $p$, where $p$ is a prime number: the prime finite field
   $\F_p$.
   
   \item an \typ{INTMOD} \kbd{Mod(a, p)}, where $p$ is a prime number: the
   prime finite field $\F_p$.
   
   \item a \typ{FFELT}, as returned by \tet{ffgen}: the corresponding finite
   field $\F_q$.
   
   \item a \typ{PADIC}, $O(p^n)$: the field $\Q_p$, where $p$-adic quantities
   will be computed to a relative accuracy of $n$ digits. We advise to input a
   model defined over $\Q$ for such curves. In any case, if you input an
   approximate model with \typ{PADIC} coefficients, it will be replaced by a lift
   to $\Q$ (an exact model ``close'' to the one that was input) and all quantities
   will then be computed in terms of this lifted model, at the given accuracy.
   
   \item a \typ{REAL} $x$: the field $\C$ of complex numbers, where floating
   point quantities are by default computed to a relative accuracy of
   \kbd{precision}$(x)$. If no such argument is given, the value of
   \kbd{realprecision} at the time \kbd{ellinit} is called will be used.
   
   This argument $D$ is indicative: the curve coefficients are checked for
   compatibility, possibly changing $D$; for instance if $D = 1$ and
   an \typ{INTMOD} is found. If inconsistencies are detected, an error is
   raised:
   \bprog
   ? ellinit([1 + O(5), 1], O(7));
    ***   at top-level: ellinit([1+O(5),1],O
    ***                 ^--------------------
    *** ellinit: inconsistent moduli in ellinit: 7 != 5
   @eprog\noindent If the curve coefficients are too general to fit any of the
   above domain categories, only basic operations, such as point addition, will
   be supported later.
   
   If the curve (seen over the domain $D$) is singular, fail and return an
   empty vector $[]$.
   \bprog
   ? E = ellinit([0,0,0,0,1]); \\ y^2 = x^3 + 1, over Q
   ? E = ellinit([0,1]);       \\ the same curve, short form
   ? E = ellinit("36a1");      \\ sill the same curve, Cremona's notations
   ? E = ellinit([0,1], 2)     \\ over F2: singular curve
   %4 = []
   ? E = ellinit(['a4,'a6] * Mod(1,5));  \\ over F_5[a4,a6], basic support !
   @eprog\noindent
   
   The result of \tet{ellinit} is an \tev{ell} structure. It contains at least
   the following information in its components:
   %
   $$ a_1,a_2,a_3,a_4,a_6,b_2,b_4,b_6,b_8,c_4,c_6,\Delta,j.$$
   %
   All are accessible via member functions. In particular, the discriminant is
   \kbd{$E$.disc}, and the $j$-invariant is \kbd{$E$.j}.
   \bprog
   ? E = ellinit([a4, a6]);
   ? E.disc
   %2 = -64*a4^3 - 432*a6^2
   ? E.j
   %3 = -6912*a4^3/(-4*a4^3 - 27*a6^2)
   @eprog
   Further components contain domain-specific data, which are in general dynamic:
   only computed when needed, and then cached in the structure.
   \bprog
   ? E = ellinit([2,3], 10^60+7);  \\ E over F_p, p large
   ? ellap(E)
   time = 4,440 ms.
   %2 = -1376268269510579884904540406082
   ? ellcard(E);  \\ now instantaneous !
   time = 0 ms.
   ? ellgenerators(E);
   time = 5,965 ms.
   ? ellgenerators(E); \\ second time instantaneous
   time = 0 ms.
   @eprog
   See the description of member functions related to elliptic curves at the
   beginning of this section.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellisoncurve(*argv):
  '''
  ellisoncurve
  Class: basic
  Section: elliptic_curves
  C-Name: ellisoncurve
  Prototype: GG
  Help: ellisoncurve(E,z): true(1) if z is on elliptic curve E, false(0) if not.
  Doc: gives 1 (i.e.~true) if the point $z$ is on the elliptic curve $E$, 0
   otherwise. If $E$ or $z$ have imprecise coefficients, an attempt is made to
   take this into account, i.e.~an imprecise equality is checked, not a precise
   one. It is allowed for $z$ to be a vector of points in which case a vector
   (of the same type) is returned.
  Variant: Also available is \fun{int}{oncurve}{GEN E, GEN z} which does not
   accept vectors of points.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellisoncurve(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elllseries(*argv):
  '''
  elllseries
  Class: basic
  Section: elliptic_curves
  C-Name: elllseries
  Prototype: GGDGp
  Help: elllseries(E,s,{A=1}): L-series at s of the elliptic curve E, where A
   a cut-off point close to 1.
  Doc: 
   $E$ being an elliptic curve, given by an arbitrary model over $\Q$ as output
   by \kbd{ellinit}, this function computes the value of the $L$-series of $E$ at
   the (complex) point $s$. This function uses an $O(N^{1/2})$ algorithm, where
   $N$ is the conductor.
   
   The optional parameter $A$ fixes a cutoff point for the integral and is best
   left omitted; the result must be independent of $A$, up to
   \kbd{realprecision}, so this allows to check the function's accuracy.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elllseries(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elllocalred(*argv):
  '''
  elllocalred
  Class: basic
  Section: elliptic_curves
  C-Name: elllocalred
  Prototype: GG
  Help: elllocalred(E,p): E being an elliptic curve, returns
   [f,kod,[u,r,s,t],c], where f is the conductor's exponent, kod is the Kodaira
   type for E at p, [u,r,s,t] is the change of variable needed to make E
   minimal at p, and c is the local Tamagawa number c_p.
  Doc: 
   calculates the \idx{Kodaira} type of the local fiber of the elliptic curve
   $E$ at the prime $p$. $E$ must be an \var{ell} structure as output by
   \kbd{ellinit}, and is assumed to have all its coefficients $a_i$ in $\Z$.
   The result is a 4-component vector $[f,kod,v,c]$. Here $f$ is the exponent of
   $p$ in the arithmetic conductor of $E$, and $kod$ is the Kodaira type which
   is coded as follows:
   
   1 means good reduction (type I$_0$), 2, 3 and 4 mean types II, III and IV
   respectively, $4+\nu$ with $\nu>0$ means type I$_\nu$;
   finally the opposite values $-1$, $-2$, etc.~refer to the starred types
   I$_0^*$, II$^*$, etc. The third component $v$ is itself a vector $[u,r,s,t]$
   giving the coordinate changes done during the local reduction;
   $u = 1$ if and only if the given equation was already minimal at $p$.
   Finally, the last component $c$ is the local \idx{Tamagawa number} $c_p$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elllocalred(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elllog(*argv):
  '''
  elllog
  Class: basic
  Section: elliptic_curves
  C-Name: elllog
  Prototype: GGGDG
  Help: elllog(E,P,G,{o}): return the discrete logarithm of the point P of
   the elliptic curve E in base G. If present, o represents the order of G.
   If not present, assume that G generates the curve.
  Doc: given two points $P$ and $G$ on the elliptic curve $E/\F_q$, returns the
   discrete logarithm of $P$ in base $G$, i.e. the smallest non-negative
   integer $n$ such that $P = [n]G$.
   See \tet{znlog} for the limitations of the underlying discrete log algorithms.
   If present, $o$ represents the order of $G$, see \secref{se:DLfun};
   the preferred format for this parameter is \kbd{[N, factor(N)]}, where $N$
   is  the order of $G$.
   
   If no $o$ is given, assume that $G$ generates the curve.
   The function also assumes that $P$ is a multiple of $G$.
   \bprog
   ? a = ffgen(ffinit(2,8),'a);
   ? E = ellinit([a,1,0,0,1]);  \\ over F_{2^8}
   ? x = a^3; y = ellordinate(E,x)[1];
   ? P = [x,y]; G = ellmul(E, P, 113);
   ? ord = [242, factor(242)]; \\ P generates a group of order 242. Initialize.
   ? ellorder(E, G, ord)
   %4 = 242
   ? e = elllog(E, P, G, ord)
   %5 = 15
   ? ellmul(E,G,e) == P
   %6 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elllog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellminimalmodel(*argv):
  '''
  ellminimalmodel
  Class: basic
  Section: elliptic_curves
  C-Name: ellminimalmodel
  Prototype: GD&
  Help: ellminimalmodel(E,{&v}): return the standard minimal integral model of
   the rational elliptic curve E. Sets v to the corresponding change of
   variables.
  Doc: return the standard minimal integral model of the rational elliptic
   curve $E$. If present, sets $v$ to the corresponding change of variables,
   which is a vector $[u,r,s,t]$ with rational components. The return value is
   identical to that of \kbd{ellchangecurve(E, v)}.
   
   The resulting model has integral coefficients, is everywhere minimal, $a_1$
   is 0 or 1, $a_2$ is 0, 1 or $-1$ and $a_3$ is 0 or 1. Such a model is
   unique, and the vector $v$ is unique if we specify that $u$ is positive,
   which we do. \sidx{minimal model}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellminimalmodel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellpow(*argv):
  '''
  ellpow
  Class: basic
  Section: elliptic_curves
  C-Name: ellmul
  Prototype: GGG
  Help: ellpow(E,z,n): deprecated alias for ellmul.
  Doc: deprecated alias for \kbd{ellmul}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellmul(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellneg(*argv):
  '''
  ellneg
  Class: basic
  Section: elliptic_curves
  C-Name: ellneg
  Prototype: GG
  Help: ellneg(E,z): opposite of the point z on elliptic curve E.
  Doc: 
   Opposite of the point $z$ on elliptic curve $E$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellneg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellorder(*argv):
  '''
  ellorder
  Class: basic
  Section: elliptic_curves
  C-Name: ellorder
  Prototype: GGDG
  Help: ellorder(E,z,{o}): order of the point z on the elliptic curve E over Q
   or a finite field, 0 if non-torsion. The parameter o, if present,
   represents a non-zero multiple of the order of z.
  Doc: gives the order of the point $z$ on the elliptic
   curve $E$, defined over $\Q$ or a finite field.
   If the curve is defined over $\Q$, return (the impossible value) zero if the
   point has infinite order.
   \bprog
   ? E = ellinit([-157^2,0]);  \\ the "157-is-congruent" curve
   ? P = [2,2]; ellorder(E, P)
   %2 = 2
   ? P = ellheegner(E); ellorder(E, P) \\ infinite order
   %3 = 0
   ? E = ellinit(ellfromj(ffgen(5^10)));
   ? ellcard(E)
   %5 = 9762580
   ? P = random(E); ellorder(E, P)
   %6 = 4881290
   ? p = 2^160+7; E = ellinit([1,2], p);
   ? N = ellcard(E)
   %8 = 1461501637330902918203686560289225285992592471152
   ? o = [N, factor(N)];
   ? for(i=1,100, ellorder(E,random(E)))
   time = 260 ms.
   @eprog
   The parameter $o$, is now mostly useless, and kept for backward
   compatibility. If present, it represents a non-zero multiple of the order
   of $z$, see \secref{se:DLfun}; the preferred format for this parameter is
   \kbd{[ord, factor(ord)]}, where \kbd{ord} is the cardinality of the curve.
   It is no longer needed since PARI is now able to compute it over large
   finite fields (was restricted to small prime fields at the time this feature
   was introduced), \emph{and} caches the result in $E$ so that it is computed
   and factored only once. Modifying the last example, we see that including
   this extra parameter provides no improvement:
   \bprog
   ? o = [N, factor(N)];
   ? for(i=1,100, ellorder(E,random(E),o))
   time = 260 ms.
   @eprog
  Variant: The obsolete form \fun{GEN}{orderell}{GEN e, GEN z} should no longer be
   used.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellorder(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellordinate(*argv):
  '''
  ellordinate
  Class: basic
  Section: elliptic_curves
  C-Name: ellordinate
  Prototype: GGp
  Help: ellordinate(E,x): y-coordinates corresponding to x-ordinate x on
   elliptic curve E.
  Doc: 
   gives a 0, 1 or 2-component vector containing
   the $y$-coordinates of the points of the curve $E$ having $x$ as
   $x$-coordinate.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellordinate(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellperiods(*argv):
  '''
  ellperiods
  Class: basic
  Section: elliptic_curves
  C-Name: ellperiods
  Prototype: GD0,L,p
  Help: ellperiods(w, {flag = 0}): w describes a complex period lattice ([w1,w2]
   or an ellinit structure). Returns normalized periods [W1,W2] generating the
   same lattice such that tau := W1/W2 satisfies Im(tau) > 0 and lies in the
   standard fundamental domain for SL2. If flag is 1, the return value is
   [[W1,W2], [eta1,eta2]], where eta1, eta2 are the quasi-periods associated to
   [W1,W2], satisfying eta1 W2 - eta2 W1 = 2 I Pi.
  Doc: Let $w$ describe a complex period lattice ($w = [w_1,w_2]$
   or an ellinit structure). Returns normalized periods $[W_1,W_2]$ generating
   the same lattice such that $\tau := W_1/W_2$ has positive imaginary part
   and lies in the standard fundamental domain for $\text{SL}_2(\Z)$.
   
   If $\fl = 1$, the function returns $[[W_1,W_2], [\eta_1,\eta_2]]$, where
   $\eta_1$ and $\eta_2$ are the quasi-periods associated to
   $[W_1,W_2]$, satisfying $\eta_1 W_2 - \eta_2 W_1 = 2 i \pi$.
   
   The output of this function is meant to be used as the first argument
   given to ellwp, ellzeta, ellsigma or elleisnum. Quasi-periods are
   needed by ellzeta and ellsigma only.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellperiods(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellrootno(*argv):
  '''
  ellrootno
  Class: basic
  Section: elliptic_curves
  C-Name: ellrootno
  Prototype: lGDG
  Help: ellrootno(E,{p}): root number for the L-function of the elliptic
   curve E/Q at a prime p (including 0, for the infinite place); global root
   number if p is omitted.
  Doc: $E$ being an \var{ell} structure over $\Q$ as output by \kbd{ellinit},
   this function computes the local root number of its $L$-series at the place
   $p$ (at the infinite place if $p = 0$). If $p$ is omitted, return the global
   root number. Note that the global root number is the sign of the functional
   equation and conjecturally is the parity of the rank of the \idx{Mordell-Weil
   group}. The equation for $E$ needs not be minimal at $p$, but if the model
   is already minimal the function will run faster.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.ellrootno(*c_arg_tuple)

def ellsigma(*argv):
  '''
  ellsigma
  Class: basic
  Section: elliptic_curves
  C-Name: ellsigma
  Prototype: GDGD0,L,p
  Help: ellsigma(L,{z='x},{flag=0}): computes the value at z of the Weierstrass
   sigma function attached to the lattice w, as given by ellperiods(,1).
   If flag = 1, returns an arbitrary determination of the logarithm of sigma.
  Doc: Computes the value at $z$ of the Weierstrass $\sigma$ function attached to
   the lattice $L$ as given by \tet{ellperiods}$(,1)$: including quasi-periods
   is useful, otherwise there are recomputed from scratch for each new $z$.
   $$ \sigma(z, L) = z \prod_{\omega\in L^*} \left(1 -
   \dfrac{z}{\omega}\right)e^{\dfrac{z}{\omega} + \dfrac{z^2}{2\omega^2}}.$$
   It is also possible to directly input $L = [\omega_1,\omega_2]$,
   or an elliptic curve $E$ as given by \kbd{ellinit} ($L = \kbd{E.omega}$).
   \bprog
   ? w = ellperiods([1,I], 1);
   ? ellsigma(w, 1/2)
   %2 = 0.47494937998792065033250463632798296855
   ? E = ellinit([1,0]);
   ? ellsigma(E) \\ at 'x, implicitly at default seriesprecision
   %4 = x + 1/60*x^5 - 1/10080*x^9 - 23/259459200*x^13 + O(x^17)
   @eprog
   
   If $\fl=1$, computes an arbitrary determination of $\log(\sigma(z))$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellsigma(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellsub(*argv):
  '''
  ellsub
  Class: basic
  Section: elliptic_curves
  C-Name: ellsub
  Prototype: GGG
  Help: ellsub(E,z1,z2): difference of the points z1 and z2 on elliptic curve E.
  Doc: 
   difference of the points $z1$ and $z2$ on the
   elliptic curve corresponding to $E$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellsub(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elltaniyama(*argv):
  '''
  elltaniyama
  Class: basic
  Section: elliptic_curves
  C-Name: elltaniyama
  Prototype: GDP
  Help: elltaniyama(E, {d = seriesprecision}): modular parametrization of
   elliptic curve E/Q.
  Doc: 
   computes the modular parametrization of the elliptic curve $E/\Q$,
   where $E$ is an \var{ell} structure as output by \kbd{ellinit}. This returns
   a two-component vector $[u,v]$ of power series, given to $d$ significant
   terms (\tet{seriesprecision} by default), characterized by the following two
   properties. First the point $(u,v)$ satisfies the equation of the elliptic
   curve. Second, let $N$ be the conductor of $E$ and $\Phi: X_0(N)\to E$
   be a modular parametrization; the pullback by $\Phi$ of the
   N\'eron differential $du/(2v+a_1u+a_3)$ is equal to $2i\pi
   f(z)dz$, a holomorphic differential form. The variable used in the power
   series for $u$ and $v$ is $x$, which is implicitly understood to be equal to
   $\exp(2i\pi z)$.
   
   The algorithm assumes that $E$ is a \emph{strong} \idx{Weil curve}
   and that the Manin constant is equal to 1: in fact, $f(x) = \sum_{n > 0}
   \kbd{ellan}(E, n) x^n$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elltaniyama(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elltatepairing(*argv):
  '''
  elltatepairing
  Class: basic
  Section: elliptic_curves
  C-Name: elltatepairing
  Prototype: GGGG
  Help: elltatepairing(E, P, Q, m): Computes the Tate pairing of the two points
   P and Q on the elliptic curve E. The point P must be of m-torsion.
  Doc: Computes the Tate pairing of the two points $P$ and $Q$ on the elliptic
   curve $E$. The point $P$ must be of $m$-torsion.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elltatepairing(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def elltors(*argv):
  '''
  elltors
  Class: basic
  Section: elliptic_curves
  C-Name: elltors0
  Prototype: GD0,L,
  Help: elltors(E,{flag=0}): torsion subgroup of elliptic curve E: order,
   structure, generators. If flag = 0, use division polynomials; if flag = 1, use
   Lutz-Nagell; if flag = 2, use Doud's algorithm.
  Doc: 
   if $E$ is an elliptic curve \emph{defined over $\Q$}, outputs the torsion
   subgroup of $E$ as a 3-component vector \kbd{[t,v1,v2]}, where \kbd{t} is the
   order of the torsion group, \kbd{v1} gives the structure of the torsion group
   as a product of cyclic groups (sorted by decreasing order), and \kbd{v2}
   gives generators for these cyclic groups. $E$ must be an \var{ell} structure
   as output by \kbd{ellinit}, defined over $\Q$.
   
   \bprog
   ?  E = ellinit([-1,0]);
   ?  elltors(E)
   %1 = [4, [2, 2], [[0, 0], [1, 0]]]
   @eprog
   Here, the torsion subgroup is isomorphic to $\Z/2\Z \times \Z/2\Z$, with
   generators $[0,0]$ and $[1,0]$.
   
   If $\fl = 0$, find rational roots of division polynomials.
   
   If $\fl = 1$, use Lutz-Nagell (\emph{much} slower).
   
   If $\fl = 2$, use Doud's algorithm: bound torsion by computing $\#E(\F_p)$
   for small primes of good reduction, then look for torsion points using
   Weierstrass $\wp$ function (and Mazur's classification). For this variant,
   $E$ must be an \var{ell}.
  Variant: Also available is \fun{GEN}{elltors}{GEN E} for \kbd{elltors(E, 0)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.elltors0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellweilpairing(*argv):
  '''
  ellweilpairing
  Class: basic
  Section: elliptic_curves
  C-Name: ellweilpairing
  Prototype: GGGG
  Help: ellweilpairing(E, P, Q, m): Computes the Weil pairing of the two points
   of m-torsion P and Q on the elliptic curve E.
  Doc: Computes the Weil pairing of the two points of $m$-torsion $P$ and $Q$
   on the elliptic curve $E$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellweilpairing(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellwp(*argv):
  '''
  ellwp
  Class: basic
  Section: elliptic_curves
  C-Name: ellwp0
  Prototype: GDGD0,L,p
  Help: ellwp(w,{z='x},{flag=0}): computes the value at z of the Weierstrass P
   function attached to the lattice w, as given by ellperiods. Optional flag
   means 0 (default), compute only P(z), 1 compute [P(z),P'(z)].
  Doc: Computes the value at $z$ of the Weierstrass $\wp$ function attached to
   the lattice $w$ as given by \tet{ellperiods}. It is also possible to
   directly input $w = [\omega_1,\omega_2]$, or an elliptic curve $E$ as given
   by \kbd{ellinit} ($w = \kbd{E.omega}$).
   \bprog
   ? w = ellperiods([1,I]);
   ? ellwp(w, 1/2)
   %2 = 6.8751858180203728274900957798105571978
   ? E = ellinit([1,1]);
   ? ellwp(E, 1/2)
   %4 = 3.9413112427016474646048282462709151389
   @eprog\noindent One can also compute the series expansion around $z = 0$:
   \bprog
   ? E = ellinit([1,0]);
   ? ellwp(E)              \\ 'x implicitly at default seriesprecision
   %5 = x^-2 - 1/5*x^2 + 1/75*x^6 - 2/4875*x^10 + O(x^14)
   ? ellwp(E, x + O(x^12)) \\ explicit precision
   %6 = x^-2 - 1/5*x^2 + 1/75*x^6 + O(x^9)
   @eprog
   
   Optional \fl\ means 0 (default): compute only $\wp(z)$, 1: compute
   $[\wp(z),\wp'(z)]$.
  Variant: For $\fl = 0$, we also have
   \fun{GEN}{ellwp}{GEN w, GEN z, long prec}, and
   \fun{GEN}{ellwpseries}{GEN E, long v, long precdl} for the power series in
   variable $v$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellwp0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellzeta(*argv):
  '''
  ellzeta
  Class: basic
  Section: elliptic_curves
  C-Name: ellzeta
  Prototype: GDGp
  Help: ellzeta(w,{z='x}): computes the value at z of the Weierstrass Zeta
   function attached to the lattice w, as given by ellperiods(,1).
  Doc: Computes the value at $z$ of the Weierstrass $\zeta$ function attached to
   the lattice $w$ as given by \tet{ellperiods}$(,1)$: including quasi-periods
   is useful, otherwise there are recomputed from scratch for each new $z$.
   $$ \zeta(z, L) = \dfrac{1}{z} + z^2\sum_{\omega\in L^*}
   \dfrac{1}{\omega^2(z-\omega)}.$$
   It is also possible to directly input $w = [\omega_1,\omega_2]$,
   or an elliptic curve $E$ as given by \kbd{ellinit} ($w = \kbd{E.omega}$).
   The quasi-periods of $\zeta$, such that
   $$\zeta(z + a\omega_1 + b\omega_2) = \zeta(z) + a\eta_1 + b\eta_2 $$
   for integers $a$ and $b$ are obtained as $\eta_i = 2\zeta(\omega_i/2)$.
   Or using directly \tet{elleta}.
   \bprog
   ? w = ellperiods([1,I],1);
   ? ellzeta(w, 1/2)
   %2 = 1.5707963267948966192313216916397514421
   ? E = ellinit([1,0]);
   ? ellzeta(E, E.omega[1]/2)
   %4 = 0.84721308479397908660649912348219163647
   @eprog\noindent One can also compute the series expansion around $z = 0$
   (the quasi-periods are useless in this case):
   \bprog
   ? E = ellinit([0,1]);
   ? ellzeta(E) \\ at 'x, implicitly at default seriesprecision
   %4 = x^-1 + 1/35*x^5 - 1/7007*x^11 + O(x^15)
   ? ellzeta(E, x + O(x^20)) \\ explicit precision
   %5 = x^-1 + 1/35*x^5 - 1/7007*x^11 + 1/1440257*x^17 + O(x^18)
   @eprog\noindent
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellzeta(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellheightmatrix(*argv):
  '''
  ellheightmatrix
  Class: basic
  Section: elliptic_curves
  C-Name: mathell
  Prototype: GGp
  Help: ellheightmatrix(E,x): gives the height matrix for vector of points x
   on elliptic curve E, assume to be a minimal model.
  Doc: 
   $x$ being a vector of points, this
   function outputs the Gram matrix of $x$ with respect to the N\'eron-Tate
   height, in other words, the $(i,j)$ component of the matrix is equal to
   \kbd{ellbil($E$,x[$i$],x[$j$])}. The rank of this matrix, at least in some
   approximate sense, gives the rank of the set of points, and if $x$ is a
   basis of the \idx{Mordell-Weil group} of $E$, its determinant is equal to
   the regulator of $E$. Note that this matrix should be divided by 2 to be in
   accordance with certain normalizations. $E$ is assumed to be integral,
   given by a minimal model.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mathell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellztopoint(*argv):
  '''
  ellztopoint
  Class: basic
  Section: elliptic_curves
  C-Name: pointell
  Prototype: GGp
  Help: ellztopoint(E,z): coordinates of point P on the curve E corresponding
   to the complex number z.
  Doc: 
   $E$ being an \var{ell} as output by
   \kbd{ellinit}, computes the coordinates $[x,y]$ on the curve $E$
   corresponding to the complex number $z$. Hence this is the inverse function
   of \kbd{ellpointtoz}. In other words, if the curve is put in Weierstrass
   form $y^2 = 4x^3 - g_2x - g_3$, $[x,y]$ represents the Weierstrass
   $\wp$-function\sidx{Weierstrass $\wp$-function} and its derivative. More
   precisely, we have
   $$x = \wp(z) - b_2/12,\quad y = \wp'(z) - (a_1 x + a_3)/2.$$
   If $z$ is in the lattice defining $E$ over $\C$, the result is the point at
   infinity $[0]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.pointell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellpointtoz(*argv):
  '''
  ellpointtoz
  Class: basic
  Section: elliptic_curves
  C-Name: zell
  Prototype: GGp
  Help: ellpointtoz(E,P): lattice point z corresponding to the point P on the
   elliptic curve E.
  Doc: 
   if $E/\C \simeq \C/\Lambda$ is a complex elliptic curve ($\Lambda =
   \kbd{E.omega}$),
   computes a complex number $z$, well-defined modulo the lattice $\Lambda$,
   corresponding to the point $P$; i.e.~such that
    $P = [\wp_\Lambda(z),\wp'_\Lambda(z)]$
   satisfies the equation
   $$y^2 = 4x^3 - g_2 x - g_3,$$
   where $g_2$, $g_3$ are the elliptic invariants.
   
   If $E$ is defined over $\R$ and $P\in E(\R)$, we have more precisely, $0 \leq
   \Re(t) < w1$ and $0 \leq \Im(t) < \Im(w2)$, where $(w1,w2)$ are the real and
   complex periods of $E$.
   \bprog
   ? E = ellinit([0,1]); P = [2,3];
   ? z = ellpointtoz(E, P)
   %2 = 3.5054552633136356529375476976257353387
   ? ellwp(E, z)
   %3 = 2.0000000000000000000000000000000000000
   ? ellztopoint(E, z) - P
   %4 = [6.372367644529809109 E-58, 7.646841173435770930 E-57]
   ? ellpointtoz(E, [0]) \\ the point at infinity
   %5 = 0
   @eprog
   
   If $E/\Q_p$ has multiplicative reduction, then $E/\bar{\Q_p}$ is analytically
   isomorphic to $\bar{\Q}_p^*/q^\Z$ (Tate curve) for some $p$-adic integer $q$.
   The behaviour is then as follows:
   
   \item If the reduction is split ($E.\kbd{tate[2]}$ is a \typ{PADIC}), we have
   an isomorphism $\phi: E(\Q_p) \simeq \Q_p^*/q^\Z$ and the function returns
   $\phi(P)\in \Q_p$.
   
   \item If the reduction is \emph{not} split ($E.\kbd{tate[2]}$ is a
   \typ{POLMOD}), we only have an isomorphism $\phi: E(K) \simeq K^*/q^\Z$ over
   the unramified quadratic extension $K/\Q_p$. In this case, the output
   $\phi(P)\in K$ is a \typ{POLMOD}.
   \bprog
   ? E = ellinit([0,-1,1,0,0], O(11^5)); P = [0,0];
   ? [u2,u,q] = E.tate; type(u) \\ split multiplicative reduction
   %2 = "t_PADIC"
   ? ellmul(E, P, 5)  \\ P has order 5
   %3 = [0]
   ? z = ellpointtoz(E, [0,0])
   %4 = 3 + 11^2 + 2*11^3 + 3*11^4 + O(11^5)
   ? z^5
   %5 = 1 + O(11^5)
   ? E = ellinit(ellfromj(1/4), O(2^6)); x=1/2; y=ellordinate(E,x)[1];
   ? z = ellpointtoz(E,[x,y]); \\ t_POLMOD of t_POL with t_PADIC coeffs
   ? liftint(z) \\ lift all p-adics
   %8 = Mod(8*u + 7, u^2 + 437)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.zell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellmodulareqn(*argv):
  '''
  ellmodulareqn
  Class: basic
  Section: elliptic_curves
  C-Name: ellmodulareqn
  Prototype: LDnDn
  Help: ellmodulareqn(N,{x},{y}): return a vector [eqn, t] where eqn is a modular
   equation of level N, for N<500, N prime. This requires the package seadata to
   be installed.  The equation is either of canonical type (t=0) or of Atkin type
   (t=1)
  Doc: return a vector [\kbd{eqn},$t$] where \kbd{eqn} is a modular equation of
   level $N$, i.e.~a bivariate polynomial with integer coefficients; $t$
   indicates the type of this equation: either \emph{canonical} ($t = 0$) or
   \emph{Atkin} ($t = 1$). This function currently requires the package
   \kbd{seadata} to be installed and is limited to $N<500$, $N$ prime.
   
   Let $j$ be the $j$-invariant function. The polynomial \kbd{eqn} satisfies
   the following functional equation, which allows to compute the values of the
   classical modular polynomial $\Phi_N$ of prime level $N$, such that
   $\Phi_N(j(\tau), j(N\tau)) = 0$, while being much smaller than the latter:
   
   \item for canonical type:
    $P(f(\tau),j(\tau)) = P(N^s/f(\tau),j(N\*\tau)) = 0$,
    where $s = 12/\gcd(12,N-1)$;
   
   \item for Atkin type:
    $P(f(\tau),j(\tau)) = P(f(\tau),j(N\*\tau)) = 0$.
   
   \noindent In both cases, $f$ is a suitable modular function (see below).
   
   The following GP function returns values of the classical modular polynomial
   by eliminating $f(\tau)$ in the above two equations, for $N\leq 31$ or
   $N\in\{41,47,59,71\}$.
   
   \bprog
   classicaleqn(N, X='X, Y='Y)=
   {
     my(E=ellmodulareqn(N), P=E[1], t=E[2], Q, d);
     if(poldegree(P,'y)>2,error("level unavailable in classicaleqn"));
     if (t == 0,
       my(s = 12/gcd(12,N-1));
       Q = 'x^(N+1) * substvec(P,['x,'y],[N^s/'x,Y]);
       d = N^(s*(2*N+1)) * (-1)^(N+1);
     ,
       Q = subst(P,'y,Y);
       d = (X-Y)^(N+1));
     polresultant(subst(P,'y,X), Q) / d;
   }
   @eprog
   
   More precisely, let $W_N(\tau)={{-1}\over{N\*\tau}}$ be the Atkin-Lehner
   involution; we have $j(W_N(\tau)) = j(N\*\tau)$ and the function $f$ also
   satisfies:
   
   \item for canonical type:
      $f(W_N(\tau)) = N^s/f(\tau)$;
   
   \item for Atkin type:
      $f(W_N(\tau)) = f(\tau)$.
   
   \noindent Furthermore, for an equation of canonical type, $f$ is the standard
   $\eta$-quotient
   $$f(\tau) = N^s \* \big(\eta(N\*\tau) / \eta(\tau) \big)^{2\*s},$$
   where $\eta$ is Dedekind's eta function, which is invariant under
   $\Gamma_0(N)$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ellmodulareqn(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Str(*argv):
  '''
  Str
  Class: basic
  Section: conversions
  C-Name: Str
  Prototype: s*
  Help: Str({x}*): concatenates its (string) argument into a single string.
  Description: 
   (gen):genstr:copy:parens      $genstr:1
   (gen,gen):genstr              Str(mkvec2($1, $2))
   (gen,gen,gen):genstr          Str(mkvec3($1, $2, $3))
   (gen,gen,gen,gen):genstr      Str(mkvec4($1, $2, $3, $4))
   (gen,...):genstr              Str(mkvecn($#, $2))
  Doc: 
   converts its argument list into a
   single character string (type \typ{STR}, the empty string if $x$ is omitted).
   To recover an ordinary \kbd{GEN} from a string, apply \kbd{eval} to it. The
   arguments of \kbd{Str} are evaluated in string context, see \secref{se:strings}.
   
   \bprog
   ? x2 = 0; i = 2; Str(x, i)
   %1 = "x2"
   ? eval(%)
   %2 = 0
   @eprog\noindent
   This function is mostly useless in library mode. Use the pair
   \tet{strtoGEN}/\tet{GENtostr} to convert between \kbd{GEN} and \kbd{char*}.
   The latter returns a malloced string, which should be freed after usage.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Str(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Strchr(*argv):
  '''
  Strchr
  Class: basic
  Section: conversions
  C-Name: Strchr
  Prototype: G
  Help: Strchr(x): converts x to a string, translating each integer into a
   character.
  Doc: 
   converts $x$ to a string, translating each integer
   into a character.
   \bprog
   ? Strchr(97)
   %1 = "a"
   ? Vecsmall("hello world")
   %2 = Vecsmall([104, 101, 108, 108, 111, 32, 119, 111, 114, 108, 100])
   ? Strchr(%)
   %3 = "hello world"
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Strchr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Strexpand(*argv):
  '''
  Strexpand
  Class: basic
  Section: conversions
  C-Name: Strexpand
  Prototype: s*
  Help: Strexpand({x}*): concatenates its (string) argument into a single
   string, performing tilde expansion.
  Doc: 
   converts its argument list into a
   single character string (type \typ{STR}, the empty string if $x$ is omitted).
   Then perform \idx{environment expansion}, see \secref{se:envir}.
   This feature can be used to read \idx{environment variable} values.
   \bprog
   ? Strexpand("$HOME/doc")
   %1 = "/home/pari/doc"
   @eprog
   
   The individual arguments are read in string context, see \secref{se:strings}.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Strexpand(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Strtex(*argv):
  '''
  Strtex
  Class: basic
  Section: conversions
  C-Name: Strtex
  Prototype: s*
  Help: Strtex({x}*): translates its (string) arguments to TeX format and
   returns the resulting string.
  Doc: 
   translates its arguments to TeX
   format, and concatenates the results into a single character string (type
   \typ{STR}, the empty string if $x$ is omitted).
   
   The individual arguments are read in string context, see \secref{se:strings}.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Strtex(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def dbg_x(*argv):
  '''
  dbg_x
  Class: basic
  Section: programming/control
  C-Name: dbgGEN
  Prototype: vGD-1,L,
  Help: dbg_x(A{,n}): print inner structure of A, complete if n is omitted, up to
   level n otherwise. Intended for debugging.
  Doc: Print the inner structure of \kbd{A}, complete if \kbd{n} is omitted, up
   to level \kbd{n} otherwise. This is useful for debugging. This is similar to
   \b{x} but does not require \kbd{A} to be an history entry. In particular,
   it can be used in the break loop.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.dbgGEN(*c_arg_tuple)

def error(*argv):
  '''
  error
  Class: basic
  Section: programming/specific
  C-Name: error0
  Prototype: vs*
  Help: error({str}*): abort script with error message str.
  Description: 
   (error):void  pari_err(0, $1)
   (?gen,...):void  pari_err(e_MISC, "${2 format_string}"${2 format_args})
  Doc: outputs its argument list (each of
   them interpreted as a string), then interrupts the running \kbd{gp} program,
   returning to the input prompt. For instance
   \bprog
   error("n = ", n, " is not squarefree!")
   @eprog\noindent
    % \syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.error0(*c_arg_tuple)

def getenv(*argv):
  '''
  getenv
  Class: basic
  Section: programming/specific
  C-Name: gp_getenv
  Prototype: s
  Help: getenv(s): value of the environment variable s, 0 if it is not defined.
  Doc: return the value of the environment variable \kbd{s} if it is defined, otherwise return 0.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gp_getenv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def readvec(*argv):
  '''
  readvec
  Class: basic
  Section: programming/specific
  C-Name: gp_readvec_file
  Prototype: D"",s,
  Help: readvec({filename}): create a vector whose components are the evaluation
   of all the expressions found in the input file filename.
  Description: 
   (str):gen      gp_readvec_file($1)
  Doc: reads in the file
   \var{filename} (subject to string expansion). If \var{filename} is
   omitted, re-reads the last file that was fed into \kbd{gp}. The return
   value is a vector whose components are the evaluation of all sequences
   of instructions contained in the file. For instance, if \var{file} contains
   \bprog
   1
   2
   3
   @eprog\noindent
   then we will get:
   \bprog
   ? \r a
   %1 = 1
   %2 = 2
   %3 = 3
   ? read(a)
   %4 = 3
   ? readvec(a)
   %5 = [1, 2, 3]
   @eprog
   In general a sequence is just a single line, but as usual braces and
   \kbd{\bs} may be used to enter multiline sequences.
  Variant: The underlying library function
   \fun{GEN}{gp_readvec_stream}{FILE *f} is usually more flexible.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gp_readvec_file(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def install(*argv):
  '''
  install
  Class: basic
  Section: programming/specific
  C-Name: gpinstall
  Prototype: vrrD"",r,D"",s,
  Help: install(name,code,{gpname},{lib}): load from dynamic library 'lib' the
   function 'name'. Assign to it the name 'gpname' in this GP session, with
   prototype 'code'. If 'lib' is omitted, all symbols known to gp
   (includes the whole 'libpari.so' and possibly others) are available.
   If 'gpname' is omitted, use 'name'.
  Doc: loads from dynamic library \var{lib} the function \var{name}. Assigns to it
   the name \var{gpname} in this \kbd{gp} session, with \emph{prototype}
   \var{code} (see below). If \var{gpname} is omitted, uses \var{name}.
   If \var{lib} is omitted, all symbols known to \kbd{gp} are available: this
   includes the whole of \kbd{libpari.so} and possibly others (such as
   \kbd{libc.so}).
   
   Most importantly, \kbd{install} gives you access to all non-static functions
   defined in the PARI library. For instance, the function \kbd{GEN addii(GEN
   x, GEN y)} adds two PARI integers, and is not directly accessible under
   \kbd{gp} (it is eventually called by the \kbd{+} operator of course):
   \bprog
   ? install("addii", "GG")
   ? addii(1, 2)
   %1 = 3
   @eprog\noindent
   It also allows to add external functions to the \kbd{gp} interpreter.
   For instance, it makes the function \tet{system} obsolete:
   \bprog
   ? install(system, vs, sys,/*omitted*/)
   ? sys("ls gp*")
   gp.c            gp.h            gp_rl.c
   @eprog\noindent This works because \kbd{system} is part of \kbd{libc.so},
   which is linked to \kbd{gp}. It is also possible to compile a shared library
   yourself and provide it to gp in this way: use \kbd{gp2c}, or do it manually
   (see the \kbd{modules\_build} variable in \kbd{pari.cfg} for hints).
   
   Re-installing a function will print a warning and update the prototype code
   if needed. However, it will not reload a symbol from the library, even if the
   latter has been recompiled.
   
   \misctitle{Prototype} We only give a simplified description here, covering
   most functions, but there are many more possibilities. The full documentation
   is available in \kbd{libpari.dvi}, see
   \bprog
     ??prototype
   @eprog
   
   \item First character \kbd{i}, \kbd{l}, \kbd{v} : return type int / long /
   void. (Default: \kbd{GEN})
   
   \item One letter for each mandatory argument, in the same order as they appear
   in the argument list: \kbd{G} (\kbd{GEN}), \kbd{\&}
   (\kbd{GEN*}), \kbd{L} (\kbd{long}), \kbd{s} (\kbd{char *}), \kbd{n}
   (variable).
   
    \item \kbd{p} to supply \kbd{realprecision} (usually \kbd{long prec} in the
    argument list), \kbd{P} to supply \kbd{seriesprecision} (usually \kbd{long
    precdl}).
   
    \noindent We also have special constructs for optional arguments and default
    values:
   
    \item \kbd{DG} (optional \kbd{GEN}, \kbd{NULL} if omitted),
   
    \item \kbd{D\&} (optional \kbd{GEN*}, \kbd{NULL} if omitted),
   
    \item \kbd{Dn} (optional variable, $-1$ if omitted),
   
   For instance the prototype corresponding to
   \bprog
     long issquareall(GEN x, GEN *n = NULL)
   @eprog\noindent is \kbd{lGD\&}.
   
   \misctitle{Caution} This function may not work on all systems, especially
   when \kbd{gp} has been compiled statically. In that case, the first use of an
   installed function will provoke a Segmentation Fault (this should never
   happen with a dynamically linked executable). If you intend to use this
   function, please check first on some harmless example such as the one above
   that it works properly on your machine.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.gpinstall(*c_arg_tuple)

def print1(*argv):
  '''
  print1
  Class: basic
  Section: programming/specific
  C-Name: print1
  Prototype: vs*
  Help: print1({str}*): outputs its string arguments (in raw format) without
   ending with newline.
  Description: 
   (?gen,...):void  pari_printf("${2 format_string}"${2 format_args})
  Doc: outputs its (string) arguments in raw
   format, without ending with a newline. Note that you can still embed newlines
   within your strings, using the \b{n} notation~!
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.print1(*c_arg_tuple)

def printf(*argv):
  '''
  printf
  Class: basic
  Section: programming/specific
  C-Name: printf0
  Prototype: vss*
  Help: printf(fmt,{x}*): prints its arguments according to the format fmt.
  Doc: This function is based on the C library command of the same name.
   It prints its arguments according to the format \var{fmt}, which specifies how
   subsequent arguments are converted for output. The format is a
   character string composed of zero or more directives:
   
   \item ordinary characters (not \kbd{\%}), printed unchanged,
   
   \item conversions specifications (\kbd{\%} followed by some characters)
   which fetch one argument from the list and prints it according to the
   specification.
   
   More precisely, a conversion specification consists in a \kbd{\%}, one or more
   optional flags (among \kbd{\#}, \kbd{0}, \kbd{-}, \kbd{+}, ` '), an optional
   decimal digit string specifying a minimal field width, an optional precision
   in the form of a period (`\kbd{.}') followed by a decimal digit string, and
   the conversion specifier (among \kbd{d},\kbd{i}, \kbd{o}, \kbd{u},
   \kbd{x},\kbd{X}, \kbd{p}, \kbd{e},\kbd{E}, \kbd{f}, \kbd{g},\kbd{G}, \kbd{s}).
   
   \misctitle{The flag characters} The character \kbd{\%} is followed by zero or
   more of the following flags:
   
   \item \kbd{\#}: The value is converted to an ``alternate form''. For
   \kbd{o} conversion (octal), a \kbd{0} is prefixed to the string. For \kbd{x}
   and \kbd{X} conversions (hexa), respectively \kbd{0x} and \kbd{0X} are
   prepended. For other conversions, the flag is ignored.
   
   \item \kbd{0}: The value should be zero padded. For
   \kbd{d},
   \kbd{i},
   \kbd{o},
   \kbd{u},
   \kbd{x},
   \kbd{X}
   \kbd{e},
   \kbd{E},
   \kbd{f},
   \kbd{F},
   \kbd{g}, and
   \kbd{G} conversions, the value is padded on the left with zeros rather than
   blanks. (If the \kbd{0} and \kbd{-} flags both appear, the \kbd{0} flag is
   ignored.)
   
   \item \kbd{-}: The value is left adjusted on the field boundary. (The
   default is right justification.) The value is padded on the right with
   blanks, rather than on the left with blanks or zeros. A \kbd{-} overrides a
   \kbd{0} if both are given.
   
   \item \kbd{` '} (a space): A blank is left before a positive number
   produced by a signed conversion.
   
   \item \kbd{+}: A sign (+ or -) is placed before a number produced by a
   signed conversion. A \kbd{+} overrides a space if both are used.
   
   \misctitle{The field width} An optional decimal digit string (whose first
   digit is non-zero) specifying a \emph{minimum} field width. If the value has
   fewer characters than the field width, it is padded with spaces on the left
   (or right, if the left-adjustment flag has been given). In no case does a
   small field width cause truncation of a field; if the value is wider than
   the field width, the field is expanded to contain the conversion result.
   Instead of a decimal digit string, one may write \kbd{*} to specify that the
   field width is given in the next argument.
   
   \misctitle{The precision} An optional precision in the form of a period
   (`\kbd{.}') followed by a decimal digit string. This gives
   the number of digits to appear after the radix character for \kbd{e},
   \kbd{E}, \kbd{f}, and \kbd{F} conversions, the maximum number of significant
   digits for \kbd{g} and \kbd{G} conversions, and the maximum number of
   characters to be printed from an \kbd{s} conversion.
   Instead of a decimal digit string, one may write \kbd{*} to specify that the
   field width is given in the next argument.
   
   \misctitle{The length modifier} This is ignored under \kbd{gp}, but
   necessary for \kbd{libpari} programming. Description given here for
   completeness:
   
   \item \kbd{l}: argument is a \kbd{long} integer.
   
   \item \kbd{P}: argument is a \kbd{GEN}.
   
   \misctitle{The conversion specifier} A character that specifies the type of
   conversion to be applied.
   
   \item \kbd{d}, \kbd{i}: A signed integer.
   
   \item \kbd{o}, \kbd{u}, \kbd{x}, \kbd{X}: An unsigned integer, converted
   to unsigned octal (\kbd{o}), decimal (\kbd{u}) or hexadecimal (\kbd{x} or
   \kbd{X}) notation. The letters \kbd{abcdef} are used for \kbd{x}
   conversions;  the letters \kbd{ABCDEF} are used for \kbd{X} conversions.
   
   \item \kbd{e}, \kbd{E}: The (real) argument is converted in the style
   \kbd{[ -]d.ddd e[ -]dd}, where there is one digit before the decimal point,
   and the number of digits after it is equal to the precision; if the
   precision is missing, use the current \kbd{realprecision} for the total
   number of printed digits. If the precision is explicitly 0, no decimal-point
   character appears. An \kbd{E} conversion uses the letter \kbd{E} rather
   than \kbd{e} to introduce the exponent.
   
   \item \kbd{f}, \kbd{F}: The (real) argument is converted in the style
   \kbd{[ -]ddd.ddd}, where the number of digits after the decimal point
   is equal to the precision; if the precision is missing, use the current
   \kbd{realprecision} for the total number of printed digits. If the precision
   is explicitly 0, no decimal-point character appears. If a decimal point
   appears, at least one digit appears before it.
   
   \item \kbd{g}, \kbd{G}: The (real) argument is converted in style
   \kbd{e} or \kbd{f} (or \kbd{E} or \kbd{F} for \kbd{G} conversions)
   \kbd{[ -]ddd.ddd}, where the total number of digits printed
   is equal to the precision; if the precision is missing, use the current
   \kbd{realprecision}. If the precision is explicitly 0, it is treated as 1.
   Style \kbd{e} is used when
   the decimal exponent is $< -4$, to print \kbd{0.}, or when the integer
   part cannot be decided given the known significant digits, and the \kbd{f}
   format otherwise.
   
   \item \kbd{c}: The integer argument is converted to an unsigned char, and the
   resulting character is written.
   
   \item \kbd{s}: Convert to a character string. If a precision is given, no
   more than the specified number of characters are written.
   
   \item \kbd{p}: Print the address of the argument in hexadecimal (as if by
   \kbd{\%\#x}).
   
   \item \kbd{\%}: A \kbd{\%} is written. No argument is converted. The complete
   conversion specification is \kbd{\%\%}.
   
   \noindent Examples:
   
   \bprog
   ? printf("floor: %d, field width 3: %3d, with sign: %+3d\n", Pi, 1, 2);
   floor: 3, field width 3:   1, with sign:  +2
   
   ? printf("%.5g %.5g %.5g\n",123,123/456,123456789);
   123.00 0.26974 1.2346 e8
   
   ? printf("%-2.5s:%2.5s:%2.5s\n", "P", "PARI", "PARIGP");
   P :PARI:PARIG
   
   \\ min field width and precision given by arguments
   ? x = 23; y=-1/x; printf("x=%+06.2f y=%+0*.*f\n", x, 6, 2, y);
   x=+23.00 y=-00.04
   
   \\ minimum fields width 5, pad left with zeroes
   ? for (i = 2, 5, printf("%05d\n", 10^i))
   00100
   01000
   10000
   100000  \\@com don't truncate fields whose length is larger than the minimum width
   ? printf("%.2f  |%06.2f|", Pi,Pi)
   3.14  |  3.14|
   @eprog\noindent All numerical conversions apply recursively to the entries
   of vectors and matrices:
   \bprog
   ? printf("%4d", [1,2,3]);
   [   1,   2,   3]
   ? printf("%5.2f", mathilbert(3));
   [ 1.00  0.50  0.33]
   
   [ 0.50  0.33  0.25]
   
   [ 0.33  0.25  0.20]
   @eprog
   \misctitle{Technical note} Our implementation of \tet{printf}
   deviates from the C89 and C99 standards in a few places:
   
   \item whenever a precision is missing, the current \kbd{realprecision} is
   used to determine the number of printed digits (C89: use 6 decimals after
   the radix character).
   
   \item in conversion style \kbd{e}, we do not impose that the
   exponent has at least two digits; we never write a \kbd{+} sign in the
   exponent; 0 is printed in a special way, always as \kbd{0.E\var{exp}}.
   
   \item in conversion style \kbd{f}, we switch to style \kbd{e} if the
   exponent is greater or equal to the precision.
   
   \item in conversion \kbd{g} and \kbd{G}, we do not remove trailing zeros
    from the fractional part of the result; nor a trailing decimal point;
    0 is printed in a special way, always as \kbd{0.E\var{exp}}.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.printf0(*c_arg_tuple)

def printsep(*argv):
  '''
  printsep
  Class: basic
  Section: programming/specific
  C-Name: printsep
  Prototype: vss*
  Help: printsep(sep,{str}*): outputs its string arguments (in raw format),
   separated by 'sep', ending with a newline.
  Doc: outputs its (string) arguments in raw format, ending with a newline.
   Successive entries are separated by \var{sep}:
   \bprog
   ? printsep(":", 1,2,3,4)
   1:2:3:4
   @eprog
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.printsep(*c_arg_tuple)

def printsep1(*argv):
  '''
  printsep1
  Class: basic
  Section: programming/specific
  C-Name: printsep1
  Prototype: vss*
  Help: printsep(sep,{str}*): outputs its string arguments (in raw format),
   separated by 'sep', without ending with a newline.
  Doc: outputs its (string) arguments in raw format, without ending with a
   newline.  Successive entries are separated by \var{sep}:
   \bprog
   ? printsep1(":", 1,2,3,4);print("|")
   1:2:3:4
   @eprog
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.printsep1(*c_arg_tuple)

def printtex(*argv):
  '''
  printtex
  Class: basic
  Section: programming/specific
  C-Name: printtex
  Prototype: vs*
  Help: printtex({str}*): outputs its string arguments in TeX format.
  Doc: outputs its (string) arguments in \TeX\ format. This output can then be
   used in a \TeX\ manuscript.
   The printing is done on the standard output. If you want to print it to a
   file you should use \kbd{writetex} (see there).
   
   Another possibility is to enable the \tet{log} default
   (see~\secref{se:defaults}).
   You could for instance do:\sidx{logfile}
   %
   \bprog
   default(logfile, "new.tex");
   default(log, 1);
   printtex(result);
   @eprog
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.printtex(*c_arg_tuple)

def Strprintf(*argv):
  '''
  Strprintf
  Class: basic
  Section: programming/specific
  C-Name: Strprintf
  Prototype: ss*
  Help: Strprintf(fmt,{x}*): returns a string built from the remaining
   arguments according to the format fmt.
  Doc: returns a string built from the remaining arguments according to the
   format fmt. The format consists of ordinary characters (not \%), printed
   unchanged, and conversions specifications. See \kbd{printf}.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.Strprintf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def warning(*argv):
  '''
  warning
  Class: basic
  Section: programming/specific
  C-Name: warning0
  Prototype: vs*
  Help: warning({str}*): display warning message str
  Description: 
   (?gen,...):void  pari_warn(warnuser, "${2 format_string}"${2 format_args})
  Doc: outputs the message ``user warning''
   and the argument list (each of them interpreted as a string).
   If colors are enabled, this warning will be in a different color,
   making it easy to distinguish.
   \bprog
   warning(n, " is very large, this might take a while.")
   @eprog
   % \syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.warning0(*c_arg_tuple)

def write(*argv):
  '''
  write
  Class: basic
  Section: programming/specific
  C-Name: write0
  Prototype: vss*
  Help: write(filename,{str}*): appends the remaining arguments (same output as
   print) to filename.
  Doc: writes (appends) to \var{filename} the remaining arguments, and appends a
   newline (same output as \kbd{print}).
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.write0(*c_arg_tuple)

def write1(*argv):
  '''
  write1
  Class: basic
  Section: programming/specific
  C-Name: write1
  Prototype: vss*
  Help: write1(filename,{str}*): appends the remaining arguments (same output as
   print1) to filename.
  Doc: writes (appends) to \var{filename} the remaining arguments without a
   trailing newline (same output as \kbd{print1}).
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.write1(*c_arg_tuple)

def writetex(*argv):
  '''
  writetex
  Class: basic
  Section: programming/specific
  C-Name: writetex
  Prototype: vss*
  Help: writetex(filename,{str}*): appends the remaining arguments (same format as
   print) to filename, in TeX format.
  Doc: as \kbd{write}, in \TeX\ format.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.writetex(*c_arg_tuple)

def next(*argv):
  '''
  next
  Class: basic
  Section: programming/control
  C-Name: next0
  Prototype: D1,L,
  Help: next({n=1}): interrupt execution of current instruction sequence, and
   start another iteration from the n-th innermost enclosing loops.
  Doc: interrupts execution of current $seq$,
   resume the next iteration of the innermost enclosing loop, within the
   current function call (or top level loop). If $n$ is specified, resume at
   the $n$-th enclosing loop. If $n$ is bigger than the number of enclosing
   loops, all enclosing loops are exited.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.next0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def pareval(*argv):
  '''
  pareval
  Class: basic
  Section: programming/parallel
  C-Name: pareval
  Prototype: G
  Help: pareval(x): parallel evaluation of the elements of the vector of
   closures x.
  Doc: parallel evaluation of the elements of \kbd{x}, where \kbd{x} is a
   vector of closures. The closures must be of arity $0$, must not access
   global variables or variables declared with \kbd{local} and must be
   free of side effects.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.pareval(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def parsum(*argv):
  '''
  parsum
  Class: basic
  Section: programming/parallel
  C-Name: parsum
  Prototype: V=GGJDG
  Help: parsum(i=a,b,expr,{x}): x plus the sum (X goes from a to b) of
   expression expr, evaluated in parallel (in random order)
  Description: 
   (gen,gen,closure,?gen):gen parsum($1, $2, $3, $4)
  Doc: sum of expression \var{expr}, initialized at $x$, the formal parameter
   going from $a$ to $b$, evaluated in parallel in random order.
   The expression \kbd{expr} must not access global variables or
   variables declared with \kbd{local()}, and must be free of side effects.
   \bprog
   parsum(i=1,1000,ispseudoprime(2^prime(i)-1))
   @eprog
   returns the numbers of prime numbers among the first $1000$ Mersenne numbers.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.parsum(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def parvector(*argv):
  '''
  parvector
  Class: basic
  Section: programming/parallel
  C-Name: parvector
  Prototype: LVJ
  Help: parvector(N,i,expr): as vector(N,i,expr) but the evaluations of expr are
   done in parallel.
  Description: 
    (small,,closure):vec    parvector($1, $3)
  Doc: As \kbd{vector(N,i,expr)} but the evaluations of \kbd{expr} are done in
   parallel. The expression \kbd{expr} must not access global variables or
   variables declared with \kbd{local()}, and must be free of side effects.
   \bprog
   parvector(10,i,quadclassunit(2^(100+i)+1).no)
   @eprog\noindent
   computes the class numbers in parallel.
   %\syn{NO}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.parvector(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ffgen(*argv):
  '''
  ffgen
  Class: basic
  Section: number_theoretical
  C-Name: ffgen
  Prototype: GDn
  Help: ffgen(q,{v}): return a generator X mod P(X) for the finite field with
   q elements. If v is given, the variable name is used to display g, else the
   variable 'x' is used. Alternative syntax, q = P(X) an irreducible
   polynomial with t_INTMOD
   coefficients, return the generator X mod P(X) of the finite field defined
   by P. If v is given, the variable name is used to display g, else the
   variable of the polynomial P is used.
  Doc: return a \typ{FFELT} generator for the finite field with $q$ elements;
   $q = p^f$ must be a prime power. This functions computes an irreducible
   monic polynomial $P\in\F_p[X]$ of degree~$f$ (via \tet{ffinit}) and
   returns $g = X \pmod{P(X)}$. If \kbd{v} is given, the variable name is used
   to display $g$, else the variable $x$ is used.
   \bprog
   ? g = ffgen(8, 't);
   ? g.mod
   %2 = t^3 + t^2 + 1
   ? g.p
   %3 = 2
   ? g.f
   %4 = 3
   ? ffgen(6)
    ***   at top-level: ffgen(6)
    ***                 ^--------
    *** ffgen: not a prime number in ffgen: 6.
   @eprog\noindent Alternative syntax: instead of a prime power $q$, one may
   input directly the polynomial $P$ (monic, irreducible, with \typ{INTMOD}
   coefficients), and the function returns the generator $g = X \pmod{P(X)}$,
   inferring $p$ from the coefficients of $P$. If \kbd{v} is given, the
   variable name is used to display $g$, else the variable of the polynomial
   $P$ is used. If $P$ is not irreducible, we create an invalid object and
   behaviour of functions dealing with the resulting \typ{FFELT}
   is undefined; in fact, it is much more costly to test $P$ for
   irreducibility than it would be to produce it via \kbd{ffinit}.
  Variant: 
   To create a generator for a prime finite field, the function
   \fun{GEN}{p_to_GEN}{GEN p, long v} returns \kbd{1+ffgen(x*Mod(1,p),v)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ffgen(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def fflog(*argv):
  '''
  fflog
  Class: basic
  Section: number_theoretical
  C-Name: fflog
  Prototype: GGDG
  Help: fflog(x,g,{o}): return the discrete logarithm of the finite field
   element x in base g. If present, o must represents the multiplicative
   order of g. If no o is given, assume that g is a primitive root.
  Doc: discrete logarithm of the finite field element $x$ in base $g$, i.e.~
   an $e$ in $\Z$ such that $g^e = o$. If
   present, $o$ represents the multiplicative order of $g$, see
   \secref{se:DLfun}; the preferred format for
   this parameter is \kbd{[ord, factor(ord)]}, where \kbd{ord} is the
   order of $g$. It may be set as a side effect of calling \tet{ffprimroot}.
   
   If no $o$ is given, assume that $g$ is a primitive root. The result is
   undefined if $e$ does not exist. This function uses
   
   \item a combination of generic discrete log algorithms (see \tet{znlog})
   
   \item a cubic sieve index calculus algorithm for large fields of degree at
   least $5$.
   
   \item Coppersmith's algorithm for fields of characteristic at most $5$.
   
   \bprog
   ? t = ffgen(ffinit(7,5));
   ? o = fforder(t)
   %2 = 5602   \\@com \emph{not} a primitive root.
   ? fflog(t^10,t)
   %3 = 10
   ? fflog(t^10,t, o)
   %4 = 10
   ? g = ffprimroot(t, &o);
   ? o   \\ order is 16806, bundled with its factorization matrix
   %6 = [16806, [2, 1; 3, 1; 2801, 1]]
   ? fforder(g, o)
   %7 = 16806
   ? fflog(g^10000, g, o)
   %8 = 10000
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.fflog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def fforder(*argv):
  '''
  fforder
  Class: basic
  Section: number_theoretical
  C-Name: fforder
  Prototype: GDG
  Help: fforder(x,{o}): multiplicative order of the finite field element x.
   Optional o represents a multiple of the order of the element.
  Doc: multiplicative order of the finite field element $x$.  If $o$ is
   present, it represents a multiple of the order of the element,
   see \secref{se:DLfun}; the preferred format for
   this parameter is \kbd{[N, factor(N)]}, where \kbd{N} is the cardinality
   of the multiplicative group of the underlying finite field.
   \bprog
   ? t = ffgen(ffinit(nextprime(10^8), 5));
   ? g = ffprimroot(t, &o);  \\@com o will be useful!
   ? fforder(g^1000000, o)
   time = 0 ms.
   %5 = 5000001750000245000017150000600250008403
   ? fforder(g^1000000)
   time = 16 ms. \\@com noticeably slower, same result of course
   %6 = 5000001750000245000017150000600250008403
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.fforder(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ffprimroot(*argv):
  '''
  ffprimroot
  Class: basic
  Section: number_theoretical
  C-Name: ffprimroot
  Prototype: GD&
  Help: ffprimroot(x, {&o}): return a primitive root of the multiplicative group
   of the definition field of the finite field element x (not necessarily the
   same as the field generated by x). If present, o is set to [ord, fa], where
   ord is the order of the group, and fa its factorization
   (useful in fflog and fforder).
  Doc: return a primitive root of the multiplicative
   group of the definition field of the finite field element $x$ (not necessarily
   the same as the field generated by $x$). If present, $o$ is set to
   a vector \kbd{[ord, fa]}, where \kbd{ord} is the order of the group
   and \kbd{fa} its factorisation \kbd{factor(ord)}. This last parameter is
   useful in \tet{fflog} and \tet{fforder}, see \secref{se:DLfun}.
   \bprog
   ? t = ffgen(ffinit(nextprime(10^7), 5));
   ? g = ffprimroot(t, &o);
   ? o[1]
   %3 = 100000950003610006859006516052476098
   ? o[2]
   %4 =
   [2 1]
   
   [7 2]
   
   [31 1]
   
   [41 1]
   
   [67 1]
   
   [1523 1]
   
   [10498781 1]
   
   [15992881 1]
   
   [46858913131 1]
   
   ? fflog(g^1000000, g, o)
   time = 1,312 ms.
   %5 = 1000000
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ffprimroot(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfgaloisconj(*argv):
  '''
  nfgaloisconj
  Class: basic
  Section: number_fields
  C-Name: galoisconj0
  Prototype: GD0,L,DGp
  Help: nfgaloisconj(nf,{flag=0},{d}): list of conjugates of a root of the
   polynomial x=nf.pol in the same number field. flag is optional (set to 0 by
   default), meaning 0: use combination of flag 4 and 1, always complete; 1:
   use nfroots; 2 : use complex numbers, LLL on integral basis (not always
   complete); 4: use Allombert's algorithm, complete if the field is Galois of
   degree <= 35 (see manual for details). nf can be simply a polynomial.
  Doc: $\var{nf}$ being a number field as output by \kbd{nfinit}, computes the
   conjugates of a root $r$ of the non-constant polynomial $x=\var{nf}[1]$
   expressed as polynomials in $r$. This also makes sense when the number field
   is not \idx{Galois} since some conjugates may lie in the field.
   $\var{nf}$ can simply be a polynomial.
   
   If no flags or $\fl=0$, use a combination of flag $4$ and $1$ and the result
   is always complete. There is no point whatsoever in using the other flags.
   
   If $\fl=1$, use \kbd{nfroots}: a little slow, but guaranteed to work in
   polynomial time.
   
   If $\fl=2$ (OBSOLETE), use complex approximations to the roots and an integral
   \idx{LLL}. The result is not guaranteed to be complete: some
   conjugates may be missing (a warning is issued if the result is not proved
   complete), especially so if the corresponding polynomial has a huge index,
   and increasing the default precision may help. This variant is slow and
   unreliable: don't use it.
   
   If $\fl=4$, use \kbd{galoisinit}: very fast, but only applies to (most) Galois
   fields. If the field is Galois with weakly
   super-solvable Galois group (see \tet{galoisinit}), return the complete list
   of automorphisms, else only the identity element. If present, $d$ is assumed to
   be a multiple of the least common denominator of the conjugates expressed as
   polynomial in a root of \var{pol}.
   
   This routine can only compute $\Q$-automorphisms, but it may be used to get
   $K$-automorphism for any base field $K$ as follows:
   \bprog
   rnfgaloisconj(nfK, R) = \\ K-automorphisms of L = K[X] / (R)
   { my(polabs, N);
     R *= Mod(1, nfK.pol);             \\ convert coeffs to polmod elts of K
     polabs = rnfequation(nfK, R);
     N = nfgaloisconj(polabs) % R;     \\ Q-automorphisms of L
     \\ select the ones that fix K
     select(s->subst(R, variable(R), Mod(s,R)) == 0, N);
   }
   K  = nfinit(y^2 + 7);
   rnfgaloisconj(K, x^4 - y*x^3 - 3*x^2 + y*x + 1)  \\ K-automorphisms of L
   @eprog
  Variant: Use directly
   \fun{GEN}{galoisconj}{GEN nf, GEN d}, corresponding to $\fl = 0$, the others
   only have historical interest.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisconj0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisexport(*argv):
  '''
  galoisexport
  Class: basic
  Section: number_fields
  C-Name: galoisexport
  Prototype: GD0,L,
  Help: galoisexport(gal,{flag}): gal being a Galois group as output by
   galoisinit, output a string representing the underlying permutation group in
   GAP notation (default) or Magma notation (flag = 1).
  Doc: \var{gal} being be a Galois group as output by \tet{galoisinit},
   export the underlying permutation group as a string suitable
   for (no flags or $\fl=0$) GAP or ($\fl=1$) Magma. The following example
   compute the index of the underlying abstract group in the GAP library:
   \bprog
   ? G = galoisinit(x^6+108);
   ? s = galoisexport(G)
   %2 = "Group((1, 2, 3)(4, 5, 6), (1, 4)(2, 6)(3, 5))"
   ? extern("echo \"IdGroup("s");\" | gap -q")
   %3 = [6, 1]
   ? galoisidentify(G)
   %4 = [6, 1]
   @eprog\noindent
   This command also accepts subgroups returned by \kbd{galoissubgroups}.
   
   To \emph{import} a GAP permutation into gp (for \tet{galoissubfields} for
   instance), the following GAP function may be useful:
   \bprog
   PermToGP := function(p, n)
     return Permuted([1..n],p);
   end;
   
   gap> p:= (1,26)(2,5)(3,17)(4,32)(6,9)(7,11)(8,24)(10,13)(12,15)(14,27)
     (16,22)(18,28)(19,20)(21,29)(23,31)(25,30)
   gap> PermToGP(p,32);
   [ 26, 5, 17, 32, 2, 9, 11, 24, 6, 13, 7, 15, 10, 27, 12, 22, 3, 28, 20, 19,
     29, 16, 31, 8, 30, 1, 14, 18, 21, 25, 23, 4 ]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisexport(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisfixedfield(*argv):
  '''
  galoisfixedfield
  Class: basic
  Section: number_fields
  C-Name: galoisfixedfield
  Prototype: GGD0,L,Dn
  Help: galoisfixedfield(gal,perm,{flag},{v=y}): gal being a Galois group as
   output by galoisinit and perm a subgroup, an element of gal.group or a vector
   of such elements, return [P,x] such that P is a polynomial defining the fixed
   field of gal[1] by the subgroup generated by perm, and x is a root of P in gal
   expressed as a polmod in gal.pol. If flag is 1 return only P. If flag is 2
   return [P,x,F] where F is the factorization of gal.pol over the field
   defined by P, where the variable v stands for a root of P.
  Description: 
   (gen, gen, ?small, ?var):vec        galoisfixedfield($1, $2, $3, $4)
  Doc: \var{gal} being be a Galois group as output by \tet{galoisinit} and
   \var{perm} an element of $\var{gal}.group$, a vector of such elements
   or a subgroup of \var{gal} as returned by galoissubgroups,
   computes the fixed field of \var{gal} by the automorphism defined by the
   permutations \var{perm} of the roots $\var{gal}.roots$. $P$ is guaranteed to
   be squarefree modulo $\var{gal}.p$.
   
   If no flags or $\fl=0$, output format is the same as for \tet{nfsubfield},
   returning $[P,x]$ such that $P$ is a polynomial defining the fixed field, and
   $x$ is a root of $P$ expressed as a polmod in $\var{gal}.pol$.
   
   If $\fl=1$ return only the polynomial $P$.
   
   If $\fl=2$ return $[P,x,F]$ where $P$ and $x$ are as above and $F$ is the
   factorization of $\var{gal}.pol$ over the field defined by $P$, where
   variable $v$ ($y$ by default) stands for a root of $P$. The priority of $v$
   must be less than the priority of the variable of $\var{gal}.pol$ (see
   \secref{se:priority}). Example:
   
   \bprog
   ? G = galoisinit(x^4+1);
   ? galoisfixedfield(G,G.group[2],2)
   %2 = [x^2 + 2, Mod(x^3 + x, x^4 + 1), [x^2 - y*x - 1, x^2 + y*x - 1]]
   @eprog\noindent
   computes the factorization  $x^4+1=(x^2-\sqrt{-2}x-1)(x^2+\sqrt{-2}x-1)$
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisfixedfield(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisidentify(*argv):
  '''
  galoisidentify
  Class: basic
  Section: number_fields
  C-Name: galoisidentify
  Prototype: G
  Help: galoisidentify(gal): gal being a Galois group as output by galoisinit,
   output the isomorphism class of the underlying abstract group as a
   two-components vector [o,i], where o is the group order, and i is the group
   index in the GAP4 small group library.
  Doc: \var{gal} being be a Galois group as output by \tet{galoisinit},
   output the isomorphism class of the underlying abstract group as a
   two-components vector $[o,i]$, where $o$ is the group order, and $i$ is the
   group index in the GAP4 Small Group library, by Hans Ulrich Besche, Bettina
   Eick and Eamonn O'Brien.
   
   This command also accepts subgroups returned by \kbd{galoissubgroups}.
   
   The current implementation is limited to degree less or equal to $127$.
   Some larger ``easy'' orders are also supported.
   
   The output is similar to the output of the function \kbd{IdGroup} in GAP4.
   Note that GAP4 \kbd{IdGroup} handles all groups of order less than $2000$
   except $1024$, so you can use \tet{galoisexport} and GAP4 to identify large
   Galois groups.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisidentify(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisinit(*argv):
  '''
  galoisinit
  Class: basic
  Section: number_fields
  C-Name: galoisinit
  Prototype: GDG
  Help: galoisinit(pol,{den}): pol being a polynomial or a number field as
   output by nfinit defining a Galois extension of Q, compute the Galois group
   and all necessary information for computing fixed fields. den is optional
   and has the same meaning as in nfgaloisconj(,4)(see manual).
  Description: 
   (gen, ?int):gal        galoisinit($1, $2)
  Doc: computes the Galois group
   and all necessary information for computing the fixed fields of the
   Galois extension $K/\Q$ where $K$ is the number field defined by
   $\var{pol}$ (monic irreducible polynomial in $\Z[X]$ or
   a number field as output by \tet{nfinit}). The extension $K/\Q$ must be
   Galois with Galois group ``weakly'' super-solvable, see below;
   returns 0 otherwise. Hence this permits to quickly check whether a polynomial
   of order strictly less than $36$ is Galois or not.
   
   The algorithm used is an improved version of the paper
   ``An efficient algorithm for the computation of Galois automorphisms'',
   Bill Allombert, Math.~Comp, vol.~73, 245, 2001, pp.~359--375.
   
   A group $G$ is said to be ``weakly'' super-solvable if there exists a
   normal series
   
   $\{1\} = H_0 \triangleleft H_1 \triangleleft \cdots \triangleleft H_{n-1}
   \triangleleft H_n$
   
   such that each $H_i$ is normal in $G$ and for $i<n$, each quotient group
   $H_{i+1}/H_i$ is cyclic, and either $H_n=G$ (then $G$ is super-solvable) or
   $G/H_n$ is isomorphic to either $A_4$ or $S_4$.
   
   In practice, almost all small groups are WKSS, the exceptions having order
   36(1 exception), 48(2), 56(1), 60(1), 72(5), 75(1), 80(1), 96(10) and $\geq
   108$.
   
   This function is a prerequisite for most of the \kbd{galois}$xxx$ routines.
   For instance:
   
   \bprog
   P = x^6 + 108;
   G = galoisinit(P);
   L = galoissubgroups(G);
   vector(#L, i, galoisisabelian(L[i],1))
   vector(#L, i, galoisidentify(L[i]))
   @eprog
   
   The output is an 8-component vector \var{gal}.
   
   $\var{gal}[1]$ contains the polynomial \var{pol}
   (\kbd{\var{gal}.pol}).
   
   $\var{gal}[2]$ is a three-components vector $[p,e,q]$ where $p$ is a
   prime number (\kbd{\var{gal}.p}) such that \var{pol} totally split
   modulo $p$ , $e$ is an integer and $q=p^e$ (\kbd{\var{gal}.mod}) is the
   modulus of the roots in \kbd{\var{gal}.roots}.
   
   $\var{gal}[3]$ is a vector $L$ containing the $p$-adic roots of
   \var{pol} as integers implicitly modulo \kbd{\var{gal}.mod}.
   (\kbd{\var{gal}.roots}).
   
   $\var{gal}[4]$ is the inverse of the Vandermonde matrix of the
   $p$-adic roots of \var{pol}, multiplied by $\var{gal}[5]$.
   
   $\var{gal}[5]$ is a multiple of the least common denominator of the
   automorphisms expressed as polynomial in a root of \var{pol}.
   
   $\var{gal}[6]$ is the Galois group $G$ expressed as a vector of
   permutations of $L$ (\kbd{\var{gal}.group}).
   
   $\var{gal}[7]$ is a generating subset $S=[s_1,\ldots,s_g]$ of $G$
   expressed as a vector of permutations of $L$ (\kbd{\var{gal}.gen}).
   
   $\var{gal}[8]$ contains the relative orders $[o_1,\ldots,o_g]$ of
   the generators of $S$ (\kbd{\var{gal}.orders}).
   
   Let $H_n$ be as above, we have the following properties:
   
   \quad\item if $G/H_n\simeq A_4$ then $[o_1,\ldots,o_g]$ ends by
   $[2,2,3]$.
   
   \quad\item if $G/H_n\simeq S_4$ then $[o_1,\ldots,o_g]$ ends by
   $[2,2,3,2]$.
   
   \quad\item for $1\leq i \leq g$ the subgroup of $G$ generated by
   $[s_1,\ldots,s_g]$ is normal, with the exception of $i=g-2$ in the
   $A_4$ case and of $i=g-3$ in the $S_A$ case.
   
   \quad\item the relative order $o_i$ of $s_i$ is its order in the
   quotient group $G/\langle s_1,\ldots,s_{i-1}\rangle$, with the same
   exceptions.
   
   \quad\item for any $x\in G$ there exists a unique family
   $[e_1,\ldots,e_g]$ such that (no exceptions):
   
   -- for $1\leq i \leq g$ we have $0\leq e_i<o_i$
   
   -- $x=g_1^{e_1}g_2^{e_2}\ldots g_n^{e_n}$
   
   If present $den$ must be a suitable value for $\var{gal}[5]$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisisabelian(*argv):
  '''
  galoisisabelian
  Class: basic
  Section: number_fields
  C-Name: galoisisabelian
  Prototype: GD0,L,
  Help: galoisisabelian(gal,{flag=0}): gal being as output by galoisinit,
   return 0 if gal is not abelian, the HNF matrix of gal over gal.gen if
   flag=0, 1 if flag is 1, and the SNF of gal is flag=2.
  Doc: \var{gal} being as output by \kbd{galoisinit}, return $0$ if
   \var{gal} is not an abelian group, and the HNF matrix of \var{gal} over
   \kbd{gal.gen} if $fl=0$, $1$ if $fl=1$.
   
   This command also accepts subgroups returned by \kbd{galoissubgroups}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisisabelian(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisisnormal(*argv):
  '''
  galoisisnormal
  Class: basic
  Section: number_fields
  C-Name: galoisisnormal
  Prototype: lGG
  Help: galoisisnormal(gal,subgrp): gal being as output by galoisinit,
   and subgrp a subgroup of gal as output by galoissubgroups,
   return 1 if subgrp is a normal subgroup of gal, else return 0.
  Doc: \var{gal} being as output by \kbd{galoisinit}, and \var{subgrp} a subgroup
   of \var{gal} as output by \kbd{galoissubgroups},return $1$ if \var{subgrp} is a
   normal subgroup of \var{gal}, else return 0.
   
   This command also accepts subgroups returned by \kbd{galoissubgroups}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.galoisisnormal(*c_arg_tuple)

def galoispermtopol(*argv):
  '''
  galoispermtopol
  Class: basic
  Section: number_fields
  C-Name: galoispermtopol
  Prototype: GG
  Help: galoispermtopol(gal,perm): gal being a Galois group as output by
   galoisinit and perm a element of gal.group, return the polynomial defining
   the corresponding Galois automorphism.
  Doc: \var{gal} being a
   Galois group as output by \kbd{galoisinit} and \var{perm} a element of
   $\var{gal}.group$, return the polynomial defining the Galois
   automorphism, as output by \kbd{nfgaloisconj}, associated with the
   permutation \var{perm} of the roots $\var{gal}.roots$. \var{perm} can
   also be a vector or matrix, in this case, \kbd{galoispermtopol} is
   applied to all components recursively.
   
   \noindent Note that
   \bprog
   G = galoisinit(pol);
   galoispermtopol(G, G[6])~
   @eprog\noindent
   is equivalent to \kbd{nfgaloisconj(pol)}, if degree of \var{pol} is greater
   or equal to $2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoispermtopol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoissubgroups(*argv):
  '''
  galoissubgroups
  Class: basic
  Section: number_fields
  C-Name: galoissubgroups
  Prototype: G
  Help: galoissubgroups(G):Output all the subgroups of G.
  Doc: outputs all the subgroups of the Galois group \kbd{gal}. A subgroup is a
   vector [\var{gen}, \var{orders}], with the same meaning
   as for $\var{gal}.gen$ and $\var{gal}.orders$. Hence \var{gen} is a vector of
   permutations generating the subgroup, and \var{orders} is the relatives
   orders of the generators. The cardinal of a subgroup is the product of the
   relative orders. Such subgroup can be used instead of a Galois group in the
   following command: \kbd{galoisisabelian}, \kbd{galoissubgroups},
   \kbd{galoisexport} and \kbd{galoisidentify}.
   
   To get the subfield fixed by a subgroup \var{sub} of \var{gal}, use
   \bprog
   galoisfixedfield(gal,sub[1])
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoissubgroups(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoissubfields(*argv):
  '''
  galoissubfields
  Class: basic
  Section: number_fields
  C-Name: galoissubfields
  Prototype: GD0,L,Dn
  Help: galoissubfields(G,{flags=0},{v}):Output all the subfields of G. flags
   have the same meaning as for galoisfixedfield.
  Doc: outputs all the subfields of the Galois group \var{G}, as a vector.
   This works by applying \kbd{galoisfixedfield} to all subgroups. The meaning of
   the flag \var{fl} is the same as for \kbd{galoisfixedfield}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoissubfields(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoisgetpol(*argv):
  '''
  galoisgetpol
  Class: basic
  Section: number_fields
  C-Name: galoisgetpol
  Prototype: LD0,L,D1,L,
  Help: galoisgetpol(a,{b},{s}): Query the galpol package for a polynomial with
   Galois group isomorphic to GAP4(a,b), totally real if s=1 (default) and
   totally complex if s=2.  The output is a vector [pol, den] where pol is the
   polynomial and den is the common denominator of the conjugates expressed
   as a polynomial in a root of pol. If b and s are omitted, return the number of
   isomorphism classes of groups of order a.
  Description: 
   (small):int               galoisnbpol($1)
   (small,):int              galoisnbpol($1)
   (small,,):int             galoisnbpol($1)
   (small,small,small):vec   galoisgetpol($1, $2 ,$3)
  Doc: Query the galpol package for a polynomial with Galois group isomorphic to
   GAP4(a,b), totally real if $s=1$ (default) and totally complex if $s=2$. The
   output is a vector [\kbd{pol}, \kbd{den}] where
   
   \item  \kbd{pol} is the polynomial of degree $a$
   
   \item \kbd{den} is the denominator of \kbd{nfgaloisconj(pol)}.
   Pass it as an optional argument to \tet{galoisinit} or \tet{nfgaloisconj} to
   speed them up:
   \bprog
   ? [pol,den] = galoisgetpol(64,4,1);
   ? G = galoisinit(pol);
   time = 352ms
   ? galoisinit(pol, den);  \\ passing 'den' speeds up the computation
   time = 264ms
   ? % == %`
   %4 = 1  \\ same answer
   @eprog
   If $b$ and $s$ are omitted, return the number of isomorphism classes of
   groups of order $a$.
  Variant: Also available is \fun{GEN}{galoisnbpol}{long a} when $b$ and $s$
   are omitted.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoisgetpol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def conjvec(*argv):
  '''
  conjvec
  Class: basic
  Section: conversions
  C-Name: conjvec
  Prototype: Gp
  Help: conjvec(z): conjugate vector of the algebraic number z.
  Doc: 
   conjugate vector representation of $z$. If $z$ is a
   polmod, equal to \kbd{Mod}$(a,T)$, this gives a vector of length
   $\text{degree}(T)$ containing:
   
   \item the complex embeddings of $z$ if $T$ has rational coefficients,
   i.e.~the $a(r[i])$ where $r = \kbd{polroots}(T)$;
   
   \item the conjugates of $z$ if $T$ has some intmod coefficients;
   
   \noindent if $z$ is a finite field element, the result is the vector of
   conjugates $[z,z^p,z^{p^2},\ldots,z^{p^{n-1}}]$ where $n=\text{degree}(T)$.
   
   \noindent If $z$ is an integer or a rational number, the result is~$z$. If
   $z$ is a (row or column) vector, the result is a matrix whose columns are
   the conjugate vectors of the individual elements of $z$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.conjvec(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def conj(*argv):
  '''
  conj
  Class: basic
  Section: conversions
  C-Name: gconj
  Prototype: G
  Help: conj(x): the algebraic conjugate of x.
  Doc: 
   conjugate of $x$. The meaning of this
   is clear, except that for real quadratic numbers, it means conjugation in the
   real quadratic field. This function has no effect on integers, reals,
   intmods, fractions or $p$-adics. The only forbidden type is polmod
   (see \kbd{conjvec} for this).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gconj(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def shiftmul(*argv):
  '''
  shiftmul
  Class: basic
  Section: operators
  C-Name: gmul2n
  Prototype: GL
  Help: shiftmul(x,n): multiply x by 2^n (n>=0 or n<0)
  Doc: multiplies $x$ by $2^n$. The difference with
   \kbd{shift} is that when $n<0$, ordinary division takes place, hence for
   example if $x$ is an integer the result may be a fraction, while for shifts
   Euclidean division takes place when $n<0$ hence if $x$ is an integer the result
   is still an integer.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gmul2n(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sqr(*argv):
  '''
  sqr
  Class: basic
  Section: transcendental
  C-Name: gsqr
  Prototype: G
  Help: sqr(x): square of x. NOT identical to x*x.
  Description: 
   (int):int        sqri($1)
   (mp):mp          gsqr($1)
   (gen):gen        gsqr($1)
  Doc: square of $x$. This operation is not completely
   straightforward, i.e.~identical to $x * x$, since it can usually be
   computed more efficiently (roughly one-half of the elementary
   multiplications can be saved). Also, squaring a $2$-adic number increases
   its precision. For example,
   \bprog
   ? (1 + O(2^4))^2
   %1 = 1 + O(2^5)
   ? (1 + O(2^4)) * (1 + O(2^4))
   %2 = 1 + O(2^4)
   @eprog\noindent
   Note that this function is also called whenever one multiplies two objects
   which are known to be \emph{identical}, e.g.~they are the value of the same
   variable, or we are computing a power.
   \bprog
   ? x = (1 + O(2^4)); x * x
   %3 = 1 + O(2^5)
   ? (1 + O(2^4))^4
   %4 = 1 + O(2^6)
   @eprog\noindent
   (note the difference between \kbd{\%2} and \kbd{\%3} above).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsqr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def abs(*argv):
  '''
  abs
  Class: basic
  Section: transcendental
  C-Name: gabs
  Prototype: Gp
  Help: abs(x): absolute value (or modulus) of x.
  Description: 
   (small):small    labs($1)
   (int):int        mpabs($1)
   (real):real      mpabs($1)
   (mp):mp          mpabs($1)
   (gen):gen:prec        gabs($1, prec)
  Doc: absolute value of $x$ (modulus if $x$ is complex).
   Rational functions are not allowed. Contrary to most transcendental
   functions, an exact argument is \emph{not} converted to a real number before
   applying \kbd{abs} and an exact result is returned if possible.
   \bprog
   ? abs(-1)
   %1 = 1
   ? abs(3/7 + 4/7*I)
   %2 = 5/7
   ? abs(1 + I)
   %3 = 1.414213562373095048801688724
   @eprog\noindent
   If $x$ is a polynomial, returns $-x$ if the leading coefficient is
   real and negative else returns $x$. For a power series, the constant
   coefficient is considered instead.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gabs(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def valuation(*argv):
  '''
  valuation
  Class: basic
  Section: conversions
  C-Name: gvaluation
  Prototype: lGG
  Help: valuation(x,p): valuation of x with respect to p.
  Doc: 
   computes the highest
   exponent of $p$ dividing $x$. If $p$ is of type integer, $x$ must be an
   integer, an intmod whose modulus is divisible by $p$, a fraction, a
   $q$-adic number with $q=p$, or a polynomial or power series in which case the
   valuation is the minimum of the valuation of the coefficients.
   
   If $p$ is of type polynomial, $x$ must be of type polynomial or rational
   function, and also a power series if $x$ is a monomial. Finally, the
   valuation of a vector, complex or quadratic number is the minimum of the
   component valuations.
   
   If $x=0$, the result is \tet{LONG_MAX} ($2^{31}-1$ for 32-bit machines or
   $2^{63}-1$ for 64-bit machines) if $x$ is an exact object. If $x$ is a
   $p$-adic numbers or power series, the result is the exponent of the zero.
   Any other type combinations gives an error.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.gvaluation(*c_arg_tuple)

def length(*argv):
  '''
  length
  Class: basic
  Section: conversions
  C-Name: glength
  Prototype: lG
  Help: length(x): number of non code words in x, number of characters for a
   string.
  Description: 
   (vecsmall):lg      lg($1)
   (vec):lg           lg($1)
   (pol):small        lgpol($1)
   (gen):small        glength($1)
  Doc: length of $x$; \kbd{\#}$x$ is a shortcut for \kbd{length}$(x)$.
   This is mostly useful for
   
   \item vectors: dimension (0 for empty vectors),
   
   \item lists: number of entries (0 for empty lists),
   
   \item matrices: number of columns,
   
   \item character strings: number of actual characters (without
   trailing \kbd{\bs 0}, should you expect it from $C$ \kbd{char*}).
   \bprog
    ? #"a string"
    %1 = 8
    ? #[3,2,1]
    %2 = 3
    ? #[]
    %3 = 0
    ? #matrix(2,5)
    %4 = 5
    ? L = List([1,2,3,4]); #L
    %5 = 4
   @eprog
   
   The routine is in fact defined for arbitrary GP types, but is awkward and
   useless in other cases: it returns the number of non-code words in $x$, e.g.
   the effective length minus 2 for integers since the \typ{INT} type has two code
   words.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.glength(*c_arg_tuple)

def max(*argv):
  '''
  max
  Class: basic
  Section: operators
  C-Name: gmax
  Prototype: GG
  Help: max(x,y): maximum of x and y
  Description: 
   (small, small):small  maxss($1, $2)
   (small, int):int      gmaxsg($1, $2)
   (int, small):int      gmaxgs($1, $2)
   (int, int):int        gmax($1, $2)
   (small, mp):mp        gmaxsg($1, $2)
   (mp, small):mp        gmaxgs($1, $2)
   (mp, mp):mp           gmax($1, $2)
   (small, gen):gen      gmaxsg($1, $2)
   (gen, small):gen      gmaxgs($1, $2)
   (gen, gen):gen        gmax($1, $2)
  Doc: creates the maximum of $x$ and $y$ when they can be compared.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gmax(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def min(*argv):
  '''
  min
  Class: basic
  Section: operators
  C-Name: gmin
  Prototype: GG
  Help: min(x,y): minimum of x and y
  Description: 
   (small, small):small  minss($1, $2)
   (small, int):int      gminsg($1, $2)
   (int, small):int      gmings($1, $2)
   (int, int):int        gmin($1, $2)
   (small, mp):mp        gminsg($1, $2)
   (mp, small):mp        gmings($1, $2)
   (mp, mp):mp           gmin($1, $2)
   (small, gen):gen      gminsg($1, $2)
   (gen, small):gen      gmings($1, $2)
   (gen, gen):gen        gmin($1, $2)
  Doc: creates the minimum of $x$ and $y$ when they can be compared.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gmin(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sign(*argv):
  '''
  sign
  Class: basic
  Section: operators
  C-Name: gsigne
  Prototype: iG
  Help: sign(x): sign of x, of type integer, real or fraction
  Description: 
   (mp):small          signe($1)
   (gen):small        gsigne($1)
  Doc: \idx{sign} ($0$, $1$ or $-1$) of $x$, which must be of
   type integer, real or fraction.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.gsigne(*c_arg_tuple)

def List(*argv):
  '''
  List
  Class: basic
  Section: conversions
  C-Name: gtolist
  Prototype: DG
  Help: List({x=[]}): transforms the vector or list x into a list. Empty list
   if x is omitted.
  Description: 
   ():list           listcreate()
   (gen):list        gtolist($1)
  Doc: 
   transforms a (row or column) vector $x$ into a list, whose components are
   the entries of $x$. Similarly for a list, but rather useless in this case.
   For other types, creates a list with the single element $x$. Note that,
   except when $x$ is omitted, this function creates a small memory leak; so,
   either initialize all lists to the empty list, or use them sparingly.
  Variant: The variant \fun{GEN}{listcreate}{void} creates an empty list.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtolist(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def lex(*argv):
  '''
  lex
  Class: basic
  Section: operators
  C-Name: lexcmp
  Prototype: iGG
  Help: lex(x,y): compare x and y lexicographically (1 if x>y, 0 if x=y, -1 if
   x<y)
  Doc: gives the result of a lexicographic comparison
   between $x$ and $y$ (as $-1$, $0$ or $1$). This is to be interpreted in quite
   a wide sense: It is admissible to compare objects of different types
   (scalars, vectors, matrices), provided the scalars can be compared, as well
   as vectors/matrices of different lengths. The comparison is recursive.
   
   In case all components are equal up to the smallest length of the operands,
   the more complex is considered to be larger. More precisely, the longest is
   the largest; when lengths are equal, we have matrix $>$ vector $>$ scalar.
   For example:
   \bprog
   ? lex([1,3], [1,2,5])
   %1 = 1
   ? lex([1,3], [1,3,-1])
   %2 = -1
   ? lex([1], [[1]])
   %3 = -1
   ? lex([1], [1]~)
   %4 = 0
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.lexcmp(*c_arg_tuple)

def listinsert(*argv):
  '''
  listinsert
  Class: basic
  Section: linear_algebra
  C-Name: listinsert
  Prototype: WGL
  Help: listinsert(L,x,n): insert x at index n in list L, shifting the
   remaining elements to the right.
  Description: 
   (list, gen, small):gen        listinsert($1, $2, $3)
  Doc: inserts the object $x$ at
   position $n$ in $L$ (which must be of type \typ{LIST}). This has
   complexity $O(\#L - n + 1)$: all the
   remaining elements of \var{list} (from position $n+1$ onwards) are shifted
   to the right.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.listinsert(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def listpop(*argv):
  '''
  listpop
  Class: basic
  Section: linear_algebra
  C-Name: listpop
  Prototype: vWD0,L,
  Help: listpop(list,{n}): removes n-th element from list. If n is
   omitted or greater than the current list length, removes last element.
  Description: 
   (list, small):void     listpop($1, $2)
  Doc: 
   removes the $n$-th element of the list
   \var{list} (which must be of type \typ{LIST}). If $n$ is omitted,
   or greater than the list current length, removes the last element.
   If the list is already empty, do nothing. This runs in time $O(\#L - n + 1)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.listpop(*c_arg_tuple)

def listput(*argv):
  '''
  listput
  Class: basic
  Section: linear_algebra
  C-Name: listput
  Prototype: WGD0,L,
  Help: listput(list,x,{n}): sets n-th element of list equal to x. If n is
   omitted or greater than the current list length, appends x.
  Description: 
   (list, gen, small):gen        listput($1, $2, $3)
  Doc: 
   sets the $n$-th element of the list
   \var{list} (which must be of type \typ{LIST}) equal to $x$. If $n$ is omitted,
   or greater than the list length, appends $x$.
   You may put an element into an occupied cell (not changing the
   list length), but it is easier to use the standard \kbd{list[n] = x}
   construct. This runs in time $O(\#L)$ in the worst case (when the list must
   be reallocated), but in time $O(1)$ on average: any number of successive
   \kbd{listput}s run in time $O(\#L)$, where $\#L$ denotes the list
   \emph{final} length.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.listput(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def listsort(*argv):
  '''
  listsort
  Class: basic
  Section: linear_algebra
  C-Name: listsort
  Prototype: vWD0,L,
  Help: listsort(L,{flag=0}): sort the list L in place. If flag is non-zero,
   suppress all but one occurence of each element in list.
  Doc: sorts the \typ{LIST} \var{list} in place, with respect to the (somewhat
   arbitrary) universal comparison function \tet{cmp}. In particular, the
   ordering is the same as for sets and \tet{setsearch} can be used on a sorted
   list.
   \bprog
   ? L = List([1,2,4,1,3,-1]); listsort(L); L
   %1 = List([-1, 1, 1, 2, 3, 4])
   ? setsearch(L, 4)
   %2 = 6
   ? setsearch(L, -2)
   %3 = 0
   @eprog\noindent This is faster than the \kbd{vecsort} command since the list
   is sorted in place: no copy is made. No value returned.
   
   If $\fl$ is non-zero, suppresses all repeated coefficients.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  PyPari.pari.listsort(*c_arg_tuple)

def matsize(*argv):
  '''
  matsize
  Class: basic
  Section: linear_algebra
  C-Name: matsize
  Prototype: G
  Help: matsize(x): number of rows and columns of the vector/matrix x as a
   2-vector.
  Doc: $x$ being a vector or matrix, returns a row vector
   with two components, the first being the number of rows (1 for a row vector),
   the second the number of columns (1 for a column vector).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matsize(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sizedigit(*argv):
  '''
  sizedigit
  Class: basic
  Section: conversions
  C-Name: sizedigit
  Prototype: lG
  Help: sizedigit(x): maximum number of decimal digits minus one of (the
   coefficients of) x.
  Doc: 
   outputs a quick bound for the number of decimal
   digits of (the components of) $x$, off by at most $1$. If you want the
   exact value, you can use \kbd{\#Str(x)}, which is slower.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.sizedigit(*c_arg_tuple)

def vecmax(*argv):
  '''
  vecmax
  Class: basic
  Section: operators
  C-Name: vecmax0
  Prototype: GD&
  Help: vecmax(x,{&v}): largest entry in the vector/matrix x. If v
   is present, set it to the index of a largest entry (indirect max).
  Description: 
    (gen):gen            vecmax($1)
    (gen, &gen):gen      vecmax0($1, &$2)
  Doc: if $x$ is a vector or a matrix, returns the largest entry of $x$,
   otherwise returns a copy of $x$. Error if $x$ is empty.
   
   If $v$ is given, set it to the index of a largest entry (indirect maximum),
   when $x$ is a vector. If $x$ is a matrix, set $v$ to coordinates $[i,j]$
   such that $x[i,j]$ is a largest entry. This flag is ignored if $x$ is not a
   vector or matrix.
   
   \bprog
   ? vecmax([10, 20, -30, 40])
   %1 = 40
   ? vecmax([10, 20, -30, 40], &v); v
   %2 = 4
   ? vecmax([10, 20; -30, 40], &v); v
   %3 = [2, 2]
   @eprog
  Variant: Also available is \fun{GEN}{vecmax}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.vecmax0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def vecmin(*argv):
  '''
  vecmin
  Class: basic
  Section: operators
  C-Name: vecmin0
  Prototype: GD&
  Help: vecmin(x,{&v}): smallest entry in the vector/matrix x. If v is
   present, set it to the index of a smallest
   entry (indirect min).
  Description: 
    (gen):gen            vecmin($1)
    (gen, &gen):gen      vecmin0($1, &$2)
  Doc: if $x$ is a vector or a matrix, returns the smallest entry of $x$,
   otherwise returns a copy of $x$. Error if $x$ is empty.
   
   If $v$ is given, set it to the index of a smallest entry (indirect minimum),
   when $x$ is a vector. If $x$ is a matrix, set $v$ to coordinates $[i,j]$ such
   that $x[i,j]$ is a smallest entry. This is ignored if $x$ is not a vector or
   matrix.
   
   \bprog
   ? vecmin([10, 20, -30, 40])
   %1 = -30
   ? vecmin([10, 20, -30, 40], &v); v
   %2 = 3
   ? vecmin([10, 20; -30, 40], &v); v
   %3 = [2, 1]
   @eprog
  Variant: Also available is \fun{GEN}{vecmin}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.vecmin0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def centerlift(*argv):
  '''
  centerlift
  Class: basic
  Section: conversions
  C-Name: centerlift0
  Prototype: GDn
  Help: centerlift(x,{v}): centered lift of x. Same as lift except for
   intmod and padic components.
  Description: 
   (pol):pol        centerlift($1)
   (vec):vec        centerlift($1)
   (gen):gen        centerlift($1)
   (pol, var):pol        centerlift0($1, $2)
   (vec, var):vec        centerlift0($1, $2)
   (gen, var):gen        centerlift0($1, $2)
  Doc: Same as \tet{lift}, except that \typ{INTMOD} and \typ{PADIC} components
   are lifted using centered residues:
   
   \item for a \typ{INTMOD} $x\in \Z/n\Z$, the lift $y$ is such that
   $-n/2<y\le n/2$.
   
   \item  a \typ{PADIC} $x$ is lifted in the same way as above (modulo
   $p^\kbd{padicprec(x)}$) if its valuation $v$ is non-negative; if not, returns
   the fraction $p^v$ \kbd{centerlift}$(x p^{-v})$; in particular, rational
   reconstruction is not attempted. Use \tet{bestappr} for this.
   
   For backward compatibility, \kbd{centerlift(x,'v)} is allowed as an alias
   for \kbd{lift(x,'v)}.
   
   \synt{centerlift}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.centerlift0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def component(*argv):
  '''
  component
  Class: basic
  Section: conversions
  C-Name: compo
  Prototype: GL
  Help: component(x,n): the n'th component of the internal representation of
   x. For vectors or matrices, it is simpler to use x[]. For list objects such
   as nf, bnf, bnr or ell, it is much easier to use member functions starting
   with ".".
  Description: 
   (error,small):gen     err_get_compo($1, $2)
   (gen,small):gen       compo($1,$2)
  Doc: extracts the $n^{\text{th}}$-component of $x$. This is to be understood
   as follows: every PARI type has one or two initial \idx{code words}. The
   components are counted, starting at 1, after these code words. In particular
   if $x$ is a vector, this is indeed the $n^{\text{th}}$-component of $x$, if
   $x$ is a matrix, the $n^{\text{th}}$ column, if $x$ is a polynomial, the
   $n^{\text{th}}$ coefficient (i.e.~of degree $n-1$), and for power series,
   the $n^{\text{th}}$ significant coefficient.
   
   For polynomials and power series, one should rather use \tet{polcoeff}, and
   for vectors and matrices, the \kbd{[$\,$]} operator. Namely, if $x$ is a
   vector, then \tet{x[n]} represents the $n^{\text{th}}$ component of $x$. If
   $x$ is a matrix, \tet{x[m,n]} represents the coefficient of row \kbd{m} and
   column \kbd{n} of the matrix, \tet{x[m,]} represents the $m^{\text{th}}$
   \emph{row} of $x$, and \tet{x[,n]} represents the $n^{\text{th}}$
   \emph{column} of $x$.
   
   Using of this function requires detailed knowledge of the structure of the
   different PARI types, and thus it should almost never be used directly.
   Some useful exceptions:
   \bprog
       ? x = 3 + O(3^5);
       ? component(x, 2)
       %2 = 81   \\ p^(p-adic accuracy)
       ? component(x, 1)
       %3 = 3    \\ p
       ? q = Qfb(1,2,3);
       ? component(q, 1)
       %5 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.compo(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def denominator(*argv):
  '''
  denominator
  Class: basic
  Section: conversions
  C-Name: denom
  Prototype: G
  Help: denominator(x): denominator of x (or lowest common denominator in case
   of an array).
  Doc: 
   denominator of $x$. The meaning of this
   is clear when $x$ is a rational number or function. If $x$ is an integer
   or a polynomial, it is treated as a rational number or function,
   respectively, and the result is equal to $1$. For polynomials, you
   probably want to use
   \bprog
   denominator( content(x) )
   @eprog\noindent
   instead. As for modular objects, \typ{INTMOD} and \typ{PADIC} have
   denominator $1$, and the denominator of a \typ{POLMOD} is the denominator
   of its (minimal degree) polynomial representative.
   
   If $x$ is a recursive structure, for instance a vector or matrix, the lcm
   of the denominators of its components (a common denominator) is computed.
   This also applies for \typ{COMPLEX}s and \typ{QUAD}s.
   
   \misctitle{Warning} Multivariate objects are created according to variable
   priorities, with possibly surprising side effects ($x/y$ is a polynomial, but
   $y/x$ is a rational function). See \secref{se:priority}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.denom(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def deriv(*argv):
  '''
  deriv
  Class: basic
  Section: polynomials
  C-Name: deriv
  Prototype: GDn
  Help: deriv(x,{v}): derivative of x with respect to v, or to the main
   variable of x if v is omitted.
  Doc: 
   derivative of $x$ with respect to the main
   variable if $v$ is omitted, and with respect to $v$ otherwise. The derivative
   of a scalar type is zero, and the derivative of a vector or matrix is done
   componentwise. One can use $x'$ as a shortcut if the derivative is with
   respect to the main variable of $x$.
   
   By definition, the main variable of a \typ{POLMOD} is the main variable among
   the coefficients from its two polynomial components (representative and
   modulus); in other words, assuming a polmod represents an element of
   $R[X]/(T(X))$, the variable $X$ is a mute variable and the derivative is
   taken with respect to the main variable used in the base ring $R$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.deriv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def diffop(*argv):
  '''
  diffop
  Class: basic
  Section: polynomials
  C-Name: diffop0
  Prototype: GGGD1,L,
  Help: diffop(x,v,d,{n=1}): apply the differential operator D to x, where D is defined
   by D(v[i])=d[i], where v is a vector of variable names. D is 0 for variables
   outside of v unless they appear as modulus of a POLMOD. If the optional parameter n
   is given, return D^n(x) instead.
  Description: 
   (gen,gen,gen,?1):gen    diffop($1, $2, $3)
   (gen,gen,gen,small):gen diffop0($1, $2, $3, $4)
  Doc: 
   Let $v$ be a vector of variables, and $d$ a vector of the same length,
   return the image of $x$ by the $n$-power ($1$ if n is not given) of the differential
   operator $D$ that assumes the value \kbd{d[i]} on the variable \kbd{v[i]}.
   The value of $D$ on a scalar type is zero, and $D$ applies componentwise to a vector
   or matrix. When applied to a \typ{POLMOD}, if no value is provided for the variable
   of the modulus, such value is derived using the implicit function theorem.
   
   Some examples:
   This function can be used to differentiate formal expressions:
   If $E=\exp(X^2)$ then we have $E'=2*X*E$. We can derivate $X*exp(X^2)$ as follow:
   \bprog
   ? diffop(E*X,[X,E],[1,2*X*E])
   %1 = (2*X^2 + 1)*E
   @eprog
   Let \kbd{Sin} and \kbd{Cos} be two function such that $\kbd{Sin}^2+\kbd{Cos}^2=1$
   and $\kbd{Cos}'=-\kbd{Sin}$. We can differentiate $\kbd{Sin}/\kbd{Cos}$ as follow,
   PARI inferring the value of $\kbd{Sin}'$ from the equation:
   \bprog
   ? diffop(Mod('Sin/'Cos,'Sin^2+'Cos^2-1),['Cos],[-'Sin])
   %1 = Mod(1/Cos^2, Sin^2 + (Cos^2 - 1))
   
   @eprog
   Compute the Bell polynomials (both complete and partial) via the Faa di Bruno formula:
   \bprog
   Bell(k,n=-1)=
   {
     my(var(i)=eval(Str("X",i)));
     my(x,v,dv);
     v=vector(k,i,if(i==1,'E,var(i-1)));
     dv=vector(k,i,if(i==1,'X*var(1)*'E,var(i)));
     x=diffop('E,v,dv,k)/'E;
     if(n<0,subst(x,'X,1),polcoeff(x,n,'X))
   }
   @eprog
  Variant: 
   For $n=1$, the function \fun{GEN}{diffop}{GEN x, GEN v, GEN d} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.diffop0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def divrem(*argv):
  '''
  divrem
  Class: basic
  Section: operators
  C-Name: divrem
  Prototype: GGDn
  Help: divrem(x,y,{v}): euclidean division of x by y giving as a
   2-dimensional column vector the quotient and the remainder, with respect to
   v (to main variable if v is omitted)
  Doc: creates a column vector with two components, the first being the Euclidean
   quotient (\kbd{$x$ \bs\ $y$}), the second the Euclidean remainder
   (\kbd{$x$ - ($x$\bs$y$)*$y$}), of the division of $x$ by $y$. This avoids the
   need to do two divisions if one needs both the quotient and the remainder.
   If $v$ is present, and $x$, $y$ are multivariate
   polynomials, divide with respect to the variable $v$.
   
   Beware that \kbd{divrem($x$,$y$)[2]} is in general not the same as
   \kbd{$x$ \% $y$}; no GP operator corresponds to it:
   \bprog
   ? divrem(1/2, 3)[2]
   %1 = 1/2
   ? (1/2) % 3
   %2 = 2
   ? divrem(Mod(2,9), 3)[2]
    ***   at top-level: divrem(Mod(2,9),3)[2
    ***                 ^--------------------
    ***   forbidden division t_INTMOD \ t_INT.
   ? Mod(2,9) % 6
   %3 = Mod(2,3)
   @eprog
  Variant: Also available is \fun{GEN}{gdiventres}{GEN x, GEN y} when $v$ is
   not needed.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.divrem(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ceil(*argv):
  '''
  ceil
  Class: basic
  Section: conversions
  C-Name: gceil
  Prototype: G
  Help: ceil(x): ceiling of x = smallest integer >= x.
  Description: 
   (small):small:parens   $1
   (int):int:copy:parens  $1
   (real):int             ceilr($1)
   (mp):int               mpceil($1)
   (gen):gen              gceil($1)
  Doc: 
   ceiling of $x$. When $x$ is in $\R$, the result is the
   smallest integer greater than or equal to $x$. Applied to a rational
   function, $\kbd{ceil}(x)$ returns the Euclidean quotient of the numerator by
   the denominator.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gceil(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def floor(*argv):
  '''
  floor
  Class: basic
  Section: conversions
  C-Name: gfloor
  Prototype: G
  Help: floor(x): floor of x = largest integer <= x.
  Description: 
   (small):small:parens   $1
   (int):int:copy:parens  $1
   (real):int             floorr($1)
   (mp):int               mpfloor($1)
   (gen):gen              gfloor($1)
  Doc: 
   floor of $x$. When $x$ is in $\R$, the result is the
   largest integer smaller than or equal to $x$. Applied to a rational function,
   $\kbd{floor}(x)$ returns the Euclidean quotient of the numerator by the
   denominator.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gfloor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def frac(*argv):
  '''
  frac
  Class: basic
  Section: conversions
  C-Name: gfrac
  Prototype: G
  Help: frac(x): fractional part of x = x-floor(x).
  Doc: 
   fractional part of $x$. Identical to
   $x-\text{floor}(x)$. If $x$ is real, the result is in $[0,1[$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gfrac(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def O(*argv):
  '''
  O
  Class: basic
  Section: polynomials
  C-Name: ggrando
  Prototype: 
  Help: O(p^e): p-adic or power series zero with precision given by e
  Doc: if $p$ is an integer
   greater than $2$, returns a $p$-adic $0$ of precision $e$. In all other
   cases, returns a power series zero with precision given by $e v$, where $v$
   is the $X$-adic valuation of $p$ with respect to its main variable.
  Variant: \fun{GEN}{zeropadic}{GEN p, long e} for a $p$-adic and
   \fun{GEN}{zeroser}{long v, long e} for a power series zero in variable $v$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ggrando(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def imag(*argv):
  '''
  imag
  Class: basic
  Section: conversions
  C-Name: gimag
  Prototype: G
  Help: imag(x): imaginary part of x.
  Doc: imaginary part of $x$. When $x$ is a quadratic number, this is the
   coefficient of $\omega$ in the ``canonical'' integral basis $(1,\omega)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gimag(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Mod(*argv):
  '''
  Mod
  Class: basic
  Section: conversions
  C-Name: gmodulo
  Prototype: GG
  Help: Mod(a,b): creates 'a modulo b'.
  Description: 
   (small, small):gen         gmodulss($1, $2)
   (small, gen):gen           gmodulsg($1, $2)
   (gen, gen):gen             gmodulo($1, $2)
  Doc: in its basic form, creates an intmod or a polmod $(a \mod b)$; $b$ must
   be an integer or a polynomial. We then obtain a \typ{INTMOD} and a
   \typ{POLMOD} respectively:
   \bprog
   ? t = Mod(2,17); t^8
   %1 = Mod(1, 17)
   ? t = Mod(x,x^2+1); t^2
   %2 = Mod(-1, x^2+1)
   @eprog\noindent If $a \% b$ makes sense and yields a result of the
   appropriate type (\typ{INT} or scalar/\typ{POL}), the operation succeeds as
   well:
   \bprog
   ? Mod(1/2, 5)
   %3 = Mod(3, 5)
   ? Mod(7 + O(3^6), 3)
   %4 = Mod(1, 3)
   ? Mod(Mod(1,12), 9)
   %5 = Mod(1, 3)
   ? Mod(1/x, x^2+1)
   %6 = Mod(-1, x^2+1)
   ? Mod(exp(x), x^4)
   %7 = Mod(1/6*x^3 + 1/2*x^2 + x + 1, x^4)
   @eprog
   If $a$ is a complex object, ``base change'' it to $\Z/b\Z$ or $K[x]/(b)$,
   which is equivalent to, but faster than, multiplying it by \kbd{Mod(1,b)}:
   \bprog
   ? Mod([1,2;3,4], 2)
   %8 =
   [Mod(1, 2) Mod(0, 2)]
   
   [Mod(1, 2) Mod(0, 2)]
   ? Mod(3*x+5, 2)
   %9 = Mod(1, 2)*x + Mod(1, 2)
   ? Mod(x^2 + y*x + y^3, y^2+1)
   %10 = Mod(1, y^2 + 1)*x^2 + Mod(y, y^2 + 1)*x + Mod(-y, y^2 + 1)
   @eprog
   
   This function is not the same as $x$ \kbd{\%} $y$, the result of which
   has no knowledge of the indended modulus $y$. Compare
   \bprog
   ? x = 4 % 5; x + 1
   %1 = 5
   ? x = Mod(4,5); x + 1
   %2 = Mod(0,5)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gmodulo(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def variable(*argv):
  '''
  variable
  Class: basic
  Section: conversions
  C-Name: gpolvar
  Prototype: DG
  Help: variable({x}): main variable of object x. Gives p for p-adic x, 0
   if no variable can be associated to x. Returns the list of user variables if
   x is omitted.
  Description: 
   (pol):var:parens:copy        $var:1
   (gen):gen        gpolvar($1)
  Doc: 
   gives the main variable of the object $x$ (the variable with the highest
   priority used in $x$), and $p$ if $x$ is a $p$-adic number. Return $0$ if
   $x$ has no variable associated to it.
   \bprog
   ? variable(x^2 + y)
   %1 = x
   ? variable(1 + O(5^2))
   %2 = 5
   ? variable([x,y,z,t])
   %3 = x
   ? variable(1)
   %4 = 0
   @eprog\noindent The construction
   \bprog
      if (!variable(x),...)
   @eprog\noindent can be used to test whether a variable is attached to $x$.
   
   If $x$ is omitted, returns the list of user variables known to the
   interpreter, by order of decreasing priority. (Highest priority is $x$,
   which always come first.)
  Variant: However, in library mode, this function should not be used for $x$
   non-\kbd{NULL}, since \tet{gvar} is more appropriate. Instead, for
   $x$ a $p$-adic (type \typ{PADIC}), $p$ is $gel(x,2)$; otherwise, use
   \fun{long}{gvar}{GEN x} which returns the variable number of $x$ if
   it exists, \kbd{NO\_VARIABLE} otherwise, which satisfies the property
   $\kbd{varncmp}(\kbd{NO\_VARIABLE}, v) > 0$ for all valid variable number
   $v$, i.e. it has lower priority than any variable.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gpolvar(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def real(*argv):
  '''
  real
  Class: basic
  Section: conversions
  C-Name: greal
  Prototype: G
  Help: real(x): real part of x.
  Doc: real part of $x$. In the case where $x$ is a quadratic number, this is the
   coefficient of $1$ in the ``canonical'' integral basis $(1,\omega)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.greal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def shift(*argv):
  '''
  shift
  Class: basic
  Section: operators
  C-Name: gshift
  Prototype: GL
  Help: shift(x,n): shift x left n bits if n>=0, right -n bits if
   n<0.
  Doc: shifts $x$ componentwise left by $n$ bits if $n\ge0$ and right by $|n|$
   bits if $n<0$. May be abbreviated as $x$ \kbd{<<} $n$ or $x$ \kbd{>>} $(-n)$.
   A left shift by $n$ corresponds to multiplication by $2^n$. A right shift of an
   integer $x$ by $|n|$ corresponds to a Euclidean division of $x$ by $2^{|n|}$
   with a remainder of the same sign as $x$, hence is not the same (in general) as
   $x \kbd{\bs} 2^n$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gshift(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def subst(*argv):
  '''
  subst
  Class: basic
  Section: polynomials
  C-Name: gsubst
  Prototype: GnG
  Help: subst(x,y,z): in expression x, replace the variable y by the
   expression z.
  Doc: replace the simple variable $y$ by the argument $z$ in the ``polynomial''
   expression $x$. Every type is allowed for $x$, but if it is not a genuine
   polynomial (or power series, or rational function), the substitution will be
   done as if the scalar components were polynomials of degree zero. In
   particular, beware that:
   
   \bprog
   ? subst(1, x, [1,2; 3,4])
   %1 =
   [1 0]
   
   [0 1]
   
   ? subst(1, x, Mat([0,1]))
     ***   at top-level: subst(1,x,Mat([0,1])
     ***                 ^--------------------
     *** subst: forbidden substitution by a non square matrix.
   @eprog\noindent
   If $x$ is a power series, $z$ must be either a polynomial, a power
   series, or a rational function. Finally, if $x$ is a vector,
   matrix or list, the substitution is applied to each individual entry.
   
   Use the function \kbd{substvec} to replace several variables at once,
   or the function \kbd{substpol} to replace a polynomial expression.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsubst(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def substpol(*argv):
  '''
  substpol
  Class: basic
  Section: polynomials
  C-Name: gsubstpol
  Prototype: GGG
  Help: substpol(x,y,z): in expression x, replace the polynomial y by the
   expression z, using remainder decomposition of x.
  Doc: replace the ``variable'' $y$ by the argument $z$ in the ``polynomial''
   expression $x$. Every type is allowed for $x$, but the same behavior
   as \kbd{subst} above apply.
   
   The difference with \kbd{subst} is that $y$ is allowed to be any polynomial
   here. The substitution is done moding out all components of $x$
   (recursively) by $y - t$, where $t$ is a new free variable of lowest
   priority. Then substituting $t$ by $z$ in the resulting expression. For
   instance
   \bprog
   ? substpol(x^4 + x^2 + 1, x^2, y)
   %1 = y^2 + y + 1
   ? substpol(x^4 + x^2 + 1, x^3, y)
   %2 = x^2 + y*x + 1
   ? substpol(x^4 + x^2 + 1, (x+1)^2, y)
   %3 = (-4*y - 6)*x + (y^2 + 3*y - 3)
   @eprog
  Variant: Further, \fun{GEN}{gdeflate}{GEN T, long v, long d} attempts to
   write $T(x)$ in the form $t(x^d)$, where $x=$\kbd{pol\_x}$(v)$, and returns
   \kbd{NULL} if the substitution fails (for instance in the example \kbd{\%2}
   above).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsubstpol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def substvec(*argv):
  '''
  substvec
  Class: basic
  Section: polynomials
  C-Name: gsubstvec
  Prototype: GGG
  Help: substvec(x,v,w): in expression x, make a best effort to replace the
   variables v1,...,vn by the expression w1,...,wn.
  Doc: $v$ being a vector of monomials of degree 1 (variables),
   $w$ a vector of expressions of the same length, replace in the expression
   $x$ all occurrences of $v_i$ by $w_i$. The substitutions are done
   simultaneously; more precisely, the $v_i$ are first replaced by new
   variables in $x$, then these are replaced by the $w_i$:
   
   \bprog
   ? substvec([x,y], [x,y], [y,x])
   %1 = [y, x]
   ? substvec([x,y], [x,y], [y,x+y])
   %2 = [y, x + y]     \\ not [y, 2*y]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsubstvec(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Col(*argv):
  '''
  Col
  Class: basic
  Section: conversions
  C-Name: gtocol0
  Prototype: GD0,L,
  Help: Col(x, {n}): transforms the object x into a column vector of dimension n.
  Description: 
   (gen):vec     gtocol($1)
  Doc: 
   transforms the object $x$ into a column vector. The dimension of the
   resulting vector can be optionally specified via the extra parameter $n$.
   
   If $n$ is omitted or $0$, the dimension depends on the type of $x$; the
   vector has a single component, except when $x$ is
   
   \item a vector or a quadratic form (in which case the resulting vector
   is simply the initial object considered as a row vector),
   
   \item a polynomial or a power series. In the case of a polynomial, the
   coefficients of the vector start with the leading coefficient of the
   polynomial, while for power series only the significant coefficients are
   taken into account, but this time by increasing order of degree.
   In this last case, \kbd{Vec} is the reciprocal function of \kbd{Pol} and
   \kbd{Ser} respectively,
   
   \item a matrix (the column of row vector comprising the matrix is returned),
   
   \item a character string (a vector of individual characters is returned).
   
   In the last two cases (matrix and character string), $n$ is meaningless and
   must be omitted or an error is raised. Otherwise, if $n$ is given, $0$
   entries are appended at the end of the vector if $n > 0$, and prepended at
   the beginning if $n < 0$. The dimension of the resulting vector is $|n|$.
   
   Note that the function \kbd{Colrev} does not exist, use \kbd{Vecrev}.
  Variant: \fun{GEN}{gtocol}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtocol0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Colrev(*argv):
  '''
  Colrev
  Class: basic
  Section: conversions
  C-Name: gtocolrev0
  Prototype: GD0,L,
  Help: Colrev(x, {n}): transforms the object x into a column vector of
   dimension n in reverse order with respect to Col(x, {n}). Empty vector if x
   is omitted.
  Description: 
   (gen):vec     gtocolrev($1)
  Doc: 
   as $\kbd{Col}(x, n)$, then reverse the result. In particular
  Variant: \fun{GEN}{gtocolrev}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtocolrev0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Pol(*argv):
  '''
  Pol
  Class: basic
  Section: conversions
  C-Name: gtopoly
  Prototype: GDn
  Help: Pol(t,{v='x}): convert t (usually a vector or a power series) into a
   polynomial with variable v, starting with the leading coefficient.
  Description: 
   (gen,?var):pol  gtopoly($1, $2)
  Doc: 
   transforms the object $t$ into a polynomial with main variable $v$. If $t$
   is a scalar, this gives a constant polynomial. If $t$ is a power series with
   non-negative valuation or a rational function, the effect is similar to
   \kbd{truncate}, i.e.~we chop off the $O(X^k)$ or compute the Euclidean
   quotient of the numerator by the denominator, then change the main variable
   of the result to $v$.
   
   The main use of this function is when $t$ is a vector: it creates the
   polynomial whose coefficients are given by $t$, with $t[1]$ being the leading
   coefficient (which can be zero). It is much faster to evaluate
   \kbd{Pol} on a vector of coefficients in this way, than the corresponding
   formal expression $a_n X^n + \dots + a_0$, which is evaluated naively exactly
   as written (linear versus quadratic time in $n$). \tet{Polrev} can be used if
   one wants $x[1]$ to be the constant coefficient:
   \bprog
   ? Pol([1,2,3])
   %1 = x^2 + 2*x + 3
   ? Polrev([1,2,3])
   %2 = 3*x^2 + 2*x + 1
   @eprog\noindent
   The reciprocal function of \kbd{Pol} (resp.~\kbd{Polrev}) is \kbd{Vec} (resp.~
   \kbd{Vecrev}).
   \bprog
   ? Vec(Pol([1,2,3]))
   %1 = [1, 2, 3]
   ? Vecrev( Polrev([1,2,3]) )
   %2 = [1, 2, 3]
   @eprog\noindent
   
   \misctitle{Warning} This is \emph{not} a substitution function. It will not
   transform an object containing variables of higher priority than~$v$.
   \bprog
   ? Pol(x + y, y)
     ***   at top-level: Pol(x+y,y)
     ***                 ^----------
     *** Pol: variable must have higher priority in gtopoly.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtopoly(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Polrev(*argv):
  '''
  Polrev
  Class: basic
  Section: conversions
  C-Name: gtopolyrev
  Prototype: GDn
  Help: Polrev(t,{v='x}): convert t (usually a vector or a power series) into a
   polynomial with variable v, starting with the constant term.
  Description: 
   (gen,?var):pol  gtopolyrev($1, $2)
  Doc: 
   transform the object $t$ into a polynomial
   with main variable $v$. If $t$ is a scalar, this gives a constant polynomial.
   If $t$ is a power series, the effect is identical to \kbd{truncate}, i.e.~it
   chops off the $O(X^k)$.
   
   The main use of this function is when $t$ is a vector: it creates the
   polynomial whose coefficients are given by $t$, with $t[1]$ being the
   constant term. \tet{Pol} can be used if one wants $t[1]$ to be the leading
   coefficient:
   \bprog
   ? Polrev([1,2,3])
   %1 = 3*x^2 + 2*x + 1
   ? Pol([1,2,3])
   %2 = x^2 + 2*x + 3
   @eprog
   The reciprocal function of \kbd{Pol} (resp.~\kbd{Polrev}) is \kbd{Vec} (resp.~
   \kbd{Vecrev}).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtopolyrev(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Ser(*argv):
  '''
  Ser
  Class: basic
  Section: conversions
  C-Name: gtoser
  Prototype: GDnDP
  Help: Ser(s,{v='x},{d=seriesprecision}): convert s into a power series with
   variable v and precision d, starting with the constant coefficient.
  Doc: transforms the object $s$ into a power series with main variable $v$
   ($x$ by default) and precision (number of significant terms) equal to
   $d$ (= the default \kbd{seriesprecision} by default). If $s$ is a
   scalar, this gives a constant power series in $v$ with precision \kbd{d}.
   If $s$ is a polynomial, the polynomial is truncated to $d$ terms if needed
   \bprog
   ? Ser(1, 'y, 5)
   %1 = 1 + O(y^5)
   ? Ser(x^2,, 5)
   %2 = x^2 + O(x^7)
   ? T = polcyclo(100)
   %3 = x^40 - x^30 + x^20 - x^10 + 1
   ? Ser(T, 'x, 11)
   %4 = 1 - x^10 + O(x^11)
   @eprog\noindent The function is more or less equivalent with multiplication by
   $1 + O(v^d)$ in theses cases, only faster.
   
   If $s$ is a vector, on the other hand, the coefficients of the vector are
   understood to be the coefficients of the power series starting from the
   constant term (as in \tet{Polrev}$(x)$), and the precision $d$ is ignored:
   in other words, in this case, we convert \typ{VEC} / \typ{COL} to the power
   series whose significant terms are exactly given by the vector entries.
   Finally, if $s$ is already a power series in $v$, we return it verbatim,
   ignoring $d$ again. If $d$ significant terms are desired in the last two
   cases, convert/truncate to \typ{POL} first.
   \bprog
   ? v = [1,2,3]; Ser(v, t, 7)
   %5 = 1 + 2*t + 3*t^2 + O(t^3)  \\ 3 terms: 7 is ignored!
   ? Ser(Polrev(v,t), t, 7)
   %6 = 1 + 2*t + 3*t^2 + O(t^7)
   ? s = 1+x+O(x^2); Ser(s, x, 7)
   %7 = 1 + x + O(x^2)  \\ 2 terms: 7 ignored
   ? Ser(truncate(s), x, 7)
   %8 = 1 + x + O(x^7)
   @eprog\noindent
   The warning given for \kbd{Pol} also applies here: this is not a substitution
   function.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtoser(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Vec(*argv):
  '''
  Vec
  Class: basic
  Section: conversions
  C-Name: gtovec0
  Prototype: GD0,L,
  Help: Vec(x, {n}): transforms the object x into a vector of dimension n.
  Description: 
   (gen):vec        gtovec($1)
  Doc: 
   transforms the object $x$ into a row vector. The dimension of the
   resulting vector can be optionally specified via the extra parameter $n$.
   
   If $n$ is omitted or $0$, the dimension depends on the type of $x$; the
   vector has a single component, except when $x$ is
   
   \item a vector or a quadratic form (in which case the resulting vector
   is simply the initial object considered as a row vector),
   
   \item a polynomial or a power series. In the case of a polynomial, the
   coefficients of the vector start with the leading coefficient of the
   polynomial, while for power series only the significant coefficients are
   taken into account, but this time by increasing order of degree.
   In this last case, \kbd{Vec} is the reciprocal function of \kbd{Pol} and
   \kbd{Ser} respectively,
   
   \item a matrix: return the vector of columns comprising the matrix.
   
   \item a character string: return the vector of individual characters.
   
   \item an error context (\typ{ERROR}): return the error components, see
   \tet{iferr}.
   
   In the last three cases (matrix, character string, error), $n$ is
   meaningless and must be omitted or an error is raised. Otherwise, if $n$ is
   given, $0$ entries are appended at the end of the vector if $n > 0$, and
   prepended at the beginning if $n < 0$. The dimension of the resulting vector
   is $|n|$. Variant: \fun{GEN}{gtovec}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtovec0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Vecrev(*argv):
  '''
  Vecrev
  Class: basic
  Section: conversions
  C-Name: gtovecrev0
  Prototype: GD0,L,
  Help: Vecrev(x, {n}): transforms the object x into a vector of dimension n
   in reverse order with respect to Vec(x, {n}). Empty vector if x is omitted.
  Description: 
   (gen):vec     gtovecrev($1)
  Doc: 
   as $\kbd{Vec}(x, n)$, then reverse the result. In particular
   In this case, \kbd{Vecrev} is the reciprocal function of \kbd{Polrev}: the
   coefficients of the vector start with the constant coefficient of the
   polynomial and the others follow by increasing degree.
  Variant: \fun{GEN}{gtovecrev}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtovecrev0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Vecsmall(*argv):
  '''
  Vecsmall
  Class: basic
  Section: conversions
  C-Name: gtovecsmall0
  Prototype: GD0,L,
  Help: Vecsmall(x, {n}): transforms the object x into a VECSMALL of dimension n.
  Description: 
   (gen):vecsmall                gtovecsmall($1)
  Doc: 
   transforms the object $x$ into a row vector of type \typ{VECSMALL}. The
   dimension of the resulting vector can be optionally specified via the extra
   parameter $n$.
   
   This acts as \kbd{Vec}$(x,n)$, but only on a limited set of objects:
   the result must be representable as a vector of small integers.
   If $x$ is a character string, a vector of individual characters in ASCII
   encoding is returned (\tet{Strchr} yields back the character string).
  Variant: \fun{GEN}{gtovecsmall}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtovecsmall0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def intformal(*argv):
  '''
  intformal
  Class: basic
  Section: polynomials
  C-Name: integ
  Prototype: GDn
  Help: intformal(x,{v}): formal integration of x with respect to v, or to the
   main variable of x if v is omitted.
  Doc: \idx{formal integration} of $x$ with respect to the variable $v$ (wrt.
   the main variable if $v$ is omitted). Since PARI cannot represent
   logarithmic or arctangent terms, any such term in the result will yield an
   error:
   \bprog
    ? intformal(x^2)
    %1 = 1/3*x^3
    ? intformal(x^2, y)
    %2 = y*x^2
    ? intformal(1/x)
      ***   at top-level: intformal(1/x)
      ***                 ^--------------
      *** intformal: domain error in intformal: residue(series, pole) != 0
   @eprog
   The argument $x$ can be of any type. When $x$ is a rational function, we
   assume that the base ring is an integral domain of characteristic zero.
   
     By  definition,   the main variable of a \typ{POLMOD} is the main variable
   among the  coefficients  from  its  two  polynomial  components
   (representative and modulus); in other words, assuming a polmod represents an
   element of $R[X]/(T(X))$, the variable $X$ is a mute variable and the
   integral is taken with respect to the main variable used in the base ring $R$.
   In particular, it is meaningless to integrate with respect to the main
   variable of \kbd{x.mod}:
   \bprog
   ? intformal(Mod(1,x^2+1), 'x)
   *** intformal: incorrect priority in intformal: variable x = x
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.integ(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def lift(*argv):
  '''
  lift
  Class: basic
  Section: conversions
  C-Name: lift0
  Prototype: GDn
  Help: lift(x,{v}):
   if v is omitted, lifts elements of Z/nZ to Z, of Qp to Q, and of K[x]/(P) to
   K[x]. Otherwise lift only polmods with main variable v.
  Description: 
   (pol):pol        lift($1)
   (vec):vec        lift($1)
   (gen):gen        lift($1)
   (pol, var):pol        lift0($1, $2)
   (vec, var):vec        lift0($1, $2)
   (gen, var):gen        lift0($1, $2)
  Doc: 
   if $v$ is omitted, lifts intmods from $\Z/n\Z$ in $\Z$,
   $p$-adics from $\Q_p$ to $\Q$ (as \tet{truncate}), and polmods to
   polynomials. Otherwise, lifts only polmods whose modulus has main
   variable~$v$. \typ{FFELT} are not lifted, nor are List elements: you may
   convert the latter to vectors first, or use \kbd{apply(lift,L)}. More
   generally, components for which such lifts are meaningless (e.g. character
   strings) are copied verbatim.
   \bprog
   ? lift(Mod(5,3))
   %1 = 2
   ? lift(3 + O(3^9))
   %2 = 3
   ? lift(Mod(x,x^2+1))
   %3 = x
   ? lift(Mod(x,x^2+1))
   %4 = x
   @eprog
   Lifts are performed recursively on an object components, but only
   by \emph{one level}: once a \typ{POLMOD} is lifted, the components of
   the result are \emph{not} lifted further.
   \bprog
   ? lift(x * Mod(1,3) + Mod(2,3))
   %4 = x + 2
   ? lift(x * Mod(y,y^2+1) + Mod(2,3))
   %5 = y*x + Mod(2, 3)   \\@com do you understand this one?
   ? lift(x * Mod(y,y^2+1) + Mod(2,3), 'x)
   %6 = Mod(y, y^2 + 1)*x + Mod(Mod(2, 3), y^2 + 1)
   ? lift(%, y)
   %7 = y*x + Mod(2, 3)
   @eprog\noindent To recursively lift all components not only by one level,
   but as long as possible, use \kbd{liftall}. To lift only \typ{INTMOD}s and
   \typ{PADIC}s components, use \tet{liftint}. To lift only \typ{POLMOD}s
   components, use \tet{liftpol}. Finally, \tet{centerlift} allows to lift
   \typ{INTMOD}s and \typ{PADIC}s using centered residues (lift of smallest
   absolute value).
  Variant: Also available is \fun{GEN}{lift}{GEN x} corresponding to
   \kbd{lift0(x,-1)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.lift0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def liftall(*argv):
  '''
  liftall
  Class: basic
  Section: conversions
  C-Name: liftall
  Prototype: G
  Help: liftall(x): lifts every element of Z/nZ to Z, of Qp to Q, and of
   K[x]/(P) to K[x].
  Description: 
   (pol):pol        liftall($1)
   (vec):vec        liftall($1)
   (gen):gen        liftall($1)
  Doc: 
   recursively lift all components of $x$ from $\Z/n\Z$ to $\Z$,
   from $\Q_p$ to $\Q$ (as \tet{truncate}), and polmods to
   polynomials. \typ{FFELT} are not lifted, nor are List elements: you may
   convert the latter to vectors first, or use \kbd{apply(liftall,L)}. More
   generally, components for which such lifts are meaningless (e.g. character
   strings) are copied verbatim.
   \bprog
   ? liftall(x * (1 + O(3)) + Mod(2,3))
   %1 = x + 2
   ? liftall(x * Mod(y,y^2+1) + Mod(2,3)*Mod(z,z^2))
   %2 = y*x + 2*z
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.liftall(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def liftint(*argv):
  '''
  liftint
  Class: basic
  Section: conversions
  C-Name: liftint
  Prototype: G
  Help: liftint(x): lifts every element of Z/nZ to Z, of Qp to Q, and of
   K[x]/(P) to K[x].
  Description: 
   (pol):pol        liftint($1)
   (vec):vec        liftint($1)
   (gen):gen        liftint($1)
  Doc: recursively lift all components of $x$ from $\Z/n\Z$ to $\Z$ and
   from $\Q_p$ to $\Q$ (as \tet{truncate}).
   \typ{FFELT} are not lifted, nor are List elements: you may
   convert the latter to vectors first, or use \kbd{apply(liftint,L)}. More
   generally, components for which such lifts are meaningless (e.g. character
   strings) are copied verbatim.
   \bprog
   ? liftint(x * (1 + O(3)) + Mod(2,3))
   %1 = x + 2
   ? liftint(x * Mod(y,y^2+1) + Mod(2,3)*Mod(z,z^2))
   %2 = Mod(y, y^2 + 1)*x + Mod(Mod(2*z, z^2), y^2 + 1)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.liftint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def liftpol(*argv):
  '''
  liftpol
  Class: basic
  Section: conversions
  C-Name: liftpol
  Prototype: G
  Help: liftpol(x): lifts every polmod component of x to polynomials
  Description: 
   (pol):pol        liftpol($1)
   (vec):vec        liftpol($1)
   (gen):gen        liftpol($1)
  Doc: recursively lift all components of $x$ which are polmods to
   polynomials. \typ{FFELT} are not lifted, nor are List elements: you may
   convert the latter to vectors first, or use \kbd{apply(liftpol,L)}. More
   generally, components for which such lifts are meaningless (e.g. character
   strings) are copied verbatim.
   \bprog
   ? liftpol(x * (1 + O(3)) + Mod(2,3))
   %1 = (1 + O(3))*x + Mod(2, 3)
   ? liftpol(x * Mod(y,y^2+1) + Mod(2,3)*Mod(z,z^2))
   %2 = y*x + Mod(2, 3)*z
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.liftpol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def numerator(*argv):
  '''
  numerator
  Class: basic
  Section: conversions
  C-Name: numer
  Prototype: G
  Help: numerator(x): numerator of x.
  Doc: 
   numerator of $x$. The meaning of this
   is clear when $x$ is a rational number or function. If $x$ is an integer
   or a polynomial, it is treated as a rational number or function,
   respectively, and the result is $x$ itself. For polynomials, you
   probably want to use
   \bprog
   numerator( content(x) )
   @eprog\noindent
   instead.
   
   In other cases, \kbd{numerator(x)} is defined to be
   \kbd{denominator(x)*x}. This is the case when $x$ is a vector or a
   matrix, but also for \typ{COMPLEX} or \typ{QUAD}. In particular since a
   \typ{PADIC} or \typ{INTMOD} has  denominator $1$, its numerator is
   itself.
   
   \misctitle{Warning} Multivariate objects are created according to variable
   priorities, with possibly surprising side effects ($x/y$ is a polynomial, but
   $y/x$ is a rational function). See \secref{se:priority}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.numer(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def padicprec(*argv):
  '''
  padicprec
  Class: basic
  Section: conversions
  C-Name: padicprec
  Prototype: lGG
  Help: padicprec(x,p): absolute p-adic precision of object x.
  Doc: absolute $p$-adic precision of the object $x$. This is the minimum
   precision of the components of $x$. The result is \tet{LONG_MAX}
   ($2^{31}-1$ for 32-bit machines or $2^{63}-1$ for 64-bit machines) if $x$ is
   an exact object.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.padicprec(*c_arg_tuple)

def polcoeff(*argv):
  '''
  polcoeff
  Class: basic
  Section: polynomials
  C-Name: polcoeff0
  Prototype: GLDn
  Help: polcoeff(x,n,{v}): coefficient of degree n of x, or the n-th component
   for vectors or matrices (for which it is simpler to use x[]). With respect
   to the main variable if v is omitted, with respect to the variable v
   otherwise.
  Description: 
   (pol, 0):gen:copy      constant_term($1)
   (gen, small, ?var):gen polcoeff0($1, $2, $3)
  Doc: coefficient of degree $n$ of the polynomial $x$, with respect to the
   main variable if $v$ is omitted, with respect to $v$ otherwise.  If $n$
   is greater than the degree, the result is zero.
   
   Naturally applies to scalars (polynomial of degree $0$), as well as to
   rational functions whose denominator is a monomial.
   It also applies to power series: if $n$ is less than the valuation, the result
   is zero. If it is greater than the largest significant degree, then an error
   message is issued.
   
    For greater flexibility, $x$ can be a vector or matrix type and the
    function then returns \kbd{component(x,n)}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polcoeff0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poldegree(*argv):
  '''
  poldegree
  Class: basic
  Section: polynomials
  C-Name: poldegree
  Prototype: lGDn
  Help: poldegree(x,{v}): degree of the polynomial or rational function x with
   respect to main variable if v is omitted, with respect to v otherwise.
   For scalar x, return 0 is x is non-zero and a negative number otherwise.
  Description: 
   (pol):small                degpol($1)
   (gen):small                degree($1)
   (gen, var):small           poldegree($1, $2)
  Doc: degree of the polynomial $x$ in the main variable if $v$ is omitted, in
   the variable $v$ otherwise.
   
   The degree of $0$ is a fixed negative number, whose exact value should not
   be used. The degree of a non-zero scalar is $0$. Finally, when $x$ is a
   non-zero polynomial or rational function, returns the ordinary degree of
   $x$. Raise an error otherwise.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.poldegree(*c_arg_tuple)

def pollead(*argv):
  '''
  pollead
  Class: basic
  Section: polynomials
  C-Name: pollead
  Prototype: GDn
  Help: pollead(x,{v}): leading coefficient of polynomial or series x, or x
   itself if x is a scalar. Error otherwise. With respect to the main variable
   of x if v is omitted, with respect to the variable v otherwise.
  Description: 
   (pol):gen:copy         leading_term($1)
   (gen):gen              pollead($1, -1)
   (gen, var):gen         pollead($1, $2)
  Doc: leading coefficient of the polynomial or power series $x$. This is
    computed with respect to the main variable of $x$ if $v$ is omitted, with
    respect to the variable $v$ otherwise.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.pollead(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def precision(*argv):
  '''
  precision
  Class: basic
  Section: conversions
  C-Name: precision0
  Prototype: GD0,L,
  Help: precision(x,{n}): if n is present, return x at precision n. If n is omitted, return real precision of object x.
  Description: 
   (real):small          prec2ndec(gprecision($1))
   (gen):int             precision0($1, 0)
   (real,0):small        prec2ndec(gprecision($1))
   (gen,0):int           precision0($1, 0)
   (real,#small):real    rtor($1, ndec2prec($2))
   (gen,#small):gen      gprec($1, $2)
   (real,small):real     precision0($1, $2)
   (gen,small):gen       precision0($1, $2)
  Doc: the function has two different behaviors according to whether $n$ is present or not.
   
   If $n$ is missing, the function returns the precision in decimal digits of the
   PARI object $x$. If $x$ is
   an exact object, the largest single precision integer is returned.
   \bprog
   ? precision(exp(1e-100))
   %1 = 134                \\ 134 significant decimal digits
   ? precision(2 + x)
   %2 = 2147483647         \\ exact object
   ? precision(0.5 + O(x))
   %3 = 28                 \\ floating point accuracy, NOT series precision
   ? precision( [ exp(1e-100), 0.5 ] )
   %4 = 28                 \\ minimal accuracy among components
   @eprog\noindent
   The return value for exact objects is meaningless since it is not even the
   same on 32 and 64-bit machines. The proper way to test whether an object is
   exact is
   \bprog
   ? isexact(x) = precision(x) == precision(0)
   @eprog
   
   If $n$ is present, the function creates a new object equal to $x$ with a new
   ``precision'' $n$. (This never changes the type of the result. In particular
   it is not possible to use it to obtain a polynomial from a power series; for
   that, see \tet{truncate}.) Now the meaning of precision is different from the
   above (floating point accuracy), and depends on the type of $x$:
   
   For exact types, no change. For $x$ a vector or a matrix, the operation is
   done componentwise.
   
   For real $x$, $n$ is the number of desired significant \emph{decimal}
   digits. If $n$ is smaller than the precision of $x$, $x$ is truncated,
   otherwise $x$ is extended with zeros.
   
   For $x$ a $p$-adic or a power series, $n$ is the desired number of
   \emph{significant} $p$-adic or $X$-adic digits, where $X$ is the main
   variable of $x$. (Note: yes, this is inconsistent.)
   Note that the precision is a priori distinct from the exponent $k$ appearing
   in $O(*^k)$; it is indeed equal to $k$ if and only if $x$ is a $p$-adic
   or $X$-adic \emph{unit}.
   \bprog
   ? precision(1 + O(x), 10)
   %1 = 1 + O(x^10)
   ? precision(x^2 + O(x^10), 3)
   %2 = x^2 + O(x^5)
   ? precision(7^2 + O(7^10), 3)
   %3 = 7^2 + O(7^5)
   @eprog\noindent
   For the last two examples, note that $x^2 + O(x^5) = x^2(1 + O(x^3))$
   indeed has 3 significant coefficients
  Variant: Also available are \fun{GEN}{gprec}{GEN x, long n} and
   \fun{long}{precision}{GEN x}. In both, the accuracy is expressed in
   \emph{words} (32-bit or 64-bit depending on the architecture).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.precision0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfbil(*argv):
  '''
  qfbil
  Class: basic
  Section: linear_algebra
  C-Name: qfbil
  Prototype: GGDG
  Help: qfbil(x,y,{q}): evaluate the bilinear form q (symmetric matrix)
   at (x,y); if q omitted, use the standard Euclidean scalar product.
  Doc: evaluate the bilinear form $q$ (symmetric matrix)
   at the vectors $(x,y)$; if $q$ omitted, use the standard Euclidean scalar
   product, corresponding to the identity matrix.
   
   Roughly equivalent to \kbd{x\til * q * y}, but a little faster and
   more convenient (does not distinguish between column and row vectors):
   \bprog
   ? x = [1,2,3]~; y = [-1,0,1]~; qfbil(x,y)
   %1 = 2
   ? q = [1,2,3;2,2,-1;3,-1,0]; qfbil(x,y, q)
   %2 = -13
   ? for(i=1,10^6, qfbil(x,y,q))
   %3 = 568ms
   ? for(i=1,10^6, x~*q*y)
   %4 = 717ms
   @eprog\noindent The associated quadratic form is also available, as
   \tet{qfnorm}, slightly faster:
   \bprog
   ? for(i=1,10^6, qfnorm(x,q))
   time = 444ms
   ? for(i=1,10^6, qfnorm(x))
   time = 176 ms.
   ? for(i=1,10^6, qfbil(x,y))
   time = 208 ms.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfbil(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfnorm(*argv):
  '''
  qfnorm
  Class: basic
  Section: linear_algebra
  C-Name: qfnorm
  Prototype: GDG
  Help: qfnorm(x,{q}): evaluate the binary quadratic form q (symmetric matrix)
   at x; if q omitted, use the standard Euclidean form.
  Doc: evaluate the binary quadratic form $q$ (symmetric matrix)
   at the vector $x$. If $q$ omitted, use the standard Euclidean form,
   corresponding to the identity matrix.
   
   Equivalent to \kbd{x\til * q * x}, but about twice faster and
   more convenient (does not distinguish between column and row vectors):
   \bprog
   ? x = [1,2,3]~; qfnorm(x)
   %1 = 14
   ? q = [1,2,3;2,2,-1;3,-1,0]; qfnorm(x, q)
   %2 = 23
   ? for(i=1,10^6, qfnorm(x,q))
   time = 384ms.
   ? for(i=1,10^6, x~*q*x)
   time = 729ms.
   @eprog\noindent We also allow \typ{MAT}s of compatible dimensions for $x$,
   and return \kbd{x\til * q * x} in this case as well:
   \bprog
   ? M = [1,2,3;4,5,6;7,8,9]; qfnorm(M) \\ Gram matrix
   %5 =
   [66  78  90]
   
   [78  93 108]
   
   [90 108 126]
   
   ? for(i=1,10^6, qfnorm(M,q))
   time = 2,144 ms.
   ? for(i=1,10^6, M~*q*M)
   time = 2,793 ms.
   @eprog
   \noindent The polar form is also available, as \tet{qfbil}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def round(*argv):
  '''
  round
  Class: basic
  Section: conversions
  C-Name: round0
  Prototype: GD&
  Help: round(x,{&e}): take the nearest integer to all the coefficients of x.
   If e is present, do not take into account loss of integer part precision,
   and set e = error estimate in bits.
  Description: 
   (small):small:parens   $1
   (int):int:copy:parens  $1
   (real):int             roundr($1)
   (mp):int               mpround($1)
   (mp, &small):int       grndtoi($1, &$2)
   (mp, &int):int         round0($1, &$2)
   (gen):gen              ground($1)
   (gen, &small):gen      grndtoi($1, &$2)
   (gen, &int):gen        round0($1, &$2)
  Doc: If $x$ is in $\R$, rounds $x$ to the nearest integer (rounding to
   $+\infty$ in case of ties), then and sets $e$ to the number of error bits,
   that is the binary exponent of the difference between the original and the
   rounded value (the ``fractional part''). If the exponent of $x$ is too large
   compared to its precision (i.e.~$e>0$), the result is undefined and an error
   occurs if $e$ was not given.
   
   \misctitle{Important remark} Contrary to the other truncation functions,
   this function operates on every coefficient at every level of a PARI object.
   For example
   $$\text{truncate}\left(\dfrac{2.4*X^2-1.7}{X}\right)=2.4*X,$$
   whereas
   $$\text{round}\left(\dfrac{2.4*X^2-1.7}{X}\right)=\dfrac{2*X^2-2}{X}.$$
   An important use of \kbd{round} is to get exact results after an approximate
   computation, when theory tells you that the coefficients must be integers.
  Variant: Also available are \fun{GEN}{grndtoi}{GEN x, long *e} and
   \fun{GEN}{ground}{GEN x}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.round0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def serreverse(*argv):
  '''
  serreverse
  Class: basic
  Section: polynomials
  C-Name: serreverse
  Prototype: G
  Help: serreverse(s): reversion of the power series s.
  Doc: reverse power series of $s$, i.e. the series $t$ such that $t(s) = x$;
   $s$ must be a power series whose valuation is exactly equal to one.
   \bprog
   ? \ps 8
   ? t = serreverse(tan(x))
   %2 = x - 1/3*x^3 + 1/5*x^5 - 1/7*x^7 + O(x^8)
   ? tan(t)
   %3 = x + O(x^8)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.serreverse(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def simplify(*argv):
  '''
  simplify
  Class: basic
  Section: conversions
  C-Name: simplify
  Prototype: G
  Help: simplify(x): simplify the object x as much as possible.
  Doc: 
   this function simplifies $x$ as much as it can. Specifically, a complex or
   quadratic number whose imaginary part is the integer 0 (i.e.~not \kbd{Mod(0,2)}
   or \kbd{0.E-28}) is converted to its real part, and a polynomial of degree $0$
   is converted to its constant term. Simplifications occur recursively.
   
   This function is especially useful before using arithmetic functions,
   which expect integer arguments:
   \bprog
   ? x = 2 + y - y
   %1 = 2
   ? isprime(x)
     ***   at top-level: isprime(x)
     ***                 ^----------
     *** isprime: not an integer argument in an arithmetic function
   ? type(x)
   %2 = "t_POL"
   ? type(simplify(x))
   %3 = "t_INT"
   @eprog
   Note that GP results are simplified as above before they are stored in the
   history. (Unless you disable automatic simplification with \b{y}, that is.)
   In particular
   \bprog
   ? type(%1)
   %4 = "t_INT"
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.simplify(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def taylor(*argv):
  '''
  taylor
  Class: basic
  Section: polynomials
  C-Name: tayl
  Prototype: GnDP
  Help: taylor(x,t,{d=seriesprecision}): taylor expansion of x with respect to
   t, adding O(t^d) to all components of x.
  Doc: Taylor expansion around $0$ of $x$ with respect to
   the simple variable $t$. $x$ can be of any reasonable type, for example a
   rational function. Contrary to \tet{Ser}, which takes the valuation into
   account, this function adds $O(t^d)$ to all components of $x$.
   \bprog
   ? taylor(x/(1+y), y, 5)
   %1 = (y^4 - y^3 + y^2 - y + 1)*x + O(y^5)
   ? Ser(x/(1+y), y, 5)
    ***   at top-level: Ser(x/(1+y),y,5)
    ***                 ^----------------
    *** Ser: main variable must have higher priority in gtoser.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.tayl(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def truncate(*argv):
  '''
  truncate
  Class: basic
  Section: conversions
  C-Name: trunc0
  Prototype: GD&
  Help: truncate(x,{&e}): truncation of x; when x is a power series,take away
   the O(X^). If e is present, do not take into account loss of integer part
   precision, and set e = error estimate in bits.
  Description: 
   (small):small:parens   $1
   (int):int:copy:parens  $1
   (real):int             truncr($1)
   (mp):int               mptrunc($1)
   (mp, &small):int       gcvtoi($1, &$2)
   (mp, &int):int         trunc0($1, &$2)
   (gen):gen              gtrunc($1)
   (gen, &small):gen      gcvtoi($1, &$2)
   (gen, &int):gen        trunc0($1, &$2)
  Doc: truncates $x$ and sets $e$ to the number of
   error bits. When $x$ is in $\R$, this means that the part after the decimal
   point is chopped away, $e$ is the binary exponent of the difference between
   the original and the truncated value (the ``fractional part''). If the
   exponent of $x$ is too large compared to its precision (i.e.~$e>0$), the
   result is undefined and an error occurs if $e$ was not given. The function
   applies componentwise on vector / matrices; $e$ is then the maximal number of
   error bits. If $x$ is a rational function, the result is the ``integer part''
   (Euclidean quotient of numerator by denominator) and $e$ is not set.
   
   Note a very special use of \kbd{truncate}: when applied to a power series, it
   transforms it into a polynomial or a rational function with denominator
   a power of $X$, by chopping away the $O(X^k)$. Similarly, when applied to
   a $p$-adic number, it transforms it into an integer or a rational number
   by chopping away the $O(p^k)$.
  Variant: The following functions are also available: \fun{GEN}{gtrunc}{GEN x}
   and \fun{GEN}{gcvtoi}{GEN x, long *e}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.trunc0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def genus2red(*argv):
  '''
  genus2red
  Class: basic
  Section: elliptic_curves
  C-Name: genus2red
  Prototype: GGDG
  Help: genus2red(Q,P,{p}): let Q,P be polynomials with integer coefficients.
   Determines the reduction at p > 2 of the
   (proper, smooth) hyperelliptic curve C/Q: y^2+Qy = P, of genus 2.
   (The special fiber X_p of the minimal regular model X of C over Z.)
  Doc: Let $Q,P$ be polynomials with integer coefficients.
   Determines the reduction at $p > 2$ of the (proper, smooth) genus~2
   curve $C/\Q$, defined by the hyperelliptic equation $y^2+Qy = P$. (The
   special fiber $X_p$ of the minimal regular model $X$ of $C$ over $\Z$.)
   If $p$ is omitted, determines the reduction type for all (odd) prime
   divisors of the discriminant.
   
   \noindent This function rewritten from an implementation of Liu's algorithm by
   Cohen and Liu (1994), \kbd{genus2reduction-0.3}, see
   \kbd{http://www.math.u-bordeaux1.fr/\til liu/G2R/}.
   
   \misctitle{CAVEAT} The function interface may change: for the
   time being, it returns $[N,\var{FaN}, T, V]$
   where $N$ is either the local conductor at $p$ or the
   global conductor, \var{FaN} is its factorization, $y^2 = T$ defines a
   minimal model over $\Z[1/2]$ and $V$ describes the reduction type at the
   various considered~$p$. Unfortunately, the program is not complete for
   $p = 2$, and we may return the odd part of the conductor only: this is the
   case if the factorization includes the (impossible) term $2^{-1}$; if the
   factorization contains another power of $2$, then this is the exact local
   conductor at $2$ and $N$ is the global conductor.
   
   \bprog
   ? default(debuglevel, 1);
   ? genus2red(0,x^6 + 3*x^3 + 63, 3)
   (potential) stable reduction: [1, []]
   reduction at p: [III{9}] page 184, [3, 3], f = 10
   %1 = [59049, Mat([3, 10]), x^6 + 3*x^3 + 63, [3, [1, []],
          ["[III{9}] page 184", [3, 3]]]]
   ? [N, FaN, T, V] = genus2red(x^3-x^2-1, x^2-x);  \\ X_1(13), global reduction
   p = 13
   (potential) stable reduction: [5, [Mod(0, 13), Mod(0, 13)]]
   reduction at p: [I{0}-II-0] page 159, [], f = 2
   ? N
   %3 = 169
   ? FaN
   %4 = Mat([13, 2])   \\ in particular, good reduction at 2 !
   ? T
   %5 = x^6 + 58*x^5 + 1401*x^4 + 18038*x^3 + 130546*x^2 + 503516*x + 808561
   ? V
   %6 = [[13, [5, [Mod(0, 13), Mod(0, 13)]], ["[I{0}-II-0] page 159", []]]]
   @eprog\noindent
   We now first describe the format of the vector $V = V_p$ in the case where
   $p$ was specified (local reduction at~$p$): it is a triple $[p, \var{stable},
   \var{red}]$. The component $\var{stable} = [\var{type}, \var{vecj}]$ contains
   information about the stable reduction after a field extension;
   depending on \var{type}s, the stable reduction is
   
   \item 1: smooth (i.e. the curve has potentially good reduction). The
         Jacobian $J(C)$ has potentially good reduction.
   
   \item 2: an elliptic curve $E$ with an ordinary double point; \var{vecj}
   contains $j$ mod $p$, the modular invariant of $E$. The (potential)
   semi-abelian reduction of $J(C)$ is the extension of an elliptic curve (with
   modular invariant $j$ mod $p$) by a torus.
   
   \item 3: a projective line with two ordinary double points. The Jacobian
   $J(C)$ has potentially multiplicative reduction.
   
   \item 4: the union of two projective lines crossing transversally at three
   points. The Jacobian $J(C)$ has potentially multiplicative reduction.
   
   \item 5: the union of two elliptic curves $E_1$ and $E_2$ intersecting
   transversally at one point; \var{vecj} contains their modular invariants
   $j_1$ and $j_2$, which may live in a quadratic extension of $\F_p$ are need
   not be distinct. The Jacobian $J(C)$ has potentially good reduction,
   isomorphic to the product of the reductions of $E_1$ and $E_2$.
   
   \item 6: the union of an elliptic curve $E$ and a projective line which has
   an ordinary double point, and these two components intersect transversally
   at one point; \var{vecj} contains $j$ mod $p$, the modular invariant of $E$.
   The (potential) semi-abelian reduction of $J(C)$ is the extension of an
   elliptic curve (with modular invariant $j$ mod $p$) by a torus.
   
   \item 7: as in type 6, but the two components are both singular. The
   Jacobian $J(C)$ has potentially multiplicative reduction.
   
   The component $\var{red} = [\var{NUtype}, \var{neron}]$ contains two data
   concerning the reduction at $p$ without any ramified field extension.
   
   The \var{NUtype} is a \typ{STR} describing the reduction at $p$ of $C$,
   following Namikawa-Ueno, \emph{The complete classification of fibers in
   pencils of curves of genus two}, Manuscripta Math., vol. 9, (1973), pages
   143-186. The reduction symbol is followed by the corresponding page number in
   this article.
   
   The second datum \var{neron} is the group of connected components (over an
   algebraic closure of $\F_p$) of the N\'eron model of $J(C)$, given as a
   finite abelian group (vector of elementary divisors).
   \smallskip
   If $p = 2$, the \var{red} component may be omitted altogether (and
   replaced by \kbd{[]}, in the case where the program could not compute it.
   When $p$ was not specified, $V$ is the vector of all $V_p$, for all
   considered $p$.
   
   \misctitle{Notes about Namikawa-Ueno types}
   
   \item A lower index is denoted between braces: for instance, \kbd{[I\obr
    2\cbr-II-5]} means \kbd{[I\_2-II-5]}.
   
   \item If $K$ and $K'$ are Kodaira symbols for singular fibers of elliptic
   curves, \kbd{[$K$-$K'$-m]} and \kbd{[$K'$-$K$-m]} are the same.
   
   \item \kbd{[$K$-$K'$-$-1$]}  is \kbd{[$K'$-$K$-$\alpha$]} in the notation of
   Namikawa-Ueno.
   
   \item The figure \kbd{[2I\_0-m]} in Namikawa-Ueno, page 159, must be denoted
   by \kbd{[2I\_0-(m+1)]}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.genus2red(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathnfmod(*argv):
  '''
  mathnfmod
  Class: basic
  Section: linear_algebra
  C-Name: hnfmod
  Prototype: GG
  Help: mathnfmod(x,d): (upper triangular) Hermite normal form of x, basis for
   the lattice formed by the columns of x, where d is a multiple of the
   non-zero determinant of this lattice.
  Doc: if $x$ is a (not necessarily square) matrix of
   maximal rank with integer entries, and $d$ is a multiple of the (non-zero)
   determinant of the lattice spanned by the columns of $x$, finds the
   \emph{upper triangular} \idx{Hermite normal form} of $x$.
   
   If the rank of $x$ is equal to its number of rows, the result is a square
   matrix. In general, the columns of the result form a basis of the lattice
   spanned by the columns of $x$. Even when $d$ is known, this is in general
   slower than \kbd{mathnf} but uses much less memory.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hnfmod(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathnfmodid(*argv):
  '''
  mathnfmodid
  Class: basic
  Section: linear_algebra
  C-Name: hnfmodid
  Prototype: GG
  Help: mathnfmodid(x,d): (upper triangular) Hermite normal form of x
   concatenated with matdiagonal(d)
  Doc: outputs the (upper triangular)
   \idx{Hermite normal form} of $x$ concatenated with the diagonal
   matrix with diagonal $d$. Assumes that $x$ has integer entries.
   Variant: if $d$ is an integer instead of a vector, concatenate $d$ times the
   identity matrix.
   \bprog
   ? m=[0,7;-1,0;-1,-1]
   %1 =
   [ 0  7]
   
   [-1  0]
   
   [-1 -1]
   ? mathnfmodid(m, [6,2,2])
   %2 =
   [2 1 1]
   
   [0 1 0]
   
   [0 0 1]
   ? mathnfmodid(m, 10)
   %3 =
   [10 7 3]
   
   [ 0 1 0]
   
   [ 0 0 1]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hnfmodid(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matfrobenius(*argv):
  '''
  matfrobenius
  Class: basic
  Section: linear_algebra
  C-Name: matfrobenius
  Prototype: GD0,L,Dn
  Help: matfrobenius(M,{flag},{v='x}): Return the Frobenius form of the square
   matrix M. If flag is 1, return only the elementary divisors as a vector of
   polynomials in the variable v. If flag is 2, return a two-components vector
   [F,B] where F is the Frobenius form and B is the basis change so that
   M=B^-1*F*B.
  Doc: returns the Frobenius form of
   the square matrix \kbd{M}. If $\fl=1$, returns only the elementary divisors as
   a vector of polynomials in the variable \kbd{v}.  If $\fl=2$, returns a
   two-components vector [F,B] where \kbd{F} is the Frobenius form and \kbd{B} is
   the basis change so that $M=B^{-1}FB$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matfrobenius(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def mathnf(*argv):
  '''
  mathnf
  Class: basic
  Section: linear_algebra
  C-Name: mathnf0
  Prototype: GD0,L,
  Help: mathnf(M,{flag=0}): (upper triangular) Hermite normal form of M, basis
   for the lattice formed by the columns of M. flag is optional whose value
   range from 0 to 3 have a binary meaning. Bit 1: complete output, returns
   a 2-component vector [H,U] such that H is the HNF of M, and U is an
   invertible matrix such that MU=H. Bit 2: allow polynomial entries, otherwise
   assume that M is integral. These use a naive algorithm; larger values
   correspond to more involved algorithms and are restricted to integer
   matrices; flag = 4: returns [H,U] using LLL reduction along the way;
   flag = 5: return [H,U,P] where P is a permutation of row indices such that
   P applied to M U is H.
  Doc: let $R$ be a Euclidean ring, equal to $\Z$ or to $K[X]$ for some field
   $K$. If $M$ is a (not necessarily square) matrix with entries in $R$, this
   routine finds the \emph{upper triangular} \idx{Hermite normal form} of $M$.
   If the rank of $M$ is equal to its number of rows, this is a square
   matrix. In general, the columns of the result form a basis of the $R$-module
   spanned by the columns of $M$.
   
   The values $0,1,2,3$ of $\fl$ have a binary meaning, analogous to the one
   in \tet{matsnf}; in this case, binary digits of $\fl$ mean:
   
   \item 1 (complete output): if set, outputs $[H,U]$, where $H$ is the Hermite
   normal form of $M$, and $U$ is a transformation matrix such that $MU=[0|H]$.
   The matrix $U$ belongs to $\text{GL}(R)$. When $M$ has a large kernel, the
   entries of $U$ are in general huge.
   
   \item 2 (generic input): \emph{Deprecated}. If set, assume that $R = K[X]$ is
   a polynomial ring; otherwise, assume that $R = \Z$. This flag is now useless
   since the routine always checks whether the matrix has integral entries.
   
   \noindent For these 4 values, we use a naive algorithm, which behaves well
   in small dimension only. Larger values correspond to different algorithms,
   are restricted to \emph{integer} matrices, and all output the unimodular
   matrix $U$. From now on all matrices have integral entries.
   
   \item $\fl=4$, returns $[H,U]$ as in ``complete output'' above, using a
   variant of \idx{LLL} reduction along the way. The matrix $U$ is provably
   small in the $L_2$ sense, and in general close to optimal; but the
   reduction is in general slow, although provably polynomial-time.
   
   If $\fl=5$, uses Batut's algorithm and output $[H,U,P]$, such that $H$ and
   $U$ are as before and $P$ is a permutation of the rows such that $P$ applied
   to $MU$ gives $H$. This is in general faster than $\fl=4$ but the matrix $U$
   is usually worse; it is heuristically smaller than with the default algorithm.
   
   When the matrix is dense and the dimension is large (bigger than 100, say),
   $\fl = 4$ will be fastest. When $M$ has maximal rank, then
   \bprog
     H = mathnfmod(M, matdetint(M))
   @eprog\noindent will be even faster. You can then recover $U$ as $M^{-1}H$.
   
   \bprog
   ? M = matrix(3,4,i,j,random([-5,5]))
   %1 =
   [ 0 2  3  0]
   
   [-5 3 -5 -5]
   
   [ 4 3 -5  4]
   
   ? [H,U] = mathnf(M, 1);
   ? U
   %3 =
   [-1 0 -1 0]
   
   [ 0 5  3 2]
   
   [ 0 3  1 1]
   
   [ 1 0  0 0]
   
   ? H
   %5 =
   [19 9 7]
   
   [ 0 9 1]
   
   [ 0 0 1]
   
   ? M*U
   %6 =
   [0 19 9 7]
   
   [0  0 9 1]
   
   [0  0 0 1]
   @eprog
   
   For convenience, $M$ is allowed to be a \typ{VEC}, which is then
   automatically converted to a \typ{MAT}, as per the \tet{Mat} function.
   For instance to solve the generalized extended gcd problem, one may use
   \bprog
   ? v = [116085838, 181081878, 314252913,10346840];
   ? [H,U] = mathnf(v, 1);
   ? U
   %2 =
   [ 103 -603    15  -88]
   
   [-146   13 -1208  352]
   
   [  58  220   678 -167]
   
   [-362 -144   381 -101]
   ? v*U
   %3 = [0, 0, 0, 1]
   @eprog\noindent This also allows to input a matrix as a \typ{VEC} of
   \typ{COL}s of the same length (which \kbd{Mat} would concatenate to
   the \typ{MAT} having those columns):
   \bprog
   ? v = [[1,0,4]~, [3,3,4]~, [0,-4,-5]~]; mathnf(v)
   %1 =
   [47 32 12]
   
   [ 0  1  0]
   
   [ 0  0  1]
   @eprog
  Variant: Also available are \fun{GEN}{hnf}{GEN M} ($\fl=0$) and
   \fun{GEN}{hnfall}{GEN M} ($\fl=1$). To reduce \emph{huge} relation matrices
   (sparse with small entries, say dimension $400$ or more), you can use the
   pair \kbd{hnfspec} / \kbd{hnfadd}. Since this is quite technical and the
   calling interface may change, they are not documented yet. Look at the code
   in \kbd{basemath/hnf\_snf.c}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mathnf0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matsnf(*argv):
  '''
  matsnf
  Class: basic
  Section: linear_algebra
  C-Name: matsnf0
  Prototype: GD0,L,
  Help: matsnf(X,{flag=0}): Smith normal form (i.e. elementary divisors) of
   the matrix X, expressed as a vector d. Binary digits of flag mean 1: returns
   [u,v,d] where d=u*X*v, otherwise only the diagonal d is returned, 2: allow
   polynomial entries, otherwise assume X is integral, 4: removes all
   information corresponding to entries equal to 1 in d.
  Doc: if $X$ is a (singular or non-singular) matrix outputs the vector of
   \idx{elementary divisors} of $X$, i.e.~the diagonal of the
   \idx{Smith normal form} of $X$, normalized so that $d_n \mid d_{n-1} \mid
   \ldots \mid d_1$.
   
   The binary digits of \fl\ mean:
   
   1 (complete output): if set, outputs $[U,V,D]$, where $U$ and $V$ are two
   unimodular matrices such that $UXV$ is the diagonal matrix $D$. Otherwise
   output only the diagonal of $D$. If $X$ is not a square matrix, then $D$
   will be a square diagonal matrix padded with zeros on the left or the top.
   
   2 (generic input): if set, allows polynomial entries, in which case the
   input matrix must be square. Otherwise, assume that $X$ has integer
   coefficients with arbitrary shape.
   
   4 (cleanup): if set, cleans up the output. This means that elementary
   divisors equal to $1$ will be deleted, i.e.~outputs a shortened vector $D'$
   instead of $D$. If complete output was required, returns $[U',V',D']$ so
   that $U'XV' = D'$ holds. If this flag is set, $X$ is allowed to be of the
   form `vector of elementary divisors' or $[U,V,D]$ as would normally be output with the cleanup flag
   unset.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matsnf0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bigomega(*argv):
  '''
  bigomega
  Class: basic
  Section: number_theoretical
  C-Name: bigomega
  Prototype: lG
  Help: bigomega(x): number of prime divisors of x, counted with multiplicity.
  Doc: number of prime divisors of the integer $|x|$ counted with
   multiplicity:
   \bprog
   ? factor(392)
   %1 =
   [2 3]
   
   [7 2]
   
   ? bigomega(392)
   %2 = 5;  \\ = 3+2
   ? omega(392)
   %3 = 2;  \\ without multiplicity
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.bigomega(*c_arg_tuple)

def eulerphi(*argv):
  '''
  eulerphi
  Class: basic
  Section: number_theoretical
  C-Name: eulerphi
  Prototype: G
  Help: eulerphi(x): Euler's totient function of x.
  Description: 
   (gen):int        eulerphi($1)
  Doc: Euler's $\phi$ (totient)\sidx{Euler totient function} function of the
   integer $|x|$, in other words $|(\Z/x\Z)^*|$.
   \bprog
   ? eulerphi(40)
   %1 = 16
   @eprog\noindent
   According to this definition we let $\phi(0) := 2$, since $\Z^* = \{-1,1\}$;
   this is consistant with \kbd{znstar(0)}: we have \kbd{znstar$(n)$.no =
   eulerphi(n)} for all $n\in\Z$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.eulerphi(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factorint(*argv):
  '''
  factorint
  Class: basic
  Section: number_theoretical
  C-Name: factorint
  Prototype: GD0,L,
  Help: factorint(x,{flag=0}): factor the integer x. flag is optional, whose
   binary digits mean 1: avoid MPQS, 2: avoid first-stage ECM (may fall back on
   it later), 4: avoid Pollard-Brent Rho and Shanks SQUFOF, 8: skip final ECM
   (huge composites will be declared prime).
  Doc: factors the integer $n$ into a product of
   pseudoprimes (see \kbd{ispseudoprime}), using a combination of the
   \idx{Shanks SQUFOF} and \idx{Pollard Rho} method (with modifications due to
   Brent), \idx{Lenstra}'s \idx{ECM} (with modifications by Montgomery), and
   \idx{MPQS} (the latter adapted from the \idx{LiDIA} code with the kind
   permission of the LiDIA maintainers), as well as a search for pure powers.
   The output is a two-column matrix as for \kbd{factor}: the first column
   contains the ``prime'' divisors of $n$, the second one contains the
   (positive) exponents.
   
   By convention $0$ is factored as $0^1$, and $1$ as the empty factorization;
   also the divisors are by default not proven primes is they are larger than
   $2^{64}$, they only failed the BPSW compositeness test (see
   \tet{ispseudoprime}). Use \kbd{isprime} on the result if you want to
   guarantee primality or set the \tet{factor_proven} default to $1$.
   Entries of the private prime tables (see \tet{addprimes}) are also included
   as is.
   
   This gives direct access to the integer factoring engine called by most
   arithmetical functions. \fl\ is optional; its binary digits mean 1: avoid
   MPQS, 2: skip first stage ECM (we may still fall back to it later), 4: avoid
   Rho and SQUFOF, 8: don't run final ECM (as a result, a huge composite may be
   declared to be prime). Note that a (strong) probabilistic primality test is
   used; thus composites might not be detected, although no example is known.
   
   You are invited to play with the flag settings and watch the internals at
   work by using \kbd{gp}'s \tet{debug} default parameter (level 3 shows
   just the outline, 4 turns on time keeping, 5 and above show an increasing
   amount of internal details).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factorint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ispowerful(*argv):
  '''
  ispowerful
  Class: basic
  Section: number_theoretical
  C-Name: ispowerful
  Prototype: lG
  Help: ispowerful(x): true(1) if x is a powerful integer (valuation at all
   primes is greater than 1), false(0) if not.
  Doc: true (1) if $x$ is a powerful integer, false (0) if not;
   an integer is powerful if and only if its valuation at all primes is
   greater than 1.
   \bprog
   ? ispowerful(50)
   %1 = 0
   ? ispowerful(100)
   %2 = 1
   ? ispowerful(5^3*(10^1000+1)^2)
   %3 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.ispowerful(*c_arg_tuple)

def issquarefree(*argv):
  '''
  issquarefree
  Class: basic
  Section: number_theoretical
  C-Name: issquarefree
  Prototype: lG
  Help: issquarefree(x): true(1) if x is squarefree, false(0) if not.
  Description: 
   (gen):bool       issquarefree($1)
  Doc: true (1) if $x$ is squarefree, false (0) if not. Here $x$ can be an
   integer or a polynomial.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.issquarefree(*c_arg_tuple)

def istotient(*argv):
  '''
  istotient
  Class: basic
  Section: number_theoretical
  C-Name: istotient
  Prototype: lGD&
  Help: istotient(x,{&N}): true(1) if x = eulerphi(n) for some integer n,
   false(0) if not. If N is given, set N = n as well.
  Doc: true (1) if $x = \phi(n)$ for some integer $n$, false (0)
   if not.
   \bprog
   ? istotient(14)
   %1 = 0
   ? istotient(100)
   %2 = 0
   @eprog
   If $N$ is given, set $N = n$ as well.
   \bprog
   ? istotient(4, &n)
   %1 = 1
   ? n
   %2 = 10
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.istotient(*c_arg_tuple)

def moebius(*argv):
  '''
  moebius
  Class: basic
  Section: number_theoretical
  C-Name: moebius
  Prototype: lG
  Help: moebius(x): Moebius function of x.
  Doc: \idx{Moebius} $\mu$-function of $|x|$. $x$ must be of type integer.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.moebius(*c_arg_tuple)

def nextprime(*argv):
  '''
  nextprime
  Class: basic
  Section: number_theoretical
  C-Name: nextprime
  Prototype: G
  Help: nextprime(x): smallest pseudoprime >= x.
  Description: 
   (gen):int        nextprime($1)
  Doc: finds the smallest pseudoprime (see
   \tet{ispseudoprime}) greater than or equal to $x$. $x$ can be of any real
   type. Note that if $x$ is a pseudoprime, this function returns $x$ and not
   the smallest pseudoprime strictly larger than $x$. To rigorously prove that
   the result is prime, use \kbd{isprime}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nextprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def numdiv(*argv):
  '''
  numdiv
  Class: basic
  Section: number_theoretical
  C-Name: numdiv
  Prototype: G
  Help: numdiv(x): number of divisors of x.
  Description: 
   (gen):int        numdiv($1)
  Doc: number of divisors of $|x|$. $x$ must be of type integer.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.numdiv(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def omega(*argv):
  '''
  omega
  Class: basic
  Section: number_theoretical
  C-Name: omega
  Prototype: lG
  Help: omega(x): number of distinct prime divisors of x.
  Doc: number of distinct prime divisors of $|x|$. $x$ must be of type integer.
   \bprog
   ? factor(392)
   %1 =
   [2 3]
   
   [7 2]
   
   ? omega(392)
   %2 = 2;  \\ without multiplicity
   ? bigomega(392)
   %3 = 5;  \\ = 3+2, with multiplicity
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.omega(*c_arg_tuple)

def precprime(*argv):
  '''
  precprime
  Class: basic
  Section: number_theoretical
  C-Name: precprime
  Prototype: G
  Help: precprime(x): largest pseudoprime <= x, 0 if x<=1.
  Description: 
   (gen):int        precprime($1)
  Doc: finds the largest pseudoprime (see
   \tet{ispseudoprime}) less than or equal to $x$. $x$ can be of any real type.
   Returns 0 if $x\le1$. Note that if $x$ is a prime, this function returns $x$
   and not the largest prime strictly smaller than $x$. To rigorously prove that
   the result is prime, use \kbd{isprime}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.precprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sigma(*argv):
  '''
  sigma
  Class: basic
  Section: number_theoretical
  C-Name: sumdivk
  Prototype: GD1,L,
  Help: sigma(x,{k=1}): sum of the k-th powers of the divisors of x. k is
   optional and if omitted is assumed to be equal to 1.
  Description: 
   (gen, ?1):int           sumdiv($1)
   (gen, 0):int            numdiv($1)
  Doc: sum of the $k^{\text{th}}$ powers of the positive divisors of $|x|$. $x$
   and $k$ must be of type integer.
  Variant: Also available is \fun{GEN}{sumdiv}{GEN n}, for $k = 1$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sumdivk(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def errname(*argv):
  '''
  errname
  Class: basic
  Section: programming/specific
  C-Name: errname
  Prototype: G
  Help: errname(E): returns the type of the error message E.
  Description: 
   (gen):errtyp err_get_num($1)
  Doc: returns the type of the error message \kbd{E} as a string.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.errname(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def getheap(*argv):
  '''
  getheap
  Class: basic
  Section: programming/specific
  C-Name: getheap
  Prototype: 
  Help: getheap(): 2-component vector giving the current number of objects in
   the heap and the space they occupy.
  Doc: returns a two-component row vector giving the
   number of objects on the heap and the amount of memory they occupy in long
   words. Useful mainly for debugging purposes.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.getheap(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sizebyte(*argv):
  '''
  sizebyte
  Class: basic
  Section: conversions
  C-Name: gsizebyte
  Prototype: lG
  Help: sizebyte(x): number of bytes occupied by the complete tree of the
   object x.
  Doc: outputs the total number of bytes occupied by the tree representing the
   PARI object $x$.
  Variant: Also available is \fun{long}{gsizeword}{GEN x} returning a
   number of \emph{words}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.gsizebyte(*c_arg_tuple)

def version(*argv):
  '''
  version
  Class: basic
  Section: programming/specific
  C-Name: pari_version
  Prototype: 
  Help: version(): returns the PARI version as [major,minor,patch] or [major,minor,patch,VCSversion].
  Doc: returns the current version number as a \typ{VEC} with three integer
   components (major version number, minor version number and patchlevel);
   if your sources were obtained through our version control system, this will
   be followed by further more precise arguments, including
   e.g.~a~\kbd{git} \emph{commit hash}.
   
   This function is present in all versions of PARI following releases 2.3.4
   (stable) and 2.4.3 (testing).
   
   Unless you are working with multiple development versions, you probably only
   care about the 3 first numeric components. In any case, the \kbd{lex} function
   offers a clever way to check against a particular version number, since it will
   compare each successive vector entry, numerically or as strings, and will not
   mind if the vectors it compares have different lengths:
   \bprog
      if (lex(version(), [2,3,5]) >= 0,
        \\ code to be executed if we are running 2.3.5 or more recent.
      ,
        \\ compatibility code
      );
   @eprog\noindent On a number of different machines, \kbd{version()} could return either of
   \bprog
    %1 = [2, 3, 4]    \\ released version, stable branch
    %1 = [2, 4, 3]    \\ released version, testing branch
    %1 = [2, 6, 1, 15174, ""505ab9b"] \\ development
   @eprog
   
   In particular, if you are only working with released versions, the first
   line of the gp introductory message can be emulated by
   \bprog
      [M,m,p] = version();
      printf("GP/PARI CALCULATOR Version %s.%s.%s", M,m,p);
    @eprog\noindent If you \emph{are} working with many development versions of
    PARI/GP, the 4th and/or 5th components can be profitably included in the
    name of your logfiles, for instance.
   
    \misctitle{Technical note} For development versions obtained via \kbd{git},
    the 4th and 5th components are liable to change eventually, but we document
    their current meaning for completeness. The 4th component counts the number
    of reachable commits in the branch (analogous to \kbd{svn}'s revision
    number), and the 5th is the \kbd{git} commit hash. In particular, \kbd{lex}
    comparison still orders correctly development versions with respect to each
    others or to released versions (provided we stay within a given branch,
    e.g. \kbd{master})!
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.pari_version(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def trap(*argv):
  '''
  trap
  Class: basic
  Section: programming/specific
  C-Name: trap0
  Prototype: DrDEDE
  Help: trap({e}, {rec}, seq): try to execute seq, trapping runtime error e (all
   of them if e omitted); sequence rec is executed if the error occurs and
   is the result of the command. THIS FUNCTION IS OBSOLETE: use "IFERR"
  Wrapper: (,_,_)
  Description: 
   (?str,?closure,?closure):gen trap0($1, $2, $3)
  Doc: THIS FUNCTION IS OBSOLETE: use \tet{iferr}, which has a nicer and much
   more powerful interface. For compatibility's sake we now describe the
   \emph{obsolete} function \tet{trap}.
   
   This function tries to
   evaluate \var{seq}, trapping runtime error $e$, that is effectively preventing
   it from aborting computations in the usual way; the recovery sequence
   \var{rec} is executed if the error occurs and the evaluation of \var{rec}
   becomes the result of the command. If $e$ is omitted, all exceptions are
   trapped. See \secref{se:errorrec} for an introduction to error recovery
   under \kbd{gp}.
   
   \bprog
   ? \\@com trap division by 0
   ? inv(x) = trap (e_INV, INFINITY, 1/x)
   ? inv(2)
   %1 = 1/2
   ? inv(0)
   %2 = INFINITY
   @eprog\noindent
   Note that \var{seq} is effectively evaluated up to the point that produced
   the error, and the recovery sequence is evaluated starting from that same
   context, it does not "undo" whatever happened in the other branch (restore
   the evaluation context):
   \bprog
   ? x = 1; trap (, /* recover: */ x, /* try: */ x = 0; 1/x)
   %1 = 0
   @eprog
   
   \misctitle{Note} The interface is currently not adequate for trapping
   individual exceptions. In the current version \vers, the following keywords
   are recognized, but the name list will be expanded and changed in the
   future (all library mode errors can be trapped: it's a matter of defining
   the keywords to \kbd{gp}):
   
   \kbd{e\_ALARM}: alarm time-out
   
   \kbd{e\_ARCH}: not available on this architecture or operating system
   
   \kbd{e\_STACK}: the PARI stack overflows
   
   \kbd{e\_INV}: impossible inverse
   
   \kbd{e\_IMPL}: not yet implemented
   
   \kbd{e\_OVERFLOW}: all forms of arithmetic overflow, including length
   or exponent overflow (when a larger value is supplied than the
   implementation can handle).
   
   \kbd{e\_SYNTAX}: syntax error
   
   \kbd{e\_MISC}: miscellaneous error
   
   \kbd{e\_TYPE}: wrong type
   
   \kbd{e\_USER}: user error (from the \kbd{error} function)
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.trap0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def intmellininvshort(*argv):
  '''
  intmellininvshort
  Class: basic
  Section: sums
  C-Name: intmellininvshort
  Prototype: GGGp
  Help: intmellininvshort(sig,z,tab): numerical integration on the
   line real(X) = sig (or sig[1]) of s(X)z^(-X)dX/(2*I*Pi), i.e. inverse Mellin
   transform of s at z. sig is coded as follows: either it is real, and then
   by default assume s(X) decreases like exp(-X). Or sig = [sigR, al], sigR is
   the abscissa of integration, and al = 0 for slowly decreasing functions, or
   al > 0 if s(X) decreases like exp(-al*X). Compulsory table tab has been
   precomputed using the command intfuncinit(t=[[-1],sig[2]],[[1],sig[2]],s)
   (with possibly its two optional additional parameters), where sig[2] = 1
   if not given. Orders of magnitude faster than intmellininv.
  Doc: numerical integration
   of $(2i\pi)^{-1}s(X)z^{-X}$ with respect to $X$ on the line $\Re(X)=sig$.
   In other words, inverse Mellin transform of $s(X)$ at the value $z$.
   Here $s(X)$ is implicitly contained in \var{tab} in \kbd{intfuncinit} format,
   typically
   \bprog
   tab = intfuncinit(T = [-1], [1], s(sig + I*T))
   @eprog\noindent
   or similar commands. Take the example of the inverse Mellin transform of
   $\Gamma(s)^3$ given in \kbd{intmellininv}:
   
   \bprog
   ? \p 105
   ? oo = [1]; \\@com for clarity
   ? A = intmellininv(s=2,4, gamma(s)^3);
   time = 2,500 ms. \\@com not too fast because of $\Gamma(s)^3$.
   \\ @com function of real type, decreasing as $\exp(-3\pi/2\cdot |t|)$
   ? tab = intfuncinit(t=[-oo, 3*Pi/2],[oo, 3*Pi/2], gamma(2+I*t)^3, 1);
   time = 1,370 ms.
   ? intmellininvshort(2,4, tab) - A
   time = 50 ms.
   %4 = -1.26... - 3.25...E-109*I \\@com 50 times faster than \kbd{A} and perfect.
   ? tab2 = intfuncinit(t=-oo, oo, gamma(2+I*t)^3, 1);
   ? intmellininvshort(2,4, tab2)
   %6 = -1.2...E-42 - 3.2...E-109*I  \\@com 63 digits lost
   @eprog\noindent
   In the computation of \var{tab}, it was not essential to include the
   \emph{exact} exponential decrease of $\Gamma(2+it)^3$. But as the last
   example shows, a rough indication \emph{must} be given, otherwise slow
   decrease is assumed, resulting in catastrophic loss of accuracy.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.intmellininvshort(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def intnuminit(*argv):
  '''
  intnuminit
  Class: basic
  Section: sums
  C-Name: intnuminit
  Prototype: GGD0,L,p
  Help: intnuminit(a,b,{m=0}): initialize tables for integrations from a to b.
   See help for intnum for coding of a and b. Possible types: compact interval,
   semi-compact (one extremity at + or - infinity) or R, and very slowly, slowly
   or exponentially decreasing, or sine or cosine oscillating at infinities.
  Doc: initialize tables for integration from
   $a$ to $b$, where $a$ and $b$ are coded as in \kbd{intnum}. Only the
   compactness, the possible existence of singularities, the speed of decrease
   or the oscillations at infinity are taken into account, and not the values.
   For instance {\tt intnuminit(-1,1)} is equivalent to {\tt intnuminit(0,Pi)},
   and {\tt intnuminit([0,-1/2],[1])} is equivalent to {\tt
   intnuminit([-1],[-1,-1/2])}. If $m$ is not given, it is computed according to
   the current precision. Otherwise the integration step is $1/2^m$. Reasonable
   values of $m$ are $m=6$ or $m=7$ for $100$ decimal digits, and $m=9$ for
   $1000$ decimal digits.
   
   The result is technical, but in some cases it is useful to know the output.
   Let $x=\phi(t)$ be the change of variable which is used. \var{tab}[1] contains
   the integer $m$ as above, either given by the user or computed from the default
   precision, and can be recomputed directly using the function \kbd{intnumstep}.
   \var{tab}[2] and \var{tab}[3] contain respectively the abscissa and weight
   corresponding to $t=0$ ($\phi(0)$ and $\phi'(0)$). \var{tab}[4] and
   \var{tab}[5] contain the abscissas and weights corresponding to positive
   $t=nh$ for $1\le n\le N$ and $h=1/2^m$ ($\phi(nh)$ and $\phi'(nh)$). Finally
   \var{tab}[6] and \var{tab}[7] contain either the abscissas and weights
   corresponding to negative $t=nh$ for $-N\le n\le -1$, or may be empty (but
   not always) if $\phi(t)$ is an odd function (implicitly we would have
   $\var{tab}[6]=-\var{tab}[4]$ and $\var{tab}[7]=\var{tab}[5]$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.intnuminit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def intnumstep(*argv):
  '''
  intnumstep
  Class: basic
  Section: sums
  C-Name: intnumstep
  Prototype: lp
  Help: intnumstep(): gives the default value of m used by all intnum and sumnum
   routines, such that the integration step is 1/2^m.
  Doc: give the value of $m$ used in all the
   \kbd{intnum} and \kbd{sumnum} programs, hence such that the integration
   step is equal to $1/2^m$.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.intnumstep(*c_arg_tuple)

def sumnuminit(*argv):
  '''
  sumnuminit
  Class: basic
  Section: sums
  C-Name: sumnuminit
  Prototype: GD0,L,D1,L,p
  Help: sumnuminit(sig, {m=0}, {sgn=1}): initialize tables for numerical
   summation. sgn is 1 (in fact >= 0), the default, for sumnum (ordinary sums)
   or -1 (in fact < 0) for sumnumalt (alternating sums). sig is as in sumnum and
   m is as in intnuminit.
  Doc: initialize tables for numerical summation using \kbd{sumnum} (with
   $\var{sgn}=1$) or \kbd{sumnumalt} (with $\var{sgn}=-1$), $sig$ is the
   abscissa of integration coded as in \kbd{sumnum}, and $m$ is as in
   \kbd{intnuminit}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sumnuminit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def padicfields(*argv):
  '''
  padicfields
  Class: basic
  Section: polynomials
  C-Name: padicfields0
  Prototype: GGD0,L,
  Help: padicfields(p, N, {flag=0}): returns polynomials generating all
   the extensions of degree N of the field of p-adic rational numbers; N is
   allowed to be a 2-component vector [n,d], in which case, returns the
   extensions of degree n and discriminant p^d. flag is optional,
   and can be 0: default, 1: return also the ramification index, the residual
   degree, the valuation of the discriminant and the number of conjugate fields,
   or 2: return only the number of extensions in a fixed algebraic closure.
  Doc: returns a vector of polynomials generating all the extensions of degree
   $N$ of the field $\Q_p$ of $p$-adic rational numbers; $N$ is
   allowed to be a 2-component vector $[n,d]$, in which case we return the
   extensions of degree $n$ and discriminant $p^d$.
   
   The list is minimal in the sense that two different polynomials generate
   non-isomorphic extensions; in particular, the number of polynomials is the
   number of classes of non-isomorphic extensions. If $P$ is a polynomial in this
   list, $\alpha$ is any root of $P$ and $K = \Q_p(\alpha)$, then $\alpha$
   is the sum of a uniformizer and a (lift of a) generator of the residue field
   of $K$; in particular, the powers of $\alpha$ generate the ring of $p$-adic
   integers of $K$.
   
   If $\fl = 1$, replace each polynomial $P$ by a vector $[P, e, f, d, c]$
   where $e$ is the ramification index, $f$ the residual degree, $d$ the
   valuation of the discriminant, and $c$ the number of conjugate fields.
   If $\fl = 2$, only return the \emph{number} of extensions in a fixed
   algebraic closure (Krasner's formula), which is much faster.
  Variant: Also available is
   \fun{GEN}{padicfields}{GEN p, long n, long d, long flag}, which computes
   extensions of $\Q_p$ of degree $n$ and discriminant $p^d$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.padicfields0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfkummer(*argv):
  '''
  rnfkummer
  Class: basic
  Section: number_fields
  C-Name: rnfkummer
  Prototype: GDGD0,L,p
  Help: rnfkummer(bnr,{subgp},{d=0}): bnr being as output by bnrinit,
   finds a relative equation for the class field corresponding to the module in
   bnr and the given congruence subgroup (the ray class field if subgp is
   omitted). d can be zero (default), or positive, and in this case the
   output is the list of all relative equations of degree d for the given bnr,
   with the same conductor as (bnr, subgp).
  Doc: \var{bnr}
   being as output by \kbd{bnrinit}, finds a relative equation for the
   class field corresponding to the module in \var{bnr} and the given
   congruence subgroup (the full ray class field if \var{subgp} is omitted).
   If $d$ is positive, outputs the list of all relative equations of
   degree $d$ contained in the ray class field defined by \var{bnr}, with
   the \emph{same} conductor as $(\var{bnr}, \var{subgp})$.
   
   \misctitle{Warning} This routine only works for subgroups of prime index. It
   uses Kummer theory, adjoining necessary roots of unity (it needs to compute a
   tough \kbd{bnfinit} here), and finds a generator via Hecke's characterization
   of ramification in Kummer extensions of prime degree. If your extension does
   not have prime degree, for the time being, you have to split it by hand as a
   tower / compositum of such extensions.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfkummer(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def matkerint(*argv):
  '''
  matkerint
  Class: basic
  Section: linear_algebra
  C-Name: matkerint0
  Prototype: GD0,L,
  Help: matkerint(x,{flag=0}): LLL-reduced Z-basis of the kernel of the matrix
   x with integral entries. flag is optional, and may be set to 0: default,
   uses LLL, 1: uses matrixqz (much slower).
  Doc: gives an \idx{LLL}-reduced $\Z$-basis
   for the lattice equal to the kernel of the matrix $x$ as columns of the
   matrix $x$ with integer entries (rational entries are not permitted).
   
   If $\fl=0$, uses an integer LLL algorithm.
   
   If $\fl=1$, uses $\kbd{matrixqz}(x,-2)$. Many orders of magnitude slower
   than the default: never use this.
  Variant: See also \fun{GEN}{kerint}{GEN x} ($\fl=0$), which is a trivial
   wrapper around
   \bprog
   ZM_lll(ZM_lll(x, 0.99, LLL_KER), 0.99, LLL_INPLACE);
   @eprog\noindent Remove the outermost \kbd{ZM\_lll} if LLL-reduction is not
   desired (saves time).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.matkerint0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qflll(*argv):
  '''
  qflll
  Class: basic
  Section: linear_algebra
  C-Name: qflll0
  Prototype: GD0,L,
  Help: qflll(x,{flag=0}): LLL reduction of the vectors forming the matrix x
   (gives the unimodular transformation matrix T such that x*T is LLL-reduced). flag is
   optional, and can be 0: default, 1: assumes x is integral, 2: assumes x is
   integral, returns a partially reduced basis,
   4: assumes x is integral, returns [K,T] where K is the integer kernel of x
   and T the LLL reduced image, 5: same as 4 but x may have polynomial
   coefficients, 8: same as 0 but x may have polynomial coefficients.
  Description: 
   (vec, ?0):vec       lll($1)
   (vec, 1):vec        lllint($1)
   (vec, 2):vec        lllintpartial($1)
   (vec, 4):vec        lllkerim($1)
   (vec, 5):vec        lllkerimgen($1)
   (vec, 8):vec        lllgen($1)
   (vec, #small):vec   $"Bad flag in qflll"
   (vec, small):vec    qflll0($1, $2)
  Doc: \idx{LLL} algorithm applied to the
   \emph{columns} of the matrix $x$. The columns of $x$ may be linearly
   dependent. The result is a unimodular transformation matrix $T$ such that $x
   \cdot T$ is an LLL-reduced basis of the lattice generated by the column
   vectors of $x$. Note that if $x$ is not of maximal rank $T$ will not be
   square. The LLL parameters are $(0.51,0.99)$, meaning that the Gram-Schmidt
   coefficients for the final basis satisfy $\mu_{i,j} \leq |0.51|$, and the
   Lov\'{a}sz's constant is $0.99$.
   
   If $\fl=0$ (default), assume that $x$ has either exact (integral or
   rational) or real floating point entries. The matrix is rescaled, converted
   to integers and the behavior is then as in $\fl = 1$.
   
   If $\fl=1$, assume that $x$ is integral. Computations involving Gram-Schmidt
   vectors are approximate, with precision varying as needed (Lehmer's trick,
   as generalized by Schnorr). Adapted from Nguyen and Stehl\'e's algorithm
   and Stehl\'e's code (\kbd{fplll-1.3}).
   
   If $\fl=2$, $x$ should be an integer matrix whose columns are linearly
   independent. Returns a partially reduced basis for $x$, using an unpublished
   algorithm by Peter Montgomery: a basis is said to be \emph{partially reduced}
   if $|v_i \pm v_j| \geq |v_i|$ for any two distinct basis vectors $v_i, \,
   v_j$.
   
   This is faster than $\fl=1$, esp. when one row is huge compared
   to the other rows (knapsack-style), and should quickly produce relatively
   short vectors. The resulting basis is \emph{not} LLL-reduced in general.
   If LLL reduction is eventually desired, avoid this partial reduction:
   applying LLL to the partially reduced matrix is significantly \emph{slower}
   than starting from a knapsack-type lattice.
   
   If $\fl=4$, as $\fl=1$, returning a vector $[K, T]$ of matrices: the
   columns of $K$ represent a basis of the integer kernel of $x$
   (not LLL-reduced in general) and $T$ is the transformation
   matrix such that $x\cdot T$ is an LLL-reduced $\Z$-basis of the image
   of the matrix $x$.
   
   If $\fl=5$, case as case $4$, but $x$ may have polynomial coefficients.
   
   If $\fl=8$, same as case $0$, but $x$ may have polynomial coefficients.
  Variant: Also available are \fun{GEN}{lll}{GEN x} ($\fl=0$),
   \fun{GEN}{lllint}{GEN x} ($\fl=1$), and \fun{GEN}{lllkerim}{GEN x} ($\fl=4$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qflll0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qflllgram(*argv):
  '''
  qflllgram
  Class: basic
  Section: linear_algebra
  C-Name: qflllgram0
  Prototype: GD0,L,
  Help: qflllgram(G,{flag=0}): LLL reduction of the lattice whose gram matrix
   is G (gives the unimodular transformation matrix). flag is optional and can
   be 0: default,1: assumes x is integral, 4: assumes x is integral,
   returns [K,T],  where K is the integer kernel of x
   and T the LLL reduced image, 5: same as 4 but x may have polynomial
   coefficients, 8: same as 0 but x may have polynomial coefficients.
  Doc: same as \kbd{qflll}, except that the
   matrix $G = \kbd{x\til * x}$ is the Gram matrix of some lattice vectors $x$,
   and not the coordinates of the vectors themselves. In particular, $G$ must
   now be a square symmetric real matrix, corresponding to a positive
   quadratic form (not necessarily definite: $x$ needs not have maximal rank).
   The result is a unimodular
   transformation matrix $T$ such that $x \cdot T$ is an LLL-reduced basis of
   the lattice generated by the column vectors of $x$. See \tet{qflll} for
   further details about the LLL implementation.
   
   If $\fl=0$ (default), assume that $G$ has either exact (integral or
   rational) or real floating point entries. The matrix is rescaled, converted
   to integers and the behavior is then as in $\fl = 1$.
   
   If $\fl=1$, assume that $G$ is integral. Computations involving Gram-Schmidt
   vectors are approximate, with precision varying as needed (Lehmer's trick,
   as generalized by Schnorr). Adapted from Nguyen and Stehl\'e's algorithm
   and Stehl\'e's code (\kbd{fplll-1.3}).
   
   $\fl=4$: $G$ has integer entries, gives the kernel and reduced image of $x$.
   
   $\fl=5$: same as $4$, but $G$ may have polynomial coefficients.
  Variant: Also available are \fun{GEN}{lllgram}{GEN G} ($\fl=0$),
   \fun{GEN}{lllgramint}{GEN G} ($\fl=1$), and \fun{GEN}{lllgramkerim}{GEN G}
   ($\fl=4$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qflllgram0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nffactor(*argv):
  '''
  nffactor
  Class: basic
  Section: number_fields
  C-Name: nffactor
  Prototype: GG
  Help: nffactor(nf,T): factor polynomial T in number field nf.
  Doc: factorization of the univariate
   polynomial $T$ over the number field $\var{nf}$ given by \kbd{nfinit}; $T$
   has coefficients in $\var{nf}$ (i.e.~either scalar, polmod, polynomial or
   column vector). The factors are sorted by increasing degree.
   
   The main variable of $\var{nf}$ must be of \emph{lower}
   priority than that of $T$, see \secref{se:priority}. However if
   the polynomial defining the number field occurs explicitly  in the
   coefficients of $T$ as modulus of a \typ{POLMOD} or as a \typ{POL}
   coefficient, its main variable must be \emph{the same} as the main variable
   of $T$. For example,
   \bprog
   ? nf = nfinit(y^2 + 1);
   ? nffactor(nf, x^2 + y); \\@com OK
   ? nffactor(nf, x^2 + Mod(y, y^2+1)); \\ @com OK
   ? nffactor(nf, x^2 + Mod(z, z^2+1)); \\ @com WRONG
   @eprog
   
   It is possible to input a defining polynomial for \var{nf}
   instead, but this is in general less efficient since parts of an \kbd{nf}
   structure will then be computed internally. This is useful in two
   situations: when you do not need the \kbd{nf} elsewhere, or when you cannot
   compute the field discriminant due to integer factorization difficulties. In
   the latter case, if you must use a partial discriminant factorization (as
   allowed by both \tet{nfdisc} or \tet{nfbasis}) to build a partially correct
   \var{nf} structure, always input \kbd{nf.pol} to \kbd{nffactor}, and not your
   makeshift \var{nf}: otherwise factors could be missed.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nffactor(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nffactormod(*argv):
  '''
  nffactormod
  Class: basic
  Section: number_fields
  C-Name: nffactormod
  Prototype: GGG
  Help: nffactormod(nf,Q,pr): factor polynomial Q modulo prime ideal pr
   in number field nf.
  Doc: factors the univariate polynomial $Q$ modulo the prime ideal \var{pr} in
   the number field $\var{nf}$. The coefficients of $Q$ belong to the number
   field (scalar, polmod, polynomial, even column vector) and the main variable
   of $\var{nf}$ must be of lower priority than that of $Q$ (see
   \secref{se:priority}). The prime ideal \var{pr} is either in
   \tet{idealprimedec} or (preferred) \tet{modprinit} format. The coefficients
   of the polynomial factors are lifted to elements of \var{nf}:
   \bprog
   ? K = nfinit(y^2+1);
   ? P = idealprimedec(K, 3)[1];
   ? nffactormod(K, x^2 + y*x + 18*y+1, P)
   %3 =
   [x + (2*y + 1) 1]
   
   [x + (2*y + 2) 1]
   ? P = nfmodprinit(K, P);  \\ convert to nfmodprinit format
   ? nffactormod(K, x^2 + y*x + 18*y+1)
   [x + (2*y + 1) 1]
   
   [x + (2*y + 2) 1]
   @eprog\noindent Same result, of course, here about 10\% faster due to the
   precomputation.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nffactormod(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfroots(*argv):
  '''
  nfroots
  Class: basic
  Section: number_fields
  C-Name: nfroots
  Prototype: DGG
  Help: nfroots({nf},x): roots of polynomial x belonging to nf (Q if
   omitted) without multiplicity.
  Doc: roots of the polynomial $x$ in the
   number field $\var{nf}$ given by \kbd{nfinit} without multiplicity (in $\Q$
   if $\var{nf}$ is omitted). $x$ has coefficients in the number field (scalar,
   polmod, polynomial, column vector). The main variable of $\var{nf}$ must be
   of lower priority than that of $x$ (see \secref{se:priority}). However if the
   coefficients of the number field occur explicitly (as polmods) as
   coefficients of $x$, the variable of these polmods \emph{must} be the same as
   the main variable of $t$ (see \kbd{nffactor}).
   
   It is possible to input a defining polynomial for \var{nf}
   instead, but this is in general less efficient since parts of an \kbd{nf}
   structure will be computed internally. This is useful in two situations: when
   you don't need the \kbd{nf}, or when you can't compute its discriminant due
   to integer factorization difficulties. In the latter case, \tet{addprimes} is
   a possibility but a dangerous one: roots will probably be missed if the
   (true) field discriminant and an \kbd{addprimes} entry are strictly divisible
   by some prime. If you have such an unsafe \var{nf}, it is safer to input
   \kbd{nf.pol}.
  Variant: See also \fun{GEN}{nfrootsQ}{GEN x},
   corresponding to $\kbd{nf} = \kbd{NULL}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfroots(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factornf(*argv):
  '''
  factornf
  Class: basic
  Section: number_fields
  C-Name: polfnf
  Prototype: GG
  Help: factornf(x,t): factorization of the polynomial x over the number field
   defined by the polynomial t.
  Doc: factorization of the univariate polynomial $x$
   over the number field defined by the (univariate) polynomial $t$. $x$ may
   have coefficients in $\Q$ or in the number field. The algorithm reduces to
   factorization over $\Q$ (\idx{Trager}'s trick). The direct approach of
   \tet{nffactor}, which uses \idx{van Hoeij}'s method in a relative setting, is
   in general faster.
   
   The main variable of $t$ must be of \emph{lower} priority than that of $x$
   (see \secref{se:priority}). However if non-rational number field elements
   occur (as polmods or polynomials) as coefficients of $x$, the variable of
   these polmods \emph{must} be the same as the main variable of $t$. For
   example
   
   \bprog
   ? factornf(x^2 + Mod(y, y^2+1), y^2+1);
   ? factornf(x^2 + y, y^2+1); \\@com these two are OK
   ? factornf(x^2 + Mod(z,z^2+1), y^2+1)
     ***   at top-level: factornf(x^2+Mod(z,z
     ***                 ^--------------------
     *** factornf: inconsistent data in rnf function.
   ? factornf(x^2 + z, y^2+1)
     ***   at top-level: factornf(x^2+z,y^2+1
     ***                 ^--------------------
     *** factornf: incorrect variable in rnf function.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polfnf(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfrootsof1(*argv):
  '''
  nfrootsof1
  Class: basic
  Section: number_fields
  C-Name: rootsof1
  Prototype: G
  Help: nfrootsof1(nf): number of roots of unity and primitive root of unity
   in the number field nf.
  Doc: Returns a two-component vector $[w,z]$ where $w$ is the number of roots of
   unity in the number field \var{nf}, and $z$ is a primitive $w$-th root
   of unity.
   \bprog
   ? K = nfinit(polcyclo(11));
   ? nfrootsof1(K)
   %2 = [22, [0, 0, 0, 0, 0, -1, 0, 0, 0, 0]~]
   ? z = nfbasistoalg(K, %[2])   \\ in algebraic form
   %3 = Mod(-x^5, x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
   ? [lift(z^11), lift(z^2)]     \\ proves that the order of z is 22
   %4 = [-1, -x^9 - x^8 - x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - x - 1]
   @eprog
   This function guesses the number $w$ as the gcd of the $\#k(v)^*$ for
   unramified $v$ above odd primes, then computes the roots in \var{nf}
   of the $w$-th cyclotomic polynomial: the algorithm is polynomial time with
   respect to the field degree and the bitsize of the multiplication table in
   \var{nf} (both of them polynomially bounded in terms of the size of the
   discriminant). Fields of degree up to $100$ or so should require less than
   one minute.
  Variant: Also available is \fun{GEN}{rootsof1_kannan}{GEN nf}, that computes
   all algebraic integers of $T_2$ norm equal to the field degree
   (all roots of $1$, by Kronecker's theorem). This is in general a little
   faster than the default when there \emph{are} roots of $1$ in the field
   (say twice faster), but can be much slower (say, \emph{days} slower), since
   the algorithm is a priori exponential in the field degree.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rootsof1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def numbpart(*argv):
  '''
  numbpart
  Class: basic
  Section: number_theoretical
  C-Name: numbpart
  Prototype: G
  Help: numbpart(n): number of partitions of n.
  Doc: gives the number of unrestricted partitions of
   $n$, usually called $p(n)$ in the literature; in other words the number of
   nonnegative integer solutions to $a+2b+3c+\cdots=n$. $n$ must be of type
   integer and $n<10^{15}$ (with trivial values $p(n) = 0$ for $n < 0$ and
   $p(0) = 1$). The algorithm uses the Hardy-Ramanujan-Rademacher formula.
   To explicitly enumerate them, see \tet{partitions}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.numbpart(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def partitions(*argv):
  '''
  partitions
  Class: basic
  Section: number_theoretical
  C-Name: partitions
  Prototype: LDGDG
  Help: partitions(k,{a=k},{n=k})): vector of partitions of the integer k.
   You can restrict the length of the partitions with parameter n (n=nmax or
   n=[nmin,nmax]), or the range of the parts with parameter a (a=amax
   or a=[amin,amax]). By default remove zeros, but one can set amin=0 to get X of
   fixed length nmax (=k by default).
  Doc: returns the vector of partitions of the integer $k$ as a sum of positive
   integers (parts); for $k < 0$, it returns the empty set \kbd{[]}, and for $k
   = 0$ the trivial partition (no parts). A partition is given by a
   \typ{VECSMALL}, where parts are sorted in nondecreasing order:
   \bprog
   ? partitions(3)
   %1 = [Vecsmall([3]), Vecsmall([1, 2]), Vecsmall([1, 1, 1])]
   @eprog\noindent correspond to $3$, $1+2$ and $1+1+1$. The number
   of (unrestricted) partitions of $k$ is given
   by \tet{numbpart}:
   \bprog
   ? #partitions(50)
   %1 = 204226
   ? numbpart(50)
   %2 = 204226
   @eprog
   
   \noindent Optional parameters $n$ and $a$ are as follows:
   
   \item $n=\var{nmax}$ (resp. $n=[\var{nmin},\var{nmax}]$) restricts
   partitions to length less than $\var{nmax}$ (resp. length between
   $\var{nmin}$ and $nmax$), where the \emph{length} is the number of nonzero
   entries.
   
   \item $a=\var{amax}$ (resp. $a=[\var{amin},\var{amax}]$) restricts the parts
   to integers less than $\var{amax}$ (resp. between $\var{amin}$ and
   $\var{amax}$).
   \bprog
   ? partitions(4, 2)  \\ parts bounded by 2
   %1 = [Vecsmall([2, 2]), Vecsmall([1, 1, 2]), Vecsmall([1, 1, 1, 1])]
   ? partitions(4,, 2) \\ at most 2 parts
   %2 = [Vecsmall([4]), Vecsmall([1, 3]), Vecsmall([2, 2])]
   ? partitions(4,[0,3], 2) \\ at most 2 parts
   %3 = [Vecsmall([4]), Vecsmall([1, 3]), Vecsmall([2, 2])]
   @eprog\noindent
   By default, parts are positive and we remove zero entries unless
   $amin\leq0$, in which case $nmin$ is ignored and $X$ is of constant length
   $\var{nmax}$:
   \bprog
   ? partitions(4, [0,3])  \\ parts between 0 and 3
   %1 = [Vecsmall([0, 0, 1, 3]), Vecsmall([0, 0, 2, 2]),\
         Vecsmall([0, 1, 1, 2]), Vecsmall([1, 1, 1, 1])]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.partitions(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factorpadic(*argv):
  '''
  factorpadic
  Class: basic
  Section: polynomials
  C-Name: factorpadic0
  Prototype: GGLD0,L,
  Help: factorpadic(pol,p,r): p-adic factorization of the polynomial pol
   to precision r.
  Doc: $p$-adic factorization
   of the polynomial \var{pol} to precision $r$, the result being a
   two-column matrix as in \kbd{factor}. Note that this is not the same
   as a factorization over $\Z/p^r\Z$ (polynomials over that ring do not form a
   unique factorization domain, anyway), but approximations in $\Q/p^r\Z$ of
   the true factorization in $\Q_p[X]$.
   \bprog
   ? factorpadic(x^2 + 9, 3,5)
   %1 =
   [(1 + O(3^5))*x^2 + O(3^5)*x + (3^2 + O(3^5)) 1]
   ? factorpadic(x^2 + 1, 5,3)
   %2 =
   [  (1 + O(5^3))*x + (2 + 5 + 2*5^2 + O(5^3)) 1]
   
   [(1 + O(5^3))*x + (3 + 3*5 + 2*5^2 + O(5^3)) 1]
   @eprog\noindent
   The factors are normalized so that their leading coefficient is a power of
   $p$. The method used is a modified version of the \idx{round 4} algorithm of
   \idx{Zassenhaus}.
   
   If \var{pol} has inexact \typ{PADIC} coefficients, this is not always
   well-defined; in this case, the polynomial is first made integral by dividing
   out the $p$-adic content,  then lifted to $\Z$ using \tet{truncate}
   coefficientwise.
   Hence we actually factor exactly a polynomial which is only $p$-adically
   close to the input. To avoid pitfalls, we advise to only factor polynomials
   with exact rational coefficients.
   
   \synt{factorpadic}{GEN f,GEN p, long r} . The function \kbd{factorpadic0} is
   deprecated, provided for backward compatibility.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factorpadic0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def padicappr(*argv):
  '''
  padicappr
  Class: basic
  Section: polynomials
  C-Name: padicappr
  Prototype: GG
  Help: padicappr(pol,a): p-adic roots of the polynomial pol congruent to a mod p.
  Doc: vector of $p$-adic roots of the
   polynomial $pol$ congruent to the $p$-adic number $a$ modulo $p$, and with
   the same $p$-adic precision as $a$. The number $a$ can be an ordinary
   $p$-adic number (type \typ{PADIC}, i.e.~an element of $\Z_p$) or can be an
   integral element of a finite extension of $\Q_p$, given as a \typ{POLMOD}
   at least one of whose coefficients is a \typ{PADIC}. In this case, the result
   is the vector of roots belonging to the same extension of $\Q_p$ as $a$.
  Variant: Also available is \fun{GEN}{Zp_appr}{GEN f, GEN a} when $a$ is a
   \typ{PADIC}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.padicappr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polrootspadic(*argv):
  '''
  polrootspadic
  Class: basic
  Section: polynomials
  C-Name: rootpadic
  Prototype: GGL
  Help: polrootspadic(x,p,r): p-adic roots of the polynomial x to precision r.
  Doc: vector of $p$-adic roots of the polynomial \var{pol}, given to
   $p$-adic precision $r$ $p$ is assumed to be a prime. Multiple roots are
   \emph{not} repeated. Note that this is not the same as the roots in
   $\Z/p^r\Z$, rather it gives approximations in $\Z/p^r\Z$ of the true roots
   living in $\Q_p$.
   \bprog
   ? polrootspadic(x^3 - x^2 + 64, 2, 5)
   %1 = [2^3 + O(2^5), 2^3 + 2^4 + O(2^5), 1 + O(2^5)]~
   @eprog
   If \var{pol} has inexact \typ{PADIC} coefficients, this is not always
   well-defined; in this case, the polynomial is first made integral by dividing
   out the $p$-adic content, then lifted
   to $\Z$ using \tet{truncate} coefficientwise. Hence the roots given are
   approximations of the roots of an exact polynomial which is $p$-adically
   close to the input. To avoid pitfalls, we advise to only factor polynomials
   with eact rational coefficients.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rootpadic(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def content(*argv):
  '''
  content
  Class: basic
  Section: number_theoretical
  C-Name: content
  Prototype: G
  Help: content(x): gcd of all the components of x, when this makes sense.
  Doc: computes the gcd of all the coefficients of $x$,
   when this gcd makes sense. This is the natural definition
   if $x$ is a polynomial (and by extension a power series) or a
   vector/matrix. This is in general a weaker notion than the \emph{ideal}
   generated by the coefficients:
   \bprog
   ? content(2*x+y)
   %1 = 1            \\ = gcd(2,y) over Q[y]
   @eprog
   
   If $x$ is a scalar, this simply returns the absolute value of $x$ if $x$ is
   rational (\typ{INT} or \typ{FRAC}), and either $1$ (inexact input) or $x$
   (exact input) otherwise; the result should be identical to \kbd{gcd(x, 0)}.
   
   The content of a rational function is the ratio of the contents of the
   numerator and the denominator. In recursive structures, if a
   matrix or vector \emph{coefficient} $x$ appears, the gcd is taken
   not with $x$, but with its content:
   \bprog
   ? content([ [2], 4*matid(3) ])
   %1 = 2
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.content(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factorback(*argv):
  '''
  factorback
  Class: basic
  Section: number_theoretical
  C-Name: factorback2
  Prototype: GDG
  Help: factorback(f,{e}): given a factorisation f, gives the factored
   object back. If this is a prime ideal factorisation you must supply the
   corresponding number field as last argument. If e is present, f has to be a
   vector of the same length, and we return the product of the f[i]^e[i].
  Description: 
   (gen):gen      factorback($1)
   (gen,):gen     factorback($1)
   (gen,gen):gen  factorback2($1, $2)
  Doc: gives back the factored object
   corresponding to a factorization. The integer $1$ corresponds to the empty
   factorization.
   
   If $e$ is present, $e$ and $f$ must be vectors of the same length ($e$ being
   integral), and the corresponding factorization is the product of the
   $f[i]^{e[i]}$.
   
   If not, and $f$ is vector, it is understood as in the preceding case with $e$
   a vector of 1s: we return the product of the $f[i]$. Finally, $f$ can be a
   regular factorization, as produced with any \kbd{factor} command. A few
   examples:
   \bprog
   ? factor(12)
   %1 =
   [2 2]
   
   [3 1]
   
   ? factorback(%)
   %2 = 12
   ? factorback([2,3], [2,1])   \\ 2^3 * 3^1
   %3 = 12
   ? factorback([5,2,3])
   %4 = 30
   @eprog
  Variant: Also available is \fun{GEN}{factorback}{GEN f} (case $e = \kbd{NULL}$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.factorback2(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def gcd(*argv):
  '''
  gcd
  Class: basic
  Section: number_theoretical
  C-Name: ggcd0
  Prototype: GDG
  Help: gcd(x,{y}): greatest common divisor of x and y.
  Description: 
   (small, small):small   cgcd($1, $2)
   (int, int):int         gcdii($1, $2)
   (gen):gen              content($1)
   (gen, gen):gen         ggcd($1, $2)
  Doc: creates the greatest common divisor of $x$ and $y$.
   If you also need the $u$ and $v$ such that $x*u + y*v = \gcd(x,y)$,
   use the \tet{bezout} function. $x$ and $y$ can have rather quite general
   types, for instance both rational numbers. If $y$ is omitted and $x$ is a
   vector, returns the $\text{gcd}$ of all components of $x$, i.e.~this is
   equivalent to \kbd{content(x)}.
   
   When $x$ and $y$ are both given and one of them is a vector/matrix type,
   the GCD is again taken recursively on each component, but in a different way.
   If $y$ is a vector, resp.~matrix, then the result has the same type as $y$,
   and components equal to \kbd{gcd(x, y[i])}, resp.~\kbd{gcd(x, y[,i])}. Else
   if $x$ is a vector/matrix the result has the same type as $x$ and an
   analogous definition. Note that for these types, \kbd{gcd} is not
   commutative.
   
   The algorithm used is a naive \idx{Euclid} except for the following inputs:
   
   \item integers: use modified right-shift binary (``plus-minus''
   variant).
   
   \item univariate polynomials with coefficients in the same number
   field (in particular rational): use modular gcd algorithm.
   
   \item general polynomials: use the \idx{subresultant algorithm} if
   coefficient explosion is likely (non modular coefficients).
   
   If $u$ and $v$ are polynomials in the same variable with \emph{inexact}
   coefficients, their gcd is defined to be scalar, so that
   \bprog
   ? a = x + 0.0; gcd(a,a)
   %1 = 1
   ? b = y*x + O(y); gcd(b,b)
   %2 = y
   ? c = 4*x + O(2^3); gcd(c,c)
   %2 = 4
   @eprog\noindent A good quantitative check to decide whether such a
   gcd ``should be'' non-trivial, is to use \tet{polresultant}: a value
   close to $0$ means that a small deformation of the inputs has non-trivial gcd.
   You may also use \tet{bezout}, which does try to compute an approximate gcd
   $d$ and provides $u$, $v$ to check whether $u x + v y$ is close to $d$.
  Variant: Also available are \fun{GEN}{ggcd}{GEN x, GEN y}, if \kbd{y} is not
   \kbd{NULL}, and \fun{GEN}{content}{GEN x}, if $\kbd{y} = \kbd{NULL}$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ggcd0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def lcm(*argv):
  '''
  lcm
  Class: basic
  Section: number_theoretical
  C-Name: glcm0
  Prototype: GDG
  Help: lcm(x,{y}): least common multiple of x and y, i.e. x*y / gcd(x,y).
  Description: 
   (int, int):int lcmii($1, $2)
   (gen):gen      glcm0($1, NULL)
   (gen, gen):gen glcm($1, $2)
  Doc: least common multiple of $x$ and $y$, i.e.~such
   that $\lcm(x,y)*\gcd(x,y) = \text{abs}(x*y)$. If $y$ is omitted and $x$
   is a vector, returns the $\text{lcm}$ of all components of $x$.
   
   When $x$ and $y$ are both given and one of them is a vector/matrix type,
   the LCM is again taken recursively on each component, but in a different way.
   If $y$ is a vector, resp.~matrix, then the result has the same type as $y$,
   and components equal to \kbd{lcm(x, y[i])}, resp.~\kbd{lcm(x, y[,i])}. Else
   if $x$ is a vector/matrix the result has the same type as $x$ and an
   analogous definition. Note that for these types, \kbd{lcm} is not
   commutative.
   
   Note that \kbd{lcm(v)} is quite different from
   \bprog
   l = v[1]; for (i = 1, #v, l = lcm(l, v[i]))
   @eprog\noindent
   Indeed, \kbd{lcm(v)} is a scalar, but \kbd{l} may not be (if one of
   the \kbd{v[i]} is a vector/matrix). The computation uses a divide-conquer tree
   and should be much more efficient, especially when using the GMP
   multiprecision kernel (and more subquadratic algorithms become available):
   \bprog
   ? v = vector(10^4, i, random);
   ? lcm(v);
   time = 323 ms.
   ? l = v[1]; for (i = 1, #v, l = lcm(l, v[i]))
   time = 833 ms.
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.glcm0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factor(*argv):
  '''
  factor
  Class: basic
  Section: number_theoretical
  C-Name: gp_factor0
  Prototype: GDG
  Help: factor(x,{lim}): factorization of x. lim is optional and can be set
   whenever x is of (possibly recursive) rational type. If lim is set return
   partial factorization, using primes < lim.
  Description: 
   (int, ?-1):vec        Z_factor($1)
   (gen, ?-1):vec        factor($1)
   (gen, small):vec      factor0($1, $2)
  Doc: general factorization function, where $x$ is a
   rational (including integers), a complex number with rational
   real and imaginary parts, or a rational function (including polynomials).
   The result is a two-column matrix: the first contains the irreducibles
   dividing $x$ (rational or Gaussian primes, irreducible polynomials),
   and the second the exponents. By convention, $0$ is factored as $0^1$.
   
   \misctitle{$\Q$ and $\Q(i)$}
   See \tet{factorint} for more information about the algorithms used.
   The rational or Gaussian primes are in fact \var{pseudoprimes}
   (see \kbd{ispseudoprime}), a priori not rigorously proven primes. In fact,
   any factor which is $\leq 10^{15}$ (whose norm is $\leq 10^{15}$ for an
   irrational Gaussian prime) is a genuine prime. Use
   \kbd{isprime} to prove primality of other factors, as in
   \bprog
   ? fa = factor(2^2^7 + 1)
   %1 =
   [59649589127497217 1]
   
   [5704689200685129054721 1]
   
   ? isprime( fa[,1] )
   %2 = [1, 1]~   \\ both entries are proven primes
   @eprog\noindent
   Another possibility is to set the global default \tet{factor_proven}, which
   will perform a rigorous primality proof for each pseudoprime factor.
   
   A \typ{INT} argument \var{lim} can be added, meaning that we look only for
   prime factors $p < \var{lim}$. The limit \var{lim} must be non-negative.
   In this case, all but the last factor are proven primes, but the remaining
   factor may actually be a proven composite! If the remaining factor is less
   than $\var{lim}^2$, then it is prime.
   \bprog
   ? factor(2^2^7 +1, 10^5)
   %3 =
   [340282366920938463463374607431768211457 1]
   @eprog\noindent
   \misctitle{Deprecated feature} Setting $\var{lim}=0$ is the same
   as setting it to $\kbd{primelimit} + 1$. Don't use this: it is unwise to
   rely on global variables when you can specify an explicit argument.
   \smallskip
   
   This routine uses trial division and perfect power tests, and should not be
   used for huge values of \var{lim} (at most $10^9$, say):
   \kbd{factorint(, 1 + 8)} will in general be faster. The latter does not
   guarantee that all small
   prime factors are found, but it also finds larger factors, and in a much more
   efficient way.
   \bprog
   ? F = (2^2^7 + 1) * 1009 * 100003; factor(F, 10^5)  \\ fast, incomplete
   time = 0 ms.
   %4 =
   [1009 1]
   
   [34029257539194609161727850866999116450334371 1]
   
   ? factor(F, 10^9)    \\ very slow
   time = 6,892 ms.
   %6 =
   [1009 1]
   
   [100003 1]
   
   [340282366920938463463374607431768211457 1]
   
   ? factorint(F, 1+8)  \\ much faster, all small primes were found
   time = 12 ms.
   %7 =
   [1009 1]
   
   [100003 1]
   
   [340282366920938463463374607431768211457 1]
   
   ? factor(F)   \\ complete factorisation
   time = 112 ms.
   %8 =
   [1009 1]
   
   [100003 1]
   
   [59649589127497217 1]
   
   [5704689200685129054721 1]
   @eprog\noindent Over $\Q$, the prime factors are sorted by increasing size.
   
   \misctitle{Rational functions}
   The polynomials or rational functions to be factored must have scalar
   coefficients. In particular PARI does not know how to factor
   \emph{multivariate} polynomials. The following domains are currently
   supported: $\Q$, $\R$, $\C$, $\Q_p$, finite fields and number fields.
   See \tet{factormod} and \tet{factorff} for
   the algorithms used over finite fields, \tet{factornf} for the algorithms
   over number fields. Over $\Q$, \idx{van Hoeij}'s method is used, which is
   able to cope with hundreds of modular factors.
   
   The routine guesses a sensible ring over which to factor: the
   smallest ring containing all coefficients, taking into account quotient
   structures induced by \typ{INTMOD}s and \typ{POLMOD}s (e.g.~if a coefficient
   in $\Z/n\Z$ is known, all rational numbers encountered are first mapped to
   $\Z/n\Z$; different moduli will produce an error). Factoring modulo a
   non-prime number is not supported; to factor in $\Q_p$, use \typ{PADIC}
   coefficients not \typ{INTMOD} modulo $p^n$.
   \bprog
   ? T = x^2+1;
   ? factor(T);                         \\ over Q
   ? factor(T*Mod(1,3))                 \\ over F_3
   ? factor(T*ffgen(ffinit(3,2,'t))^0)  \\ over F_{3^2}
   ? factor(T*Mod(Mod(1,3), t^2+t+2))   \\ over F_{3^2}, again
   ? factor(T*(1 + O(3^6))              \\ over Q_3, precision 6
   ? factor(T*1.)                       \\ over R, current precision
   ? factor(T*(1.+0.*I))                \\ over C
   ? factor(T*Mod(1, y^3-2))            \\ over Q(2^{1/3})
   @eprog\noindent In most cases, it is clearer and simpler to call an
   explicit variant than to rely on the generic \kbd{factor} function and
   the above detection mechanism:
   \bprog
   ? factormod(T, 3)           \\ over F_3
   ? factorff(T, 3, t^2+t+2))  \\ over F_{3^2}
   ? factorpadic(T, 3,6)       \\ over Q_3, precision 6
   ? nffactor(y^3-2, T)        \\ over Q(2^{1/3})
   ? polroots(T)               \\ over C
   @eprog
   
   Note that factorization of polynomials is done up to
   multiplication by a constant. In particular, the factors of rational
   polynomials will have integer coefficients, and the content of a polynomial
   or rational function is discarded and not included in the factorization. If
   needed, you can always ask for the content explicitly:
   \bprog
   ? factor(t^2 + 5/2*t + 1)
   %1 =
   [2*t + 1 1]
   
   [t + 2 1]
   
   ? content(t^2 + 5/2*t + 1)
   %2 = 1/2
   @eprog\noindent
   The irreducible factors are sorted by increasing degree.
   See also \tet{nffactor}.
  Variant: This function should only be used by the \kbd{gp} interface. Use
   directly \fun{GEN}{factor}{GEN x} or \fun{GEN}{boundfact}{GEN x, ulong lim}.
   The obsolete function \fun{GEN}{factor0}{GEN x, long lim} is kept for
   backward compatibility.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gp_factor0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def idealfactorback(*argv):
  '''
  idealfactorback
  Class: basic
  Section: number_fields
  C-Name: idealfactorback
  Prototype: GGDGD0,L,
  Help: idealfactorback(nf,f,{e},{flag = 0}): given a factorisation f, gives the
   ideal product back. If e is present, f has to be a
   vector of the same length, and we return the product of the f[i]^e[i]. If
   flag is non-zero, perform idealred along the way.
  Doc: gives back the ideal corresponding to a factorization. The integer $1$
   corresponds to the empty factorization.
   If $e$ is present, $e$ and $f$ must be vectors of the same length ($e$ being
   integral), and the corresponding factorization is the product of the
   $f[i]^{e[i]}$.
   
   If not, and $f$ is vector, it is understood as in the preceding case with $e$
   a vector of 1s: we return the product of the $f[i]$. Finally, $f$ can be a
   regular factorization, as produced by \kbd{idealfactor}.
   \bprog
   ? nf = nfinit(y^2+1); idealfactor(nf, 4 + 2*y)
   %1 =
   [[2, [1, 1]~, 2, 1, [1, 1]~] 2]
   
   [[5, [2, 1]~, 1, 1, [-2, 1]~] 1]
   
   ? idealfactorback(nf, %)
   %2 =
   [10 4]
   
   [0  2]
   
   ? f = %1[,1]; e = %1[,2]; idealfactorback(nf, f, e)
   %3 =
   [10 4]
   
   [0  2]
   
   ? % == idealhnf(nf, 4 + 2*y)
   %4 = 1
   @eprog
   If \kbd{flag} is non-zero, perform ideal reductions (\tet{idealred}) along the
   way. This is most useful if the ideals involved are all \emph{extended}
   ideals (for instance with trivial principal part), so that the principal parts
   extracted by \kbd{idealred} are not lost. Here is an example:
   \bprog
   ? f = vector(#f, i, [f[i], [;]]);  \\ transform to extended ideals
   ? idealfactorback(nf, f, e, 1)
   %6 = [[1, 0; 0, 1], [2, 1; [2, 1]~, 1]]
   ? nffactorback(nf, %[2])
   %7 = [4, 2]~
   @eprog
   The extended ideal returned in \kbd{\%6} is the trivial ideal $1$, extended
   with a principal generator given in factored form. We use \tet{nffactorback}
   to recover it in standard form.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.idealfactorback(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polisirreducible(*argv):
  '''
  polisirreducible
  Class: basic
  Section: polynomials
  C-Name: isirreducible
  Prototype: lG
  Help: polisirreducible(pol): true(1) if pol is an irreducible non-constant
   polynomial, false(0) if pol is reducible or constant.
  Doc: \var{pol} being a polynomial (univariate in the present version \vers),
   returns 1 if \var{pol} is non-constant and irreducible, 0 otherwise.
   Irreducibility is checked over the smallest base field over which \var{pol}
   seems to be defined.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.isirreducible(*c_arg_tuple)

def newtonpoly(*argv):
  '''
  newtonpoly
  Class: basic
  Section: number_fields
  C-Name: newtonpoly
  Prototype: GG
  Help: newtonpoly(x,p): Newton polygon of polynomial x with respect to the
   prime p.
  Doc: gives the vector of the slopes of the Newton
   polygon of the polynomial $x$ with respect to the prime number $p$. The $n$
   components of the vector are in decreasing order, where $n$ is equal to the
   degree of $x$. Vertical slopes occur iff the constant coefficient of $x$ is
   zero and are denoted by \tet{LONG_MAX}, the biggest single precision
   integer representable on the machine ($2^{31}-1$ (resp.~$2^{63}-1$) on 32-bit
   (resp.~64-bit) machines), see \secref{se:valuation}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.newtonpoly(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nffactorback(*argv):
  '''
  nffactorback
  Class: basic
  Section: number_fields
  C-Name: nffactorback
  Prototype: GGDG
  Help: nffactorback(nf,f,{e}): given a factorisation f, returns
   the factored object back as an nf element.
  Doc: gives back the \kbd{nf} element corresponding to a factorization.
   The integer $1$ corresponds to the empty factorization.
   
   If $e$ is present, $e$ and $f$ must be vectors of the same length ($e$ being
   integral), and the corresponding factorization is the product of the
   $f[i]^{e[i]}$.
   
   If not, and $f$ is vector, it is understood as in the preceding case with $e$
   a vector of 1s: we return the product of the $f[i]$. Finally, $f$ can be a
   regular factorization matrix.
   \bprog
   ? nf = nfinit(y^2+1);
   ? nffactorback(nf, [3, y+1, [1,2]~], [1, 2, 3])
   %2 = [12, -66]~
   ? 3 * (I+1)^2 * (1+2*I)^3
   %3 = 12 - 66*I
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nffactorback(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poldisc(*argv):
  '''
  poldisc
  Class: basic
  Section: polynomials
  C-Name: poldisc0
  Prototype: GDn
  Help: poldisc(pol,{v}): discriminant of the polynomial pol, with respect to main
   variable if v is omitted, with respect to v otherwise.
  Description: 
   (pol):gen        discsr($1)
   (gen):gen        poldisc0($1, -1)
   (gen, var):gen   poldisc0($1, $2)
  Doc: discriminant of the polynomial
   \var{pol} in the main variable if $v$ is omitted, in $v$ otherwise. The
   algorithm used is the \idx{subresultant algorithm}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.poldisc0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polresultant(*argv):
  '''
  polresultant
  Class: basic
  Section: polynomials
  C-Name: polresultant0
  Prototype: GGDnD0,L,
  Help: polresultant(x,y,{v},{flag=0}): resultant of the polynomials x and y,
   with respect to the main variables of x and y if v is omitted, with respect
   to the variable v otherwise. flag is optional, and can be 0: default,
   uses either the subresultant algorithm, a modular algorithm or Sylvester's
   matrix, depending on the inputs; 1 uses Sylvester's matrix (should always be
   slower than the default).
  Doc: resultant of the two
   polynomials $x$ and $y$ with exact entries, with respect to the main
   variables of $x$ and $y$ if $v$ is omitted, with respect to the variable $v$
   otherwise. The algorithm assumes the base ring is a domain. If you also need
   the $u$ and $v$ such that $x*u + y*v = \text{Res}(x,y)$, use the
   \tet{polresultantext} function.
   
   If $\fl=0$ (default), uses the the algorithm best suited to the inputs,
   either the \idx{subresultant algorithm} (Lazard/Ducos variant, generic case),
   a modular algorithm (inputs in $\Q[X]$) or Sylvester's matrix (inexact
   inputs).
   
   If $\fl=1$, uses the determinant of Sylvester's matrix instead; this should
   always be slower than the default.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polresultant0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polsym(*argv):
  '''
  polsym
  Class: basic
  Section: polynomials
  C-Name: polsym
  Prototype: GL
  Help: polsym(x,n): column vector of symmetric powers of the roots of x up to n.
  Doc: creates the column vector of the \idx{symmetric powers} of the roots of the
   polynomial $x$ up to power $n$, using Newton's formula.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polsym(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def poldiscreduced(*argv):
  '''
  poldiscreduced
  Class: basic
  Section: polynomials
  C-Name: reduceddiscsmith
  Prototype: G
  Help: poldiscreduced(f): vector of elementary divisors of Z[a]/f'(a)Z[a],
   where a is a root of the polynomial f.
  Doc: reduced discriminant vector of the
   (integral, monic) polynomial $f$. This is the vector of elementary divisors
   of $\Z[\alpha]/f'(\alpha)\Z[\alpha]$, where $\alpha$ is a root of the
   polynomial $f$. The components of the result are all positive, and their
   product is equal to the absolute value of the discriminant of~$f$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.reduceddiscsmith(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def rnfcharpoly(*argv):
  '''
  rnfcharpoly
  Class: basic
  Section: number_fields
  C-Name: rnfcharpoly
  Prototype: GGGDn
  Help: rnfcharpoly(nf,T,a,{var='x}): characteristic polynomial of a
   over nf, where a belongs to the algebra defined by T over nf. Returns a
   polynomial in variable var (x by default).
  Doc: characteristic polynomial of
   $a$ over $\var{nf}$, where $a$ belongs to the algebra defined by $T$ over
   $\var{nf}$, i.e.~$\var{nf}[X]/(T)$. Returns a polynomial in variable $v$
   ($x$ by default).
   \bprog
   ? nf = nfinit(y^2+1);
   ? rnfcharpoly(nf, x^2+y*x+1, x+y)
   %2 = x^2 + Mod(-y, y^2 + 1)*x + 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.rnfcharpoly(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polsturm(*argv):
  '''
  polsturm
  Class: basic
  Section: polynomials
  C-Name: sturmpart
  Prototype: lGDGDG
  Help: polsturm(pol,{a},{b}): number of real roots of the squarefree polynomial
   pol in the interval ]a,b] (which are respectively taken to be -oo or +oo when
   omitted).
  Doc: number of real roots of the real squarefree polynomial \var{pol} in the
   interval $]a,b]$, using Sturm's algorithm. $a$ (resp.~$b$) is taken to be
   $-\infty$ (resp.~$+\infty$) if omitted.
  Variant: Also available is \fun{long}{sturm}{GEN pol} (total number of real
   roots).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  return PyPari.pari.sturmpart(*c_arg_tuple)

def polsylvestermatrix(*argv):
  '''
  polsylvestermatrix
  Class: basic
  Section: polynomials
  C-Name: sylvestermatrix
  Prototype: GG
  Help: polsylvestermatrix(x,y): forms the sylvester matrix associated to the
   two polynomials x and y. Warning: the polynomial coefficients are in
   columns, not in rows.
  Doc: forms the Sylvester matrix
   corresponding to the two polynomials $x$ and $y$, where the coefficients of
   the polynomials are put in the columns of the matrix (which is the natural
   direction for solving equations afterwards). The use of this matrix can be
   essential when dealing with polynomials with inexact entries, since
   polynomial Euclidean division doesn't make much sense in this case.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sylvestermatrix(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def gcdext(*argv):
  '''
  gcdext
  Class: basic
  Section: number_theoretical
  C-Name: gcdext0
  Prototype: GG
  Help: gcdext(x,y): returns [u,v,d] such that d=gcd(x,y) and u*x+v*y=d.
  Doc: Returns $[u,v,d]$ such that $d$ is the gcd of $x,y$,
   $x*u+y*v=\gcd(x,y)$, and $u$ and $v$ minimal in a natural sense.
   The arguments must be integers or polynomials. \sidx{extended gcd}
   \sidx{Bezout relation}
   \bprog
   ? [u, v, d] = gcdext(32,102)
   %1 = [16, -5, 2]
   ? d
   %2 = 2
   ? gcdext(x^2-x, x^2+x-2)
   %3 = [-1/2, 1/2, x - 1]
   @eprog
   
   If $x,y$ are polynomials in the same variable and \emph{inexact}
   coefficients, then compute $u,v,d$ such that $x*u+y*v = d$, where $d$
   approximately divides both and $x$ and $y$; in particular, we do not obtain
   \kbd{gcd(x,y)} which is \emph{defined} to be a scalar in this case:
   \bprog
   ? a = x + 0.0; gcd(a,a)
   %1 = 1
   
   ? gcdext(a,a)
   %2 = [0, 1, x + 0.E-28]
   
   ? gcdext(x-Pi, 6*x^2-zeta(2))
   %3 = [-6*x - 18.8495559, 1, 57.5726923]
   @eprog\noindent For inexact inputs, the output is thus not well defined
   mathematically, but you obtain explicit polynomials to check whether the
   approximation is close enough for your needs.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gcdext0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bezoutres(*argv):
  '''
  bezoutres
  Class: basic
  Section: polynomials
  C-Name: polresultantext0
  Prototype: GGDn
  Help: bezoutre(A,B,{v}): deprecated alias for polresultantext
  Doc: deprecated alias for \kbd{polresultantext}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polresultantext0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def characteristic(*argv):
  '''
  characteristic
  Class: basic
  Section: conversions
  C-Name: characteristic
  Prototype: mG
  Help: characteristic(x): characteristic of the base ring over which x is
   defined
  Doc: 
   returns the characteristic of the base ring over which $x$ is defined (as
   defined by \typ{INTMOD} and \typ{FFELT} components). The function raises an
   exception if incompatible primes arise from \typ{FFELT} and \typ{PADIC}
   components.
   \bprog
   ? characteristic(Mod(1,24)*x + Mod(1,18)*y)
   %1 = 6
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.characteristic(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ffinit(*argv):
  '''
  ffinit
  Class: basic
  Section: number_theoretical
  C-Name: ffinit
  Prototype: GLDn
  Help: ffinit(p,n,{v='x}): monic irreducible polynomial of degree n over F_p[v].
  Description: 
   (int, small, ?var):pol        ffinit($1, $2, $3)
  Doc: computes a monic polynomial of degree $n$ which is irreducible over
    $\F_p$, where $p$ is assumed to be prime. This function uses a fast variant
    of Adleman and Lenstra's algorithm.
   
   It is useful in conjunction with \tet{ffgen}; for instance if
   \kbd{P = ffinit(3,2)}, you can represent elements in $\F_{3^2}$ in term of
   \kbd{g = ffgen(P,'t)}. This can be abbreviated as
   \kbd{g = ffgen(3\pow2, 't)}, where the defining polynomial $P$ can be later
   recovered as \kbd{g.mod}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ffinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ffnbirred(*argv):
  '''
  ffnbirred
  Class: basic
  Section: number_theoretical
  C-Name: ffnbirred0
  Prototype: GLD0,L,
  Help: ffnbirred(q,n{,fl=0}): number of monic irreducible polynomials over F_q, of
   degree n (fl=0, default) or at most n (fl=1).
  Description: 
   (int, small, ?0):int      ffnbirred($1, $2)
   (int, small, 1):int       ffsumnbirred($1, $2)
   (int, small, ?small):int  ffnbirred0($1, $2, $3)
  Doc: computes the number of monic irreducible polynomials over $\F_q$ of degree exactly $n$,
   ($\fl=0$ or omited) or at most $n$ ($\fl=1$).
  Variant: Also available are
    \fun{GEN}{ffnbirred}{GEN q, long n} (for $\fl=0$)
    and \fun{GEN}{ffsumnbirred}{GEN q, long n} (for $\fl=1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ffnbirred0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def addprimes(*argv):
  '''
  addprimes
  Class: basic
  Section: number_theoretical
  C-Name: addprimes
  Prototype: DG
  Help: addprimes({x=[]}): add primes in the vector x to the prime table to
   be used in trial division. x may also be a single integer. Composite
   "primes" are NOT allowed!
  Doc: adds the integers contained in the
   vector $x$ (or the single integer $x$) to a special table of
   ``user-defined primes'', and returns that table. Whenever \kbd{factor} is
   subsequently called, it will trial divide by the elements in this table.
   If $x$ is empty or omitted, just returns the current list of extra
   primes.
   
   The entries in $x$ must be primes: there is no internal check, even if
   the \tet{factor_proven} default is set. To remove primes from the list use
   \kbd{removeprimes}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.addprimes(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def isprime(*argv):
  '''
  isprime
  Class: basic
  Section: number_theoretical
  C-Name: gisprime
  Prototype: GD0,L,
  Help: isprime(x,{flag=0}): true(1) if x is a (proven) prime number, false(0)
   if not. If flag is 0 or omitted, use a combination of algorithms. If flag is
   1, the primality is certified by the Pocklington-Lehmer Test. If flag is 2,
   the primality is certified using the APRCL test.
  Description: 
   (int, ?0):bool        isprime($1)
   (int, 1):bool         plisprime($1, 0)
   (int, 2):gen          plisprime($1, 1)
   (gen, ?small):gen     gisprime($1, $2)
  Doc: true (1) if $x$ is a prime
   number, false (0) otherwise. A prime number is a positive integer having
   exactly two distinct divisors among the natural numbers, namely 1 and
   itself.
   
   This routine proves or disproves rigorously that a number is prime, which can
   be very slow when $x$ is indeed prime and has more than $1000$ digits, say.
   Use \tet{ispseudoprime} to quickly check for compositeness. See also
   \kbd{factor}. It accepts vector/matrices arguments, and is then applied
   componentwise.
   
   If $\fl=0$, use a combination of Baillie-PSW pseudo primality test (see
   \tet{ispseudoprime}), Selfridge ``$p-1$'' test if $x-1$ is smooth enough, and
   Adleman-Pomerance-Rumely-Cohen-Lenstra (APRCL) for general $x$.
   
   If $\fl=1$, use Selfridge-Pocklington-Lehmer ``$p-1$'' test and output a
   primality certificate as follows: return
   
   \item 0 if $x$ is composite,
   
   \item 1 if $x$ is small enough that passing Baillie-PSW test guarantees
   its primality (currently $x < 2^{64}$, as checked by Jan Feitsma),
   
   \item $2$ if $x$ is a large prime whose primality could only sensibly be
   proven (given the algorithms implemented in PARI) using the APRCL test.
   
   \item Otherwise ($x$ is large and $x-1$ is smooth) output a three column
   matrix as a primality certificate. The first column contains prime
   divisors $p$ of $x-1$ (such that $\prod p^{v_p(x-1)} > x^{1/3}$), the second
   the corresponding elements $a_p$ as in Proposition~8.3.1 in GTM~138 , and the
   third the output of isprime(p,1).
   
   The algorithm fails if one of the pseudo-prime factors is not prime, which is
   exceedingly unlikely and well worth a bug report. Note that if you monitor
   \kbd{isprime} at a high enough debug level, you may see warnings about
   untested integers being declared primes. This is normal: we ask for partial
   factorisations (sufficient to prove primality if the unfactored part is not
   too large), and \kbd{factor} warns us that the cofactor hasn't been tested.
   It may or may not be tested later, and may or may not be prime. This does
   not affect the validity of the whole \kbd{isprime} procedure.
   
   If $\fl=2$, use APRCL.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gisprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ispseudoprime(*argv):
  '''
  ispseudoprime
  Class: basic
  Section: number_theoretical
  C-Name: gispseudoprime
  Prototype: GD0,L,
  Help: ispseudoprime(x,{flag}): true(1) if x is a strong pseudoprime, false(0)
   if not. If flag is 0 or omitted, use BPSW test, otherwise use strong
   Rabin-Miller test for flag randomly chosen bases.
  Description: 
   (int,?0):bool      BPSW_psp($1)
   (int,#small):bool  millerrabin($1,$2)
   (int,small):bool   ispseudoprime($1, $2)
   (gen,?small):bool  gispseudoprime($1, $2)
  Doc: true (1) if $x$ is a strong pseudo
   prime (see below), false (0) otherwise. If this function returns false, $x$
   is not prime; if, on the other hand it returns true, it is only highly likely
   that $x$ is a prime number. Use \tet{isprime} (which is of course much
   slower) to prove that $x$ is indeed prime.
   The function accepts vector/matrices arguments, and is then applied
   componentwise.
   
   If $\fl = 0$, checks whether $x$ is a Baillie-Pomerance-Selfridge-Wagstaff
   pseudo prime (strong Rabin-Miller pseudo prime for base $2$, followed by
   strong Lucas test for the sequence $(P,-1)$, $P$ smallest positive integer
   such that $P^2 - 4$ is not a square mod $x$).
   
   There are no known composite numbers passing this test, although it is
   expected that infinitely many such numbers exist. In particular, all
   composites $\leq 2^{64}$ are correctly detected (checked using
   \kbd{http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html}).
   
   If $\fl > 0$, checks whether $x$ is a strong Miller-Rabin pseudo prime  for
   $\fl$ randomly chosen bases (with end-matching to catch square roots of $-1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gispseudoprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def prime(*argv):
  '''
  prime
  Class: basic
  Section: number_theoretical
  C-Name: prime
  Prototype: L
  Help: prime(n): returns the n-th prime (n C-integer).
  Doc: the $n^{\text{th}}$ prime number
   \bprog
   ? prime(10^9)
   %1 = 22801763489
   @eprog\noindent Uses checkpointing and a naive $O(n)$ algorithm.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.prime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def primepi(*argv):
  '''
  primepi
  Class: basic
  Section: number_theoretical
  C-Name: primepi
  Prototype: G
  Help: primepi(x): the prime counting function pi(x) = #{p <= x, p prime}.
  Description: 
   (gen):int        primepi($1)
  Doc: the prime counting function. Returns the number of
   primes $p$, $p \leq x$.
   \bprog
   ? primepi(10)
   %1 = 4;
   ? primes(5)
   %2 = [2, 3, 5, 7, 11]
   ? primepi(10^11)
   %3 = 4118054813
   @eprog\noindent Uses checkpointing and a naive $O(x)$ algorithm.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.primepi(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def primes(*argv):
  '''
  primes
  Class: basic
  Section: number_theoretical
  C-Name: primes0
  Prototype: G
  Help: primes(n): returns the vector of the first n primes (integer), or the
   primes in interval n = [a,b].
  Doc: creates a row vector whose components are the first $n$ prime numbers.
   (Returns the empty vector for $n \leq 0$.) A \typ{VEC} $n = [a,b]$ is also
   allowed, in which case the primes in $[a,b]$ are returned
   \bprog
   ? primes(10)     \\ the first 10 primes
   %1 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
   ? primes([0,29])  \\ the primes up to 29
   %2 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
   ? primes([15,30])
   %3 = [17, 19, 23, 29]
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.primes0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def randomprime(*argv):
  '''
  randomprime
  Class: basic
  Section: number_theoretical
  C-Name: randomprime
  Prototype: DG
  Help: randomprime({N = 2^31}): returns a strong pseudo prime in [2, N-1].
  Doc: returns a strong pseudo prime (see \tet{ispseudoprime}) in $[2,N-1]$.
   A \typ{VEC} $N = [a,b]$ is also allowed, with $a \leq b$ in which case a
   pseudo prime $a \leq p \leq b$ is returned; if no prime exists in the
   interval, the function will run into an infinite loop. If the upper bound
   is less than $2^{64}$ the pseudo prime returned is a proven prime.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.randomprime(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def removeprimes(*argv):
  '''
  removeprimes
  Class: basic
  Section: number_theoretical
  C-Name: removeprimes
  Prototype: DG
  Help: removeprimes({x=[]}): remove primes in the vector x from the prime table.
   x can also be a single integer. List the current extra primes if x is omitted.
  Doc: removes the primes listed in $x$ from
   the prime number table. In particular \kbd{removeprimes(addprimes())} empties
   the extra prime table. $x$ can also be a single integer. List the current
   extra primes if $x$ is omitted.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.removeprimes(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfauto(*argv):
  '''
  qfauto
  Class: basic
  Section: linear_algebra
  C-Name: qfauto0
  Prototype: GDG
  Help: qfauto(G,{fl}): automorphism group of the positive definite quadratic form
   G.
  Doc: 
   $G$ being a square and symmetric matrix with integer entries representing a
   positive definite quadratic form, outputs the automorphism group of the
   associate lattice.
   Since this requires computing the minimal vectors, the computations can
   become very lengthy as the dimension grows. $G$ can also be given by an
   \kbd{qfisominit} structure.
   See \kbd{qfisominit} for the meaning of \var{fl}.
   
   The output is a two-components vector $[o,g]$ where $o$ is the group order
   and $g$ is the list of generators (as a vector). For each generators $H$,
   the equality $G={^t}H\*G\*H$ holds.
   
   The interface of this function is experimental and will likely change in the
   future.
   
   This function implements an algorithm of Plesken and Souvignier, following
   Souvignier's implementation.
  Variant: Also available is \fun{GEN}{qfauto}{GEN G, GEN fl}
   where $G$ is a vector of \kbd{zm}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfauto0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfautoexport(*argv):
  '''
  qfautoexport
  Class: basic
  Section: linear_algebra
  C-Name: qfautoexport
  Prototype: GD0,L,
  Help: qfautoexport(qfa,{flag}): qfa being an automorphism group as output by
   qfauto, output a string representing the underlying matrix group in
   GAP notation (default) or Magma notation (flag = 1).
  Doc: \var{qfa} being an automorphism group as output by
   \tet{qfauto}, export the underlying matrix group as a string suitable
   for (no flags or $\fl=0$) GAP or ($\fl=1$) Magma. The following example
   computes the size of the matrix group using GAP:
   \bprog
   ? G = qfauto([2,1;1,2])
   %1 = [12, [[-1, 0; 0, -1], [0, -1; 1, 1], [1, 1; 0, -1]]]
   ? s = qfautoexport(G)
   %2 = "Group([[-1, 0], [0, -1]], [[0, -1], [1, 1]], [[1, 1], [0, -1]])"
   ? extern("echo \"Order("s");\" | gap -q")
   %3 = 12
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfautoexport(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfisom(*argv):
  '''
  qfisom
  Class: basic
  Section: linear_algebra
  C-Name: qfisom0
  Prototype: GGDG
  Help: qfisom(G,H,{fl}): find an isomorphism between the integral positive
   definite quadratic forms G and H if it exists. G can also be given by a
   qfisominit structure which is preferable if several forms need to be compared
   to G.
  Doc: 
   $G$, $H$ being square and symmetric matrices with integer entries representing
   positive definite quadratic forms, return an invertible matrix $S$ such that
   $G={^t}S\*H\*S$. This defines a isomorphism between the corresponding lattices.
   Since this requires computing the minimal vectors, the computations can
   become very lengthy as the dimension grows.
   See \kbd{qfisominit} for the meaning of \var{fl}.
   
   $G$ can also be given by an \kbd{qfisominit} structure which is preferable if
   several forms $H$ need to be compared to $G$.
   
   This function implements an algorithm of Plesken and Souvignier, following
   Souvignier's implementation.
  Variant: Also available is \fun{GEN}{qfisom}{GEN G, GEN H, GEN fl}
   where $G$ is a vector of \kbd{zm}, and $H$ is a \kbd{zm}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfisom0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def qfisominit(*argv):
  '''
  qfisominit
  Class: basic
  Section: linear_algebra
  C-Name: qfisominit0
  Prototype: GDG
  Help: qfisominit(G,{fl}): G being a square and symmetric matrix representing an
   integral positive definite quadratic form, this function return a structure
   allowing to compute isomorphisms between G and other quadratic form faster.
  Doc: 
   $G$ being a square and symmetric matrix with integer entries representing a
   positive definite quadratic form, return an \kbd{isom} structure allowing to
   compute isomorphisms between $G$ and other quadratic forms faster.
   
   The interface of this function is experimental and will likely change in future
   release.
   
   If present, the optional parameter \var{fl} must be a \typ{VEC} with two
   components. It allows to specify the invariants used, which can make the
   computation faster or slower. The components are
   
   \item \kbd{fl[1]} Depth of scalar product combination to use.
   
   \item \kbd{fl[2]} Maximum level of Bacher polynomials to use.
   
   Since this function computes the minimal vectors, it can become very lengthy
   as the dimension of $G$ grows.
  Variant: Also available is
   \fun{GEN}{qfisominit}{GEN F, GEN fl}
   where $F$ is a vector of \kbd{zm}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.qfisominit0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def random(*argv):
  '''
  random
  Class: basic
  Section: conversions
  C-Name: genrand
  Prototype: DG
  Help: random({N=2^31}): random object, depending on the type of N.
   Integer between 0 and N-1 (t_INT), int mod N (t_INTMOD), element in a finite
   field (t_FFELT), point on an elliptic curve (ellinit mod p or over a finite
   field).
  Doc: 
   returns a random element in various natural sets depending on the
   argument $N$.
   
   \item \typ{INT}: returns an integer
   uniformly distributed between $0$ and $N-1$. Omitting the argument
   is equivalent to \kbd{random(2\pow31)}.
   
   \item \typ{REAL}: returns a real number in $[0,1[$ with the same accuracy as
   $N$ (whose mantissa has the same number of significant words).
   
   \item \typ{INTMOD}: returns a random intmod for the same modulus.
   
   \item \typ{FFELT}: returns a random element in the same finite field.
   
   \item \typ{VEC} of length $2$, $N = [a,b]$: returns an integer uniformly
   distributed between $a$ and $b$.
   
   \item \typ{VEC} generated by \kbd{ellinit} over a finite field $k$
   (coefficients are \typ{INTMOD}s modulo a prime or \typ{FFELT}s): returns a
   ``random'' $k$-rational \emph{affine} point on the curve. More precisely
   if the curve has a single point (at infinity!) we return it; otherwise
   we return an affine point by drawing an abscissa uniformly at
   random until \tet{ellordinate} succeeds. Note that this is definitely not a
   uniform distribution over $E(k)$, but it should be good enough for
   applications.
   
   \item \typ{POL} return a random polynomial of degree at most the degree of $N$.
   The coefficients are drawn by applying \kbd{random} to the leading
   coefficient of $N$.
   
   \bprog
   ? random(10)
   %1 = 9
   ? random(Mod(0,7))
   %2 = Mod(1, 7)
   ? a = ffgen(ffinit(3,7), 'a); random(a)
   %3 = a^6 + 2*a^5 + a^4 + a^3 + a^2 + 2*a
   ? E = ellinit([3,7]*Mod(1,109)); random(E)
   %4 = [Mod(103, 109), Mod(10, 109)]
   ? E = ellinit([1,7]*a^0); random(E)
   %5 = [a^6 + a^5 + 2*a^4 + 2*a^2, 2*a^6 + 2*a^4 + 2*a^3 + a^2 + 2*a]
   ? random(Mod(1,7)*x^4)
   %6 = Mod(5, 7)*x^4 + Mod(6, 7)*x^3 + Mod(2, 7)*x^2 + Mod(2, 7)*x + Mod(5, 7)
   
   @eprog
   These variants all depend on a single internal generator, and are
   independent from your operating system's random number generators.
   A random seed may be obtained via \tet{getrand}, and reset
   using \tet{setrand}: from a given seed, and given sequence of \kbd{random}s,
   the exact same values will be generated. The same seed is used at each
   startup, reseed the generator yourself if this is a problem. Note that
   internal functions also call the random number generator; adding such a
   function call in the middle of your code will change the numbers produced.
   
   \misctitle{Technical note}
   Up to
   version 2.4 included, the internal generator produced pseudo-random numbers
   by means of linear congruences, which were not well distributed in arithmetic
   progressions. We now
   use Brent's XORGEN algorithm, based on Feedback Shift Registers, see
   \kbd{http://wwwmaths.anu.edu.au/\til{}brent/random.html}. The generator has period
   $2^{4096}-1$, passes the Crush battery of statistical tests of L'Ecuyer and
   Simard, but is not suitable for cryptographic purposes: one can reconstruct
   the state vector from a small sample of consecutive values, thus predicting
   the entire sequence.
  Variant: 
    Also available: \fun{GEN}{ellrandom}{GEN E} and \fun{GEN}{ffrandom}{GEN a}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.genrand(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def getrand(*argv):
  '''
  getrand
  Class: basic
  Section: programming/specific
  C-Name: getrand
  Prototype: 
  Help: getrand(): current value of random number seed.
  Doc: returns the current value of the seed used by the
   pseudo-random number generator \tet{random}. Useful mainly for debugging
   purposes, to reproduce a specific chain of computations. The returned value
   is technical (reproduces an internal state array), and can only be used as an
   argument to \tet{setrand}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.getrand(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def setrand(*argv):
  '''
  setrand
  Class: basic
  Section: programming/specific
  C-Name: setrand
  Prototype: vG
  Help: setrand(n): reset the seed of the random number generator to n.
  Doc: reseeds the random number generator using the seed $n$. No value is
   returned. The seed is either a technical array output by \kbd{getrand}, or a
   small positive integer, used to generate deterministically a suitable state
   array. For instance, running a randomized computation starting by
   \kbd{setrand(1)} twice will generate the exact same output.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  PyPari.pari.setrand(*c_arg_tuple)

def polgraeffe(*argv):
  '''
  polgraeffe
  Class: basic
  Section: polynomials
  C-Name: polgraeffe
  Prototype: G
  Help: polgraeffe(f): returns the Graeffe transform g of f, such that
   g(x^2) = f(x)f(-x)
  Doc: returns the \idx{Graeffe} transform $g$ of $f$, such that $g(x^2) = f(x)
   f(-x)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polgraeffe(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polroots(*argv):
  '''
  polroots
  Class: basic
  Section: polynomials
  C-Name: roots
  Prototype: Gp
  Help: polroots(x): complex roots of the polynomial x using
   Schonhage's method, as modified by Gourdon.
  Doc: complex roots of the polynomial
   \var{pol}, given as a column vector where each root is repeated according to
   its multiplicity. The precision is given as for transcendental functions: in
   GP it is kept in the variable \kbd{realprecision} and is transparent to the
   user, but it must be explicitly given as a second argument in library mode.
   
   The algorithm used is a modification of A.~Sch\"onhage\sidx{Sch\"onage}'s
   root-finding algorithm, due to and originally implemented by X.~Gourdon.
   Barring bugs, it is guaranteed to converge and to give the roots to the
   required accuracy.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.roots(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def galoissubcyclo(*argv):
  '''
  galoissubcyclo
  Class: basic
  Section: number_fields
  C-Name: galoissubcyclo
  Prototype: GDGD0,L,Dn
  Help: galoissubcyclo(N,H,{fl=0},{v}):Compute a polynomial (in variable v)
   defining the subfield of Q(zeta_n) fixed by the subgroup H of (Z/nZ)*. N can
   be an integer n, znstar(n) or bnrinit(bnfinit(y),[n,[1]],1). H can be given
   by a generator, a set of generator given by a vector or a HNF matrix (see
   manual). If flag is 1, output only the conductor of the abelian extension.
   If flag is 2 output [pol,f] where pol is the polynomial and f the conductor.
  Doc: computes the subextension
   of $\Q(\zeta_n)$ fixed by the subgroup $H \subset (\Z/n\Z)^*$. By the
   Kronecker-Weber theorem, all abelian number fields can be generated in this
   way (uniquely if $n$ is taken to be minimal).
   
   \noindent The pair $(n, H)$ is deduced from the parameters $(N, H)$ as follows
   
   \item $N$ an integer: then $n = N$; $H$ is a generator, i.e. an
   integer or an integer modulo $n$; or a vector of generators.
   
   \item $N$ the output of \kbd{znstar($n$)}. $H$ as in the first case
   above, or a matrix, taken to be a HNF left divisor of the SNF for $(\Z/n\Z)^*$
   (of type \kbd{$N$.cyc}), giving the generators of $H$ in terms of \kbd{$N$.gen}.
   
   \item $N$ the output of \kbd{bnrinit(bnfinit(y), $m$, 1)} where $m$ is a
   module. $H$ as in the first case, or a matrix taken to be a HNF left
   divisor of the SNF for the ray class group modulo $m$
   (of type \kbd{$N$.cyc}), giving the generators of $H$ in terms of \kbd{$N$.gen}.
   
   In this last case, beware that $H$ is understood relatively to $N$; in
   particular, if the infinite place does not divide the module, e.g if $m$ is
   an integer, then it is not a subgroup of $(\Z/n\Z)^*$, but of its quotient by
   $\{\pm 1\}$.
   
   If $fl=0$, compute a polynomial (in the variable \var{v}) defining the
   the subfield of $\Q(\zeta_n)$ fixed by the subgroup \var{H} of $(\Z/n\Z)^*$.
   
   If $fl=1$, compute only the conductor of the abelian extension, as a module.
   
   If $fl=2$, output $[pol, N]$, where $pol$ is the polynomial as output when
   $fl=0$ and $N$ the conductor as output when $fl=1$.
   
   The following function can be used to compute all subfields of
   $\Q(\zeta_n)$ (of exact degree \kbd{d}, if \kbd{d} is set):
   \bprog
   polsubcyclo(n, d = -1)=
   { my(bnr,L,IndexBound);
     IndexBound = if (d < 0, n, [d]);
     bnr = bnrinit(bnfinit(y), [n,[1]], 1);
     L = subgrouplist(bnr, IndexBound, 1);
     vector(#L,i, galoissubcyclo(bnr,L[i]));
   }
   @eprog\noindent
   Setting \kbd{L = subgrouplist(bnr, IndexBound)} would produce subfields of exact
   conductor $n\infty$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.galoissubcyclo(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polsubcyclo(*argv):
  '''
  polsubcyclo
  Class: basic
  Section: polynomials
  C-Name: polsubcyclo
  Prototype: LLDn
  Help: polsubcyclo(n,d,{v='x}): finds an equation (in variable v) for the d-th
   degree subfields of Q(zeta_n). Output is a polynomial, or a vector of
   polynomials if there are several such fields or none.
  Doc: gives polynomials (in variable $v$) defining the sub-Abelian extensions
   of degree $d$ of the cyclotomic field $\Q(\zeta_n)$, where $d\mid \phi(n)$.
   
   If there is exactly one such extension the output is a polynomial, else it is
   a vector of polynomials, possibly empty. To get a vector in all cases,
   use \kbd{concat([], polsubcyclo(n,d))}.
   
   The function \tet{galoissubcyclo} allows to specify exactly which
   sub-Abelian extension should be computed.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polsubcyclo(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def nfsubfields(*argv):
  '''
  nfsubfields
  Class: basic
  Section: number_fields
  C-Name: nfsubfields
  Prototype: GD0,L,
  Help: nfsubfields(pol,{d=0}): find all subfields of degree d of number field
   defined by pol (all subfields if d is null or omitted). Result is a vector of
   subfields, each being given by [g,h], where g is an absolute equation and h
   expresses one of the roots of g in terms of the root x of the polynomial
   defining nf.
  Doc: finds all subfields of degree
   $d$ of the number field defined by the (monic, integral) polynomial
   \var{pol} (all subfields if $d$ is null or omitted). The result is a vector
   of subfields, each being given by $[g,h]$, where $g$ is an absolute equation
   and $h$ expresses one of the roots of $g$ in terms of the root $x$ of the
   polynomial defining $\var{nf}$. This routine uses J.~Kl\"uners's algorithm
   in the general case, and B.~Allombert's \tet{galoissubfields} when \var{nf}
   is Galois (with weakly supersolvable Galois group).\sidx{Galois}\sidx{subfield}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nfsubfields(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrL1(*argv):
  '''
  bnrL1
  Class: basic
  Section: number_fields
  C-Name: bnrL1
  Prototype: GDGD0,L,p
  Help: bnrL1(bnr, {H}, {flag=0}): bnr being output by bnrinit(,,1) and
   H being a square matrix defining a congruence subgroup of bnr (the
   trivial subgroup if omitted), for each character of bnr trivial on this
   subgroup, compute L(1, chi) (or equivalently the first non-zero term c(chi)
   of the expansion at s = 0). The binary digits of flag mean 1: if 0 then
   compute the term c(chi) and return [r(chi), c(chi)] where r(chi) is the
   order of L(s, chi) at s = 0, or if 1 then compute the value at s = 1 (and in
   this case, only for non-trivial characters), 2: if 0 then compute the value
   of the primitive L-function associated to chi, if 1 then compute the value
   of the L-function L_S(s, chi) where S is the set of places dividing the
   modulus of bnr (and the infinite places), 3: return also the characters.
  Doc: let \var{bnr} be the number field data output by \kbd{bnrinit(,,1)} and
   \var{H} be a square matrix defining a congruence subgroup of the
   ray class group corresponding to \var{bnr} (the trivial congruence subgroup
   if omitted). This function returns, for each \idx{character} $\chi$ of the ray
   class group which is trivial on $H$, the value at $s = 1$ (or $s = 0$) of the
   abelian $L$-function associated to $\chi$. For the value at $s = 0$, the
   function returns in fact for each $\chi$ a vector $[r_\chi, c_\chi]$ where
   $$L(s, \chi) = c \cdot s^r + O(s^{r + 1})$$
   \noindent near $0$.
   
   The argument \fl\ is optional, its binary digits
   mean 1: compute at $s = 0$ if unset or $s = 1$ if set, 2: compute the
   primitive $L$-function associated to $\chi$ if unset or the $L$-function
   with Euler factors at prime ideals dividing the modulus of \var{bnr} removed
   if set (that is $L_S(s, \chi)$, where $S$ is the
   set of infinite places of the number field together with the finite prime
   ideals dividing the modulus of \var{bnr}), 3: return also the character if
   set.
   \bprog
   K = bnfinit(x^2-229);
   bnr = bnrinit(K,1,1);
   bnrL1(bnr)
   @eprog\noindent
   returns the order and the first non-zero term of $L(s, \chi)$ at $s = 0$
   where $\chi$ runs through the characters of the class group of
   $K = \Q(\sqrt{229})$. Then
   \bprog
   bnr2 = bnrinit(K,2,1);
   bnrL1(bnr2,,2)
   @eprog\noindent
   returns the order and the first non-zero terms of $L_S(s, \chi)$ at $s = 0$
   where $\chi$ runs through the characters of the class group of $K$ and $S$ is
   the set of infinite places of $K$ together with the finite prime $2$. Note
   that the ray class group modulo $2$ is in fact the class group, so
   \kbd{bnrL1(bnr2,0)} returns the same answer as \kbd{bnrL1(bnr,0)}.
   
   This function will fail with the message
   \bprog
    *** bnrL1: overflow in zeta_get_N0 [need too many primes].
   @eprog\noindent if the approximate functional equation requires us to sum
   too many terms (if the discriminant of $K$ is too large).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrL1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrrootnumber(*argv):
  '''
  bnrrootnumber
  Class: basic
  Section: number_fields
  C-Name: bnrrootnumber
  Prototype: GGD0,L,p
  Help: bnrrootnumber(bnr,chi,{flag=0}): returns the so-called Artin Root
   Number, i.e. the constant W appearing in the functional equation of the
   Hecke L-function associated to chi. Set flag = 1 if the character is known
   to be primitive.
  Doc: if $\chi=\var{chi}$ is a
   \idx{character} over \var{bnr}, not necessarily primitive, let
   $L(s,\chi) = \sum_{id} \chi(id) N(id)^{-s}$ be the associated
   \idx{Artin L-function}. Returns the so-called \idx{Artin root number}, i.e.~the
   complex number $W(\chi)$ of modulus 1 such that
   %
   $$\Lambda(1-s,\chi) = W(\chi) \Lambda(s,\overline{\chi})$$
   %
   \noindent where $\Lambda(s,\chi) = A(\chi)^{s/2}\gamma_\chi(s) L(s,\chi)$ is
   the enlarged L-function associated to $L$.
   
   The generators of the ray class group are needed, and you can set $\fl=1$ if
   the character is known to be primitive. Example:
   
   \bprog
   bnf = bnfinit(x^2 - x - 57);
   bnr = bnrinit(bnf, [7,[1,1]], 1);
   bnrrootnumber(bnr, [2,1])
   @eprog\noindent
   returns the root number of the character $\chi$ of
   $\Cl_{7\infty_1\infty_2}(\Q(\sqrt{229}))$ defined by $\chi(g_1^ag_2^b)
   = \zeta_1^{2a}\zeta_2^b$. Here $g_1, g_2$ are the generators of the
   ray-class group given by \kbd{bnr.gen} and $\zeta_1 = e^{2i\pi/N_1},
   \zeta_2 = e^{2i\pi/N_2}$ where $N_1, N_2$ are the orders of $g_1$ and
   $g_2$ respectively ($N_1=6$ and $N_2=3$ as \kbd{bnr.cyc} readily tells us).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrrootnumber(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnrstark(*argv):
  '''
  bnrstark
  Class: basic
  Section: number_fields
  C-Name: bnrstark
  Prototype: GDGp
  Help: bnrstark(bnr,{subgroup}): bnr being as output by
   bnrinit(,,1), finds a relative equation for the class field corresponding to
   the module in bnr and the given congruence subgroup (the trivial subgroup if
   omitted) using Stark's units. The ground field and the class field must be
   totally real.
  Doc: \var{bnr} being as output by \kbd{bnrinit(,,1)}, finds a relative equation
   for the class field corresponding to the modulus in \var{bnr} and the given
   congruence subgroup (as usual, omit $\var{subgroup}$ if you want the whole ray
   class group).
   
   The main variable of \var{bnr} must not be $x$, and the ground field and the
   class field must be totally real. When the base field is $\Q$, the vastly
   simpler \tet{galoissubcyclo} is used instead. Here is an example:
   \bprog
   bnf = bnfinit(y^2 - 3);
   bnr = bnrinit(bnf, 5, 1);
   bnrstark(bnr)
   @eprog\noindent
   returns the ray class field of $\Q(\sqrt{3})$ modulo $5$. Usually, one wants
   to apply to the result one of
   \bprog
   rnfpolredabs(bnf, pol, 16)     \\@com compute a reduced relative polynomial
   rnfpolredabs(bnf, pol, 16 + 2) \\@com compute a reduced absolute polynomial
   @eprog
   
   The routine uses \idx{Stark units} and needs to find a suitable auxiliary
   conductor, which may not exist when the class field is not cyclic over the
   base. In this case \kbd{bnrstark} is allowed to return a vector of
   polynomials defining \emph{independent} relative extensions, whose compositum
   is the requested class field. It was decided that it was more useful
   to keep the extra information thus made available, hence the user has to take
   the compositum herself.
   
   Even if it exists, the auxiliary conductor may be so large that later
   computations become unfeasible. (And of course, Stark's conjecture may simply
   be wrong.) In case of difficulties, try \tet{rnfkummer}:
   \bprog
   ? bnr = bnrinit(bnfinit(y^8-12*y^6+36*y^4-36*y^2+9,1), 2, 1);
   ? bnrstark(bnr)
     ***   at top-level: bnrstark(bnr)
     ***                 ^-------------
     *** bnrstark: need 3919350809720744 coefficients in initzeta.
     *** Computation impossible.
   ? lift( rnfkummer(bnr) )
   time = 24 ms.
   %2 = x^2 + (1/3*y^6 - 11/3*y^4 + 8*y^2 - 5)
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnrstark(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polzagier(*argv):
  '''
  polzagier
  Class: basic
  Section: polynomials
  C-Name: polzag
  Prototype: LL
  Help: polzagier(n,m): Zagier's polynomials of index n,m.
  Doc: creates Zagier's polynomial $P_n^{(m)}$ used in
   the functions \kbd{sumalt} and \kbd{sumpos} (with $\fl=1$). One must have $m\le
   n$. The exact definition can be found in ``Convergence acceleration of
   alternating series'', Cohen et al., Experiment.~Math., vol.~9, 2000, pp.~3--12.
   
   %@article {MR2001m:11222,
   %    AUTHOR = {Cohen, Henri and Rodriguez Villegas, Fernando and Zagier, Don},
   %     TITLE = {Convergence acceleration of alternating series},
   %   JOURNAL = {Experiment. Math.},
   %    VOLUME = {9},
   %      YEAR = {2000},
   %    NUMBER = {1},
   %     PAGES = {3--12},
   %}
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polzag(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bnfisintnorm(*argv):
  '''
  bnfisintnorm
  Class: basic
  Section: number_fields
  C-Name: bnfisintnorm
  Prototype: GG
  Help: bnfisintnorm(bnf,x): compute a complete system of solutions (modulo
   units of positive norm) of the absolute norm equation N(a)=x, where a
   belongs to the maximal order of big number field bnf (if bnf is not
   certified, this depends on GRH).
  Doc: computes a complete system of
   solutions (modulo units of positive norm) of the absolute norm equation
   $\Norm(a)=x$,
   where $a$ is an integer in $\var{bnf}$. If $\var{bnf}$ has not been certified,
   the correctness of the result depends on the validity of \idx{GRH}.
   
   See also \tet{bnfisnorm}.
  Variant: The function \fun{GEN}{bnfisintnormabs}{GEN bnf, GEN a}
   returns a complete system of solutions modulo units of the absolute norm
   equation $|\Norm(x)| = |a|$. As fast as \kbd{bnfisintnorm}, but solves
   the two equations $\Norm(x) = \pm a$ simultaneously.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bnfisintnorm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def thue(*argv):
  '''
  thue
  Class: basic
  Section: polynomials
  C-Name: thue
  Prototype: GGDG
  Help: thue(tnf,a,{sol}): solve the equation P(x,y)=a, where tnf was created
   with thueinit(P), and sol, if present, contains the solutions of Norm(x)=a
   modulo units in the number field defined by P. If tnf was computed without
   assuming GRH (flag 1 in thueinit), the result is unconditional. If tnf is a
   polynomial, compute thue(thueinit(P,0), a).
  Doc: returns all solutions of the equation
   $P(x,y)=a$ in integers $x$ and $y$, where \var{tnf} was created with
   $\kbd{thueinit}(P)$. If present, \var{sol} must contain the solutions of
   $\Norm(x)=a$ modulo units of positive norm in the number field
   defined by $P$ (as computed by \kbd{bnfisintnorm}). If there are infinitely
   many solutions, an error will be issued.
   
   It is allowed to input directly the polynomial $P$ instead of a \var{tnf},
   in which case, the function first performs \kbd{thueinit(P,0)}. This is
   very wasteful if more than one value of $a$ is required.
   
   If \var{tnf} was computed without assuming GRH (flag $1$ in \tet{thueinit}),
   then the result is unconditional. Otherwise, it depends in principle of the
   truth of the GRH, but may still be unconditionally correct in some
   favourable cases. The result is conditional on the GRH if
   $a\neq \pm 1$ and, $P$ has a single irreducible rational factor, whose
   associated tentative class number $h$ and regulator $R$ (as computed
   assuming the GRH) satisfy
   
   \item $h > 1$,
   
   \item $R/0.2 > 1.5$.
   
   Here's how to solve the Thue equation $x^{13} - 5y^{13} = - 4$:
   \bprog
   ? tnf = thueinit(x^13 - 5);
   ? thue(tnf, -4)
   %1 = [[1, 1]]
   @eprog\noindent In this case, one checks that \kbd{bnfinit(x\pow13 -5).no}
   is $1$. Hence, the only solution is $(x,y) = (1,1)$, and the result is
   unconditional. On the other hand:
   \bprog
   ? P = x^3-2*x^2+3*x-17; tnf = thueinit(P);
   ? thue(tnf, -15)
   %2 = [[1, 1]]  \\ a priori conditional on the GRH.
   ? K = bnfinit(P); K.no
   %3 = 3
   ? K.reg
   %4 = 2.8682185139262873674706034475498755834
   @eprog
   This time the result is conditional. All results computed using this
   particular \var{tnf} are likewise conditional, \emph{except} for a right-hand
   side of $\pm 1$.
   The above result is in fact correct, so we did not just disprove the GRH:
   \bprog
   ? tnf = thueinit(x^3-2*x^2+3*x-17, 1 /*unconditional*/);
   ? thue(tnf, -15)
   %4 = [[1, 1]]
   @eprog
   Note that reducible or non-monic polynomials are allowed:
   \bprog
   ? tnf = thueinit((2*x+1)^5 * (4*x^3-2*x^2+3*x-17), 1);
   ? thue(tnf, 128)
   %2 = [[-1, 0], [1, 0]]
   @eprog\noindent Reducible polynomials are in fact much easier to handle.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.thue(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def thueinit(*argv):
  '''
  thueinit
  Class: basic
  Section: polynomials
  C-Name: thueinit
  Prototype: GD0,L,p
  Help: thueinit(P,{flag=0}): initialize the tnf corresponding to P, that will
   be used to solve Thue equations P(x,y) = some-integer. If flag is non-zero,
   certify the result unconditionnaly. Otherwise, assume GRH (much faster of
   course).
  Doc: initializes the \var{tnf} corresponding to $P$, a univariate polynomial
   with integer coefficients. The result is meant to be used in conjunction with
   \tet{thue} to solve Thue equations $P(X / Y)Y^{\deg P} = a$, where $a$ is an
   integer.
   
   If $\fl$ is non-zero, certify results unconditionally. Otherwise, assume
   \idx{GRH}, this being much faster of course. In the latter case, the result
   may still be unconditionally correct, see \tet{thue}. For instance in most
   cases where $P$ is reducible (not a pure power of an irreducible), \emph{or}
   conditional computed class groups are trivial \emph{or} the right hand side
   is $\pm1$, then results are always unconditional.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.thueinit(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def agm(*argv):
  '''
  agm
  Class: basic
  Section: transcendental
  C-Name: agm
  Prototype: GGp
  Help: agm(x,y): arithmetic-geometric mean of x and y.
  Doc: arithmetic-geometric mean of $x$ and $y$. In the
   case of complex or negative numbers, the optimal AGM is returned
   (the largest in absolute value over all choices of the signs of the square
   roots).  $p$-adic or power series arguments are also allowed. Note that
   a $p$-adic agm exists only if $x/y$ is congruent to 1 modulo $p$ (modulo
   16 for $p=2$). $x$ and $y$ cannot both be vectors or matrices.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.agm(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def cos(*argv):
  '''
  cos
  Class: basic
  Section: transcendental
  C-Name: gcos
  Prototype: Gp
  Help: cos(x): cosine of x.
  Doc: cosine of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gcos(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def cotan(*argv):
  '''
  cotan
  Class: basic
  Section: transcendental
  C-Name: gcotan
  Prototype: Gp
  Help: cotan(x): cotangent of x.
  Doc: cotangent of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gcotan(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def exp(*argv):
  '''
  exp
  Class: basic
  Section: transcendental
  C-Name: gexp
  Prototype: Gp
  Help: exp(x): exponential of x.
  Description: 
   (real):real         mpexp($1)
   (mp):mp:prec        gexp($1, prec)
   (gen):gen:prec      gexp($1, prec)
  Doc: exponential of $x$.
   $p$-adic arguments with positive valuation are accepted.
  Variant: For a \typ{PADIC} $x$, the function
   \fun{GEN}{Qp_exp}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gexp(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def expm1(*argv):
  '''
  expm1
  Class: basic
  Section: transcendental
  C-Name: gexpm1
  Prototype: Gp
  Help: expm1(x): exp(x)-1.
  Description: 
   (real):real         mpexpm1($1)
  Doc: return $\exp(x)-1$, computed in a way that is also accurate
   when the real part of $x$ is near $0$. Only accept real or complex arguments.
   A naive direct computation would suffer from catastrophic cancellation;
   PARI's direct computation of $\exp(x)$ alleviates this well known problem at
   the expense of computing $\exp(x)$ to a higher accuracy when $x$ is small.
   Using \kbd{expm1} is recommanded instead:
   \bprog
   ? default(realprecision, 10000); x = 1e-100;
   ? a = expm1(x);
   time = 4 ms.
   ? b = exp(x)-1;
   time = 28 ms.
   ? default(realprecision, 10040); x = 1e-100;
   ? c = expm1(x);  \\ reference point
   ? abs(a-c)/c  \\ relative error in expm1(x)
   %7 = 0.E-10017
   ? abs(b-c)/c  \\ relative error in exp(x)-1
   %8 = 1.7907031188259675794 E-9919
   @eprog\noindent As the example above shows, when $x$ is near $0$,
   \kbd{expm1} is both faster and more accurate than \kbd{exp(x)-1}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gexpm1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def log(*argv):
  '''
  log
  Class: basic
  Section: transcendental
  C-Name: glog
  Prototype: Gp
  Help: log(x): natural logarithm of x.
  Description: 
   (gen):gen:prec        glog($1, prec)
  Doc: principal branch of the natural logarithm of
   $x \in \C^*$, i.e.~such that $\text{Im(log}(x))\in{} ]-\pi,\pi]$.
   The branch cut lies
   along the negative real axis, continuous with quadrant 2, i.e.~such that
   $\lim_{b\to 0^+} \log (a+bi) = \log a$ for $a \in\R^*$. The result is complex
   (with imaginary part equal to $\pi$) if $x\in \R$ and $x < 0$. In general,
   the algorithm uses the formula
   $$\log(x) \approx {\pi\over 2\text{agm}(1, 4/s)} - m \log 2, $$
   if $s = x 2^m$ is large enough. (The result is exact to $B$ bits provided
   $s > 2^{B/2}$.) At low accuracies, the series expansion near $1$ is used.
   
   $p$-adic arguments are also accepted for $x$, with the convention that
   $\log(p)=0$. Hence in particular $\exp(\log(x))/x$ is not in general equal to
   1 but to a $(p-1)$-th root of unity (or $\pm1$ if $p=2$) times a power of $p$.
  Variant: For a \typ{PADIC} $x$, the function
   \fun{GEN}{Qp_log}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.glog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sin(*argv):
  '''
  sin
  Class: basic
  Section: transcendental
  C-Name: gsin
  Prototype: Gp
  Help: sin(x): sine of x.
  Doc: sine of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsin(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sqrt(*argv):
  '''
  sqrt
  Class: basic
  Section: transcendental
  C-Name: gsqrt
  Prototype: Gp
  Help: sqrt(x): square root of x.
  Description: 
   (real):gen           sqrtr($1)
   (gen):gen:prec       gsqrt($1, prec)
  Doc: principal branch of the square root of $x$, defined as $\sqrt{x} =
   \exp(\log x / 2)$. In particular, we have
   $\text{Arg}(\text{sqrt}(x))\in{} ]-\pi/2, \pi/2]$, and if $x\in \R$ and $x<0$,
   then the result is complex with positive imaginary part.
   
   Intmod a prime $p$, \typ{PADIC} and \typ{FFELT} are allowed as arguments. In
   the first 2 cases (\typ{INTMOD}, \typ{PADIC}), the square root (if it
   exists) which is returned is the one whose first $p$-adic digit is in the
   interval $[0,p/2]$. For other arguments, the result is undefined.
  Variant: For a \typ{PADIC} $x$, the function
   \fun{GEN}{Qp_sqrt}{GEN x} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsqrt(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sqrtn(*argv):
  '''
  sqrtn
  Class: basic
  Section: transcendental
  C-Name: gsqrtn
  Prototype: GGD&p
  Help: sqrtn(x,n,{&z}): nth-root of x, n must be integer. If present, z is
   set to a suitable root of unity to recover all solutions. If it was not
   possible, z is set to zero.
  Doc: principal branch of the $n$th root of $x$,
   i.e.~such that $\text{Arg}(\text{sqrt}(x))\in{} ]-\pi/n, \pi/n]$. Intmod
   a prime and $p$-adics are allowed as arguments.
   
   If $z$ is present, it is set to a suitable root of unity allowing to
   recover all the other roots. If it was not possible, z is
   set to zero. In the case this argument is present and no square root exist,
   $0$ is returned instead or raising an error.
   \bprog
   ? sqrtn(Mod(2,7), 2)
   %1 = Mod(4, 7)
   ? sqrtn(Mod(2,7), 2, &z); z
   %2 = Mod(6, 7)
   ? sqrtn(Mod(2,7), 3)
     ***   at top-level: sqrtn(Mod(2,7),3)
     ***                 ^-----------------
     *** sqrtn: nth-root does not exist in gsqrtn.
   ? sqrtn(Mod(2,7), 3,  &z)
   %2 = 0
   ? z
   %3 = 0
   @eprog
   
   The following script computes all roots in all possible cases:
   \bprog
   sqrtnall(x,n)=
   { my(V,r,z,r2);
     r = sqrtn(x,n, &z);
     if (!z, error("Impossible case in sqrtn"));
     if (type(x) == "t_INTMOD" || type(x)=="t_PADIC",
       r2 = r*z; n = 1;
       while (r2!=r, r2*=z;n++));
     V = vector(n); V[1] = r;
     for(i=2, n, V[i] = V[i-1]*z);
     V
   }
   addhelp(sqrtnall,"sqrtnall(x,n):compute the vector of nth-roots of x");
   @eprog\noindent
  Variant: If $x$ is a \typ{PADIC}, the function
   \fun{GEN}{Qp_sqrt}{GEN x, GEN n, GEN *z} is also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsqrtn(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def tan(*argv):
  '''
  tan
  Class: basic
  Section: transcendental
  C-Name: gtan
  Prototype: Gp
  Help: tan(x): tangent of x.
  Doc: tangent of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtan(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Euler(*argv):
  '''
  Euler
  Class: basic
  Section: transcendental
  C-Name: mpeuler
  Prototype: p
  Help: Euler=Euler(): Euler's constant with current precision.
  Description: 
   ():real:prec        mpeuler(prec)
  Doc: Euler's constant $\gamma=0.57721\cdots$. Note that
   \kbd{Euler} is one of the few reserved names which cannot be used for
   user variables.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mpeuler(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Catalan(*argv):
  '''
  Catalan
  Class: basic
  Section: transcendental
  C-Name: mpcatalan
  Prototype: p
  Help: Catalan=Catalan(): Catalan's number with current precision.
  Description: 
   ():real:prec        mpcatalan(prec)
  Doc: Catalan's constant $G = \sum_{n>=0}\dfrac{(-1)^n}{(2n+1)^2}=0.91596\cdots$.
   Note that \kbd{Catalan} is one of the few reserved names which cannot be
   used for user variables.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mpcatalan(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def Pi(*argv):
  '''
  Pi
  Class: basic
  Section: transcendental
  C-Name: mppi
  Prototype: p
  Help: Pi=Pi(): the constant pi, with current precision.
  Description: 
   ():real:prec        mppi(prec)
  Doc: the constant $\pi$ ($3.14159\cdots$). Note that \kbd{Pi} is one of the few
   reserved names which cannot be used for user variables.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mppi(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sqrtnint(*argv):
  '''
  sqrtnint
  Class: basic
  Section: number_theoretical
  C-Name: sqrtnint
  Prototype: GL
  Help: sqrtnint(x,n): integer n-th root of x, where x is non-negative integer.
  Description: 
   (gen,small):int sqrtnint($1, $2)
  Doc: returns the integer $n$-th root of $x$, i.e. the largest integer $y$ such
   that $y^n \leq x$, where $x$ is a non-negative integer.
   \bprog
   ? N = 120938191237; sqrtnint(N, 5)
   %1 = 164
   ? N^(1/5)
   %2 = 164.63140849829660842958614676939677391
   @eprog\noindent The special case $n = 2$ is \tet{sqrtint}
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sqrtnint(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def teichmuller(*argv):
  '''
  teichmuller
  Class: basic
  Section: transcendental
  C-Name: teich
  Prototype: G
  Help: teichmuller(x): teichmuller character of p-adic number x.
  Doc: Teichm\"uller character of the $p$-adic number $x$, i.e. the unique
   $(p-1)$-th root of unity congruent to $x / p^{v_p(x)}$ modulo $p$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.teich(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bernfrac(*argv):
  '''
  bernfrac
  Class: basic
  Section: transcendental
  C-Name: bernfrac
  Prototype: L
  Help: bernfrac(x): Bernoulli number B_x, as a rational number.
  Doc: Bernoulli number\sidx{Bernoulli numbers} $B_x$,
   where $B_0=1$, $B_1=-1/2$, $B_2=1/6$,\dots, expressed as a rational number.
   The argument $x$ should be of type integer.
  '''
  c_params = []
  c_params.append(argv[0])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bernfrac(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bernpol(*argv):
  '''
  bernpol
  Class: basic
  Section: transcendental
  C-Name: bernpol
  Prototype: LDn
  Help: bernpol(n, {v = 'x}): Bernoulli polynomial B_n, in variable v.
  Doc: \idx{Bernoulli polynomial} $B_n$ in variable $v$.
   \bprog
   ? bernpol(1)
   %1 = x - 1/2
   ? bernpol(3)
   %2 = x^3 - 3/2*x^2 + 1/2*x
   @eprog
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bernpol(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def bernreal(*argv):
  '''
  bernreal
  Class: basic
  Section: transcendental
  C-Name: bernreal
  Prototype: Lp
  Help: bernreal(x): Bernoulli number B_x, as a real number with the current
   precision.
  Doc: Bernoulli number\sidx{Bernoulli numbers}
   $B_x$, as \kbd{bernfrac}, but $B_x$ is returned as a real number
   (with the current precision).
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.bernreal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def acosh(*argv):
  '''
  acosh
  Class: basic
  Section: transcendental
  C-Name: gacosh
  Prototype: Gp
  Help: acosh(x): inverse hyperbolic cosine of x.
  Doc: principal branch of $\text{cosh}^{-1}(x) = 2
    \log(\sqrt{(x+1)/2} + \sqrt{(x-1)/2})$. In particular,
   $\text{Re}(\text{acosh}(x))\geq 0$ and
   $\text{In}(\text{acosh}(x))\in ]-\pi,\pi]0$; if $x\in \R$ and $x<1$, then
   $\text{acosh}(x)$ is complex.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gacosh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def acos(*argv):
  '''
  acos
  Class: basic
  Section: transcendental
  C-Name: gacos
  Prototype: Gp
  Help: acos(x): arc cosine of x.
  Doc: principal branch of $\text{cos}^{-1}(x) = -i \log (x + i\sqrt{1-x^2})$.
   In particular, $\text{Re(acos}(x))\in [0,\pi]$ and if $x\in \R$ and $|x|>1$,
   then $\text{acos}(x)$ is complex. The branch cut is in two pieces:
   $]-\infty,-1]$ , continuous with quadrant II, and $[1,+\infty[$, continuous
   with quadrant IV. We have $\text{acos}(x) = \pi/2 - \text{asin}(x)$ for all
   $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gacos(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def arg(*argv):
  '''
  arg
  Class: basic
  Section: transcendental
  C-Name: garg
  Prototype: Gp
  Help: arg(x): argument of x,such that -pi<arg(x)<=pi.
  Doc: argument of the complex number $x$, such that $-\pi<\text{arg}(x)\le\pi$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.garg(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def asinh(*argv):
  '''
  asinh
  Class: basic
  Section: transcendental
  C-Name: gasinh
  Prototype: Gp
  Help: asinh(x): inverse hyperbolic sine of x.
  Doc: principal branch of $\text{sinh}^{-1}(x) = \log(x + \sqrt{1+x^2})$. In
   particular $\text{Im(asinh}(x))\in [-\pi/2,\pi/2]$.
   The branch cut is in two pieces: [-i oo ,-i],  continuous with quadrant III
   and [i,+i oo [ continuous with quadrant I.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gasinh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def asin(*argv):
  '''
  asin
  Class: basic
  Section: transcendental
  C-Name: gasin
  Prototype: Gp
  Help: asin(x): arc sine of x.
  Doc: principal branch of $\text{sin}^{-1}(x) = -i \log(ix + \sqrt{1 - x^2})$.
   In particular, $\text{Re(asin}(x))\in [-\pi/2,\pi/2]$ and if $x\in \R$ and
   $|x|>1$ then $\text{asin}(x)$ is complex. The branch cut is in two pieces:
   $]-\infty,-1]$, continuous with quadrant II, and $[1,+\infty[$ continuous
   with quadrant IV. The function satisfies $i \text{asin}(x) =
   \text{asinh}(ix)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gasin(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def atan(*argv):
  '''
  atan
  Class: basic
  Section: transcendental
  C-Name: gatan
  Prototype: Gp
  Help: atan(x): arc tangent of x.
  Doc: principal branch of $\text{tan}^{-1}(x) = \log ((1+ix)/(1-ix)) /
   2i$. In particular the real part of $\text{atan}(x))$ belongs to
   $]-\pi/2,\pi/2[$.
   The branch cut is in two pieces:
   $]-i\infty,-i[$, continuous with quadrant IV, and $]i,+i \infty[$ continuous
   with quadrant II. The function satisfies $i \text{atan}(x) =
   -i\text{atanh}(ix)$ for all $x\neq \pm i$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gatan(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def atanh(*argv):
  '''
  atanh
  Class: basic
  Section: transcendental
  C-Name: gatanh
  Prototype: Gp
  Help: atanh(x): inverse hyperbolic tangent of x.
  Doc: principal branch of $\text{tanh}^{-1}(x) = log ((1+x)/(1-x)) / 2$. In
   particular the imaginary part of $\text{atanh}(x)$ belongs to
   $[-\pi/2,\pi/2]$; if $x\in \R$ and $|x|>1$ then $\text{atanh}(x)$ is complex.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gatanh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def cosh(*argv):
  '''
  cosh
  Class: basic
  Section: transcendental
  C-Name: gcosh
  Prototype: Gp
  Help: cosh(x): hyperbolic cosine of x.
  Doc: hyperbolic cosine of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gcosh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def gammah(*argv):
  '''
  gammah
  Class: basic
  Section: transcendental
  C-Name: ggammah
  Prototype: Gp
  Help: gammah(x): gamma of x+1/2 (x integer).
  Doc: gamma function evaluated at the argument $x+1/2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ggammah(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def gamma(*argv):
  '''
  gamma
  Class: basic
  Section: transcendental
  C-Name: ggamma
  Prototype: Gp
  Help: gamma(s): gamma function at s, a complex or p-adic number, or a series.
  Doc: For $s$ a complex number, evaluates Euler's gamma
   function \sidx{gamma-function}
   $$\Gamma(s)=\int_0^\infty t^{s-1}\exp(-t)\,dt.$$
   Error if $s$ is a non-positive integer, where $\Gamma$ has a pole.
   
   For $s$ a \typ{PADIC}, evaluates the Morita gamma function at $s$, that
   is the unique continuous $p$-adic function on the $p$-adic integers
   extending $\Gamma_p(k)=(-1)^k \prod_{j<k}'j$, where the prime means that $p$
   does not divide $j$.
   \bprog
   ? gamma(1/4 + O(5^10))
   %1= 1 + 4*5 + 3*5^4 + 5^6 + 5^7 + 4*5^9 + O(5^10)
   ? algdep(%,4)
   %2 = x^4 + 4*x^2 + 5
   @eprog
  Variant: For a \typ{PADIC} $x$, the function \fun{GEN}{Qp_gamma}{GEN x} is
   also available.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ggamma(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def lngamma(*argv):
  '''
  lngamma
  Class: basic
  Section: transcendental
  C-Name: glngamma
  Prototype: Gp
  Help: lngamma(x): logarithm of the gamma function of x.
  Doc: principal branch of the logarithm of the gamma function of $x$. This
   function is analytic on the complex plane with non-positive integers
   removed, and can have much larger arguments than \kbd{gamma} itself.
   
   For $x$ a power series such that $x(0)$ is not a pole of \kbd{gamma},
   compute the Taylor expansion. (PARI only knows about regular power series
   and can't include logarithmic terms.)
   \bprog
   ? lngamma(1+x+O(x^2))
   %1 = -0.57721566490153286060651209008240243104*x + O(x^2)
   ? lngamma(x+O(x^2))
    ***   at top-level: lngamma(x+O(x^2))
    ***                 ^-----------------
    *** lngamma: domain error in lngamma: valuation != 0
   ? lngamma(-1+x+O(x^2))
    *** lngamma: Warning: normalizing a series with 0 leading term.
    ***   at top-level: lngamma(-1+x+O(x^2))
    ***                 ^--------------------
    *** lngamma: domain error in intformal: residue(series, pole) != 0
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.glngamma(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def psi(*argv):
  '''
  psi
  Class: basic
  Section: transcendental
  C-Name: gpsi
  Prototype: Gp
  Help: psi(x): psi-function at x.
  Doc: the $\psi$-function of $x$, i.e.~the logarithmic derivative
   $\Gamma'(x)/\Gamma(x)$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gpsi(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sinh(*argv):
  '''
  sinh
  Class: basic
  Section: transcendental
  C-Name: gsinh
  Prototype: Gp
  Help: sinh(x): hyperbolic sine of x.
  Doc: hyperbolic sine of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gsinh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def tanh(*argv):
  '''
  tanh
  Class: basic
  Section: transcendental
  C-Name: gtanh
  Prototype: Gp
  Help: tanh(x): hyperbolic tangent of x.
  Doc: hyperbolic tangent of $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gtanh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def factorial(*argv):
  '''
  factorial
  Class: basic
  Section: number_theoretical
  C-Name: mpfactr
  Prototype: Lp
  Help: factorial(x): factorial of x, the result being given as a real number.
  Doc: factorial of $x$. The expression $x!$ gives a result which is an integer,
   while $\kbd{factorial}(x)$ gives a real number.
  Variant: \fun{GEN}{mpfact}{long x} returns $x!$ as a \typ{INT}.
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.mpfactr(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sumformal(*argv):
  '''
  sumformal
  Class: basic
  Section: polynomials
  C-Name: sumformal
  Prototype: GDn
  Help: sumformal(f,{v}): formal sum of f with respect to v, or to the
   main variable of f if v is omitted.
  Doc: \idx{formal sum} of the polynomial expression $f$ with respect to the
   main variable if $v$ is omitted, with respect to the variable $v$ otherwise;
   it is assumed that the base ring has characteristic zero. In other words,
   considering $f$ as a polynomial function in the variable $v$,
   returns $F$, a polynomial in $v$ vanishing at $0$, such that $F(b) - F(a)
   = sum_{v = a+1}^b f(v)$:
   \bprog
   ? sumformal(n)  \\ 1 + ... + n
   %1 = 1/2*n^2 + 1/2*n
   ? f(n) = n^3+n^2+1;
   ? F = sumformal(f(n))  \\ f(1) + ... + f(n)
   %3 = 1/4*n^4 + 5/6*n^3 + 3/4*n^2 + 7/6*n
   ? sum(n = 1, 2000, f(n)) == subst(F, n, 2000)
   %4 = 1
   ? sum(n = 1001, 2000, f(n)) == subst(F, n, 2000) - subst(F, n, 1000)
   %5 = 1
   ? sumformal(x^2 + x*y + y^2, y)
   %6 = y*x^2 + (1/2*y^2 + 1/2*y)*x + (1/3*y^3 + 1/2*y^2 + 1/6*y)
   ? x^2 * y + x * sumformal(y) + sumformal(y^2) == %
   %7 = 1
   @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sumformal(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def dilog(*argv):
  '''
  dilog
  Class: basic
  Section: transcendental
  C-Name: dilog
  Prototype: Gp
  Help: dilog(x): dilogarithm of x.
  Doc: principal branch of the dilogarithm of $x$,
   i.e.~analytic continuation of the power series $\log_2(x)=\sum_{n\ge1}x^n/n^2$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.dilog(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def eta(*argv):
  '''
  eta
  Class: basic
  Section: transcendental
  C-Name: eta0
  Prototype: GD0,L,p
  Help: eta(z,{flag=0}): if flag=0, returns prod(n=1,oo, 1-q^n), where
   q = exp(2 i Pi z) if z is a complex scalar (belonging to the upper half plane);
   q = z if z is a p-adic number or can be converted to a power series.
   If flag is non-zero, the function only applies to complex scalars and returns
   the true eta function, with the factor q^(1/24) included.
  Doc: Variants of \idx{Dedekind}'s $\eta$ function.
   If $\fl = 0$, return $\prod_{n=1}^\infty(1-q^n)$, where $q$ depends on $x$
   in the following way:
   
   \item $q = e^{2i\pi x}$ if $x$ is a \emph{complex number} (which must then
   have positive imaginary part); notice that the factor $q^{1/24}$ is
   missing!
   
   \item $q = x$ if $x$ is a \typ{PADIC}, or can be converted to a
   \emph{power series} (which must then have positive valuation).
   
   If $\fl$ is non-zero, $x$ is converted to a complex number and we return the
   true $\eta$ function, $q^{1/24}\prod_{n=1}^\infty(1-q^n)$,
   where $q = e^{2i\pi x}$.
  Variant: 
   Also available is \fun{GEN}{trueeta}{GEN x, long prec} ($\fl=1$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.eta0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def erfc(*argv):
  '''
  erfc
  Class: basic
  Section: transcendental
  C-Name: gerfc
  Prototype: Gp
  Help: erfc(x): complementary error function.
  Doc: complementary error function, analytic continuation of
   $(2/\sqrt\pi)\int_x^\infty e^{-t^2}\,dt = \kbd{incgam}(1/2,x^2)/\sqrt\pi$,
   where the latter expression extends the function definition from real $x$ to
   all complex $x \neq 0$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gerfc(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def zeta(*argv):
  '''
  zeta
  Class: basic
  Section: transcendental
  C-Name: gzeta
  Prototype: Gp
  Help: zeta(s): Riemann zeta function at s with s a complex or a p-adic number.
  Doc: For $s$ a complex number, Riemann's zeta
   function \sidx{Riemann zeta-function} $\zeta(s)=\sum_{n\ge1}n^{-s}$,
   computed using the \idx{Euler-Maclaurin} summation formula, except
   when $s$ is of type integer, in which case it is computed using
   Bernoulli numbers\sidx{Bernoulli numbers} for $s\le0$ or $s>0$ and
   even, and using modular forms for $s>0$ and odd.
   
   For $s$ a $p$-adic number, Kubota-Leopoldt zeta function at $s$, that
   is the unique continuous $p$-adic function on the $p$-adic integers
   that interpolates the values of $(1 - p^{-k}) \zeta(k)$ at negative
   integers $k$ such that $k \equiv 1 \pmod{p-1}$ (resp. $k$ is odd) if
   $p$ is odd (resp. $p = 2$).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.gzeta(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def hyperu(*argv):
  '''
  hyperu
  Class: basic
  Section: transcendental
  C-Name: hyperu
  Prototype: GGGp
  Help: hyperu(a,b,x): U-confluent hypergeometric function.
  Doc: $U$-confluent hypergeometric function with
   parameters $a$ and $b$. The parameters $a$ and $b$ can be complex but
   the present implementation requires $x$ to be positive.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hyperu(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def incgam(*argv):
  '''
  incgam
  Class: basic
  Section: transcendental
  C-Name: incgam0
  Prototype: GGDGp
  Help: incgam(s,x,{g}): incomplete gamma function. g is optional and is the
   precomputed value of gamma(s).
  Doc: incomplete gamma function $\int_x^\infty e^{-t}t^{s-1}\,dt$, extended by
   analytic continuation to all complex $x, s$ not both $0$. The relative error
   is bounded in terms of the precision of $s$ (the accuracy of $x$ is ignored
   when determining the output precision). When $g$ is given, assume that
   $g=\Gamma(s)$. For small $|x|$, this will speed up the computation.
  Variant: Also available is \fun{GEN}{incgam}{GEN s, GEN x, long prec}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2].ref)
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.incgam0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def incgamc(*argv):
  '''
  incgamc
  Class: basic
  Section: transcendental
  C-Name: incgamc
  Prototype: GGp
  Help: incgamc(s,x): complementary incomplete gamma function.
  Doc: complementary incomplete gamma function.
   The arguments $x$ and $s$ are complex numbers such that $s$ is not a pole of
   $\Gamma$ and $|x|/(|s|+1)$ is not much larger than 1 (otherwise the
   convergence is very slow). The result returned is $\int_0^x
   e^{-t}t^{s-1}\,dt$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.incgamc(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besselh1(*argv):
  '''
  besselh1
  Class: basic
  Section: transcendental
  C-Name: hbessel1
  Prototype: GGp
  Help: besselh1(nu,x): H^1-bessel function of index nu and argument x.
  Doc: $H^1$-Bessel function of index \var{nu} and argument $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hbessel1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besselh2(*argv):
  '''
  besselh2
  Class: basic
  Section: transcendental
  C-Name: hbessel2
  Prototype: GGp
  Help: besselh2(nu,x): H^2-bessel function of index nu and argument x.
  Doc: $H^2$-Bessel function of index \var{nu} and argument $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.hbessel2(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besseli(*argv):
  '''
  besseli
  Class: basic
  Section: transcendental
  C-Name: ibessel
  Prototype: GGp
  Help: besseli(nu,x): I-bessel function of index nu and argument x.
  Doc: $I$-Bessel function of index \var{nu} and
   argument $x$. If $x$ converts to a power series, the initial factor
   $(x/2)^\nu/\Gamma(\nu+1)$ is omitted (since it cannot be represented in PARI
   when $\nu$ is not integral).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.ibessel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besselj(*argv):
  '''
  besselj
  Class: basic
  Section: transcendental
  C-Name: jbessel
  Prototype: GGp
  Help: besselj(nu,x): J-bessel function of index nu and argument x.
  Doc: $J$-Bessel function of index \var{nu} and
   argument $x$. If $x$ converts to a power series, the initial factor
   $(x/2)^\nu/\Gamma(\nu+1)$ is omitted (since it cannot be represented in PARI
   when $\nu$ is not integral).
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.jbessel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besseljh(*argv):
  '''
  besseljh
  Class: basic
  Section: transcendental
  C-Name: jbesselh
  Prototype: GGp
  Help: besseljh(n,x): J-bessel function of index n+1/2 and argument x, where
   n is a non-negative integer.
  Doc: $J$-Bessel function of half integral index.
   More precisely, $\kbd{besseljh}(n,x)$ computes $J_{n+1/2}(x)$ where $n$
   must be of type integer, and $x$ is any element of $\C$. In the
   present version \vers, this function is not very accurate when $x$ is small.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.jbesselh(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besseln(*argv):
  '''
  besseln
  Class: basic
  Section: transcendental
  C-Name: nbessel
  Prototype: GGp
  Help: besseln(nu,x): N-bessel function of index nu and argument x.
  Doc: $N$-Bessel function of index \var{nu} and argument $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.nbessel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def ellj(*argv):
  '''
  ellj
  Class: basic
  Section: elliptic_curves
  C-Name: jell
  Prototype: Gp
  Help: ellj(x): elliptic j invariant of x.
  Doc: 
   elliptic $j$-invariant. $x$ must be a complex number
   with positive imaginary part, or convertible into a power series or a
   $p$-adic number with positive valuation.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.jell(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def besselk(*argv):
  '''
  besselk
  Class: basic
  Section: transcendental
  C-Name: kbessel
  Prototype: GGp
  Help: besselk(nu,x): K-bessel function of index nu and argument x.
  Doc: $K$-Bessel function of index \var{nu} and argument $x$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.kbessel(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def polylog(*argv):
  '''
  polylog
  Class: basic
  Section: transcendental
  C-Name: polylog0
  Prototype: LGD0,L,p
  Help: polylog(m,x,{flag=0}): m-th polylogarithm of x. flag is optional, and
   can be 0: default, 1: D_m~-modified m-th polylog of x, 2: D_m-modified m-th
   polylog of x, 3: P_m-modified m-th polylog of x.
  Doc: one of the different polylogarithms, depending on \fl:
   
   If $\fl=0$ or is omitted: $m^\text{th}$ polylogarithm of $x$, i.e.~analytic
   continuation of the power series $\text{Li}_m(x)=\sum_{n\ge1}x^n/n^m$
   ($x < 1$). Uses the functional equation linking the values at $x$ and $1/x$
   to restrict to the case $|x|\leq 1$, then the power series when
   $|x|^2\le1/2$, and the power series expansion in $\log(x)$ otherwise.
   
   Using $\fl$, computes a modified $m^\text{th}$ polylogarithm of $x$.
   We use Zagier's notations; let $\Re_m$ denote $\Re$ or $\Im$ depending
   on whether $m$ is odd or even:
   
   If $\fl=1$: compute $\tilde D_m(x)$, defined for $|x|\le1$ by
   $$\Re_m\left(\sum_{k=0}^{m-1} \dfrac{(-\log|x|)^k}{k!}\text{Li}_{m-k}(x)
   +\dfrac{(-\log|x|)^{m-1}}{m!}\log|1-x|\right).$$
   
   If $\fl=2$: compute $D_m(x)$, defined for $|x|\le1$ by
   $$\Re_m\left(\sum_{k=0}^{m-1}\dfrac{(-\log|x|)^k}{k!}\text{Li}_{m-k}(x)
   -\dfrac{1}{2}\dfrac{(-\log|x|)^m}{m!}\right).$$
   
   If $\fl=3$: compute $P_m(x)$, defined for $|x|\le1$ by
   $$\Re_m\left(\sum_{k=0}^{m-1}\dfrac{2^kB_k}{k!}(\log|x|)^k\text{Li}_{m-k}(x)
   -\dfrac{2^{m-1}B_m}{m!}(\log|x|)^m\right).$$
   
   These three functions satisfy the functional equation
   $f_m(1/x) = (-1)^{m-1}f_m(x)$.
  Variant: Also available is
   \fun{GEN}{gpolylog}{long m, GEN x, long prec} (\fl = 0).
  '''
  c_params = []
  c_params.append(argv[0])
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_params.append(argv[3])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.polylog0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def sumdedekind(*argv):
  '''
  sumdedekind
  Class: basic
  Section: number_theoretical
  C-Name: sumdedekind
  Prototype: GG
  Help: sumdedekind(h,k): Dedekind sum associated to h,k
  Doc: returns the \idx{Dedekind sum} associated to the integers $h$ and $k$,
    corresponding to a fast implementation of
    \bprog
     s(h,k) = sum(n = 1, k-1, (n/k)*(frac(h*n/k) - 1/2))
    @eprog
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.sumdedekind(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def theta(*argv):
  '''
  theta
  Class: basic
  Section: transcendental
  C-Name: theta
  Prototype: GGp
  Help: theta(q,z): Jacobi sine theta-function.
  Doc: Jacobi sine theta-function
   $$ \theta_1(z, q) = 2q^{1/4} \sum_{n\geq 0} (-1)^n q^{n(n+1)} \sin((2n+1)z).$$
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.theta(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def thetanullk(*argv):
  '''
  thetanullk
  Class: basic
  Section: transcendental
  C-Name: thetanullk
  Prototype: GLp
  Help: thetanullk(q,k): k-th derivative at z=0 of theta(q,z).
  Doc: $k$-th derivative at $z=0$ of $\kbd{theta}(q,z)$.
  Variant: 
   \fun{GEN}{vecthetanullk}{GEN q, long k, long prec} returns the vector
   of all $\dfrac{d^i\theta}{dz^i}(q,0)$ for all odd $i = 1, 3, \dots, 2k-1$.
   \fun{GEN}{vecthetanullk_tau}{GEN tau, long k, long prec} returns
   \kbd{vecthetanullk\_tau} at $q = \exp(2i\pi \kbd{tau})$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.thetanullk(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def eint1(*argv):
  '''
  eint1
  Class: basic
  Section: transcendental
  C-Name: veceint1
  Prototype: GDGp
  Help: eint1(x,{n}): exponential integral E1(x). If n is present and x > 0,
   computes the vector of the first n values of the exponential integral E1(n.x)
  Doc: exponential integral $\int_x^\infty \dfrac{e^{-t}}{t}\,dt =
   \kbd{incgam}(0, x)$, where the latter expression extends the function
   definition from real $x > 0$ to all complex $x \neq 0$.
   
   If $n$ is present, we must have $x > 0$; the function returns the
   $n$-dimensional vector $[\kbd{eint1}(x),\dots,\kbd{eint1}(nx)]$. Contrary to
   other transcendental functions, and to the default case ($n$ omitted), the
   values are correct up to a bounded \emph{absolute}, rather than relative,
   error $10^-n$, where $n$ is \kbd{precision}$(x)$ if $x$ is a \typ{REAL}
   and defaults to \kbd{realprecision} otherwise. (In the most important
   application, to the computation of $L$-functions via approximate functional
   equations, those values appear as weights in long sums and small individual
   relative errors are less useful than controlling the absolute error.) This is
   faster than repeatedly calling \kbd{eint1($i$ * x)}, but less precise.
  Variant: Also available is \fun{GEN}{eint1}{GEN x, long prec}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1].ref)
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.veceint1(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def weber(*argv):
  '''
  weber
  Class: basic
  Section: transcendental
  C-Name: weber0
  Prototype: GD0,L,p
  Help: weber(x,{flag=0}): One of Weber's f function of x. flag is optional,
   and can be 0: default, function f(x)=exp(-i*Pi/24)*eta((x+1)/2)/eta(x),
   1: function f1(x)=eta(x/2)/eta(x)
   2: function f2(x)=sqrt(2)*eta(2*x)/eta(x). Note that
   j = (f^24-16)^3/f^24 = (f1^24+16)^3/f1^24 = (f2^24+16)^3/f2^24.
  Doc: one of Weber's three $f$ functions.
   If $\fl=0$, returns
   $$f(x)=\exp(-i\pi/24)\cdot\eta((x+1)/2)\,/\,\eta(x) \quad\hbox{such that}\quad
   j=(f^{24}-16)^3/f^{24}\,,$$
   where $j$ is the elliptic $j$-invariant  (see the function \kbd{ellj}).
   If $\fl=1$, returns
   $$f_1(x)=\eta(x/2)\,/\,\eta(x)\quad\hbox{such that}\quad
   j=(f_1^{24}+16)^3/f_1^{24}\,.$$
   Finally, if $\fl=2$, returns
   $$f_2(x)=\sqrt{2}\eta(2x)\,/\,\eta(x)\quad\hbox{such that}\quad
   j=(f_2^{24}+16)^3/f_2^{24}.$$
   Note the identities $f^8=f_1^8+f_2^8$ and $ff_1f_2=\sqrt2$.
  Variant: Also available are \fun{GEN}{weberf}{GEN x, long prec},
   \fun{GEN}{weberf1}{GEN x, long prec} and \fun{GEN}{weberf2}{GEN x, long prec}.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_params.append(argv[2])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.weber0(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

def lambertw(*argv):
  '''
  lambertw
  Class: basic
  Section: transcendental
  C-Name: glambertW
  Prototype: Gp
  Help: lambertw(y): solution of the implicit equation x*exp(x)=y.
  Doc: Lambert $W$ function, solution of the implicit equation $xe^x=y$,
   for $y > 0$.
  '''
  c_params = []
  c_params.append(argv[0].ref)
  c_params.append(argv[1])
  c_arg_tuple = tuple(c_params)
  av = PyPari.getAVMA()
  result = PyPari.Gen(ctypes.cast(PyPari.pari.gclone(PyPari.pari.glambertW(*c_arg_tuple)), ctypes.POINTER(ctypes.c_long)))
  PyPari.setAVMA(av)
  return result

