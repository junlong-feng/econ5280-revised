# Chapter 2 Matrix Algebra              

<font size="5">Junlong Feng</font>

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/junlong-feng/econ5280/main?filepath=Lecture2_MatrixAlgebra.ipynb)

## Outline

* Motivation: We use matrix to represent data.
* Notation: Scalar, vector and matrix.
* Addition and multiplication: Corner stones of all our calculations.
* Invertibility: $\neq 0$?
* Eigenvalue and positive (semi-) definiteness: $>0(\geq0)$ in matrix sense.
* Matrix calculus: Taking derivatives.  

## 1. Motivation

Recall in the previous lecture, we 

* used an $n\times 1$ column vector $Y$ to collect $n$ observations of the dependent variable, and a $n\times k$ matrix $X$ to collect all $n$ observations of the $k$ independent variables, and
* we said we are going to propose formulae based on $Y$ and $X$ to study the dice of God, so
* it is crucial to understand basics about matrices before we understand our formulae built on them.

## 2. Notation

We will play with three types of objects:

* Scalar: a number. 
* Vector: a list of numbers arranged in a column or a row. 
  * e.g. $A=\begin{pmatrix}a_{1}\\a_{2}\end{pmatrix}$, $B=(b_{1},b_{2})$. $a_{1},a_{2},b_{1},b_{2}$ are all scalars. $A$ is a $2\times 1$ vector (reads two by one). $B$ is a $1\times 2$ vector.
  * All vectors are treated as columns in this course if not otherwise stated.
* Matrix: a rectangle array of numbers. 
  * When we say a matrix is $n\times k$, we mean there are $n$ rows and $k$ columns in it. 
  * e.g. $A=\begin{pmatrix}a_{11}& a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\end{pmatrix}$. $A$ is a $2\times 3$ matrix. $a_{ij}$ denotes the **entry** in the $i$-th row and the $j$-th column, or we say $a_{ij}$ is the $(i,j)$-th entry in $A$, denoted by $A_{ij}$.

A few remarks are in order: 

* All **numbers** in the above definitions can be real or complex, deterministic or **random**.
* A scalar is a $1\times 1$ matrix and a vector is a $?\times 1$ (column) or $1\times ?$ (row) matrix. 
* Hence, matrix provides a unified perspective to view all kinds of mathematical objects we will see.

Sometimes it is convenient to represent a matrix in terms of the column or row vectors instead of entries. 

* e.g. $A=\begin{pmatrix}a_{11}& a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\end{pmatrix}$ can also be written as $A=(A_{1},A_{2},A_{3})$ where $A_{j}=\begin{pmatrix}a_{1j}\\a_{2j}\end{pmatrix}$ for all $j=1,2,3$.

### 2.1 Transpose

**Definition**. $A'$ is the transpose of $A$ if the $(i,j)$-th entry in $A'$ is the $(j,i)$-th entry in $A$ for all $i,j$.

Remarks:

* If $A$ is $n\times k$, then $A'$ is $k\times n$.
* e.g. if $A=\begin{pmatrix}a_{11}& a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\end{pmatrix}$, then $A'=\begin{pmatrix}a_{11}& a_{21}\\a_{12}&a_{22}\\a_{13}&a_{23}\end{pmatrix}$.
* Recall that we can write $A=(A_{1},A_{2},A_{3})$ where $A_{j}$ is the $j$-th column in $A$. Then you can verify that $A'=\begin{pmatrix}A_{1}'\\A_{2}'\\A_{3}'\end{pmatrix}$.

### 2.2 Symmetry

**Definition**. Matrix $A$ is symmetric if $A=A'$.

**Definition**. An $n\times k$ matrix $A$ is square if $n=k$.

Remarks.

* A symmetric matrix must be a square matrix.
* Matrix $A$ is symmetric if and only if $A_{ij}=A_{ji}$ for all $i,j$.

**Definition**. A symmetric matrix is **diagonal** if $A_{ij}=0$ for all $i\neq j$. A diagonal matrix is an **identity** matrix if $A_{ii}=1$ for all $i$.

Remarks.

* Example of a diagonal matrix: $A=\begin{pmatrix}2&0\\0&-4\end{pmatrix}$.
* Example of an identity matrix: $A=\begin{pmatrix}1&0\\0&1\end{pmatrix}$.

## 3. Addition and Multiplication

### 3.1 Addition

You can **only add up matrices of the same order**. For instance, if $A$ is $2\times 3$ and $B=2\times 3$, you can calculate their sum by:
$$
\begin{align*}
A+B&=\begin{pmatrix}a_{11}& a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\end{pmatrix}+\begin{pmatrix}b_{11}& b_{12}&b_{13}\\b_{21}&b_{22}&b_{23}\end{pmatrix}\\
&=\begin{pmatrix}a_{11}+b_{11}& a_{12}+b_{12}&a_{13}+b_{13}\\a_{21}+b_{21}&a_{22}+b_{22}&a_{23}+b_{23}\end{pmatrix}.
\end{align*}
$$

### 3.2 Scalar Multiplication

You can multiple a matrix (including vector and scalar) by a scalar: $cA=Ac=(cA_{ij}).$

### 3.3 Multiplication of Matrices

**Definition**. Two matrices $A$ and $B$ are **conformable** if the number of **columns** of $A$ is equal to the number of **rows** of $B$.

* e.g. If $A$ is $k\times r$ and $B$ is $r\times s$, they are conformable. Their product is a $k\times s$ matrix.

How do we calculate matrix multiplication? Let's start with a simple case.

**Multiplication of vectors: Inner product**. Let $a=(a_{1},a_{2},...,a_{k})'$ and $b=(b_{1},b_{2},...,b_{k})'$.

* First check whether they are conformable. 
  * No.
* What about $a'$ and $b$? 
  * Yes.
* Define $a'b=a_{1}b_{1}+a_{2}b_{2}\cdots a_{k}b_{k}\equiv \sum_{i=1}^{k}a_{i}b_{i}$.
* Inner product of two vectors yields a scalar.

**Multiplication of matrices: Representation I**. Let $A$ be $n\times k$ and $B$ be $n\times p$. Let $C\equiv A'B$. We have known that $C$ is $k\times p$. The $(a,b)$-th entry in $C$ is defined as:
$$
C_{ab}=\sum_{i=1}^{n}A_{ia}B_{ib}.
$$
**Multiplication of matrices: Representation II**. Another representation is:
$$
C=A_{1}'B_{1}+A'_{2}B_{2}+\cdots+ A'_{n} B_{n}=\sum_{i=1}^{n}A'_{i}B_{i},
$$
where $A_{i}$ and $B_{i}$ are the $i$-th rows in $A$ and $B$ respectively. Each $A'_{i}B_{i}$ ($k\times 1\times 1\times p$) is a $k\times p$ matrix, so the sum is also $k\times p$.

```
```





~~~R
A=matrix(c(1,2,3,4,5,6),3,2)
B=matrix(c(2,3,4,5,6,7),3,2)
## 1. Calculate AB by R
t(A)%*%B
## 2. Representation 1
C1=matrix(0,2,2)
for (i in 1:2){
  for (j in 1:2){
    C1[i,j]=sum(A[,i]*B[,j])
  }
}
C1
## 3. Representation 2
C2=matrix(0,2,2)
tA1=matrix(A[1,],2,1);tA2=matrix(A[2,],2,1)
tA3=matrix(A[3,],2,1);B1=matrix(B[1,],1,2)
B2=matrix(B[2,],1,2);B3=matrix(B[3,],1,2)
C2=tA1%*%B1+tA2%*%B2+tA3%*%B3
C2
~~~

Some final remarks:

* $A$ and $B$ are conformable **DOES NOT IMPLY** that $B$ and $A$ are conformable.
* In the above example, $B$ and $A$ are not conformable if $p\neq k$.
* So, that $AB$ exists **DOES NOT IMPLY** that $BA$ also exists.
* So of course $AB$ is **not necessarily** equal to $BA$. Matrix multiplication is **not commutative**.
* However, it is associative and distributive: under conformability, we have

$$
A(BC)=ABC=(AB)C, \text{\ \ }A(B+C)=AB+AC.
$$

* Another useful result: $(AB)'=B'A'$. No need to worry about conformability of $B'$ and $A'$; they must be conformable as long as $A$ and $B$ are.

## 4. Invertibility

Invertibility is only defined for square matrices. We denote the inverse of a matrix by $A^{-1}$. Like the scalar case, $A^{-1}A=AA^{-1}=I$, where $I$ is the identity matrix of the same order of $A$. 

However, not all square matrices are invertible. Consider a real scalar, i.e., a $1\times 1$ real matrix $a$. Scalar $a$ is invertible (i.e., $1/a$ is well-defined) if and only if $a\neq 0$. Unfortunately, this does not **directly** carry to matrix, though it does provides some insights.

Let us start from a diagonal matrix $A=\begin{pmatrix}a_{11}&0\\0&a_{22}\end{pmatrix}$. You can think that this matrix represents a 2-D space: $a_{11}$ is the length of the $x$-direction and $a_{22}$ is the length of the $y$-direction. **$A$ is invertible if and only if both directions have nonzero length**.

* Roughly speaking, a $k\times k$ invertible matrix is such that none of its $k$-dimensions is *fake*; you cannot find a space of lower dimension to contain it.
* Imagine $A=\begin{pmatrix}1&0\\0&0\end{pmatrix}$. This one is NOT invertible because it has a *fake* or *redundant dimension*; a 1-D space would have been enough to contain it.

If we formalize this idea, this means none of the rows or columns of a matrix can be expressed as a linear combination of other rows or columns. Otherwise, the row or column that can be expressed by others is redundant. 

**Definition**. For a set of vectors, if none of them can be written as a linear combination of others, then they are **linearly independent**.

**Definition**. The **rank** of a matrix is the number of linearly independent columns (or equivalently, rows). 

* By definition, for a $k\times r$ matrix, it is always true that $rank(A)\leq \min(k,r)$. 
* A matrix is **full-rank** if $rank(A)= \min(k,r)$.
* A matrix $A$ is invertible **if and only if** $A$ is **square** and **full-rank**. 

Some useful results:

* $rank(AB)\leq \min(rank(A),rank(B))$.
* **IMPORTANT**. $rank(AA')=rank(A'A)=rank(A)$. Since $AA'$ and $A'A$ are square no matter whether $A$ itself is square, this means $AA'$ and $A'A$ are both invertible **if and only if** $A$ is full-rank.
* If $A$ and $B$ are both invertible, then $(AB)^{-1}=B^{-1}A^{-1}$.
* For an invertible diagonal $A=\text{diag}(a_{ii})$, $A^{-1}=\text{diag}(1/a_{ii})$.
* For a $2\times 2$ matrix $A=\begin{pmatrix}a&b\\c&d\end{pmatrix}$, its inverse has a simple analytical form: $A^{-1}=\frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}$.

## 5. Eigenvalue and Positive (Semi-)Definiteness

For two real numbers $a$ and $b$, we know $a>b$ if and only if $a-b>0$. Like invertibility, can we borrow this idea and define how to compare two matrices?

First, how do we define $a-b>0$? For any real number $c$, $a-b>0$ if and only if $c^2\cdot (a-b)>0,\forall c\neq 0$. Similarly, $a-b\geq 0$ if and only if $c^{2}\cdot (a-b)\geq 0,\forall c$.

We borrow this idea and define a *positive* matrix:

**Definition.** A **symmetric** $k\times k$ matrix $A$ is positive definite (p.d.) if for any $k\times 1$ nonzero vector $x$, $x'Ax>0$. $A$ is positive semi-definite (p.s.d.) if for any $k\times 1$ vector $x$, $x'Ax\geq 0$.

We can write $A>0$ if $A$ is p.d. and $A\geq 0$ if $A$ is p.s.d. (We can also define negative definite (n.d.) and negative semi-definite (n.s.d.) similarly.)

The definition is simple but not easy to use; it requires one to calculate $x'Ax$ for all possible $x$. Would prefer a more straightforward way to verify whether a given symmetric matrix is p.d. or p.s.d.

Again, lets' start with a diagonal matrix $A=\text{diag}(a_{ii}),i=1,...,k$. You can verify that for an arbitrary vector $x=(x_{1},...,x_{k})'$,
$$
x'Ax=\sum_{i=1}^{k}x_{i}^{2}a_{ii}.
$$
Therefore, if $x'Ax>0$ all $x\neq 0$, it has to be the case that $a_{ii}>0$ for all $i$. Similarly, $x'Ax\geq 0$ for all $x$ if it $a_{ii}\geq 0$ for all $i$.

Now how about an arbitrary non-diagonal symmetric matrix? We have the following result:

**Theorem (Spectral Decomposition)**. For a real $k\times k$ symmetric matrix $A$, there exist a $k\times k$ matrix $H$ and a diagonal matrix $\Lambda$ such that i) $A=H\Lambda H'$ and ii) $HH'=H'H=I$ where $I$ is the $k\times k$ identity matrix.

* For a square but **non-symmetric** matrix $A$, there exist an invertible $H$ and a diagonal matrix $\Lambda$ such that $A=H\Lambda H^{-1}$. The decomposition for symmetric matrices is a special case where $H^{-1}=H'$. This is why in that case $HH'=H'H=I$.
* **Eigenvectors are not unique**. If all eigenvalues are distinct, eigenvectors are unique up to column-wise permutation column sign change.

**Definition**. The diagonal elements in $\Lambda$ are called **eigenvalues** of $A$. The column vectors in $H$ are called **eigenvectors** of $A$. Matrices $\Lambda$ and $A$ are **similar**.

By similarity, $A>(\geq)0$ if and only if $\Lambda>(\geq)0$, and we have already known that $\Lambda>(\geq)0$ if and only if all its diagonal elements are $>(\geq)0$. Therefore,

* A symmetric real matrix $A$ is p.d. or p.s.d. if and only if all its eigenvalues are positive or non-negative.
* One can verify that $AA'$ (and $A'A$) must be p.s.d. for any $A$ because $AA'=H\Lambda H'H\Lambda H=H\Lambda^{2}H'$, and thus all the eigenvalues of $AA'$ are nonnegative. Moreover, if $A$ is full-rank (i.e., all diagonal elements in $\Lambda$ are nonzero), then $AA'$ is p.d.

Finally, spectral decomposition and matrix similarity can also help us to know the rank. We said a square diagonal matrix is full-rank if none of its dimensions is redundant, i.e., all the diagonal elements are nonzero. For a non-diagonal matrix, one can think about it as $rotating$ the non-redundant dimensions of a diagonal matrix, and thus it should maintain all the fundamental geometric properties. Indeed, the rank of a square matrix is the same as the rank of the diagonal matrix of its eigenvalues. Therefore:

* The rank of a square matrix $A$ is equal to the number of nonzero eigenvalues.
* A square matrix is invertible if and only if it is full-rank, or equivalently, if and only it none of its eigenvalues is zero.
* Suppose $A$ is invertible. We know $A=H\Lambda H^{-1}$. By the above point, $\Lambda$ is invertible. Therefore, by the formula we learned in [Section 4](#4 Invertibility), $A^{-1}=H\Lambda^{-1}H^{-1}$. Hence, the inverse of a matrix has the same eigenvectors as the original matrix, and each eigenvalue is the inverse of the corresponding original one.

```R
library(Matrix)
X=matrix(c(2,1,1,2),2,2) # you can try a matrix that is not full rank
ed=eigen(X) # spectral decomposition
rankMatrix(X) # calculate the rank; need package "Matrix"
ed$values # call eigenvalues
ed$vectors # call eigenvectors
Xinv=solve(X) # taking inverse of a matrix
edinv=eigen(Xinv)
edinv$values
edinv$vectors # eigen vectors are the same as X up to permutation and sign change.
```

## 6 Calculus

It's useful to know how to take derivatives of the following functions w.r.t. a vector:

* For two vectors $a$ and $x$, $\partial_{a}(x'a)=\partial_{a}(a'x)=x$.
* For vector $a$, $\partial_{a}(a'a)=2a$.
* For matrix $X$ and a conformable vector $\beta$, $\partial_{\beta}X\beta=X'$ .
* Chain rule holds but you need to be careful about conformability: 

$$
\begin{align*}
\partial_{\beta}(Y-X\beta)'(Y-X\beta)=&\left(\partial_{\beta}(Y-X\beta)\right)\cdot \partial_{(Y-X\beta)}\left[(Y-X\beta)'(Y-X\beta)\right]\\
=&X'\cdot 2(Y-X\beta)\\
=&2(X'Y-X'X\beta).
\end{align*}
$$

## Further Readings

* Singular value and singular value decomposition (SVD) for rectangle matrices.
* Matrix cookbook: <https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf>.
* Trace and determinant.
*  Idempotent matrices.
