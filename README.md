# SCuMpy

![Pond SCM](pond-scum.jpg)

Welcome to SCuMpy! The name is intended
to (1) evoke the name of the amazing 
Python library NumPy, and (2) be a 
portmanteau of SCM (Structural Causal Model) 
and python.

A linear SCM is a DAG 
whose arrows are labelled by path coefficients
(a.k.a. gains). If you want to learn 
more about linear SCM, explained 
using my notational conventions 
which are also SCuMpy’s notational conventions, 
check out my free, open 
source book <a href="https://qbnets.wordpress.com/2020/11/30/my-free-book-bayesuvius-on-bayesian-networks/">
Bayesuvius</a>. 
Look in the chapter entitled 
“Linear Deterministic Bnets with External Noise”.

In SCuMpy, we use the following simple notation. 
Let $\underline{a}, \underline{b}$ be any 2 
random variables. (I underline random variables instead
of the usual convention of capitalizing them.) Then

$$\text{Covariance}(\underline{a}, \underline{b})=
\langle\underline{a}, 
\underline{b}\rangle$$


$$\text{Standard Deviation of }\underline{a} =
\sigma_{\underline{a}} = \sqrt{
\langle\underline{a}, 
\underline{a}\rangle}$$

$$\text{Correlation}(\underline{a}, \underline{b}) =
\frac{\langle\underline{a}, 
\underline{b}\rangle}
{\sigma_{\underline{a}}\sigma_{\underline{b}}}$$

$$\frac{\partial\underline{a}}{
\partial\underline{b}}=
\frac{\langle\underline{a}, 
\underline{b}\rangle}
{\langle\underline{a}, 
\underline{a}\rangle}
$$

 Linear SCM are described by 
a system of linear equations of the form

$$\underline{x}_j = \sum_i \alpha_{j|i}\quad 
\underline{x}_i + \underline{\epsilon}_j$$

where the 
$x_j$ are the the nodes,
the $\alpha_{j|i}$
are
the path coefficients 
(a.k.a. gains), and the
$\underline{\epsilon}_j$ 
are the external noise.
The $\underline{\epsilon}_j$ are
root nodes with
zero mean,
and zero covariance among themselves.

We can view this as either

1. a linear system 
of equations with the unknowns 
$\underline{x}_j$. We can solve for these
unknowns using basic Linear Algebra.
Once we solve for the unknowns,
we can calculate $\langle\underline{x}_j, \underline{x}_k\rangle$.
2. a linear system 
of equations with the
unknowns $\alpha_{j|i}$. We can solve for these
unknowns using basic Linear Algebra.

At this point, SCuMpy can 
do (1) and (2), symbolically, using
the excellent
Python symbolic manipulator 
<a href="https://en.wikipedia.org/wiki/SymPy">SymPy</a>.
To test SCuMpy, we used the 20 DAGs 
defined in 
the "G&B-trols" paper:

<a href=https://ftp.cs.ucla.edu/pub/stat_ser/r493.pdf>
A Crash Course in Good and Bad Controls</a>,
by Carlos Cinelli, Andrew Forney and Judea Pearl

In the SCuMpy folder called "jupyter_notebooks",
you will find 20 notebooks where (1) is done
symbolically,
and another 20 notebooks where (2) is done symbolically, 
for each of the 20 DAGs in the G&B-trols paper. 





