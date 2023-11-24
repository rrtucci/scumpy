# SCuMpy

![Pond SCM](pond-scum.jpg)Pond SCM

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

$$\text{Mean Value of } \underline{a}=
\langle\underline{a}\rangle$$


$$\text{Covariance}(\underline{a}, \underline{b})=
\langle\underline{a}, 
\underline{b}\rangle$$


$$\text{Standard Deviation of }\underline{a} =
\sigma_{\underline{a}} = \sqrt{
\langle\underline{a}, 
\underline{a}\rangle}$$

$$\text{Correlation}(\underline{a}, \underline{b}) =
\rho_{
\underline{a}, \underline{b}
}=
\frac{\langle\underline{a}, 
\underline{b}\rangle}
{\sigma_{\underline{a}}\sigma_{\underline{b}}}$$

$$\frac{\partial\underline{b}}{
\partial\underline{a}}=
\frac{\langle\underline{a}, 
\underline{b}\rangle}
{\langle\underline{a}, 
\underline{a}\rangle}
$$

 Linear SCM are described by 
a system of linear equations of the form

$$\underline{x}_i = \sum_j \alpha_{i|j} \underline{x}_j + \underline {\epsilon}_i$$

where the 
$x_i$ are the internal nodes,
the $\alpha_{i|j}$
are
the path coefficients 
(a.k.a. gains), and the
$\underline{\epsilon}_i$ 
are the external nodes 
that inject noise into the system.
The $\underline{\epsilon}_i$ are
root nodes with
zero covariance
with each other.

We can view this as either

1. a linear system 
of equations with the unknowns 
$\underline{x}_i$. We can solve for these
unknowns using basic Linear Algebra.
Once we solve for the unknowns,
we can calculate $\langle\underline{x}_i, \underline{x}_k\rangle$.
2. a linear system 
of equations with the
unknowns $\alpha_{i|j}$. We can solve for these
unknowns using basic Linear Algebra.

SCuMpy
takes as input a DAG 
expressed as a dot file. A dot file
is a text file 
describing a single (usually)  DAG
in the dot language.
The dot  language is the language used to
describe DAGs by
the graph rendering software GraphViz.
SCuMpy stores all its dot files in a
folder entitled "dot_atlas".
Another term for "dot atlas"
is "DAG atlas".

At this point, 
given a DAG dot file as input,
SCuMpy can 
do (1) and (2), symbolically, using
the excellent
Python symbolic manipulator 
<a href="https://en.wikipedia.org/wiki/SymPy">SymPy</a>.
To test SCuMpy, we used the 20 DAGs 
defined in 
the "G&B-trols" paper:

> <a href=https://ftp.cs.ucla.edu/pub/stat_ser/r493.pdf>A Crash Course in Good and Bad Controls</a>,
by Carlos Cinelli, Andrew Forney and Judea Pearl


In the SCuMpy folder called "jupyter_notebooks",
you will find 20 notebooks where (1) is done
symbolically,
and another 20 notebooks where (2) is done symbolically, 
for each of the 20 DAGs in the G&B-trols paper. 

In addition, we have notebooks
for the 
[back door](https://github.com/rrtucci/scumpy/blob/master/jupyter_notebooks/back-door.ipynb), 
[front door](https://github.com/rrtucci/scumpy/blob/master/jupyter_notebooks/front-door.ipynb) 
and 
[napkin](https://github.com/rrtucci/scumpy/blob/master/jupyter_notebooks/napkin.ipynb) 
examples in
Pearl's the "Book of Why".

SCuMpy can also be used to prove rigorously,
from its symbolic output,
whether a do-query is identifiable
for a particular DAG. It can do this without
using the fairly complicated
Do Calculus rules, which are the most general
way of establishing identifiability. [See this notebook](https://github.com/rrtucci/scumpy/blob/master/jupyter_notebooks/unconfounded-children.ipynb)


SCuMpy can also do some numeric 
calculations. It can calculate
numerically the arrow
gain $\alpha_{i|j}$ for each arrow 
$\underline{x}_j\rightarrow\underline{x}_i$
of the DAG. The estimation algorithm 
requires as input a file which 
contains a dataset with the node 
names as column labels, and with 
node values as rows.
[See this notebook](https://github.com/rrtucci/scumpy/blob/master/jupyter_notebooks/estimating-gains.ipynb)

Exciting news. SCuMpy can now do linear SCM
with feedback loops. See [this blog post](https://qbnets.wordpress.com/2023/02/19/scumpy-can-now-do-linear-scm-with-feedback-loops/)
of mine for more details. This can be used to 
do Causal Inference with time series (a.k.a. panel data).

## Installation Instructions
See [this blog post](https://qbnets.wordpress.com/2023/01/26/first-version-of-scumpy-released-and-how-to-install-it-for-python-beginners/) of mine.



