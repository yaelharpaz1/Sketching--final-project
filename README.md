# Sketching--final-project
Implementation of the paper 'Fast and accurate least-mean-squares solvers', Maalouf, Alaa, Ibrahim Jubran, sand Dan Feldman. Advances in Neural Information Processing Systems. 2019.

Motivation
The algorithm in this paper is for boosting the performance of existing LMS solvers like Linear/ Ridge/ Lasso/ Elastic-Net regression, while aiming to maintain the optimal solution both accurately and fast.
This algorithm gets a ﬁnite set of n d-dimensional real vectors and returns a weighted subset of d+1 vectors. Caratheodory’s Theorem computes such subset in O(n^2 d^2) time, while here this subset is computed in O(nd) time using techniques known as sketches and coresets.
