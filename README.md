
Matrix Algebra. Sparse and dense matrices I'm using for years.

Matrix class (dense matrices) is imo rather relliable and versatile. It contains lot of useful matrix functions. I have been developing and maintaing it for a long time. I frequently face problems in regards to Matrix Algebra and I use this class almost always. I find it very good for prototyping, testing, algorithms etc.

SparseMatrix<T> template class has been implemented as a standard vector of standard maps. I have succesfully used it many times and I have been sometimes very happy for having such a tool. I've used it mostly, in problems of mathematical optimization (Internal Point Optimization in particular) and in various geodetic problems (least squares geodetic adjustments). 

Regarding the basic problem of solving hyge system of linear equations - I have tested SparseMatrix<T> against boost implementation of sparse matrix (boost::numeric::ublas::compressed_matrix<T>), Pardiso and Eigen. My implementation is much faster the


