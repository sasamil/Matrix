
Matrix Algebra. Sparse and dense matrices I'm using for years.

Matrix class (dense matrices) is imo rather relliable and versatile. It contains a lot of useful matrix functions. And, a lot of operators, as well - the complex Matrix expressions (e.g. Matrix X = inv(trans(A)*P*A) * trans(A)*P*f ) are possible with this class. I have been developing and maintaing it for a long time. I frequently face problems in regards to Matrix Algebra and I use this class almost always. I find it very good for prototyping, testing, algorithms etc.

SparseMatrix<T> template class has been implemented as a standard vector of standard maps. I have succesfully used it and I have been many times, so happy for having such a tool. I've used it mostly, in problems of mathematical optimization (Internal Point Optimization in particular) and in various geodetic problems (least squares geodetic adjustments). 

Regarding the basic problem of solving hyge system of linear equations - I have tested SparseMatrix<T> against boost implementation of sparse matrix (boost::numeric::ublas::compressed_matrix<T>), Pardiso (Institute of Computational Science, Universita della Svizzera italiana, Lugano, Switzerland) and Eigen. For some geodetic problems (it means mainly diagonal form and no need for pivoting), it turns out that SparseMatrix<T> is much faster than the boost's implementation. In additition, it has been slightly faster then Pardiso and slightly slower than Eigen.


