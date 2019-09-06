
Matrix Algebra. Sparse and dense matrices I've been using for years.

Matrix class is for dense matrices. It is (imo) reliable and versatile. It contains a lot of useful matrix functions. And, a lot of operators, as well - the complex Matrix expressions (e.g. Matrix X = inv(trans(A)&#42;P&#42;A) &#42; trans(A)&#42;P&#42;f ) are quite possible with this class. I frequently face problems in regards to Matrix Algebra and I use this class almost always. I find it very good for prototyping, testing, algorithms etc. This class has not been optimized for performancies. It has been designed for correctness, not for speed.

SparseMatrix<T> template class has been implemented as a standard vector of standard maps. I have been using it succesfully and I have been so happy, so many times, for having such a tool. (makes difference between crashes, overflows and getting job done within seconds) I have used it mostly, in mathematical optimization algorithms and specially, in various geodetic problems containing huge matrices. 

Regarding the basic problem of solving very big system of linear equations, I have tested SparseMatrix<T> against boost implementation of sparse matrix (boost::numeric::ublas::compressed_matrix<T>), Pardiso (Institute of Computational Science, Universita della Svizzera italiana, Lugano, Switzerland) and Eigen. For some geodetic problems (it means, mainly diagonal form and no need for pivoting), it turns out that SparseMatrix<T> is much faster than the boost's implementation. In additition, it has been slightly faster then Pardiso and slightly slower than Eigen. <img src="https://raw.githubusercontent.com/sasamil/WMS-TMS-Maker-Qt-GUI/master/icons/emoticons/eusa_wall.gif" alt="fail to break the wall" height="15" width="25">


