#include "MantidKernel/Matrix.h"
#include "MantidKernel/V3D.h"
#include "MantidKernel/Exception.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidKernel/MersenneTwister.h"

#include <Eigen/LU>

using Mantid::Kernel::TimeSeriesProperty;

namespace Mantid {

namespace Kernel {
//
#define fabs(x) std::fabs((x)*1.0)

namespace {
//-------------------------------------------------------------------------
// Utility methods and function objects in anom_matrix.cols()mous namespace
//-------------------------------------------------------------------------
/**
\class PIndex
\author S. Ansell
\date Aug 2005
\version 1.0
\brief Class  to fill an index with a progressive count
*/
template <typename T> struct PIndex {
private:
  int count; ///< counter
public:
  /// Constructor
  PIndex() : count(0) {}
  /// functional
  std::pair<T, int> operator()(const T &A) {
    return std::pair<T, int>(A, count++);
  }
};

/**
\class PSep
\author S. Ansell
\date Aug 2005
\version 1.0
\brief Class to access the second object in index pair.
*/
template <typename T> struct PSep {
  /// Functional to the second object
  int operator()(const std::pair<T, int> &A) { return A.second; }
};

/**
* Function to take a vector and sort the vector
* so as to produce an index. Leaves the vector unchanged.
* @param pVec :: Input vector
* @param Index :: Output vector
*/
template <typename T>
void indexSort(const std::vector<T> &pVec, std::vector<int> &Index) {
  Index.resize(pVec.size());
  std::vector<typename std::pair<T, int>> PartList;
  PartList.resize(pVec.size());

  transform(pVec.begin(), pVec.end(), PartList.begin(), PIndex<T>());
  sort(PartList.begin(), PartList.end());
  transform(PartList.begin(), PartList.end(), Index.begin(), PSep<T>());
  return;
}
template void indexSort(const std::vector<double> &, std::vector<int> &);
template void indexSort(const std::vector<float> &, std::vector<int> &);
template void indexSort(const std::vector<int> &, std::vector<int> &);
}

template <typename T> std::vector<T> Matrix<T>::getVector() const {
  return std::vector<T>(m_matrix.data(), m_matrix.data() + m_matrix.size());
}
//
template <typename T>
Matrix<T>::Matrix(const size_t nrow, const size_t ncol, const bool makeIdentity)
/**
  Constructor with pre-set sizes. Matrix is zeroed
  @param nrow :: number of rows
  @param ncol :: number of columns
  @param makeIdentity :: flag for the constructor to return an identity matrix
*/
{

  if (!makeIdentity) {
    m_matrix = InternalMatrixType(nrow, ncol);
  } else {
    m_matrix = InternalMatrixType::Identity(nrow, ncol);
  }
}

template <typename T>
Matrix<T>::Matrix(const std::vector<T> &A, const std::vector<T> &B)
/**
  Constructor to take two vectors and multiply them to
  construct a matrix. (assuming that we have columns x row
  vector.
  @param A :: Column vector to multiply
  @param B :: Row vector to multiply
*/
{
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> x(A.data(),
                                                          m_matrix.rows());
  Eigen::Map<const Eigen::Matrix<T, 1, Eigen::Dynamic>> y(B.data(),
                                                          m_matrix.cols());

  m_matrix = x * y;
}
//
template <typename T> Matrix<T>::Matrix(const std::vector<T> &data) {
  size_t numel = data.size();
  size_t nxt = static_cast<size_t>(sqrt(double(numel)));
  size_t test = nxt * nxt;
  if (test != numel) {
    throw(std::invalid_argument(
        "number of elements in input vector have to be square of some value"));
  }

  m_matrix = Eigen::Map<const InternalMatrixType>(data.data(), nxt, nxt);
}

template <typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &A)
/**
  Matrix addition THIS + A
  If the size is different then 0 is added where appropiate
  Matrix A is not expanded.
  @param A :: Matrix to add
  @return Matrix(this + A)
*/
{
  m_matrix += A.m_matrix;

  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &A)
/**
  Matrix subtractoin THIS - A
  If the size is different then 0 is added where appropiate
  Matrix A is not expanded.
  @param A :: Matrix to add
  @return Ma
*/
{
  m_matrix -= A.m_matrix;

  return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &A) const
/**
  Matrix addition THIS + A
  If the size is different then 0 is added where appropiate
  Matrix A is not expanded.
  @param A :: Matrix to add
  @return Matrix(this + A)
*/
{
  Matrix<T> X(*this);
  return X += A;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &A) const
/**
  Matrix subtraction THIS - A
  If the size is different then 0 is subtracted where
  appropiate. This matrix determines the size
  @param A :: Matrix to add
  @return Matrix(this + A)
*/
{
  Matrix<T> X(*this);
  return X -= A;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &A) const
/**
  Matrix multiplication THIS * A
  @param A :: Matrix to multiply by  (this->row must == A->columns)
  @throw MisMatch<size_t> if there is a size mismatch.
  @return Matrix(This * A)
*/
{
  if (m_matrix.cols() != A.m_matrix.rows())
    throw Kernel::Exception::MisMatch<size_t>(
        m_matrix.cols(), A.m_matrix.rows(), "Matrix::operator*(Matrix)");

  Matrix<T> X(*this);
  X.m_matrix *= A.m_matrix;

  return X;
}

template <typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T> &Vec) const
/**
  Matrix multiplication THIS * Vec to produce a vec
  @param Vec :: size of vector > this->nrows
  @throw MisMatch<size_t> if there is a size mismatch.
  @return Matrix(This * Vec)
*/
{
  if (m_matrix.cols() > Vec.size())
    throw Kernel::Exception::MisMatch<size_t>(m_matrix.cols(), Vec.size(),
                                              "Matrix::operator*(Vec)");

  Eigen::Map<const InternalColVectorType> x(Vec.data(), Vec.size());

  T *dataPointer = static_cast<InternalColVectorType>(m_matrix * x).data();

  return std::vector<T>(dataPointer, dataPointer + Vec.size());
}

template <typename T>
V3D Matrix<T>::operator*(const V3D &Vx) const
/**
  Matrix multiplication THIS * V
  @param Vx :: Colunm vector to multiply by
  @throw MisMatch<size_t> if there is a size mismatch.
  @return Matrix(This * A)
*/
{
  if (m_matrix.cols() != 3 || m_matrix.rows() > 3)
    throw Kernel::Exception::MisMatch<size_t>(m_matrix.cols(), 3,
                                              "Matrix::operator*(V3D)");

  return V3D(m_matrix.template cast<double>() * Vx.getVector());
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &Value) const
/**
  Matrix multiplication THIS * Value
  @param Value :: Scalar to multiply by
  @return V * (this)
*/
{
  Matrix<T> X(*this);
  X.m_matrix *= Value;
  return X;
}

/**
  Matrix multiplication THIS *= A
  Note that we call operator* to avoid the problem
  of changing matrix size.
 @param A :: Matrix to multiply by  (this->row must == A->columns)
 @return This *= A
*/
template <typename T> Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &A) {
  if (m_matrix.cols() != A.m_matrix.rows())
    throw Kernel::Exception::MisMatch<size_t>(
        m_matrix.cols(), A.m_matrix.rows(), "Matrix*=(Matrix<T>)");
  // This construct to avoid the problem of changing size
  m_matrix *= A.m_matrix;
  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(const T &Value)
/**
  Matrix multiplication THIS * Value
  @param Value :: Scalar to multiply matrix by
  @return *this
*/
{
  m_matrix *= Value;
  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator/=(const T &Value)
/**
  Matrix divishio THIS / Value
  @param Value :: Scalar to multiply matrix by
  @return *this
*/
{
  m_matrix /= Value;
  return *this;
}

template <typename T>
bool Matrix<T>::operator!=(const Matrix<T> &A) const
/**
Element by Element comparison
@param A :: Matrix to check
@return true :: on succees
@return false :: failure
*/
{
  return !(this->operator==(A));
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T> &A) const
/**
Element by element comparison within tolerance.
Tolerance means that the value must be > tolerance
and less than (diff/max)>tolerance

Always returns 0 if the Matrix have different sizes
@param A :: matrix to check.
@return true on success
*/
{
  return this->equals(A, 1e-8);
}

//---------------------------------------------------------------------------------------
template <typename T>
bool Matrix<T>::equals(const Matrix<T> &A, const double Tolerance) const
/**
Element by element comparison within tolerance.
Tolerance means that the value must be > tolerance
and less than (diff/max)>tolerance

Always returns 0 if the Matrix have different sizes
@param A :: matrix to check.
@param Tolerance :: tolerance in comparing elements
@return true on success
*/
{
  if (&A != this) // this == A == always true
  {
    return (m_matrix - A.m_matrix).isZero(static_cast<T>(Tolerance));
  }
  // this == this is true
  return true;
}

//---------------------------------------------------------------------------------------
/** Element by element comparison of <
Always returns false if the Matrix have different sizes
@param A :: matrix to check.
@return true if this < A
*/
template <typename T> bool Matrix<T>::operator<(const Matrix<T> &A) const {
  if (&A == this) // this < A == always false
    return false;

  if (A.m_matrix.rows() != m_matrix.rows() ||
      A.m_matrix.cols() != m_matrix.cols())
    return false;

  return !((m_matrix - A.m_matrix).array() > 0.0).any();
}

//---------------------------------------------------------------------------------------
/** Element by element comparison of >=
Always returns false if the Matrix have different sizes
@param A :: matrix to check.
@return true if this >= A
*/
template <typename T> bool Matrix<T>::operator>=(const Matrix<T> &A) const {
  if (&A == this)
    return true;

  if (A.m_matrix.rows() != m_matrix.rows() ||
      A.m_matrix.cols() != m_matrix.cols())
    return false;

  return !((m_matrix - A.m_matrix).array() < 0.0).any();
}

/**
  Swap rows I and J
  @param RowI :: row I to swap
  @param RowJ :: row J to swap
*/
template <typename T>
void Matrix<T>::swapRows(const size_t RowI, const size_t RowJ) {
  m_matrix.row(RowI).swap(m_matrix.row(RowJ));
}

/**
  Swap columns I and J
  @param colI :: col I to swap
  @param colJ :: col J to swap
*/
template <typename T>
void Matrix<T>::swapCols(const size_t colI, const size_t colJ) {
  m_matrix.col(colI).swap(m_matrix.col(colJ));
}

template <typename T>
void Matrix<T>::zeroMatrix()
/**
  Zeros all elements of the matrix
*/
{
  if (m_matrix.rows() * m_matrix.cols()) {
    m_matrix.setZero();
  }
}

template <typename T>
void Matrix<T>::identityMatrix()
/**
  Makes the matrix an idenity matrix.
  Zeros all the terms outside of the square
*/
{
  if (m_matrix.rows() * m_matrix.cols()) {
    m_matrix.setIdentity(m_matrix.rows(), m_matrix.cols());
  }
}
template <typename T>
void Matrix<T>::setColumn(const size_t nCol, const std::vector<T> &newCol) {
  if (nCol >= m_matrix.cols()) {
    throw(std::invalid_argument("nCol requested> nCol availible"));
  }

  m_matrix.col(nCol) =
      Eigen::Map<const InternalColVectorType>(newCol.data(), newCol.size());
}
template <typename T>
void Matrix<T>::setRow(const size_t nRow, const std::vector<T> &newRow) {
  if (nRow >= m_matrix.rows()) {
    throw(std::invalid_argument("nRow requested> nRow availible"));
  }

  m_matrix.row(nRow) =
      Eigen::Map<const InternalRowVectorType>(newRow.data(), newRow.size());
}

template <typename T>
Matrix<T> Matrix<T>::preMultiplyByDiagonal(const std::vector<T> &Dvec) const
/**
  Creates a diagonal matrix D from the given vector Dvec and
  PRE-multiplies the matrix by it (i.e. D * M).
  @param Dvec :: diagonal matrix (just centre points)
  @return D*this
*/
{
  if (Dvec.size() != m_matrix.rows()) {
    std::ostringstream cx;
    cx << "Matrix::preMultiplyByDiagonal Size: " << Dvec.size() << " "
       << m_matrix.rows() << " " << m_matrix.cols();
    throw std::runtime_error(cx.str());
  }

  return Matrix<T>(Eigen::Map<const InternalColVectorType>(
                       Dvec.data(), Dvec.size()).asDiagonal() *
                   m_matrix);
}

template <typename T>
Matrix<T> Matrix<T>::postMultiplyByDiagonal(const std::vector<T> &Dvec) const
/**
  Creates a diagonal matrix D from the given vector Dvec and
  POST-multiplies the matrix by it (i.e. M * D).
  @param Dvec :: diagonal matrix (just centre points)
  @return this*D
*/
{
  if (Dvec.size() != m_matrix.cols()) {
    std::ostringstream cx;
    cx << "Error Matrix::bDiaognal size:: " << Dvec.size() << " "
       << m_matrix.rows() << " " << m_matrix.cols();
    throw std::runtime_error(cx.str());
  }

  return Matrix<T>(m_matrix *
                   Eigen::Map<const InternalColVectorType>(
                       Dvec.data(), Dvec.size()).asDiagonal());
}

template <typename T>
void Matrix<T>::setMem(const size_t newRows, const size_t newCols) {
  m_matrix.resize(newRows, newCols);
}

template <typename T>
Matrix<T> Matrix<T>::Tprime() const
/**
  Transpose the matrix :
  Has transpose for a square matrix case.
  @return M^T
*/
{
  if (!m_matrix.rows() * m_matrix.cols())
    return *this;

  Matrix<T> MT(*this);
  MT.Transpose();
  return MT;
}

template <typename T>
Matrix<T> &Matrix<T>::Transpose()
/**
  Transpose the matrix :
  Has a inplace transpose for a square matrix case.
  @return this^T
*/
{
  m_matrix = m_matrix.transpose();
  return *this;
}

template <>
void Matrix<int>::GaussJordan(Kernel::Matrix<int> &)
/**
  Not valid for Integer
  @throw std::invalid_argument
*/
{
  throw std::invalid_argument(
      "Gauss-Jordan inversion not valid for integer matrix");
}

template <typename T>
void Matrix<T>::GaussJordan(Matrix<T> &B)
/**
  Invert this matrix in place using Gauss-Jordan elimination.
  Matrix will be replaced by its inverse.
  @param B :: [input, output] Must have same dimensions as A. Returned as
  identity matrix. (?)
  @throw std::invalid_argument on input error
  @throw std::runtime_error if singular
 */
{}

template <typename T>
T Matrix<T>::Invert()
/**
  If the Matrix is square then invert the matrix
  using LU decomposition
  @return Determinant (0 if the matrix is singular)
*/
{
  Eigen::FullPivLU<InternalDoubleMatrixType> lu(
      m_matrix.template cast<double>());
  m_matrix = lu.inverse().template cast<T>();
  return static_cast<T>(lu.determinant());
}

template <typename T>
T Matrix<T>::determinant() const
/**
  Calculate the derminant of the matrix
  @return Determinant of matrix.
*/
{
  if (m_matrix.rows() != m_matrix.cols())
    throw Kernel::Exception::MisMatch<size_t>(
        m_matrix.rows(), m_matrix.cols(),
        "Determinant error :: Matrix is not m_matrix.rows()N");

  return static_cast<T>(m_matrix.template cast<double>().determinant());
}

template <typename T>
T Matrix<T>::factor()
/**
   Gauss jordan diagonal factorisation
   The diagonal is left as the values,
   the lower part is zero.
   @return the factored matrix
*/
{
  if (m_matrix.rows() != m_matrix.cols() || m_matrix.rows() < 1)
    throw std::runtime_error("Matrix::factor Matrix is not m_matrix.rows()N");

  return determinant();
}

template <typename T>
void Matrix<T>::normVert()
/**
  Normalise EigenVectors
  Assumes that they have already been calculated
*/
{
  m_matrix.rowwise().normalize();
}

template <typename T>
T Matrix<T>::compSum() const
/**
  Add up each component sums for the matrix
  @return \f$ \sum_i \sum_j V_{ij}^2 \f$
 */
{
  return m_matrix.sum();
}

template <typename T>
std::vector<T> Matrix<T>::Diagonal() const
/**
  Returns the diagonal form as a vector
  @return Diagonal elements
*/
{
  auto diag = m_matrix.diagonal();
  return std::vector<T>(diag.data(), diag.data() + diag.size());
}

template <typename T>
T Matrix<T>::Trace() const
/**
  Calculates the trace of the matrix
  @return Trace of matrix
*/
{
  return m_matrix.trace();
}

template <typename T>
void Matrix<T>::sortEigen(Matrix<T> &DiagMatrix)
/**
  Sorts the eigenvalues into increasing
  size. Moves the EigenVectors correspondingly
  @param DiagMatrix :: matrix of the EigenValues
*/
{
  if (m_matrix.cols() != m_matrix.rows() ||
      m_matrix.rows() != DiagMatrix.m_matrix.rows() ||
      m_matrix.rows() != DiagMatrix.m_matrix.cols()) {
    std::cerr << "Matrix not Eigen Form" << std::endl;
    throw(std::invalid_argument(" Matrix is not in an eigenvalue format"));
  }
  std::vector<int> index;
  std::vector<T> X = DiagMatrix.Diagonal();
  indexSort(X, index);
  Matrix<T> EigenVec(*this);
  for (size_t Icol = 0; Icol < m_matrix.rows(); Icol++) {
    for (size_t j = 0; j < m_matrix.rows(); j++) {
      m_matrix(j, Icol) = EigenVec.m_matrix(j, index[Icol]);
    }
    DiagMatrix.m_matrix(Icol, Icol) = X[index[Icol]];
  }

  return;
}

template <typename T>
int Matrix<T>::Diagonalise(Matrix<T> &EigenVec, Matrix<T> &DiagMatrix) const
/**
  Attempt to diagonalise the matrix IF symmetric
  @param EigenVec :: (output) the Eigenvectors matrix
  @param DiagMatrix :: the diagonal matrix of eigenvalues
  @return :: 1  on success 0 on failure
*/
{
  if (m_matrix.rows() != m_matrix.cols() || m_matrix.rows() < 1) {
    std::cerr << "Matrix not square" << std::endl;
    return 0;
  }
  for (size_t i = 0; i < m_matrix.rows(); i++)
    for (size_t j = i + 1; j < m_matrix.rows(); j++)
      if (fabs(m_matrix(i, j) - m_matrix(j, i)) > 1e-6) {
        std::cerr << "Matrix not symmetric" << std::endl;
        std::cerr << (*this);
        return 0;
      }
  return 0;
}

template <typename T>
bool Matrix<T>::isRotation() const
/** Check if a matrix represents a proper rotation
@ return :: true/false
*/
{
  if (this->m_matrix.rows() != this->m_matrix.cols())
    throw(std::invalid_argument("matrix is not square"));
  //  std::cout << "Matrix determinant-1 is " << (this->determinant()-1) <<
  //  std::endl;
  if (fabs(this->determinant() - 1) > 1e-5) {
    return false;
  } else {
    Matrix<T> prod(m_matrix.rows(), m_matrix.cols()),
        ident(m_matrix.rows(), m_matrix.cols(), true);
    prod = this->operator*(this->Tprime());
    //    std::cout << "Matrix * Matrix' = " << std::endl << prod << std::endl;
    return prod.equals(ident, 1e-5);
  }
}

template <typename T>
bool Matrix<T>::isOrthogonal() const
/** Check if a matrix is orthogonal. Same as isRotation, but allows determinant
to be -1
@ return :: true/false
*/
{
  if (this->m_matrix.rows() != this->m_matrix.cols())
    throw(std::invalid_argument("matrix is not square"));
  if (fabs(fabs(this->determinant()) - 1.) > 1e-5) {
    return false;
  } else {
    Matrix<T> prod(m_matrix.rows(), m_matrix.cols()),
        ident(m_matrix.rows(), m_matrix.cols(), true);
    prod = this->operator*(this->Tprime());
    return prod.equals(ident, 1e-7);
  }
}

template <typename T>
std::vector<T> Matrix<T>::toRotation()
/**
  Transform the matrix to a rotation matrix, by normalizing each column to 1
  @return :: a vector of scaling factors
  @throw :: std::invalid_argument if the absolute value of the determinant is
  less then 1e-10 or not square matrix
*/
{
  if (this->m_matrix.rows() != this->m_matrix.cols())
    throw(std::invalid_argument("matrix is not square"));
  if (fabs(this->determinant()) < 1e-10)
    throw(std::invalid_argument("Determinant is too small"));
  // step 1: orthogonalize the matrix
  for (size_t i = 0; i < this->m_matrix.cols(); ++i) {
    double spself = 0.;
    for (size_t j = 0; j < this->m_matrix.rows(); ++j)
      spself += (m_matrix(j, i) * m_matrix(j, i));
    for (size_t k = i + 1; k < this->m_matrix.cols(); ++k) {
      double spother = 0;
      for (size_t j = 0; j < this->m_matrix.rows(); ++j)
        spother += (m_matrix(j, i) * m_matrix(j, k));
      for (size_t j = 0; j < this->m_matrix.rows(); ++j)
        m_matrix(j, k) -= static_cast<T>(m_matrix(j, i) * spother / spself);
    }
  }
  // step 2: get scales and rescsale the matrix
  std::vector<T> scale(this->m_matrix.rows());
  T currentScale;
  for (size_t i = 0; i < this->m_matrix.cols(); ++i) {
    currentScale = T(0.);
    for (size_t j = 0; j < this->m_matrix.rows(); ++j)
      currentScale += (m_matrix(j, i) * m_matrix(j, i));
    currentScale = static_cast<T>(sqrt(static_cast<double>(currentScale)));
    if (currentScale < 1e-10)
      throw(std::invalid_argument("Scale is too small"));
    scale[i] = currentScale;
  }
  Matrix<T> scalingMatrix(m_matrix.rows(), m_matrix.cols()),
      change(m_matrix.rows(), m_matrix.cols(), true);
  for (size_t i = 0; i < this->m_matrix.cols(); ++i)
    scalingMatrix[i][i] = static_cast<T>(1.0 / scale[i]);
  *this = this->operator*(scalingMatrix);
  if (this->determinant() < 0.) {
    scale[0] = -scale[0];
    change[0][0] = static_cast<T>(-1);
    *this = this->operator*(change);
  }
  return scale;
}

template <typename T>
void Matrix<T>::print() const
/**
  Simple print out routine
 */
{
  write(std::cout, 10);
  return;
}

/** set matrix elements ito random values  in the range from  rMin to rMax*/
template <typename T>
void Matrix<T>::setRandom(size_t seed, double rMin, double rMax) {
  MersenneTwister rng(seed, rMin, rMax);

  for (size_t i = 0; i < m_matrix.rows(); i++) {
    for (size_t j = 0; j < m_matrix.cols(); j++) {
      m_matrix(i, j) = static_cast<T>(rng.nextValue());
    }
  }
}

template <typename T>
void Matrix<T>::write(std::ostream &Fh, const int blockCnt) const
/**
  Write out function for blocks of 10 Columns
  @param Fh :: file stream for output
  @param blockCnt :: number of columns per line (0 == full)
*/
{
  std::ios::fmtflags oldFlags = Fh.flags();
  Fh.setf(std::ios::floatfield, std::ios::scientific);
  const size_t blockNumber((blockCnt > 0) ? blockCnt : m_matrix.cols());
  size_t BCnt(0);
  do {
    const size_t ACnt = BCnt;
    BCnt += blockNumber;
    if (BCnt > m_matrix.cols()) {
      BCnt = m_matrix.cols();
    }

    if (ACnt) {
      Fh << " ----- " << ACnt << " " << BCnt << " ------ " << std::endl;
    }
    for (size_t i = 0; i < m_matrix.rows(); i++) {
      for (size_t j = ACnt; j < BCnt; j++) {
        Fh << std::setw(10) << m_matrix(i, j) << "  ";
      }
      Fh << std::endl;
    }
  } while (BCnt < m_matrix.cols());

  Fh.flags(oldFlags);
  return;
}

template <typename T>
std::string Matrix<T>::str() const
/**
  Convert the matrix into a simple linear string expression
  @return String value of output
*/
{
  std::ostringstream cx;
  for (size_t i = 0; i < m_matrix.rows(); i++) {
    for (size_t j = 0; j < m_matrix.cols(); j++) {
      cx << std::setprecision(6) << m_matrix(i, j) << " ";
    }
  }
  return cx.str();
}

/**
 * Write an object to a stream. Format will be
 * Matrix(nrows,ncols)x_00,x_01...,x_10,x_11
 * @param os :: output stream
 * @param matrix :: Matrix to write out
 * @return The output stream (of)
*/
template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix) {
  dumpToStream(os, matrix, ',');
  return os;
}

/**
 * Write a Matrix to a stream. Format will be
 * Matrix(nrowsSEPncols)x_00SEPx_01...SEPx_10SEPx_11
 * @param os :: output stream
 * @param matrix :: Matrix to write out
 * @param delimiter :: A character to use as delimiter for the string
*/
template <typename T>
void dumpToStream(std::ostream &os, const Kernel::Matrix<T> &matrix,
                  const char delimiter) {
  size_t nrows(matrix.numRows()), ncols(matrix.numCols());
  os << "Matrix(" << nrows << delimiter << ncols << ")";
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      os << matrix[i][j];
      if (i < nrows - 1 || j < ncols - 1)
        os << delimiter;
    }
  }
}

/**
* Fill an object from a stream. Format should be
* Matrix(nrows,ncols)x_00,x_01...,x_10,x_11
* @param is :: A stream object
* @param in :: An object to fill
* @returns A reference to the stream
*/
template <typename T>
std::istream &operator>>(std::istream &is, Kernel::Matrix<T> &in) {
  fillFromStream(is, in, ',');
  return is;
}

/**
* Fill a Matrix from a stream using the given separator. Format should be
* Matrix(nrowsSEPncols)x_00SEPx_01...SEPx_10SEPx_11
* where SEP is replaced by the given separator
* @param is :: A stream object
* @param in :: An Matrix object to fill
* @param delimiter :: A single character separator that delimits the entries
*/
template <typename T>
void fillFromStream(std::istream &is, Kernel::Matrix<T> &in,
                    const char delimiter) {
  // Stream should start with Matrix(
  char dump;
  std::string start(7, ' ');
  for (int i = 0; i < 7; ++i) {
    is >> dump;
    start[i] = dump;
    if (!is)
      throw std::invalid_argument(
          "Unexpected character when reading Matrix from stream.");
  }
  if (start != "Matrix(")
    throw std::invalid_argument("Incorrect input format for Matrix stream.");
  // Now read a nrows,ncols and )
  size_t nrows(0), ncols(0);
  is >> nrows;
  if (!is)
    throw std::invalid_argument("Expected number of rows when reading Matrix "
                                "from stream, found something else.");
  is >> dump;
  is >> ncols;
  if (!is)
    throw std::invalid_argument("Expected number of columns when reading "
                                "Matrix from stream, found something else.");
  is >> dump;
  if (dump != ')')
    throw std::invalid_argument("Expected closing parenthesis after ncols when "
                                "reading Matrix from stream, found something "
                                "else.");

  // Resize the matrix
  in.setMem(nrows, ncols);

  // Use getline with the delimiter set to "," to read
  std::string value_str;
  size_t row(0), col(0);
  while (!is.eof() && std::getline(is, value_str, delimiter)) {
    try {
      T value = boost::lexical_cast<T>(value_str);
      in.m_matrix(row, col) = value;
    } catch (boost::bad_lexical_cast &) {
      throw std::invalid_argument(
          "Unexpected type found while reading Matrix from stream: \"" +
          value_str + "\"");
    }
    ++col;
    if (col == ncols) // New row
    {
      col = 0;
      ++row;
    }
  }
}

///\cond TEMPLATE

// Symbol definitions for common types
template class MANTID_KERNEL_DLL Matrix<double>;
template class MANTID_KERNEL_DLL Matrix<int>;
template class MANTID_KERNEL_DLL Matrix<float>;

template MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &,
                                                    const DblMatrix &);
template MANTID_KERNEL_DLL void dumpToStream(std::ostream &, const DblMatrix &,
                                             const char);
template MANTID_KERNEL_DLL std::istream &operator>>(std::istream &,
                                                    DblMatrix &);
template MANTID_KERNEL_DLL void fillFromStream(std::istream &, DblMatrix &,
                                               const char);

template MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &,
                                                    const Matrix<float> &);
template MANTID_KERNEL_DLL void dumpToStream(std::ostream &,
                                             const Matrix<float> &, const char);
template MANTID_KERNEL_DLL std::istream &operator>>(std::istream &,
                                                    Matrix<float> &);
template MANTID_KERNEL_DLL void fillFromStream(std::istream &, Matrix<float> &,
                                               const char);

template MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &,
                                                    const IntMatrix &);
template MANTID_KERNEL_DLL void dumpToStream(std::ostream &, const IntMatrix &,
                                             const char);
template MANTID_KERNEL_DLL std::istream &operator>>(std::istream &,
                                                    IntMatrix &);
template MANTID_KERNEL_DLL void fillFromStream(std::istream &, IntMatrix &,
                                               const char);
///\endcond TEMPLATE

} // namespace Kernel
} // namespace
