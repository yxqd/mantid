#ifndef MANTID_ALGORITHMS_MAXENT_H_
#define MANTID_ALGORITHMS_MAXENT_H_

#include "MantidAlgorithms/DllConfig.h"
#include "MantidAPI/Algorithm.h"
namespace Mantid {
namespace Algorithms {

/** MaxEnt : TODO: DESCRIPTION

  Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
  National Laboratory & European Spallation Source

  This file is part of Mantid.

  Mantid is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mantid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  File change history is stored at: <https://github.com/mantidproject/mantid>
  Code Documentation is available at: <http://doxygen.mantidproject.org>
*/

/// Auxiliary class to store the search directions (xi), quadratic coefficients
/// (c1, s1, c2, s2), angle and chi-sq
class SearchDirections {
  // SB eq. 21
  // SB eq. 24
public:
  SearchDirections(size_t dim, size_t points) {
    xi = Kernel::DblMatrix(dim, points);
    eta = Kernel::DblMatrix(dim, points);
    s1 = Kernel::DblMatrix(dim, 1);
    c1 = Kernel::DblMatrix(dim, 1);
    s2 = Kernel::DblMatrix(dim, dim);
    c2 = Kernel::DblMatrix(dim, dim);
  };
  Kernel::DblMatrix xi;  // Image space
  Kernel::DblMatrix eta; // Data space
  Kernel::DblMatrix s1;  // S_mu
  Kernel::DblMatrix c1;  // C_mu
  Kernel::DblMatrix s2;  // g_mu_nu
  Kernel::DblMatrix c2;  // M_mu_nu
  double chisq;
  double angle;
};

class DLLExport MaxEnt : public API::Algorithm {
public:
  /// Constructor
  MaxEnt();
  /// Destructor
  virtual ~MaxEnt();

  /// Algorithm's name
  virtual const std::string name() const;
  /// Algorithm's version
  virtual int version() const;
  /// Algorithm's category
  virtual const std::string category() const;
  /// Algorithm's summary
  virtual const std::string summary() const;

private:
  /// Initialise the algorithm's properties
  void init();
  /// Run the algorithm
  void exec();
  /// Validate the input properties
  std::map<std::string, std::string> validateInputs();
  /// Transforms from image space to data space
  std::vector<double> opus(const std::vector<double> &input);
  /// Transforms from data space to image space
  std::vector<double> tropus(const std::vector<double> &input);
  /// Calculates chi-square
  double getChiSq(const std::vector<double> &data,
                  const std::vector<double> &errors,
                  const std::vector<double> &dataCalc);
  /// Calculates the gradient of Chi
  std::vector<double> getCGrad(const std::vector<double> &data,
                               const std::vector<double> &errors,
                               const std::vector<double> &dataCalc);
  /// Calculates the gradient of S (entropy)
  std::vector<double> getSGrad(const std::vector<double> &image,
                               double background);
  /// Calculates the search directions and the quadratic coefficients
  SearchDirections calculateSearchDirections(const std::vector<double> &data,
                                             const std::vector<double> &error,
                                             const std::vector<double> &image,
                                             double background);
};

} // namespace Algorithms
} // namespace Mantid

#endif /* MANTID_ALGORITHMS_MAXENT_H_ */