#ifndef MANTID_CUSTOMINTERFACES_DATACOMPARISONHELPER_H_
#define MANTID_CUSTOMINTERFACES_DATACOMPARISONHELPER_H_

#include "MantidAPI/MatrixWorkspace.h"
#include "qwt_data.h"
#include <string>
#include <vector>

namespace MantidQt {
namespace CustomInterfaces {
namespace DataComparisonHelper {

/**
* Creates QwtData using X and Y values from the workspace spectra.
* @param ws :: Workspace with X and Y values to use
* @param wsIndex :: Workspace index to use
* @return Pointer to created QwtData
*/
QwtArrayData curveDataFromWs(Mantid::API::MatrixWorkspace_const_sptr ws,
                             size_t wsIndex) {
  const double *x = &ws->x(wsIndex)[0];
  const double *y = &ws->y(wsIndex)[0];
  size_t size = ws->blocksize();

  return QwtArrayData(x, y, size);
}

}
} // namespace CustomInterfaces
} // namespace MantidQt

#endif /* MANTID_CUSTOMINTERFACES_DATACOMPARISONHELPER_H_ */
