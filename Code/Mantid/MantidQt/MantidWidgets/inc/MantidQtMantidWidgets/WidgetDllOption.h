#ifndef MANTIDQT_MANTIDWIDGETS_DLLOPTION_H_
#define MANTIDQT_MANTIDWIDGETS_DLLOPTION_H_

#include <MantidKernel/System.h>

#ifdef IN_MANTIDQT_MANTIDWIDGETS
#define EXPORT_OPT_MANTIDQT_MANTIDWIDGETS DLLExport
#else
#define EXPORT_OPT_MANTIDQT_MANTIDWIDGETS DLLImport
#endif /* IN_MANTIDQT_MANTIDWIDGETS */

#endif //MANTIDQT_MANTIDWIDGETS_DLLOPTION_H_
