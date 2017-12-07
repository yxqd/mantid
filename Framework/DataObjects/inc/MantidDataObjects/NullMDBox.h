#ifndef MANTID_DATAOBJECTS_NULLMDBOX_H_
#define MANTID_DATAOBJECTS_NULLMDBOX_H_

#include "MantidDataObjects/DllConfig.h"
#include "MantidDataObjects/MDBoxBase.h"


namespace Mantid {
namespace DataObjects {

TMDE_CLASS
class MANTID_DATAOBJECTS_DLL NullMDBox : public MDBoxBase<MDE, nd>  {
public:
  NullMDBox(Mantid::API::BoxController *const splitter, const uint32_t depth,
  const std::vector<Mantid::Geometry::MDDimensionExtents<coord_t>> & extentsVector,
  const size_t boxID = UNDEF_SIZET) :
    MDBoxBase<MDE, nd>(splitter, depth, boxID, extentsVector) {
    if (this->m_BoxController->getNDims() != nd)
      throw std::invalid_argument(
        "NullMDBox::ctor(): controller passed has the wrong number of dimensions.");
  }


  // ----------------------------- ISaveable Methods
  // ------------------------------------------------------
  Kernel::ISaveable *getISaveable() override {
    throw std::runtime_error("Not implemented");
  }

  Kernel::ISaveable *getISaveable() const override {
    throw std::runtime_error("Not implemented");
  }

  void setFileBacked(const uint64_t /*fileLocation*/, const size_t /*fileSize*/,
                     const bool /*markSaved*/) override {
    throw std::runtime_error("Not implemented");
  }

  void setFileBacked() override {
    throw std::runtime_error("Not implemented");
  }

  void clearFileBacked(bool ) override {
    throw std::runtime_error("Not implemented");
  }

  //-----------------------------------------------------------------------------------------------
  void saveAt(API::IBoxControllerIO *const /* */,
              uint64_t /*position*/) const override {
    throw std::runtime_error("Not implemented");
  }

  void loadAndAddFrom(API::IBoxControllerIO *const /* */, uint64_t /*position*/,
                      size_t /* Size */) override {
    throw std::runtime_error("Not implemented");
  }

  void reserveMemoryForLoad(uint64_t /* Size */) override {
    throw std::runtime_error("Not implemented");
  }

  void clearDataFromMemory() override {
    throw std::runtime_error("Not implemented");
  }

  void clear() override {
    throw std::runtime_error("Not implemented");
  }

  uint64_t getNPoints() const override;
  size_t getDataInMemorySize() const override {
    throw std::runtime_error("Not implemented");
  }

  uint64_t getTotalDataSize() const override {
    throw std::runtime_error("Not implemented");
  }

  size_t getNumDims() const override {
    throw std::runtime_error("Not implemented");
  }

  size_t getNumMDBoxes() const override {
    throw std::runtime_error("Not implemented");
  }

  size_t getNumChildren() const override { return 0; }

  bool isBox() const override {
    throw std::runtime_error("Not implemented");
  }

  API::IMDNode *getChild(size_t /*index*/) override {
    throw std::runtime_error("MDBox does not have children.");
  }

  /// Sets the children from a vector of children
  void setChildren(const std::vector<API::IMDNode *> & /*boxes*/,
                   const size_t /*indexStart*/,
                   const size_t /*indexEnd*/) override {
    throw std::runtime_error("MDBox cannot have children.");
  }


  std::vector<MDE> *getEventsCopy() override{
    throw std::runtime_error("Not implemented");
  }

  void getEventsData(std::vector<coord_t> &,
                     size_t &) const override {
    throw std::runtime_error("Not implemented");
  }

  void setEventsData(const std::vector<coord_t> &) override {
    throw std::runtime_error("Not implemented");
  }

  size_t addEvent(const MDE &) override {
    throw std::runtime_error("Not implemented");
  }

  size_t addEventUnsafe(const MDE &) override {
    throw std::runtime_error("Not implemented");
  }

  size_t addEvents(const std::vector<MDE> &) override {
    throw std::runtime_error("Not implemented");
  }

  size_t addEventsUnsafe(const std::vector<MDE> &) override{
    throw std::runtime_error("Not implemented");
  }

  void buildAndAddEvent(const signal_t , const signal_t ,
                        const std::vector<coord_t> &, uint16_t ,
                        uint32_t ) override{
    throw std::runtime_error("Not implemented");
  }

  void buildAndAddEventUnsafe(const signal_t , const signal_t ,
                              const std::vector<coord_t> &,
                              uint16_t , uint32_t ) override{
    throw std::runtime_error("Not implemented");
  }

  size_t buildAndAddEvents(const std::vector<signal_t> &,
                           const std::vector<coord_t> &,
                           const std::vector<uint16_t> &,
                           const std::vector<uint32_t> &) override {
    throw std::runtime_error("Not implemented");
  }

  //---------------------------------------------------------------------------------------------------------------------------------
  void centerpointBin(MDBin<MDE, nd> &, bool *) const override {
    throw std::runtime_error("Not implemented");
  }

  void
  generalBin(MDBin<MDE, nd> &bin,
             Mantid::Geometry::MDImplicitFunction &function) const override{
    throw std::runtime_error("Not implemented");
  }

  void splitAllIfNeeded(Mantid::Kernel::ThreadScheduler * /*ts*/ = nullptr)
  override { /* Do nothing with a box default. */
    throw std::runtime_error("Not implemented");
  }

  //---------------------------------------------------------------------------------------------------------------------------------
  /** Recalculate signal and various averages dependent on signal and the signal
   * coordinates */
  void refreshCache(Kernel::ThreadScheduler * /*ts*/ = nullptr) override{
    throw std::runtime_error("Not implemented");
  }

  void calculateCentroid(coord_t *) const override {
    throw std::runtime_error("Not implemented");
  }

  void calculateCentroid(coord_t *, const int ) const override{
    throw std::runtime_error("Not implemented");
  }

  coord_t *getCentroid() const override{
    throw std::runtime_error("Not implemented");
  }

  void integrateSphere(
    Mantid::API::CoordTransform &radiusTransform, coord_t radiusSquared,
    signal_t &signal, signal_t &errorSquared,
    const coord_t innerRadiusSquared = 0.0,
    const bool useOnePercentBackgroundCorrection = true) const override{
    throw std::runtime_error("Not implemented");
  }

  void centroidSphere(Mantid::API::CoordTransform &radiusTransform,
                      coord_t radiusSquared, coord_t *centroid,
                      signal_t &signal) const override{
    throw std::runtime_error("Not implemented");
  }

  void integrateCylinder(Mantid::API::CoordTransform &,
                         const coord_t , const coord_t ,
                         signal_t &, signal_t &,
                         std::vector<signal_t> &) const override{
    throw std::runtime_error("Not implemented");
  }

  //------------------------------------------------------------------------------------------------------------------------------------
  void getBoxes(std::vector<MDBoxBase<MDE, nd> *> &, size_t /*maxDepth*/,
  bool /*leafOnly*/){
    throw std::runtime_error("Not implemented");
  }
  void getBoxes(std::vector<API::IMDNode *> &boxes, size_t /*maxDepth*/,
                bool /*leafOnly*/) override{
    throw std::runtime_error("Not implemented");
  }

  void getBoxes(std::vector<MDBoxBase<MDE, nd> *> &, size_t ,
  bool leafOnly, Mantid::Geometry::MDImplicitFunction *){
    throw std::runtime_error("Not implemented");
  }

  void getBoxes(std::vector<API::IMDNode *> &, size_t ,
                bool ,
                Mantid::Geometry::MDImplicitFunction *) override{
    throw std::runtime_error("Not implemented");
  }
  //------------------------------------------------------------------------------------------------------------------------------------
  void transformDimensions(std::vector<double> &,
                           std::vector<double> &) override{
    throw std::runtime_error("Not implemented");
  }
  //------------------------------------------------------------------------------------------------------------------------------------
  /* Getter to determine if masking is applied.
  @return true if masking is applied.   */
  bool getIsMasked() const override {
    throw std::runtime_error("Not implemented");
  }
  /// Setter for masking the box
  void mask() override{
    throw std::runtime_error("Not implemented");
  }
  /// Setter for unmasking the box
  void unmask() override{
    throw std::runtime_error("Not implemented");
  }

private:
  /// private default copy constructor as the only correct constructor is the
  /// one with the boxController;
  NullMDBox(const NullMDBox &);
};

} // namespace DataObjects
} // namespace Mantid

#endif /* MANTID_DATAOBJECTS_NULLMDBOX_H_ */