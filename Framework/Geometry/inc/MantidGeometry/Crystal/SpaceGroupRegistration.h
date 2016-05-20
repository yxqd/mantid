#ifndef MANTID_GEOMETRY_SPACEGROUPREGISTRATION_H_
#define MANTID_GEOMETRY_SPACEGROUPREGISTRATION_H_

#include "MantidGeometry/Crystal/SpaceGroupFactory.h"
#include "MantidKernel/RegistrationHelper.h"

namespace Mantid {
namespace Geometry {

// Registration helpers are in their own namespace so that names can be readable
namespace SpaceGroupRegistration {

/// Return a copy of the input string without spaces in it.
std::string copy_remove_spaces(const std::string &str) {
  std::string result;

  std::remove_copy_if(str.cbegin(), str.cend(), std::back_inserter(result),
                      [](const char &ch) { return ch == ' '; });

  return result;
}

// Using this name fits the "grammar" for registration better.
using Subscribe = Kernel::RegistrationHelper;

/**
 * @class SpaceGroupSubscriber
 *
 * This class assists registering space groups into
 * the factory at compile time. This is done using the comma operator, which
 * evaluates all statements, discards their value and returns only the last
 * value. It is used like this:
 *
 *  using namespace Mantid::Geometry::SpaceGroupRegistration;
 *
 *  GeneratedSpaceGroupSubscriber
 *      GeneratedSpaceGroup(SpaceGroupFactory::Instance());
 *
 *  Subscribe groups(
 *        (GeneratedSpaceGroup(....),
 *         GeneratedSpaceGroup(...))
 *        );
 *
 * The additional parantheses are required because the comma operator has lowest
 * precedence. Please note that very long "comma-chains" increase compile times.
 */
template <typename GeneratorType> class SpaceGroupSubscriber {
public:
  SpaceGroupSubscriber(SpaceGroupFactoryImpl &factory) : m_factory(factory) {}
  virtual ~SpaceGroupSubscriber() = default;

  /// Subscribe the specified space group to the factory. Returns 1 so that it
  /// works with RegistrationHelper
  int operator()(size_t number, const std::string &hmSymbol,
                 const std::string &generators) const {
    subscribeToFactory(number, hmSymbol, generators);

    return 1;
  }

  /// Subscribe the space group and create aliases.
  int operator()(size_t number, const std::string &hmSymbol,
                 const std::string &generators,
                 const std::string &aliases) const {
    subscribeToFactory(number, hmSymbol, generators);
    registerAliasesToFactory(hmSymbol, aliases);

    return 1;
  }

protected:
  /// Performs the actual subscription and registeres a default alias
  /// (Hermann-Mauguin symbol without spaces).
  virtual void subscribeToFactory(size_t number, const std::string &hmSymbol,
                                  const std::string &generators) const {
    m_factory.subscribeUsingGenerator<GeneratorType>(number, hmSymbol,
                                                     generators);
    registerAliasesToFactory(hmSymbol, copy_remove_spaces(hmSymbol));
  }

  /// Registers the specified aliases for the symbol to the factory.
  void registerAliasesToFactory(const std::string &hmSymbol,
                                const std::string &aliases) const {
    m_factory.registerAliases(hmSymbol,
                              aliases + "," + copy_remove_spaces(aliases));
  }

  SpaceGroupFactoryImpl &m_factory;
  SpaceGroupSubscriber() = default;
};

template <typename GeneratorType>
class OrthorhombicSpaceGroupRegistrationHelper
    : public SpaceGroupSubscriber<GeneratorType> {
public:
  using SpaceGroupSubscriber<GeneratorType>::SpaceGroupSubscriber;

protected:
  /// Subscribes default aliases for each generated space group.
  void subscribeToFactory(size_t number, const std::string &hmSymbol,
                          const std::string &generators) const override {
    std::vector<std::string> additionalSymbols =
        this->m_factory.template subscribeOrthorhombicSpaceGroup<GeneratorType>(
            number, hmSymbol, generators);

    // Base symbol alias
    this->registerAliasesToFactory(hmSymbol, copy_remove_spaces(hmSymbol));

    // All generated symbols get an alias too
    for (const auto &symbol : additionalSymbols) {
      this->registerAliasesToFactory(symbol, copy_remove_spaces(symbol));
    }
  }
};

using GeneratedSpaceGroupSubscriber =
    SpaceGroupSubscriber<AlgorithmicSpaceGroupGenerator>;

using TabulatedSpaceGroupSubscriber =
    SpaceGroupSubscriber<TabulatedSpaceGroupGenerator>;

using TransformedSpaceGroupSubscriber =
    SpaceGroupSubscriber<TransformationSpaceGroupGenerator>;

using OrthorhombicSpaceGroupSubscriber =
    OrthorhombicSpaceGroupRegistrationHelper<AlgorithmicSpaceGroupGenerator>;

using TransformedOrthorhombicSpaceGroupSubscriber =
    OrthorhombicSpaceGroupRegistrationHelper<TransformationSpaceGroupGenerator>;
}
}
}
#endif /* MANTID_GEOMETRY_SPACEGROUPREGISTRATION_H_ */
