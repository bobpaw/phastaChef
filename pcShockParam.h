#ifndef PC_SHOCKPARAM_H
#define PC_SHOCKPARAM_H

namespace pc {
  /**
   * Scoped const double "enum" values for Shock_Param field.
   */
  class ShockParam {
  public:
    static constexpr double NONE = 0; ///< No detected shock.
    static constexpr double SHOCK = 1; ///< Detected shock.
    static constexpr double FRAGMENT = 2; ///< Defragmented element.
    static constexpr double EXTEND = 3; ///< Shock extension.

  private:
    // Static only class has private constructor.
    ShockParam();
  };
} // namespace pc

#endif // PC_SHOCKPARAM_H
