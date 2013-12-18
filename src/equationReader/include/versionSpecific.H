// If your OpenFOAM version number is 1.7.x or less, comment out the next line:
// #define ThisIsFoamVersion2

#ifdef ThisIsFoamVersion2
#   ifndef MathConstantScope
#      define MathConstantScope constant::mathematical
#   endif
#else
#   ifndef MathConstantScope
#      define MathConstantScope mathematicalConstant
#   endif
#endif
