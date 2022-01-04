/** 
@file simplexGrid.hh
@brief included by mesh.hh, this header file includes all required headers to perform seriell mesh refinement.
*/

//exclude standard header files for doxygen documentation (cond, endcond)
///@cond
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <math.h>
#include <stdio.h>
///@endcond

#include "array/array_lib.hh"
#include "auxiliary/auxiliary.hh"
#include "adaption/adaption_seriell.hh"