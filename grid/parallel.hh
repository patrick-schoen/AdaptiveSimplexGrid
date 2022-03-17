/** 
@file parallel.hh
@brief this header file includes all required headers to perform parallel mesh refinement and parallel communication between elements, faces, edges and nodes.
*/

#include "adaption/refine_parallel.hh"
#include "parallel/decomposition.hh"
#include "parallel/decomposition2D.hh"
#include "parallel/globalNodeId.hh"
#include "parallel/metisDecomposition.hh"
#include "parallel/mpiCommunication.hh"
#include "parallel/sendMesh.hh"
#include "parallel/meshPartition.hh"
#include "parallel/spaceFillingCurve.hh"