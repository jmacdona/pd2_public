#ifndef __stickmatch_h
#define __stickmatch_h

#include<boost/shared_ptr.hpp>
#include "pose/chain.h"
#include "stick.h"
#include "bundle.h"
#include "utils/vector3d.h"


//  MICHAEL I. SADOWSKI 2010

using namespace PRODART::UTILS;

namespace PRODART
{
namespace STICK
{
typedef std::vector< vector3d >  vector3d_vector;
boost::shared_ptr<bundle>match_sticks(boost::shared_ptr<POSE::chain>chain, double gap_pen);
double linescore(boost::shared_ptr<POSE::chain>chain, int start, int end, double &dens,int print);
}
}

#endif
