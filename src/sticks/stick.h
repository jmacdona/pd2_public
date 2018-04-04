#ifndef __stick_h
#define __stick_h

#include <boost/shared_ptr.hpp>
#include <vector>

#include "utils/vector3d.h"

//  MICHAEL I. SADOWSKI 201

namespace PRODART
{
namespace STICK
{
class stick;

typedef boost::shared_ptr<stick> stick_shared_ptr;
typedef std::vector<stick_shared_ptr> stick_shared_ptr_vector;

class stick
	{
	friend inline boost::shared_ptr<stick>new_stick();
	friend inline boost::shared_ptr<stick>new_stick(UTILS::vector3d &start, UTILS::vector3d &end, int seqstart, int seqend, double d, boost::shared_ptr<POSE::chain> &ch);
	protected: 

	UTILS::vector3d loc;
	double dens;
	int *endis;
	UTILS::vector3d *endvs;
	boost::shared_ptr<POSE::chain> _chain;

	public: 

	inline UTILS::vector3d get_loc()
		{
		return loc;
		}
	inline double get_dens()
		{
		return dens;
		}

	inline int start()
		{
		return endis[0];
		}

	inline int end()
		{
		return endis[1];
		}

	~stick()
		{
		delete [] endis;
		delete [] endvs;
		}

	private:

	stick()
		{
		endis=new int[2];
		endvs=new UTILS::vector3d[2];
		}
	stick(UTILS::vector3d &start, UTILS::vector3d &end,int seqstart, int seqend, double d, boost::shared_ptr<POSE::chain> &ch)
		{
		endis=new int[2];
		endvs=new UTILS::vector3d[2];

		endvs[0]=start;
		endvs[1]=end;
		
		endis[0]=seqstart;
		endis[1]=seqend;

		dens=d;

		_chain = ch;
		}
	};
boost::shared_ptr<stick>new_stick()
	{
 	boost::shared_ptr<stick> sp(new stick());
	return sp;
	}
boost::shared_ptr<stick>new_stick(UTILS::vector3d &start, UTILS::vector3d &end, int seqstart, int seqend, double d, boost::shared_ptr<POSE::chain> &ch)
	{
 	boost::shared_ptr<stick> sp(new stick(start, end, seqstart, seqend, d,ch));
	return sp;
	}
}
}

#endif
