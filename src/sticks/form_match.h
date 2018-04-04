#ifndef __form_match_h
#define __form_match_h

#include <boost/shared_ptr.hpp>
#include "stick.h"
#include "utils/vector3d.h"

//  MICHAEL I. SADOWSKI 201

namespace PRODART
{
namespace STICK
{
class form_match;

typedef boost::shared_ptr<form_match> form_match_shared_ptr;
typedef std::vector<form_match_shared_ptr> form_match_shared_ptr_vector;

class form_match
	{
	friend boost::shared_ptr<form_match>new_form_match();	 // needs a different constructor
	
	public:

	

	inline double get_pack_score()
		{
		return pack_score;
		}
	inline double get_rmsd()
		{
		return rmsd;
		}

	bundle_shared_ptr_vector get_sticks_a()
		{
		return sticks_a;
		}
	bundle_shared_ptr_vector get_sticks_b()
		{
		return sticks_b;
		}
	// TODO overlaod some operators
	
	UTILS::vector3d move_to(UTILS::vector3d point)
		{
		// move all the sticks then return new centroid
		return loc;
		}
	private:
	       // this holds pretty much everything	
	       
		bundle_shared_pointer_vector sticks_a;
		bundle_shared_pointer_vector sticks_b;

		double pack_score;
		double rmsd;

		// centroid
		
		UTILS::vector3d loc;

		inline form_match()
			{
			}

	};
boost::shared_ptr<form_match>new_form_match()
	{
	boost::shared_ptr<form>nformmatch(new form_match);
	return nformmatch;
	}
}
}
#endif
