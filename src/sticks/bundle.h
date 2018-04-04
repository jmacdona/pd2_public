#ifndef __bundle_h
#define __bundle_h

#include<boost/shared_ptr.hpp>
#include<vector>
#include "pose/chain.h"
#include "pose/pose.h"
#include "pose/atom.h"
#include "stick.h"
//#include "sheet.h"
#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<fstream>

//  MICHAEL I. SADOWSKI 201

using namespace std;

using std::ofstream;

namespace PRODART
{
namespace STICK
{
class bundle;

typedef boost::shared_ptr<bundle> bundle_shared_ptr;
typedef std::vector<bundle_shared_ptr> bundle_shared_ptr_vector;

class bundle
	{
	friend inline boost::shared_ptr<bundle> new_bundle(STICK::stick_shared_ptr_vector sticks, boost::shared_ptr<POSE::chain>chain);
	//friend class sheet;

	public:
	
	stick_shared_ptr stick(int n)	
		{
		return _sticks[n];
		}



	void write_annotated_pdb(std::ofstream &pdbout)
		{
		int i,cur=0;

		for(i=0;i<=_chain->get_last_internal_residue_index();i++)
			_chain->get_ca(i)->set_b_factor(0.0);
			
		for(i=0;i<=_chain->get_last_internal_residue_index();i++)
			{
			if(i>=_sticks[cur]->start())
				{
				if(i<=_sticks[cur]->end())
					_chain->get_ca(i)->set_b_factor(_sticks[cur]->get_dens());
				else
					++cur;
				}
			if(cur==(int)_sticks.size()) break;
			}

		_chain->get_pose()->suppress_remark_output();
		//_chain->get_pose()->CAonly=true;
		_chain->get_pose()->outputPdb(pdbout, true);
		}

	private:
	
	bundle(STICK::stick_shared_ptr_vector &sticks, boost::shared_ptr<POSE::chain> &chain)
		{
		_sticks=sticks;
		_chain=chain;
		}
	STICK::stick_shared_ptr_vector _sticks;
	boost::shared_ptr<POSE::chain>_chain;
	};


inline boost::shared_ptr<bundle> new_bundle(stick_shared_ptr_vector sticks, boost::shared_ptr<POSE::chain>chain)
	{
	boost::shared_ptr<bundle>nbundle(new bundle(sticks,chain));
	return nbundle;
	}

}
}
#endif
