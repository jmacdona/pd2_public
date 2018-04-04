#include "stickmatch.h"
#include "pose/chain.h"
#include "utils/line_fit.h"

using std::cerr;
using std::cout;
using std::endl;

using namespace PRODART::UTILS;

namespace PRODART
{
namespace STICK
{
boost::shared_ptr<bundle> match_sticks(boost::shared_ptr<POSE::chain>chain, double gap_pen)
	{
	int stick_win=21;
	int stick_min=3;
	double min_dens=0.5;		 
	double min_lsc=0.1;

	double **mat=new double *[chain->length()], **dens=new double *[chain->length()];
	int **xat=new int *[chain->length()];
	int **yat=new int *[chain->length()];

	for(int k=0;k<chain->length();k++)
		{
		xat[k]=new int[stick_win];
		yat[k]=new int[stick_win];
		mat[k]=new double[stick_win];
		dens[k]=new double[stick_win];

		for(int l=0;l<stick_win;l++)
			{
			xat[k][l]=-1;
			yat[k][l]=-1;
			mat[k][l]=0;
			}	
		}

	// fill the matrix

	for(int k=0;k<chain->get_last_internal_residue_index()-(stick_min+1);k++)
		{
		double gap=gap_pen*stick_min;

		for(int w=1;w<stick_win;w++)
			{
			int pos2=w+k;

			if(pos2>chain->get_last_internal_residue_index()-1)
				break;

			double lsc=0.0;

			if(w>=stick_min)
				{
				int print=0;
				lsc=linescore(chain,k,pos2,dens[k][w],print);//-gap;
				gap+=gap_pen;
				}

			mat[k][w]=lsc*gap;


			if(lsc<min_lsc || dens[k][w]<min_dens || (dens[k][w]>2.2 && w<4))
				{
				dens[k][w]=0.0;
				mat[k][w]=0.0*gap;
				}

			cout << "MAT[" << k <<"]["<<w<<"]"<< " " <<mat[k][w] << " lsc: " <<lsc << " dens " << dens[k][w] << endl;
			}
		}

	// cumulate scoring

	double overall_best_score=-100000.0;
	int bx=0,by=0;
	for(int k=0;k<chain->get_last_internal_residue_index()-(stick_min+1);k++)
		{
		for(int l=1;l<stick_win;l++)
			{
			double curr_best_score=-100000.0;
			int lx=0, ly=0;

			// search back 
			for(int m=1;m<stick_win;m++)
				{
				int pos3=k-(m+1);
				if(pos3<0) break;

			//	cout << " k " << k << " l " << l << " pos3 " << pos3 << " m " << m << " mat[pos3][m] " << mat[pos3][m] << " mat[k][l] " << mat[k][l] << endl;
					
				if(mat[pos3][m]+mat[k][l] > curr_best_score)
					{
					lx=pos3;
					ly=m;
					curr_best_score=mat[pos3][m]+mat[k][l];
					}
				}
			if(curr_best_score != -100000.0)
				mat[k][l]=curr_best_score;
			xat[k][l]=lx;
			yat[k][l]=ly;

			//	cout << "k " << k << " l " << l << " curr best score: " << curr_best_score << " lx " << lx << " ly " << ly << endl;
			if(curr_best_score>overall_best_score)
				{
				bx=k;
				by=l;
			//	cout<<"overall best was " << overall_best_score << " now " << curr_best_score << " at " << bx << " " << by <<endl;
				overall_best_score=curr_best_score;
				}
			}
		}

	cout << "Best score : " << overall_best_score << endl;
	// traceback and fill stick array
	
	stick_shared_ptr_vector sticks;

	while(bx!=-1 && by!=-1)
		{
		// dens is only set if we have matched a stick here
		cout << bx << " " << by << " " << mat[bx][by] << " ";
		if(dens[bx][by]>0.1) 
			{
			UTILS::vector3d v3dstart=chain->get_ca_pos(bx);
			UTILS::vector3d v3dend  =chain->get_ca_pos(bx+by);
			boost::shared_ptr<stick>s=new_stick(v3dstart, v3dend, bx, bx+by, dens[bx][by],chain);	
			sticks.insert(sticks.begin(),s);
			cout << "Stick " << sticks.size() << " start: " << bx << " end " << bx+by << " dens " <<dens[bx][by];
			}
		cout <<endl;
		int tx=xat[bx][by];
		int ty=yat[bx][by];
		bx=tx;
		by=ty;
		}	

	// make bundle
	
	boost::shared_ptr<bundle> ret=new_bundle(sticks, chain);

	// free up
	
	for(int k=0;k<chain->length();k++)
		delete [] mat[k];
	delete [] mat;

	return ret;
	}

double linescore(boost::shared_ptr<POSE::chain>chain,int start, int end, double &dens, int print)
	{
	// use linefit_3d to get the basis 
	
	vector3d_vector v;
	
	for(int k=start;k<end+1;k++)
		v.push_back(chain->get_ca_pos(k));

	vector3d a,b;

	// store eigenvalues here
	
	double ev[3];
	double dvar=0.0;

	line_fit3d(v,a,b,dens,dvar,ev);


	// now calculate the mysterious "linefit" constant
	
	double mlfc=exp(-5.0*dvar);

	if(print)
		cout << "\tdens = " << dens << " dvar = " << dvar << " mlfc = " << mlfc << " evs: " << sqrt(ev[0]) << " " << sqrt(ev[1]) << " " << sqrt(ev[2]) << endl;

	// transform to the "equivalent ellipsoid"

	double ratoi=2000.0;

	if((ev[1]+ev[2])>0)
		ratoi=sqrt(ev[0])/(sqrt(ev[1])+sqrt(ev[2]));
		//ratoi=ev[0]/(ev[1]+ev[2]);

//	cout << "\tratoi: " << ratoi << endl;

	mlfc*=ratoi;

//	cout << "\tfinal: " << mlfc << endl;

	return mlfc;	
	}
}
}
