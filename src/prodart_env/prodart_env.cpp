//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * prodart_env.cpp
 *
 *  Created on: 25 Feb 2010
 *      Author: jmacdona
 */
#include "prodart_env.h"

#ifdef USING_SCONS_BUILDER 
#include "prodart_env_var.h"
#endif


using PRODART::UTILS::options_manager;
using boost::any_cast;
using boost::any;
using boost::trim;
using boost::split;
using boost::is_any_of;
using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::vector;

using namespace boost::filesystem;
using namespace boost;

typedef vector<string> string_vector;


namespace PRODART {
namespace ENV {

namespace {
prodart_env *instance = NULL;
boost::once_flag env_once_flag = BOOST_ONCE_INIT;
}


class cust_timer : public boost::timer {
public:
#ifdef PD2_MPI_MODE
	/*
	int rank_of_mpi_process;
	cust_timer() : boost::timer(){
		mpi::communicator world;
		rank_of_mpi_process = world.rank();
	}
	*/
#endif

	~cust_timer(){
		if (instance){
			if (ENV::get_option_value<bool>("print_elapsed_time")){
//#ifdef PD2_MPI_MODE
				//cerr << "\n" << "prodart_env: elapsed time (seconds): " << this->elapsed() << "\trank:" << rank_of_mpi_process << "\n" << endl;
//#else
				cerr << "\n" << "prodart_env: elapsed time (seconds): " << this->elapsed() << "\n" << endl;
				//cerr << "prodart_env: max elapsed time (seconds): " << this->elapsed_max() << endl;
//#endif
			}
		}
	}
};

namespace {
cust_timer creation_timer;
}





prodart_env::prodart_env(){
	command_line_parsed = false;
	register_common_options();
	rand_num_incr = 0;

}

//destructor never really called
prodart_env::~prodart_env(){
	delete instance;
	instance = 0;
}

void prodart_env::Init(){
	if (!instance){
		//TODO might be dodgy with threads also output messy
		signal(SIGSEGV, backtrace_handler);
		signal(SIGINT, interrupt_handler);
		instance = new prodart_env;
		creation_timer.restart();
	}
	//return instance;
}


prodart_env* prodart_env::Instance(){
	/*
	if (!instance){
		instance = new prodart_env;
		creation_timer.restart();
	}
	*/
	boost::call_once(&prodart_env::Init, env_once_flag);

	return instance;
}

/* inlined
options_manager& prodart_env::get_options_manager(){
	return this->options;
}
const options_manager& prodart_env::get_options_manager() const{
	return this->options;
}
*/

boost::mutex cmdln_mutex;
bool prodart_env::is_command_line_parsed() const{
	boost::mutex::scoped_lock lock(cmdln_mutex);
	return command_line_parsed;
}

void prodart_env::register_common_options(){


	ifstream f("/dev/urandom");
	MTRand::uint32 seed;
	f.read(reinterpret_cast<char*>(&seed), sizeof(seed));

	options.register_option("pose:io:pdb:i", string(""), "input PDB file path" );
	options.register_option("pose:io:pdb:r", string(""), "input reference PDB file path" );
	options.register_option("pose:io:pdb:o", string(""), "output PDB file path" );
	options.register_option("pose:io:pdb:list:l", string(""), "input PDB file list path" );
	options.register_option("pose:io:pdb:traj:t", string(""), "input PDB trajectory file path" );
	options.register_option("database", string("./database/"), "database directory path" );
	options.register_option("database:database_paths", string("./database/paths.dat"), "database directory path settings file", true );
	options.register_option("random_seed", (MTRand::uint32)seed, "random seed initialised by /dev/urandom by default" );
	//4567856

	options.register_option("prodart_version", string(_PRODART_VERSION_), "version" );
	options.register_option("hg_version_number", string(HG_REV), "Mercurial version number from hg id" , true);
	options.register_option("svn_version_number", string(SVN_REV), "Subversion version number from svnversion", true );
	options.register_option("compile_date", string(_COMPILE_DATE_), "compile date" );

	options.register_option("output_root", string("pd2_temp"), "output file root name to prepend misc output files" );

	options.register_option("print_elapsed_time", bool(true), "print elapsed time on destruction", true );

	options.register_option("sim:ca:ca_nb_cutoff", double(21.01), "CA non-bonded cutoff", true ); // 3/10/13: changed from 13.5 to accommodated GO potential
	options.register_option("sim:ca:ca_cell_margin", double(100), "CA sim cell margin", true );
	options.register_option("sim:ca:move_margin", double(4.5), "CA non-bonded cutoff", true ); // 3/10/13: changed from 4.5 to accommodated GO potential

	options.register_option("sim:ca:no_nb_cutoff", double(9.01), "N' O' non-bonded cutoff", true ); // 18/10/13: was 13.5 but should be sim:ca:move_margin + pseudo_hb_dist_cutoff
	options.register_option("sim:ca:no_cell_margin", double(100), "N' O' sim cell margin", true );
	options.register_option("sim:ca:pseudo_hb_dist_cutoff", double(4.5), "N' O' hbond cutoff", true );

	options.register_option("sim:bb:ho_nb_cutoff", double(13.5), "H O non-bonded cutoff", true );
	options.register_option("sim:bb:ho_cell_margin", double(100), "H O sim cell margin", true );
	options.register_option("sim:bb:hb_dist_cutoff", double(2.6), "H O hbond cutoff", true );
	options.register_option("sim:bb:move_margin", double(2.0), "H O hbond cutoff", true );

	options.register_option("sim:bb:all_bb_nb_cutoff", double(10.0), "all bb atoms non-bonded cutoff", true );
	options.register_option("sim:bb:all_bb_cell_margin", double(100), "all bb atoms sim cell margin", true );

	options.register_option("potential:bb:bb_frag3_mq:count_cutoff", int(2), "frag3_mq count cutoff", true );

	options.register_option("as_geom_pot:motif_pdb", string(""), "active site motif PDB file", true );
	options.register_option("as_geom_pot:mapping", string(""), "fixed mapping file", true );

	options.register_option("loop_model:loop_mask", string(""), "loop mask e.g. 00111000", true );

	options.register_option("bb_motif_anneal_protocol:start_as_pot_wt", double(0.2), "", true);
	options.register_option("bb_motif_anneal_protocol:final_as_pot_wt", double(7.0), "", true);
	options.register_option("bb_motif_anneal_protocol:start_min_as_pot_wt", double(0.5), "", true);
	options.register_option("bb_motif_anneal_protocol:final_min_as_pot_wt", double(10.0), "", true);
	options.register_option("bb_motif_anneal_protocol:start_ca_as_pot_wt", double(1.0), "", true);
	options.register_option("bb_motif_anneal_protocol:final_ca_as_pot_wt", double(10.0), "", true);

	options.register_option("rosetta_contraints:cst_ca_cutoff", double(14.0), "", true);
	options.register_option("rosetta_contraints:cst_all_atom_cutoff", double(14.0), "", true);

	options.register_option("minimiser:lin_min_tol", double(0.01), "", true); //0.01
	options.register_option("minimiser:tol", double(0.01), "", true);        //1e-3
	options.register_option("minimiser:step_size", double(0.00000000000001), "", true); //0.0001  //0.01

	options.register_option("restraints::rstfile", string(""), "restraint file", true );


	char *db_env = std::getenv("PRODART_DATABASE_PATH");
	if (db_env != NULL){
		path db_pth(db_env);
		path db_pth_fl(db_pth / "paths.dat");
		options.set_option_value("database", db_pth.string());
		options.set_option_value("database:database_paths", db_pth_fl.string() );
	}

}

boost::mutex new_MTR_mutex;
boost::shared_ptr<MTRand> prodart_env::get_new_random_num_gen(){

	boost::mutex::scoped_lock lock(new_MTR_mutex);

#ifndef PD2_MPI_MODE
	// non-MPI compile
	const MTRand::uint32 rn_seed = options.get_option_value<MTRand::uint32>("random_seed") + rand_num_incr++;
	boost::shared_ptr<MTRand> new_mtrand(new MTRand(rn_seed));
#else
	// if compiling in MPI mode to ensure different random seeds per process
	mpi::communicator world;
	const MTRand::uint32 rank = (MTRand::uint32)world.rank();
	const MTRand::uint32 rn_seed = options.get_option_value<MTRand::uint32>("random_seed") + rand_num_incr++ + rank;
	boost::shared_ptr<MTRand> new_mtrand(new MTRand(rn_seed));
#endif

	// run through first 1000
	for (int i = 0; i < 1000; i++){
		new_mtrand->randInt();
	}

#ifndef NDEBUG
		cerr << "prodart_env: get_new_random_num_gen: made new MTRand." << endl;
#endif

	return new_mtrand;

}

boost::mutex thread_MTR_mutex;
MTRand::MTRand_shared_ptr prodart_env::get_random_num_gen(boost::thread::id tid){

	boost::mutex::scoped_lock lock(thread_MTR_mutex);

	if (rand_gen_store.find(tid) != rand_gen_store.end()){
		return rand_gen_store[tid];
	}
	else {

#ifndef NDEBUG
		cerr << "prodart_env: get_random_num_gen: making new MTRand for thread id: " << tid << endl;
#endif


		boost::shared_ptr<MTRand> new_rg = get_new_random_num_gen();
		rand_gen_store[tid] = new_rg;
		return new_rg;
	}
}

bool prodart_env::get_database_path_option(int argc, char *argv[]){

	string_vector args(argc);

	for (int i = 0; i < argc; i++){
		args[i] = string(argv[i]);
		trim(args[i]);
	}

	for (int i = 1; i < argc; i++){
		if (args[i].compare("--database") == 0){
			if (i+1 < argc){
				cerr << "PRODART2: init: database path set: " << args[i+1] << endl;
				options.set_option_value("database", args[i+1]);
				path db_pth(args[i+1]);
				path db_pth_fl(db_pth / "paths.dat");
				options.set_option_value("database:database_paths", db_pth_fl.string() );
			}
			else {
				return false;
			}
		}
	}

	return true;
}

bool prodart_env::init(int argc, char *argv[]){

	//comment out later:
	//cerr << "prodart_env: initialising..." << endl;

	//const bool get_db_path_result =
	this->get_database_path_option(argc, argv);


	const bool paths_result = this->parse_paths_dat();
	if (paths_result == false){
		cerr << "\nPRODART2: init: database paths.dat parsing errors" << endl;
		return false;
	}


	const bool sclp_result = options.parse_cmd_line(argc, argv);
	if (sclp_result == false){
		cerr << "\nPRODART2: init: command line parsing errors" << endl;
		return false;
	}

	command_line_parsed = true;
	return true;
}


bool prodart_env::parse_paths_dat(){


	const bool hide_path_opts = true;

	path path_dat_path(options.get_option_value<string>("database:database_paths"));

	ifstream input(path_dat_path);

	if (!input.is_open()){
		cerr << "\nPRODART2: init: could not open database:database_paths file: " << path_dat_path.string() << endl;
		return false;
	}

	string lineStr;


	long length, lineNum = 0 ;


	string_vector SplitVec;
	while ( !input.eof() ) {
		getline(input, lineStr);
		string resStr;
		lineNum++;

		length = lineStr.length();

		//cout << endl << lineNum << " " << length << " ";

		if (length > 0) {
			split( SplitVec, lineStr, is_any_of("\t ") );
			if ( SplitVec[0].substr(0,1).compare("#") != 0
	                 && SplitVec[0].substr(0,1).compare("") != 0
	                 && SplitVec.size() >= 2 ){

				string default_prepend("database:path:");

				string paraName = SplitVec[0];
				trim(paraName);
				string paraVal = SplitVec[1];
				trim(paraVal);

				default_prepend.append(paraName);
				paraName = default_prepend;

				path db_pth(options.get_option_value<string>("database"));
				path file_path(db_pth / paraVal);

				options.register_option(paraName, file_path.string(), "database file path setting", hide_path_opts );

				//ifstream test_open(file_path);

				if (!exists(file_path)){
					cerr << "\nPRODART2: init: WARNING file does not exist: " << file_path.string() << endl;
				}
				else {
					//test_open.close();
				}

			}
		}
	}

	return true;
}


void backtrace_handler(int sig){
	  void *array[200];
	  size_t size;

	  // get void*'s for all entries on the stack
	  size = backtrace(array, 200);

	  // print out all the frames to stderr

	  fprintf(stderr, "\n\nError: segmentation fault: SIGSEGV: signal %d:\n", sig);
#ifdef PD2_MPI_MODE
	  mpi::communicator world;
	  fprintf(stderr, "Error: segmentation fault: SIGSEGV: MPI rank:%d:\n", world.rank());
#endif
	  backtrace_symbols_fd(array, size, STDERR_FILENO);
	  exit(1);
}

void interrupt_handler(int sig){
	  void *array[200];
	  size_t size;

	  // get void*'s for all entries on the stack
	  size = backtrace(array, 200);

	  // print out all the frames to stderr
	  fprintf(stderr, "\n\nError: interrupt detected: SIGINT: signal %d:\n", sig);
#ifdef PD2_MPI_MODE
	  mpi::communicator world;
	  fprintf(stderr, "Error: interrupt detected: SIGINT: MPI rank:%d:\n", world.rank());
#endif
	  backtrace_symbols_fd(array, size, STDERR_FILENO);
	  exit(1);
}


}
}


