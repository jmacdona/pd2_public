import sys
if '-h' in sys.argv or '--help' in sys.argv or '--help' in sys.argv:
	print '''
loop relax script

arguments <reference pdb> <pdb to relax> <start resnum> <end resnum>

James MacDonald, 2011 '''
	sys.exit(0)

ref_pdb = sys.argv[1]
input_pdb = sys.argv[2]
start_res = int(sys.argv[3])
end_res = int(sys.argv[4])
output_pdb = input_pdb



from rosetta import *
rosetta.init()
pose = Pose(input_pdb)
ref_pose = Pose(ref_pdb)
start_pose = Pose()
start_pose.assign(pose)

start_res = pose.pdb_info().pdb2pose(" ",start_res) + 1
end_res = pose.pdb_info().pdb2pose(" ",end_res) - 1

print "start_res: " + str(start_res)
print "end_res: " + str(end_res)

mm4060 = MoveMap()
mm4060.set_bb_true_range(start_res,end_res)
loop = Loop(start_res,end_res)
loop.auto_choose_cutpoint(pose)
set_single_loop_fold_tree(pose, loop)
minmover = MinMover()
minmover.movemap(mm4060)
scorefxn = create_score_function('standard')
scorefxn.apply_patch_from_file("score12")
std_scorefxn = create_score_function('standard')
std_scorefxn.apply_patch_from_file("score12")
add_cutpoint_variants(pose)
scorefxn.set_weight(chainbreak,100.0)
scorefxn.set_weight(coordinate_constraint,1.0)
minmover.score_function(scorefxn)
minmover.min_type("dfpmin")


task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
task_pack.fix_everything()
for x in range(start_res, end_res+1):
	if pose.residue(x).name() != "CYD":
		print x
		task_pack.set_pack_residue(x, True)
print task_pack
packmover = PackRotamersMover(scorefxn,task_pack)

print std_scorefxn.show(pose)

fa_rep_wt = std_scorefxn.get_weight(fa_rep)

best_score = std_scorefxn(pose)
best_pose = Pose()
best_pose.assign(pose)

print "Best score: " + str(best_score)



for x in range(1, 15+1):
	scorefxn.set_weight(fa_rep,0.02*fa_rep_wt)
	scorefxn.set_weight(chainbreak,100.0)
	packmover.score_function(scorefxn)
	packmover.apply(pose)
	print scorefxn(pose)
	minmover.tolerance(0.01)
	minmover.score_function(scorefxn)
	minmover.apply(pose)
	print scorefxn(pose)

        scorefxn.set_weight(fa_rep,0.250*fa_rep_wt)
        scorefxn.set_weight(chainbreak,100.0)
        packmover.score_function(scorefxn)
        packmover.apply(pose)
        print scorefxn(pose)
	minmover.tolerance(0.01)
        minmover.score_function(scorefxn)
        minmover.apply(pose)
        print scorefxn(pose)

        scorefxn.set_weight(fa_rep,0.550*fa_rep_wt)
        scorefxn.set_weight(chainbreak,100.0)
        packmover.score_function(scorefxn)
        packmover.apply(pose)
        print scorefxn(pose)
        minmover.tolerance(0.01)
        minmover.score_function(scorefxn)
        minmover.apply(pose)
        print scorefxn(pose)

        scorefxn.set_weight(fa_rep,1.0*fa_rep_wt)
        scorefxn.set_weight(chainbreak,200.0)
        packmover.score_function(scorefxn)
        packmover.apply(pose)
        print scorefxn(pose)
        minmover.tolerance(0.00001)
        minmover.score_function(scorefxn)
        minmover.apply(pose)
        print scorefxn(pose)

	this_score = std_scorefxn(pose)
	if this_score <= best_score:
		best_pose.assign(pose)
		best_score = this_score
		print "Best score: " + str(best_score)
	else:
		pose.assign(best_pose)

print "Best score: " + str(best_score)

print std_scorefxn.show(best_pose)
energy = std_scorefxn(best_pose)
loops = Loops()
loops.add_loop(loop)

best_pose.dump_pdb(output_pdb)

out_f = open(output_pdb, 'a')

Lrms = loop_rmsd(start_pose, ref_pose, loops, False)
out_ln = "RMSD_start_to_ref: " + str(Lrms)
out_f.writelines( out_ln + "\n" )
print out_ln

Lrms = loop_rmsd(best_pose, start_pose, loops, False)
out_ln = "RMSD_best_to_start: " + str(Lrms)
out_f.writelines( out_ln + "\n"  )
print out_ln

Lrms = loop_rmsd(best_pose, ref_pose, loops, False)
out_ln = "RMSD_best_to_ref: " + str(Lrms) + " " + str(energy)
out_f.writelines( out_ln + "\n"  )
print out_ln

out_ln = "ENERGY: " + str(energy)
out_f.writelines( out_ln + "\n"  )
print out_ln


out_f.close()






