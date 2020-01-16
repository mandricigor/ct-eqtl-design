

import os
import stat
import sys
import single_cell_calc as scc


exname = sys.argv[1]

budget = int(sys.argv[2])


low_cells = 500
high_cells = 2750
low_individuals = 40
high_individuals = 120
experiments = scc.exp_design(budget, low_cells, high_cells, low_individuals, high_individuals, multi=8) # in single-cell-calc library prep is 0





#---igor4---
#budget = 30000
#low_cells = 500
#high_cells = 2750
#low_individuals = 40
#high_individuals = 120
#experiments = scc.exp_design_fixed_lane_capacity(budget, low_cells, high_cells, low_individuals, high_individuals, max_multi=16)



#---igor3---
#budget = 35000
#low_cells = 500
#high_cells = 2750
#low_individuals = 40
#high_individuals = 120
#experiments = scc.exp_design_fixed_lane_capacity(budget, low_cells, high_cells, low_individuals, high_individuals, max_multi=16)



#---igor2---
#budget = 35000
#low_cells = 500
#high_cells = 2750
#low_individuals = 40
#high_individuals = 120
#experiments = scc.exp_design(budget, low_cells, high_cells, low_individuals, high_individuals, diff_cell=250, multi=8)




filenames = []

for indiv in experiments:
    for total_cells, singlets, multiplets, singlets_reads, multiplets_reads in experiments[indiv]:
        filename = "run.%s.%s.%s.%s.%s.%s.%s.paper.sh" % (exname, indiv, total_cells, singlets, multiplets, singlets_reads, multiplets_reads)
        filenames.append(filename)
        with open(filename, "w") as f:
            f.write("python3.7 super_igor_down_10x_super_multiplex.py %s %s %s %s %s %s\n" % (exname, indiv, singlets, multiplets, singlets_reads, multiplets_reads))
            st = os.stat(filename)
            os.chmod(filename, st.st_mode | stat.S_IEXEC)

with open("myFunc_%s.sh" % budget, "w") as f:
    for x, y in enumerate(filenames):
        f.write("if [ $1 == %s ];then ./%s ;fi\n" % (x, y))

st = os.stat("myFunc_%s.sh" % budget)
os.chmod("myFunc_%s.sh" % budget, st.st_mode | stat.S_IEXEC)


with open("myFuncFastWrapper_%s.sh" % budget, "w") as f:
    f.write("#!/bin/bash\n# myFuncFastWrapper_%s.sh\necho $SGE_TASK_ID\n./myFunc_%s.sh $SGE_TASK_ID" % (budget, budget))

st = os.stat("myFuncFastWrapper_%s.sh" % budget)
os.chmod("myFuncFastWrapper_%s.sh" % budget, st.st_mode | stat.S_IEXEC)


if not os.path.exists(exname):
    os.makedirs(exname)


#print ("Recommended launching:")
print ("qsub -cwd -V -N %s -l h_data=80G,time=24:00:00 -t 1-%s:1 myFuncFastWrapper_%s.sh" % (exname, len(filenames), budget))


