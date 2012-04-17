# please run "sh makedata0.sh" first to generate all the information file
# function combineseedscript(flag, classnum, randseed, allflag, maximgperfold, sdaflag, weightflag, cidx, gidx, nbrs, classifierid)

alias qarray='perl /export/home/install/bin/r2/qarray.pl'

qarray -m "combineseedscript(%%)" "41,20,25,1,140,1,1,11,5,0,1"
