#!/usr/bin/python3.10
import os

def printfun(rank, blocklist,itercount,mode):
     for s in blocklist:
        os.system("echo -n {},{}, >> results.txt".format(rank,s))
        os.system("grep  Elapsed ser-mpi{}_{}_{}.txt >> results.txt".format(mode,rank,s))
     os.system("sed -i \"s/Elapsed time         = //g\" results.txt")
     os.system("sed -i \"s/ (s)//g\" results.txt")

os.system("rm results.txt")
itercount=100
for mode in ["","--enzyme"]:
  printfun(1, [96],itercount,mode)
  printfun(8, [96],itercount,mode)
  printfun(27, [96],itercount,mode)
  printfun(64, [96],itercount,mode)
