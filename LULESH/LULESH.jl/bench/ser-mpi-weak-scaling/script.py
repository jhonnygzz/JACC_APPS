#!/usr/bin/python3.10
import os
import pathlib
scriptdir = pathlib.Path(__file__).parent.resolve()

def printfun(rank, blocklist,itercount):
  os.chdir(str(scriptdir)+"/../../")
  for s in blocklist:
    for mode in ["","--enzyme"]:
      os.system("./mpiexecjl  -bind-to socket --project -np {}  julia --project examples/benchmark.jl -s  --mpi {} {} {} > ser-mpi{}_{}_{}.txt".format(rank,mode,s,itercount, mode,rank,s))
      os.system("mv *.txt bench/ser-mpi-weak-scaling/")
  os.chdir(scriptdir)


itercount=10
printfun(1, [96],itercount)
printfun(8, [96],itercount)
printfun(27, [96],itercount)
printfun(64, [96],itercount)







# #!/usr/bin/python3.10
# import os
# import pathlib
# scriptdir = pathlib.Path(__file__).parent.resolve()

# def printfun(rank, blocklist,itercount):
#   os.chdir(str(scriptdir)+"/../../")
#   for s in blocklist:
#     # for mode in ["","--enzyme"]:
#     for mode in [""]:
#       # os.system("julia -p 4 examples/benchmark.jl -s --mpi {} {} {} > ser-mpi{}_{}_{}.txt".format(mode, s, itercount, mode, rank, s))
#       os.system("./mpiexecjl  -bind-to socket --project -np {}  julia --project -p 4 examples/benchmark.jl -s  --mpi {} {} {} > ser-mpi{}_{}_{}.txt".format(rank,mode,s,itercount, mode,rank,s))
#       os.system("mv *.txt bench/ser-mpi-weak-scaling/")
#   os.chdir(scriptdir)


# itercount=1
# printfun(1, [96],itercount)
# printfun(8, [96],itercount)
# printfun(27, [96],itercount)
# printfun(64, [96],itercount)

# # import os
# # import pathlib
# # import subprocess

# # scriptdir = pathlib.Path(__file__).parent.resolve()

# # def printfun(rank, blocklist, itercount):
# #     os.chdir(scriptdir / "../../")
# #     for s in blocklist:
# #         for mode in ["", "--enzyme"]:
# #             cmd = [
# #                 "./mpiexecjl", "-bind-to", "socket", "--project", "-np", str(rank),
# #                 "julia", "--project", "examples/benchmark.jl", "-s", "--mpi", mode, str(s), str(itercount)
# #             ]
# #             output_file = f"ser-mpi{mode}_{rank}_{s}.txt"
# #             with open(output_file, "w") as f:
# #                 print(f"Running command: {' '.join(cmd)}")
# #                 result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
# #                 if result.returncode != 0:
# #                     print(f"Command failed with return code {result.returncode}")
# #             subprocess.run(["mv", "*.txt", "bench/ser-mpi-weak-scaling/"])
# #     os.chdir(scriptdir)

# # itercount = 1
# # printfun(1, [96], itercount)
# # printfun(8, [96], itercount)
# # printfun(27, [96], itercount)
# # printfun(64, [96], itercount)