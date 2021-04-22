sys.setrecursionlimit(10000)

STEPSFRAGSIM=["0","10","20","30","40","50","100","150","200","300","400","500","600"]

rule all:
    input: expand("simulations/gen_{steps}_n1000_l200.fa",steps=STEPSFRAGSIM)

for p in STEPSFRAGSIM:
    rule:
        input: "simulations/gen_{param}.fa".format(param=p)
        output: "simulations/gen_{param}_n1000_l200.fa".format(param=p)
        shell: "/home/incerta/jana/Software/gargammel/src/fragSim  -n 1000 -l 200 {input} | {output}"


