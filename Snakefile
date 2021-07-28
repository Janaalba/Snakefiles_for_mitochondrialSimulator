sys.setrecursionlimit(10000)

MAXVAL=600

rule all:
    input:
        expand("simulations/gen_{maxval}.fa",maxval=MAXVAL)


rule simulations_mt:
    input: "simulations/gen_0.fa"
    output: expand("simulations/gen_{idx}.fa", idx=range(1, MAXVAL+1))
    shell: "/home/incerta/gabriel/Software/mitochondrialSimulator/mitochondrialSimulator.py -g {MAXVAL} -o simulations/gen --name=generation /home/incerta/gabriel/Software/mitochondrialSimulator/rCRS_0_005.conf\
 {input}"



