sys.setrecursionlimit(10000)

MAXVAL=600
STEPS10=range(10, MAXVAL, 10)

rule all:
    input:
        expand("simulations/gen_{maxval}.fa",maxval=MAXVAL),
        expand("simulations/gen_{steps}.nw",steps=STEPS10),


rule simulations_mt:
    input: "simulations/gen_0.fa"
    output: expand("simulations/gen_{idx}.fa", idx=range(1, MAXVAL+1))
    shell: "/home/incerta/gabriel/Software/mitochondrialSimulator/mitochondrialSimulator.py -g {MAXVAL} -o simulations/gen --name=generation /home/incerta/gabriel/Software/mitochondrialSimulator/rCRS_0_005.conf\
 {input}"


for p in STEPS10:
    rule:
        input:
            input_1="simulations/gen_0.fa",
            input_2="simulations/gen_{param}.fa".format(param=p)
        output: "simulations/gen_{param}.nw".format(param=p)
        params:
            str(p),
            out_prefix="simulations/gen_{param}".format(param=p)
        shell: "/home/people/stud007/incerta/Software/EMBOSS-6.5.7/emboss/needle -gapopen 10.0 -gapextend 0.5 -outfile {output} -asequence {input.input_1} -bsequence {input.input_2}"

