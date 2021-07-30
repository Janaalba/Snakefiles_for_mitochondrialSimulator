sys.setrecursionlimit(10000)

MAXVAL=600
STEPS10=range(10, MAXVAL, 10)
STEPSMSA=["0","10","20","30","40","50","100","150","200","300","400","500","600"]
LENGTHS=["30","35","40","45","50","55","60","65","75","85","100","150","200","400"]
RATES=["0.00010","0.00015","0.00020","0.00025","0.0003","0.00035","0.0004","0.00045","0.0005","0.0006","0.0007","0.0008","0.0009","0.001","0.003","0.005","0.007","0.009","0.01","0.03","0.05","0.07","0.09","0.1","0.3","0.5","0.7","0.9"]
DAMAGE=["dhigh","dmid","single","none"]
NUMFRAGS=1000000
FRAGNUM=175000

rule all:
    input:
        expand("simulations/gen_{maxval}.fa",maxval=MAXVAL),
        expand("simulations/gen_{steps}.nw",steps=STEPS10)


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

rule runmsacat:
    input:  expand("simulations/gen_{idx}.fa",idx=STEPSMSA)
    output: "simulations/all.fasta"
    shell: "cat {input} > {output}"

#run mafft
rule runmafft:
    input:  "simulations/all.fasta"
    output:  "simulations/all_mafft.fasta"
    shell: "/home/incerta/gabriel/Software/mafft-7.475-without-extensions/core/mafft --auto {input} > {output}"

#run prank
rule runprank:
    input:  "simulations/all.fasta"
    output: "simulations/all_prank.best.fas"
    params:
        outprefix="simulations/all_prank"
    shell: "/home/incerta/gabriel/Software/prank-msa/src/prank -d={input} -showall -o={params.outprefix} -DNA"


rule faidx:
    input: "simulations/gen_{step}.fa"
    output: "simulations/gen_{step}.fa.fai"
    shell: "/home/incerta/gabriel/Software/samtools/samtools faidx {input}"

