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
        expand("simulations/gen_{steps}.nw",steps=STEPS10),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}.fa.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}.fa.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_adpt.fa",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_s1.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_s2.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_o.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_o_s{rate}.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_o_r1_s{rate}.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
        expand("simulations/gen_{steps}_n{nfrag}_l{fraglen}_d{dam}_o_r2_s{rate}.fq.gz",steps=STEPSMSA,nfrag=NUMFRAGS,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
        expand("simulations/numtS_n{fragn}_l{fraglen}.fa.gz",fragn=FRAGNUM,fraglen=LENGTHS),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}.fa.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_adpt.fa",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s1.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s2.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_s{rate}.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r1_s{rate}.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
        expand("simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r2_s{rate}.fq.gz",fragn=FRAGNUM,fraglen=LENGTHS,dam=DAMAGE,rate=RATES),
#        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o"
#        "simulations/all_prank.best.fas"
        "simulations/all_mafft.fasta"


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

rule fragsim:
    input:
        input_1="simulations/gen_{step}.fa",
        input_2="simulations/gen_{step}.fa.fai"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/fragSim  -n {wildcards.nfrags} -l {wildcards.fraglen} --circ generation_{wildcards.step} {input.input_1} | gzip > {output}"


rule deamsim:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}.fa.gz"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/deamSim  -matfile {wildcards.dam}  {input} | gzip > {output}"

rule adptsim:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}.fa.gz"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_adpt.fa"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/adptSim   -l 125 -artp {output}  {input}"

rule art:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_adpt.fa"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s1.fq",
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s2.fq"
    wildcard_constraints:
        fraglen="\d+"
    params:
        out_prefix="simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s"
    shell:
        "/home/incerta/jana/Software/gargammel/art_src_MountRainier_Linux/art_illumina -ss HS25 -amp -na -p -l 125 -c 1   -i {input} -o {params.out_prefix}"

rule adptsimz:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_adpt.fa"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_adpt.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "gzip  {input}"

rule artz1:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s1.fq"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s1.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
       "gzip {input}"

rule artz2:
    input:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s2.fq"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s2.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "gzip {input}"

rule trimmed:
    input:
        input_1="simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s1.fq.gz",
        input_2="simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_s2.fq.gz"
    output:
        "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o.fq.gz"
    params:
        out_prefix="simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/gabriel/Software/leeHom/src/leeHom --ancientdna -fq1 {input.input_1} -fq2 {input.input_2} -fqo {params.out_prefix}"

rule subsamp_mt:
    input: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o.fq.gz"
    output: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"

rule subsamp_r1:
    input: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o_r1.fq.gz"
    output: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o_r1_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"

rule subsamp_r2:
    input: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o_r2.fq.gz"
    output: "simulations/gen_{step}_n{nfrags}_l{fraglen}_d{dam}_o_r2_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"
        

#same process for numtS

rule faidx_numtS:
    input: "simulations/numtS.fasta"
    output: "simulations/numtS.fasta.fai"
    shell: "/home/incerta/gabriel/Software/samtools/samtools faidx {input}"

rule fragsim_numtS:
    input:
        input_1="simulations/numtS.fasta",
        input_2="simulations/numtS.fasta.fai"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/fragSim  -n {wildcards.fragn} -l {wildcards.fraglen} {input.input_1} | gzip > {output}"


rule deamsim_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}.fa.gz"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/deamSim  -matfile {wildcards.dam}  {input} | gzip > {output}"


rule adptsim_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}.fa.gz"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_adpt.fa"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/jana/Software/gargammel/src/adptSim   -l 125 -artp {output}  {input}"

rule art_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_adpt.fa"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s1.fq",
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s2.fq"
    wildcard_constraints:
        fraglen="\d+"
    params:
        out_prefix="simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s"
    shell:
        "/home/incerta/jana/Software/gargammel/art_src_MountRainier_Linux/art_illumina -ss HS25 -amp -na -p -l 125 -c 1   -i {input} -o {params.out_prefix}"


rule adptsimz_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_adpt.fa"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_adpt.fa.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "gzip  {input}"


rule artz1_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s1.fq"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s1.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
       "gzip {input}"

rule artz2_numtS:
    input:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s2.fq"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s2.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "gzip {input}"


rule trimmed_numtS:
    input:
        input_1="simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s1.fq.gz",
        input_2="simulations/numtS_n{fragn}_l{fraglen}_d{dam}_s2.fq.gz"
    output:
        "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o.fq.gz"
    params:
        out_prefix="simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o"
    wildcard_constraints:
        fraglen="\d+"
    shell:
        "/home/incerta/gabriel/Software/leeHom/src/leeHom --ancientdna -fq1 {input.input_1} -fq2 {input.input_2} -fqo {params.out_prefix}"


rule subsamp_numtS:
    input: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o.fq.gz"
    output: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"

rule subsamp_numtS_r1:
    input: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r1.fq.gz"
    output: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r1_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"

rule subsamp_numtS_r2:
    input: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r2.fq.gz"
    output: "simulations/numtS_n{fragn}_l{fraglen}_d{dam}_o_r2_s{rate}.fq.gz"
    wildcard_constraints:
        fraglen="\d+"
    shell: "/home/incerta/gabriel/Software/seqtk/seqtk sample -s 321 {input} {wildcards.rate} |gzip > {output}"
