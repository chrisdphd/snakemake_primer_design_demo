SAMPLES = ['RPL23a', 'KRAS_ex1_G12', '1.4kb_nr_KRAS']

rule all:
        input:
                expand('output/{sample}.primers_fwd', sample=SAMPLES), expand('output/{sample}.primers_rev', sample=SAMPLES)

rule generate_fwd:
	input:
		"input/{sample}.seq"
	output:
		"output/{sample}.primers_fwd"
	shell:
		"""./DNA_substringer_20240223.py -f {input} > {output}"""


rule generate_rev:
        input:
                "input/{sample}.seq"
        output:
                "output/{sample}.primers_rev"
        shell:
                """./DNA_substringer_20210327.py -f {input} -r > {output}"""


#rule collate:
#        input:
#                in1 = expand('{sample}.primers_fwd', sample=SAMPLES),
#                in2 = expand('{sample}.primers_rev', sample=SAMPLES)
#        output:
#                '{sample}.primers'
#        shell:
#                'cat {sample}.primers_fwd {sample}.primers_rev > {output}'
