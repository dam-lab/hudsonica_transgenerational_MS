import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""description""")
parser.add_argument('--input', required=True, dest='Input', type=str, help="Input file in varscan format")
parser.add_argument('--output', required=True, dest='Output', type=str, help="Output file in sync format.")

args = parser.parse_args()


InputFile=open(args.Input,"r")
OutputFile=open(args.Output,"w")

base_dict = {'A':0, 'T':1, 'C':2, 'G':3}

###
# define conversion function

def convert_to_sync(ref_allele, var_allele, varscan_input):
    ct_ref=varscan_input.split(":")[2]
    ct_var=varscan_input.split(":")[3]
    sync_line = "0:0:0:0:0:0".split(":")
    sync_line[base_dict.get(ref)]=ct_ref
    sync_line[base_dict.get(var)]=ct_var
    return ":".join(sync_line)


firstline = True

for line in InputFile:

    if firstline:    #skip first line
        firstline = False
        continue

    line=line.split('\n')[0]
    cols=line.split('\t')

    ref=cols[2]
    var=cols[3]

    line_out = []

    for pop_input in cols[10:]:
        if len(line_out) == 0:
            line_out = convert_to_sync(ref, var, pop_input)

        else:
            line_out = line_out + "\t" + convert_to_sync(ref, var, pop_input)

    line_out = "\t".join(cols[0:3]) + "\t" + line_out

    OutputFile.write(line_out+'\n')

OutputFile.close()
InputFile.close()
