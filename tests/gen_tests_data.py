from textwrap import wrap

all_aaa_file = "../data/all_aaa.fa"
with open(all_aaa_file, 'w') as of:
    of.write(">AllAs"+"\n")
    all_aaa_str = "A" * 80
    all_aaa_str = all_aaa_str + "\n"
    for x in range(40):
        of.write(all_aaa_str)
    of.write("A" * 24)
    of.write("\n")

all_aaa_file = "../data/repat_4a4t.fa"
with open(all_aaa_file, 'w') as of:
    of.write(">RepeatAT"+"\n")
    all_aaa_str = "AAAATTTT" * 10
    all_aaa_str = all_aaa_str + "\n"
    for x in range(40):
        of.write(all_aaa_str)
    of.write("AAAATTTT" * 3)
    of.write("\n")

all_aaa_file = "../data/all_aaa_small.fa"
with open(all_aaa_file, 'w') as of:
    of.write(">AllAsSmall"+"\n")
    all_aaa_str = "A" * 80
    all_aaa_str = all_aaa_str + "\n"
    for x in range(10):
        of.write(all_aaa_str)
    of.write("A" * 12)
    of.write("\n")

def repeat_input(input_file, header, nas, nts, nlength):
    build_str = ""
    while nlength > 0:
        if nas > 0:
            build_str = build_str + "A"*nas
        if nts > 0:
            build_str = build_str + "T"*nts
        nlength = nlength - (nas + nts)
    line_list = wrap(build_str, 80)
    out_str = "\n".join(line_list)
    with open(input_file, 'w') as of:
        of.write(header+"\n")
        of.write(out_str+"\n")

repeat_input("../data/repeat_8ats_2048.fa", ">Repat8AT2048", 8, 1, 2048000)
repeat_input("../data/repeat_12ats_2048.fa", ">Repeat12AT2048", 12, 1, 2048000)
repeat_input("../data/repeat_16ats_2048.fa", ">Repeat16AT2048", 16, 1, 2048000)
repeat_input("../data/repeat_20ats_2048.fa", ">Repeat20AT2048", 20, 1, 2048000)
repeat_input("../data/repeat_24ats_2048.fa", ">Repeat24AT2048",  24, 1, 2048000)
repeat_input("../data/repeat_32ats_2048.fa", ">Repeat32AT2048", 32, 1, 2048000)
repeat_input("../data/repeat_64ats_2048.fa", ">Repeat32AT2048", 64, 1, 2048000)
repeat_input("../data/repeat_128ats_2048.fa", ">Repeat32AT2048", 128, 1, 2048000)
repeat_input("../data/repeat_as_2048.fa", ">RepatAs2048", 32, 0, 2048000)
