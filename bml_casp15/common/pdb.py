import os, sys

def parse_pdb_row(row, param):
    result = ""
    if param == "anum":
        result = row[6:6 + 5]
    elif param == "aname":
        result = row[12:12 + 4]
    elif param == "altloc":
        result = row[16:16 + 1]
    elif param == "rname":
        result = row[17:17 + 3]
    elif param == "rnum":
        result = row[22:22 + 5]
    elif param == "chain":
        result = row[21:21 + 1]
    elif param == "x":
        result = row[30:30 + 8]
    elif param == "y":
        result = row[38:38 + 8]
    elif param == "z":
        result = row[46:46 + 8]
    else:
        die(f"Invalid row[{row}] or parameter[{param}]")

    return result.lstrip().rstrip()


def extract_pdb(pdb, newpdb, start, end):
    if start > end:
        die(f"wrong index <start:{start}, end:{end}>\n")

    contents = open(pdb, "r").readlines()

    new_PDBlines = []

    for line in contents:

        if line.find('ATOM') < 0:
            continue

        this_rnum = int(parse_pdb_row(line, "rnum"))

        if this_rnum >= start and this_rnum <= end:
            new_PDBlines += [line]

    # (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
    resCounter = 0
    atomCounter = 0
    prevrNum = "XX"
    with open(newpdb, "w") as f:

        for line in new_PDBlines:

            if line.find('ATOM') < 0:
                continue

            this_rnum = parse_pdb_row(line, "rnum")

            if prevrNum != this_rnum:
                prevrNum = this_rnum
                resCounter = resCounter + 1

            atomCounter = atomCounter + 1

            rnum_string = format(str(resCounter), '>4')
            anum_string = format(str(resCounter), '>5')

            row = line[0:6] + anum_string + line[11:11 + 5] + " " + \
                  line[17:17 + 3] + " " + " " + rnum_string + " " + line[27:]

            f.write(row)

        f.write("END\n")