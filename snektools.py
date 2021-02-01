def read_lines(filepath):
    lines = open(filepath).readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
    return lines


def list_to_tab(lol, fout, rwa='w'):
    """write lists of lists to tab separated files
    rwa: can be write w or append a"""
    with open(fout, rwa) as outfile:
        for i in lol:
            line = str(i[0])
            for j in range(1, len(i)):
                line = line + '\t' + str(i[j])
            line = line + '\n'
            outfile.write(line)


def pretty_print(lol):
    """pretty print list of lists as tab separated"""
    for i in lol:
        line = str(i[0])
        for j in range(1, len(i)):
            line = line + '\t' + str(i[j])
        print(line)


def file_to_list(fin, out='int'):
    """from a file path which is supposed to be a list of int numbers,
    extract all the numbers and return a list of ints
    takes care of empty lines at the end if present"""
    fhandle = open(fin, 'r')
    items = fhandle.readlines()
    fhandle.close()

    outlist = list()

    for i in items:
        if i == '\n':
            continue
        else:
            if out == 'int':
                outlist.append(int(i.strip()))
            elif out == 'str':
                outlist.append(i.strip())

    # this chunk of code here doesn't work if there are empty lines in between the file instead of just at the end
    # for i in range(len(items)):
    #     if items[i] == '\n':
    #         items.pop(i)
    #     else:
    #         if out == 'int':
    #             items[i] = int(items[i].strip())
    #         elif out == 'str':
    #             items[i] = items[i].strip()q

    return outlist
