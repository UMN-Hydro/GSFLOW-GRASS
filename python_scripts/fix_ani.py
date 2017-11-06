import re

base_dir = '/home/awickert/Downloads/Shullcas_spinup/outputs/PRMS_GSFLOW/'
hru_filename = 'Shullcas.ani.nhru'
segment_filename = 'Shullcas.ani.nsegment'

for filename in [hru_filename, segment_filename]:
    infile = file(base_dir + filename, 'r')
    outfile = file(base_dir + filename + '.corrected', 'w')
    for line in infile:
        if line[:2] == '  ':
            p = re.compile("\-[0-9]{2}\-")
            for m in p.finditer(line):
                if m.start():
                    break
            _start = m.start() - 4 # space for the year
            outfile.write(line[_start:])
        else:
            outfile.write(line)


