import sys

class Log_xyz:
    def logtoxyz(filename):
        elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 
                    7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 
                    13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 
                    18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 
                    23: 'V',24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 
                    28: 'Ni', 29: 'Cu', 30: 'Zn', 33: 'As', 35: 'Br', 
                    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 49: 'In', 
                    53: 'I', 55: 'Cs', 74: 'W', 75: 'Re', 77: 'Ir', 
                    78: 'Pt', 79: 'Au', 82: 'Pb'}
        open_old = open(filename, 'r')
        rline = open_old.readlines()
        starts = []
        ends = []
        for i in reversed(range(len(rline))):
            if 'Standard orientation:' in rline[i]:
                starts.append(i)
            elif ' Cite this work as' in rline[i]:
                break
        starts.sort()
        for j in range(1,len(starts)+1):  
            for m in range(starts[j-1] + 5, len(rline)):
                if '---' in rline[m]:
                    ends.append(m)
                    break
        xyzs=[] 
        for num in range(len(starts)):  
            xyz = []
            for line in rline[starts[num] + 5: ends[num]]:
                element =elements[int(line.split()[1])]
                xyz.append(element + line[30:-1] + '\n')  
            xyzs.append(xyz)
        open_old.close()
        return xyzs
    
    def logtoxyzfile(filename):
        xyzs=Log_xyz.logtoxyz(filename)
        xyz=xyzs[len(xyzs)-1]
        with open(filename[:-4]+'.xyz','w') as w:
            w.write(''.join(xyz))
            w.close
        
filename=sys.argv[1]
Log_xyz.logtoxyzfile(filename)
