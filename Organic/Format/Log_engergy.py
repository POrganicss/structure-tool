
import sys

class Log_energy:
    def getenergy(filename):
        energy={}
        open_old = open(filename, 'r')
        lines = open_old.readlines()
        for line in lines:
            if 'SCF Done:' in line:
                energy['SP']=line.split('=')[1].split(' ')[2].replace(' ','')
                
            if 'Sum of electronic and zero-point Energies' in line:
                energy['ZPE']=line.split('=')[1].replace(' ','')
            if 'Zero-point correction' in line:
                energy['dZPE']=line.split('=')[1].replace(' ','').replace('(Hartree/Particle)','')
                   
            if 'Sum of electronic and thermal Energies' in line:
                energy['U']=line.split('=')[1].replace(' ','')
            if 'Thermal correction to Energy' in line:
                energy['dU']=line.split('=')[1].replace(' ','')
                    
            if 'Sum of electronic and thermal Enthalpies' in line:
                energy['H']=line.split('=')[1].replace(' ','')
            if 'Thermal correction to Enthalpy' in line:
                energy['dH']=line.split('=')[1].replace(' ','')
                    
            if 'Sum of electronic and thermal Free Energies' in line:
                energy['G']=line.split('=')[1].replace(' ','')
            if 'Thermal correction to Gibbs Free Energy' in line:
                energy['dG']=line.split('=')[1].replace(' ','')
            if 'Error termination' in line:
                return False
        open_old.close()
        return energy
    
    def logtoenergy(filename):
        energy=Log_energy.getenergy(filename)
        if energy==False:
            with open(filename[:-4]+'.eng','w') as f:
                f.write('Error termination')
        else:
            with open(filename[:-4]+'.eng','w') as f:
                f.write('Gibbs free energy:'+energy['G'])
                f.write('Enthalpies:'+energy['H'])
                f.write('Energies:'+energy['U'])
                f.write('zero-point Energies:'+energy['ZPE'])
                f.write('SP:'+energy['SP'])
                f.write('\n')
                f.write('\n')
                f.write('dG:'+energy['dG'])
                f.write('dH:'+energy['dH'])
                f.write('dU:'+energy['dU'])
                f.write('dZPE:'+energy['dZPE'])
                f.close
        
filename=sys.argv[1]
Log_energy.logtoenergy(filename)
