import sys

class Log_success:
    def getsuccess(filename):
        open_old = open(filename, 'r')
        lines = open_old.readlines()
        for line in lines:
            if 'Error termination' in line:
                return False
            if 'Normal termination' in line:
                return True
        return '-1'
    
    def logtosuccess(filename):
        success=Log_success.getsuccess(filename)
        if success==False:
            with open(filename[:-4]+'.res','w') as f:
                f.write('Error termination')
        elif success==True:
            with open(filename[:-4]+'.res','w') as f:
                f.write('Normal termination')
        else:
            with open(filename[:-4]+'.res','w') as f:
                f.write('Other situations')
filename=sys.argv[1]
Log_success.logtosuccess(filename)
