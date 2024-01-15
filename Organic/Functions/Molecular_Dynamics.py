import os
import sys
path = os.getcwd()
path = path.replace("\\", "/")+'/Organic'
sys.path.append(path+'/Format')
sys.path.append(path+'/Tool')
sys.path.append(path+'/Applications')
sys.path.append(path+'/Functions')
sys.path.append(path+'/resources')

from Datatransmission import Datatransmission as Dt

class Molecular_Dynamics:
    YPath='C:\\Users\\10282\\Documents\\MD'
    def run():
        code=''
        code = 'cd ' + Molecular_Dynamics.YPath+'\n'
        code=code+'dir'
        print(Dt.Command(code))

Molecular_Dynamics.run()
        