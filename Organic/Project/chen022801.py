import __init__
from Applications.Openeye import Openeye as oe
from Format.Pdf import Pdf



   
path="D:\\BOrganic\\Desktop\\cs01\\Test24022301\\report\\report.pdf"
pages=Pdf.get_pages(path)
lines=pages[1].splitlines()
print(lines[25:-1])

 
scores=oe.get_score(path)
T_scores=oe.trans_score(scores)

print(T_scores)
 
for i,Total in enumerate(T_scores["Total Score"]):
    sum=float(T_scores["Shape"][i])+float(T_scores["Hydrogen Bond"][i])+float(T_scores["Protein Desolvation"][i])+float(T_scores["Ligand Desolvation"][i])
    if sum-float(Total)<=0.01:
        print("第"+str(i+1)+'个数据正确')
    else:
        print("第"+str(i+1)+'个数据错误')
        print("第"+str(i+1)+'个Shape数据:' +str(T_scores["Shape"][i]))
        print("第"+str(i+1)+'个Hydrogen Bond数据:' +str(T_scores["Hydrogen Bond"][i]))
        print("第"+str(i+1)+'个Protein Desolvation数据:' +str(T_scores["Protein Desolvation"][i]))
        print("第"+str(i+1)+'个Ligand Desolvation数据:' +str(T_scores["Ligand Desolvation"][i]))
        print("第"+str(i+1)+'个Total Score数据:' +str(Total))
  