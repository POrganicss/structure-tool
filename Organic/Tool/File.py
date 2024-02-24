import os
from openbabel import openbabel
import glob
#提供数据的读取和保存
class File:
    #保存文件；filename:路径+文件名; content:输出的文本内容
    def tofile(filename='test', *args):
        
        if isinstance(args[0], list):
            if filename.lower().endswith(('.xlsx', '.xls')):
                import xlsxwriter as xw
                workbook = xw.Workbook(filename)
                worksheet1 = workbook.add_worksheet("Sheet1")

                if len(args) > 1 or (isinstance(args[0][0], list) and len(args[0]) > 1):
                    if isinstance(args[0][0], list) and len(args[0]) > 1:
                        args = args[0]

                    maxnum = max(len(a) for a in args)
                    for a in args:
                        while len(a) < maxnum:
                            a.append(0)

                    for i, row_data in enumerate(zip(*args), start=1):
                        worksheet1.write_row(f'A{i}', row_data)

                elif len(args) == 1:
                    for i, data in enumerate(args[0], start=1):
                        worksheet1.write_row(f'A{i}', [data])

                workbook.close()

            else:
                content = '\n'.join(args[0])
                try:
                    with open(filename, 'w') as file:
                        file.write(content)
                except Exception as e:
                    print(e)

        elif isinstance(args[0], openbabel.OBMol):  # 如果args[0]是OBMol数据
            conv = openbabel.OBConversion()
            output_format = filename.split('.')[-1] if '.' in filename else 'cdxml'
            conv.SetOutFormat(output_format)
            cdxml_string = conv.WriteString(args[0])

            if cdxml_string is not None:
                with open(filename, "w") as file:
                    file.write(cdxml_string)
            else:
                raise ValueError("Failed to convert molecule")
        elif isinstance(args[0], str):
            content = args[0]
            with open(filename,'w') as file:
                file.write(content)
                file.close
        
    def getdata(filename):
        filename = os.path.normpath(filename)
        with open(filename,'r',encoding='utf-8')as r:
            content=r.read()
            r.close
        return content

    def toexcel(filename, *args):
        import xlsxwriter as xw
        filename = os.path.normpath(filename)
        if not filename.lower().endswith(('.xlsx', '.xls')):
            filename=filename+'.xlsx'
        workbook = xw.Workbook(filename)
        worksheet1 = workbook.add_worksheet("Sheet1")

        if len(args) > 1 or (isinstance(args[0][0], list) and len(args[0]) > 1):
            if isinstance(args[0][0], list) and len(args[0]) > 1:
                args = args[0]

            maxnum = max(len(a) for a in args)
            for a in args:
                while len(a) < maxnum:
                    a.append(0)

            for i, row_data in enumerate(zip(*args), start=1):
                worksheet1.write_row(f'A{i}', row_data)

        elif len(args) == 1:
            for i, data in enumerate(args[0], start=1):
                worksheet1.write_row(f'A{i}', [data])

        workbook.close()
            
    def getexcel(filename,column_index=0):
        import pandas as pd
        filename =filename.replace("\\", "\\\\")
        
        # 使用pandas读取Excel文件
        df = pd.read_excel(filename,header=None)

        # 将指定列的数据存储在列表中
        #column_data = df[column_name].tolist()
        column_data = df.iloc[:, column_index].tolist() 
         
        return column_data
    
    def copy(Patha,Pathb,*args):
        import shutil
        # 确保目标文件夹存在
        if not os.path.exists(Pathb):
            os.makedirs(Pathb)
            
        for filename in args:
            source_file = os.path.join(Patha, filename)
            destination_file = os.path.join(Pathb, filename)
            shutil.copy(source_file, destination_file)
        
    def cp(Path,nameA,NameB,*args):
        File.copy(os.path.join(Path,nameA),os.path.join(Path,NameB),*args)
        
    def create_files(Path,*args):
        for filename in args:
            file = os.path.join(Path, filename)
            if not os.path.exists(file):
                os.makedirs(file)

    def getfile_name(path,_format):
        _format="."+_format
        file_list = []
        for file in os.listdir(path):
            if file.lower().endswith(_format.lower()):  # 使用lower()函数将文件名和目标格式都转换为小写，然后进行比较
                file_path = os.path.join(path, file)
                new_file_name = os.path.join(path, os.path.splitext(file)[0] + _format.lower())  # 将实际的文件名后缀改为小写
                os.rename(file_path, new_file_name)  # 修改实际文件名的后缀为小写
                file_list.append(os.path.splitext(os.path.basename(new_file_name))[0])  # 获取不带后缀的文件名并添加到列表中
        return file_list

