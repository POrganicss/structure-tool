
#提供数据的读取和保存
import os


class File:
    #保存文件；filename:路径+文件名; content:输出的文本内容
    def tofile(filename='test', *args):
        filename = os.path.normpath(filename)
        from openbabel import openbabel
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
    
    def getdata(filename):
        filename = os.path.normpath(filename)
        with open(filename,'r',encoding='utf-8')as r:
            content=r.read()
            r.close
        return content

    def toexcel(filename, *args):
        filename = os.path.normpath(filename)
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
            
    def getexcel(filename):
        filename = os.path.normpath(filename)
        import pandas as pd

        # 读取Excel文件
        excel_file_path = filename  # 将路径替换为你的Excel文件路径
        column_name = 'Smile'  # 将列名替换为你要读取的列名

        # 使用pandas读取Excel文件
        df = pd.read_excel(excel_file_path)

        # 将指定列的数据存储在列表中
        column_data = df[column_name].tolist()
            
        return column_data
    
    def totemp(content,filename='test'):
        print()
        

