import os
import psutil

#提供数据的读取和保存
class File:
    
    # -----------------保存和读取文件-----------------#
    def save(path,*args):
        match path.split('.')[-1].lower():  # 根据文件扩展名判断
            case 'xlsx' | 'xls':
                File.savexcel(path, *args)
            case 'cdxml':
                File.savecdxml(path, args[0])
            case _:
                list_separator=", "
                line_separator="\n"
                prefix="{"
                postfix="}"
                match args:
                    case [element] if isinstance(element, str):
                        cleaned_data =element.replace("\n", " ").strip()
                        File.savedata(path, cleaned_data)

                    case [element] if isinstance(element, list):
                        formatted_data = f"{prefix}{list_separator.join(map(str, element))}{postfix}"
                        cleaned_data =formatted_data.replace("\n", " ").strip()
                        File.savedata(path, cleaned_data)

                    case elements if all(isinstance(el, str) for el in elements):
                        cleaned_elements = [el.replace("\n", " ").strip() for el in elements]
                        data = line_separator.join(cleaned_elements)
                        File.savedata(path, data)

                    case elements if all(isinstance(el, list) for el in elements):
                        formatted_data = line_separator.join(
                            f"{prefix}{list_separator.join(map(str, el))}{postfix}" for el in elements
                        )
                        cleaned_data =formatted_data.replace("\n", " ").strip()
                        File.savedata(path, cleaned_data)

                    case _:
                        print("Unsupported argument type.")
                        return
                           
    def read(path,format_="list"):
        if os.path.isdir(path):
            return File.getname(path, format_)
        else:
            match format_:
                case "list":
                    file_extension = path.split('.')[-1].lower()
                    if file_extension in ['xlsx', 'xls']:
                        return File.getexcel(path)
                    elif file_extension == 'pdf':
                        return File.getpdf(path)
                    else:
                        with open(path, 'r') as file:
                            return file.readlines()

                case "str":
                    match path.split('.')[-1].lower():
                        case "pdf":
                            return '\n'.join(File.getpdf(path))
                        
                        case _:
                            with open(path, 'r') as file:
                                return file.read()

                case _:
                    return "未知的格式"
          
          
    # --------------------保存文件-------------------#           
    def savexcel(path,*args):
        import xlsxwriter as xw
        workbook = xw.Workbook(path)
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
      
    def savecdxml(path,content):
        from openbabel import openbabel
        conv = openbabel.OBConversion()
        output_format = path.split('.')[-1] if '.' in path else 'cdxml'
        conv.SetOutFormat(output_format)
        cdxml_string = conv.WriteString(content)
        if cdxml_string is not None:
            with open(path, "w") as file:
                file.write(cdxml_string)
        else:
            raise ValueError("Failed to convert molecule")
    
    def savedata(file_path, data):
        try:
            available_memory = psutil.virtual_memory().available
            chunk_size = max(int(available_memory / 4), 1 * 1024 * 1024)
            if len(data) <= chunk_size:
                with open(file_path, 'w') as file:
                    file.write(data)
            else:
                with open(file_path, 'w') as file:
                    for i in range(0, len(data), chunk_size):
                        file.write(data[i:i + chunk_size])
        except IOError as e:
            print(f"Error writing to file {file_path}: {e}")
            raise


    # --------------------读取文件-------------------#
    def getexcel(path, orientation='column'):
        import pandas as pd
        xls = pd.ExcelFile(path)
        non_empty_data = []

        for sheet_name in xls.sheet_names:
            df = pd.read_excel(xls, sheet_name=sheet_name, header=None)
            if df.empty:
                continue

            if df.shape[0] == 1:  # 单行数据
                non_empty_data.append(df.iloc[0].tolist())
            elif df.shape[1] == 1:  # 单列数据
                non_empty_data.append(df.iloc[:, 0].tolist())
            else:  # 多行多列数据
                data = df.values.tolist() if orientation == 'row' else df.T.values.tolist()
                non_empty_data.append(data)
        return non_empty_data[0] if len(non_empty_data) == 1 else non_empty_data
        def getname(path, _format, recursive=False):
            from pathlib import Path
        _format = "." + _format.lower()
        file_list = []
        path = Path(path)

        # 根据是否递归选择适当的搜索方法
        if recursive:
            files = list(path.rglob('*' + _format))  # 递归搜索所有子目录
        else:
            files = list(path.glob('*' + _format))  # 仅搜索当前目录

        for file_path in files:
            if file_path.is_file():  # 确保是文件
                new_file_name = file_path.with_suffix(_format)  # 修改后缀为小写
                file_path.rename(new_file_name)  # 重命名文件以更新后缀
                file_list.append(new_file_name.stem)  # 添加不带后缀的文件名到列表中
        return file_list

    def getpdf(path):
        import PyPDF2
        contents = []
        with open(path, 'rb') as file:
            reader = PyPDF2.PdfReader(file)
            num_pages = len(reader.pages)
            for page_num in range(num_pages):
                contents.append(reader.pages[page_num].extract_text())
        return contents
    
    def getname(path,_format):
        if not os.path.exists(path):
            print(f"路径 {path} 不存在。")
            return []  # 如果路径不存在，返回空列表

        _format = "." + _format.lower()  # 确保格式以点开始，并转为小写
        file_list = []
        
        for file in os.listdir(path):
            if file.lower().endswith(_format):  # 将文件名转换为小写后进行比较
                file_path = os.path.join(path, file)
                base_name, ext = os.path.splitext(file)
                
                # 如果文件扩展名已经是小写，不需要重命名
                if ext.lower() == ext:
                    new_file_name = file
                else:
                    new_file_name = base_name + _format  # 创建新文件名
                    new_file_path = os.path.join(path, new_file_name)
                    os.rename(file_path, new_file_path)  # 修改文件扩展名为小写
                file_list.append(base_name)  # 添加不带扩展名的文件名到列表
        
        return file_list
    
    
    # --------------------文件处理-------------------#
    def copy(Source_path, destination_path, *filenames):
        import shutil
        try:
            # 确保目标文件夹存在
            os.makedirs(destination_path, exist_ok=True)

            for filename in filenames:
                source_file = os.path.join(Source_path, filename)
                destination_file = os.path.join(destination_path, filename)

                # 检查源文件是否存在
                if not os.path.exists(source_file):
                    print(f"警告：源文件 {source_file} 不存在，已跳过。")
                    continue

                # 如果目标文件存在，则尝试删除它
                if os.path.exists(destination_file):
                    try:
                        os.remove(destination_file)
                    except Exception as e:
                        print(f"警告：无法删除目标文件 {destination_file}，可能是因为文件正在使用中。")
                        base, ext = os.path.splitext(destination_file)
                        counter = 1  # 用于生成新文件名的计数器
                        new_destination = f"{base}_exist{counter:02}{ext}"
                        # 如果重命名的文件也存在，则增加计数器直到找到未被使用的名称
                        while os.path.exists(new_destination):
                            counter += 1
                            new_destination = f"{base}_exist{counter:02}{ext}"
                        destination_file = new_destination  # 更新最终的目标文件路径

                # 将文件从源路径复制到目标路径
                shutil.copy(source_file, destination_file)
                print(f"文件 {filename} 已从 {Source_path} 复制到 {destination_file}。")
        
        except Exception as e:
            print(f"复制文件时发生错误：{e}")
        
    def cp(Path,folderA,folderB,*args):
        File.copy(os.path.join(Path,folderA),os.path.join(Path,folderB),*args)

    def create_directories(path, *dirnames):
        # 首先检查根路径是否存在，如果不存在，则创建
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"已创建根目录：{path}")

        for dirname in dirnames:
            # 对每个目录名称，构建完整的文件夹路径
            dir_path = os.path.join(path, dirname)

            # 检查该文件夹路径是否存在，如果不存在，则创建
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
                print(f"已创建目录：{dir_path}")
            else:
                print(f"目录 {dir_path} 已存在，跳过创建。")
