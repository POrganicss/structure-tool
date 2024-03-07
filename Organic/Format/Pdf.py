import re
import __init__
from Tool.File import File
import concurrent.futures

import fitz
class Pdf:

    def getnum(line):
        # 使用正则表达式提取数字，包括小数和负数
        numbers = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        # 将字符串数字转换为实际数字类型
        return numbers[-1]  
    
    def get_pages(path):
        pages = []

        def extract_text_from_page(page):
            return page.get_text()

        with fitz.open(path) as doc:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # 提交每一页的提取任务
                futures = [executor.submit(extract_text_from_page, doc.load_page(page_num)) for page_num in range(doc.page_count)]
                
                # 获取每一页的提取结果
                for future in concurrent.futures.as_completed(futures):
                    try:
                        page_text = future.result()
                        pages.append(page_text)
                    except Exception as e:
                        print(f"An error occurred: {e}")

        return pages
            