
#提供数据验证功能
class Verify:
    
    def Is(elemnet,format):
        if not isinstance(elemnet,format):
            print(str(elemnet)+'错误')

    def Is(elemnets, formats):
        for element,format in zip(elemnets,formats):
            if not isinstance(element, format):
                print(str(element)+'数据错误')

    # 判断单一字符是否为数字
    @staticmethod
    def isnum(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    # 判断一组字符是否为数字
    @staticmethod
    def isnums(strings):
       return all(Verify.isnum(s) for s in strings)
