
import re
import __init__


class Group:
       
    def reorder_smiles(smiles, start_num=1):
        pattern = re.compile(r'(\d+)')
        numbers = pattern.findall(smiles)
        unique_numbers = sorted(set(numbers), key=numbers.index)
        max_num = start_num
        for num in unique_numbers:
            smiles = smiles.replace(num, str(max_num), 1)
            max_num += 1
        return smiles, max_num
    
    def addgroups(result, position, Fragments, bracket):
        
        
        Results = []
        result,maxnum=Group.reorder_smiles(result)
        frg1, frg2 =result[:position], result[position:]
        
        frg1, frg2 = frg1.strip(), frg2.strip()
        
        new_Fragments=[]
        for Fragment in Fragments:
            new_Fragments.append(Group.reorder_smiles(Fragment,maxnum)[0])
        
        # for Fragment in new_Fragments:
        #     if Fragment.strip()=='':
        #         Results.append("".join(frg1 + frg2))
        #     else:
        #         Fragment = f"({Fragment})" if bracket == 1 else Fragment
        #         Results.append("".join(frg1 + Fragment + frg2))
        
        
        Results = [frg1 + (f"({Fragment})" if bracket == 1 else Fragment) 
                + frg2 if Fragment.strip() else frg1 + frg2
                    for Fragment in new_Fragments
                    ]        
                
        return Results
     
    def getnum(original_position, operations):
        """
        计算经过一系列操作后的新位置。
        
        :param original_position: 原始的位置
        :param operations: 操作的列表，每个元素是一个元组(operation_type, positionA, positionB)
                            其中 operation_type 有 'add', 'remove', 'replace', 'cyclize', 'decyclize' 五种类型
                            'replace' 操作的元组可以包含四个元素，但对于 'cyclize' 和 'decyclize' 只需要三个
        :return: 新的位置
        """
        new_position = original_position
        for operation in operations:
            operation_type = operation[0]
            match operation_type:
                case 'add':
                    start_position, new_char_count = operation[1:3]
                    if start_position <= original_position:
                        new_position += new_char_count
                case 'remove':
                    start_position, affected_char_count = operation[1:3]
                    if start_position < original_position:
                        new_position = max(start_position, new_position - min(affected_char_count, original_position - start_position))
                case 'replace':
                    start_position, affected_char_count, new_char_count = operation[1:4]
                    if start_position < original_position:
                        if start_position + affected_char_count > original_position:
                            new_position = start_position + min(original_position - start_position, new_char_count)
                        else:
                            new_position += new_char_count - affected_char_count
                case 'cyclize' | 'decyclize':
                    positionA, positionB = operation[1:3]
                    # 对于环化和去环化，如果任一位置在原始位置之前，将影响位置
                    if positionA < original_position:
                        new_position += 1 if operation_type == 'cyclize' else -1
                    if positionB < original_position:
                        new_position += 1 if operation_type == 'cyclize' else -1
        return new_position

    def addgroup(smiles,operationss,positions,group):
         # 检查输入列表长度是否相同
        if len(smiles) != len(positions):
            raise ValueError("The length of 'smiles' and 'positions' must be the same.")
        # 更新后的 SMILES 列表
        updated_smiles = []
        # 遍历 SMILES 和对应的插入位置
        for sm, pos,op in zip(smiles, positions,operationss):
            pos=Group.getnum(pos,op)
            # 检查插入位置是否在合理范围内
            if pos < 0 or pos > len(sm):
                raise ValueError("The insertion position is out of range.")
            sm,num=Group.reorder_smiles(sm)
            group_new=Group.reorder_smiles(group,num)
            # 在指定位置插入基团
            new_sm = sm[:pos] + group_new + sm[pos:]
            updated_smiles.append(new_sm)
            op.append(('add',pos,len(group_new)))
        # 返回更新后的 SMILES 列表
        return updated_smiles,operationss
        
    def removegroup(smiles,operationss,positions,group):
        # 检查输入列表长度是否相同
        if len(smiles) != len(positions):
            raise ValueError("The length of 'smiles' and 'positions' must be the same.")

        # 更新后的 SMILES 列表
        updated_smiles = []
        
        # 遍历 SMILES 和对应的插入位置
        for sm, pos, op in zip(smiles, positions, operationss):
            pos = Group.getnum(pos, op)  # 获取经过之前操作更新后的位置
            # 检查插入位置是否在合理范围内
            if pos < 0 or pos + len(group) > len(sm):
                raise ValueError("The removal position is out of range.")
            
            # 检查是否确实为指定的基团
            if sm[pos:pos+len(group)] != group:
                raise ValueError("The specified group does not match the group at the given position.")
            
            # 移除指定位置的基团
            new_sm = sm[:pos] + sm[pos+len(group):]
            updated_smiles.append(new_sm)

            # 记录这一操作
            op.append(('remove', pos, len(group)))

        # 返回更新后的 SMILES 列表
        return updated_smiles, operationss
    
    
    def replacegroup(smiles, operationss, positions, old_group, new_group):
        # 检查输入列表长度是否相同
        if len(smiles) != len(positions):
            raise ValueError("The length of 'smiles' and 'positions' must be the same.")
        
        # 更新后的 SMILES 列表
        updated_smiles = []
        
        # 遍历 SMILES 和对应的插入位置
        for sm, pos, op in zip(smiles, positions, operationss):
            pos = Group.getnum(pos, op)  # 获取经过之前操作更新后的位置
            # 检查替换位置是否在合理范围内
            if pos < 0 or pos + len(old_group) > len(sm):
                raise ValueError("The replacement position is out of range.")
            
            # 检查是否确实为指定的基团
            if sm[pos:pos+len(old_group)] != old_group:
                raise ValueError("The specified group does not match the group at the given position.")
            
            # 替换指定位置的基团
            new_sm = sm[:pos] + new_group + sm[pos+len(old_group):]
            updated_smiles.append(new_sm)

            # 记录这一操作
            op.append(('replace', pos, len(old_group), len(new_group)))

        # 返回更新后的 SMILES 列表和操作列表
        return updated_smiles, operationss
    
    
    def cyclizegroup(smiles, operationss, positionsAs, positionsBs):
        # 检查输入列表长度是否相同
        if len(smiles) != len(positionsAs) or len(smiles) != len(positionsBs):
            raise ValueError("The lengths of 'smsiles', 'positionsA' and 'positionsB' must be the same.")

        # 更新后的 SMILES 列表
        updated_smiles = []

        # 设置环的标记数字，SMILES中通常使用1至9的数字来标记环
        # 注意：在实际使用中，应确保不与字符串中现有的环数字冲突
        

        # 遍历 SMILES 和对应的插入位置
        for idx, (sm, op, posA, posB) in enumerate(zip(smiles, operationss,positionsAs,positionsBs)):
            posA = Group.getnum(posA, op)  # 获取经过之前操作更新后的位置A
            posB = Group.getnum(posB, op)  # 获取经过之前操作更新后的位置B

            # 检查插入位置是否在合理范围内
            if posA < 0 or posA >= len(sm) or posB < 0 or posB >= len(sm):
                raise ValueError("The cyclization position is out of range.")
            
            sm,num = Group.reorder_smiles(sm)
            cycle_number=num+1
            
            # 在指定位置插入环化标记
            new_sm = sm[:posA] + f"{sm[posA]}{cycle_number}" + sm[posA + 1:]
            new_sm = new_sm[:posB] + f"{new_sm[posB]}{cycle_number}" + new_sm[posB + 1:]
            updated_smiles.append(new_sm)

            # 记录这一操作
            op.append(('cyclize', posA, posB))  # 这里'1'表示添加了一个标记

            # 可能需要为下一个环使用不同的数字
            cycle_number = (cycle_number % 9) + 1  # 循环使用1到9的数字

        # 返回更新后的 SMILES 列表和操作列表
        return updated_smiles, operationss
        
    
    def decyclizegroup(smiles, operationss, positionsAs, positionsBs):
        # 检查输入列表长度是否相同
        if len(smiles) != len(positionsAs) or len(smiles) != len(positionsBs):
            raise ValueError("The lengths of 'smiles', 'positionsA', and 'positionsB' must be the same.")
        
        # 更新后的 SMILES 列表
        updated_smiles = []

        # 遍历 SMILES 和对应的去环化位置
        for idx, (sm, op, posA, posB) in enumerate(zip(smiles, operationss, positionsAs, positionsBs)):
            posA_updated = Group.getnum(posA, op)  # 获取经过之前操作更新后的位置A
            posB_updated = Group.getnum(posB, op)  # 获取经过之前操作更新后的位置B

            # 检查去环化位置是否在合理范围内
            if posA_updated < 0 or posA_updated >= len(sm) or posB_updated < 0 or posB_updated >= len(sm):
                raise ValueError(f"The decyclization position is out of range for molecule {idx+1}.")

            # 检查预期的环化数字是否匹配
            if sm[posA_updated] == sm[posB_updated] and sm[posA_updated].isdigit():
                # 执行去环化操作：移除环化数字
                new_sm = sm[:posA_updated] + sm[posA_updated + 1:]
                new_sm = new_sm[:posB_updated - 1] + new_sm[posB_updated:]  # 注意：第二次删除要调整位置
                updated_smiles.append(new_sm)

                # 记录这一操作
                op.append(('decyclize', posA_updated, posB_updated))

            else:
                # 如果没有匹配的环化数字，保留原SMILES字符串
                updated_smiles.append(sm)
                print(f"No matching cycle numbers at positions {posA} and {posB} for molecule {idx+1}.")

        # 返回更新后的 SMILES 列表和操作列表
        return updated_smiles, operationss
    
