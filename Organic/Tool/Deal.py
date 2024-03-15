
import re
import time

import __init__
from Tool.Executor import Executor
import itertools
from itertools import product
from Tool.File import File

class Deal:
    
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
        result,maxnum=Deal.reorder_smiles(result)
        frg1, frg2 =result[:position], result[position:]
        
        frg1, frg2 = frg1.strip(), frg2.strip()
        
        new_Fragments=[]
        for Fragment in Fragments:
            new_Fragments.append(Deal.reorder_smiles(Fragment,maxnum)[0])
        
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
     
    def getsmiles(skeleton:str, positions:list, Fragments:list, brackets:list):
        
        Result = [skeleton]
        Fragments = {position: Fragments[i] for i, position in enumerate(positions)}
        #brackets = [1] * len(positions) if len(brackets) == 1 and brackets[0] == -1 else brackets
        #brackets = {position: brackets[i] for i, position in enumerate(positions)}
        
        brackets = {position: (1 if len(brackets) == 1 and brackets[0] == -1 else brackets[i])
                 for i, position in enumerate(positions)}
        
        positions.sort(reverse=True)
        for position in positions:
            num=len(Result)
            Result=Executor.ThreadExecutor(Deal.addgroups, Result, [position]*num, 
                                           [Fragments[position]]*num, 
                                           [brackets[position]]*num)
            Result=[item for sublist in Result for item in sublist]
        return Result
    
    def insert_substituents(skeleton, positions, substituents_list):
        """
        插入取代基到指定的SMILES字符串中的特定位置。
        
        :param skeleton: 原始的SMILES字符串。
        :param positions: 要插入取代基的位置列表，以化学编号为准（即从1开始）。
        :param substituents_list: 所有可能的取代基列表的列表。
        :return: 所有可能的SMILES字符串列表。
        """
        # 将化学位置转换为Python的索引（从0开始）
        positions = [pos - 1 for pos in positions]  # 化学位置转换为编程位置
        # 生成所有可能的取代基组合
        all_combinations = itertools.product(*substituents_list)
        # 存储所有生成的SMILES
        generated_smiles = []
        # 遍历每个取代基组合
        for substituents in all_combinations:
            modified_skeleton = list(skeleton)  # 将SMILES字符串转换为字符列表以方便插入
            # 添加取代基
            for sub, pos in zip(substituents, positions):
                modified_skeleton.insert(pos, sub)  # 在正确的位置插入取代基
                # 由于每插入一个取代基，后面的索引就会增加，需要对后续位置进行调整
                positions = [p + len(sub) for p in positions]
            # 将修改后的列表转换回字符串，并添加到结果列表中
            generated_smiles.append(''.join(modified_skeleton))
            # 重置positions为初始状态以供下一组取代基使用
            positions = [pos - sum(len(sub) for sub in substituents) for pos in positions]
        return generated_smiles
    
    def insert_substituents_rdkit(skeleton_smiles, positions, substituents):
        """
        使用RDKit在SMILES字符串的特定位置插入取代基
        
        :param skeleton_smiles: 原始骨架分子的SMILES字符串
        :param positions: 要插入取代基的原子位置列表，基于1的索引
        :param substituents: 要插入的取代基的SMILES字符串列表
        :return: 所有可能的带取代基的SMILES字符串列表
        """
        from rdkit import Chem
        from rdkit.Chem import rdChemReactions
        mol = Chem.MolFromSmiles(skeleton_smiles)
        product_smiles = []
        
        # 每个位置和相应的取代基创建反应SMARTS
        for position, sub_smiles in zip(positions, substituents):
            # 创建一个分子碎片，其中包含一个虚拟原子（D）作为连接点
            sub_mol = Chem.MolFromSmiles(f"[*:1]{sub_smiles}")
            # 使用RDKit的AddAttachmentPoint方法添加连接点
            sub_mol_with_attachment = Chem.AddAttachmentPoint(sub_mol, 1, 0)
            
            # 创建一个反应SMARTS，以此来表示在指定位置插入取代基
            reaction_smarts = f"[C:1]>>[C:1]({Chem.MolToSmiles(sub_mol_with_attachment)})[*:1]"
            reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
            
            # 执行反应
            products = reaction.RunReactants((mol,))
            for product_tuple in products:
                for product in product_tuple:
                    # 清理反应后的分子，并将其转换为SMILES
                    cleaned_product = Chem.RemoveHs(product)  # 移除多余的氢
                    product_smiles.append(Chem.MolToSmiles(cleaned_product))
        
        return product_smiles

    def preprocess_combinations(substitution_positions, substituents_lists, bracket_flags,link,binding_map):
        all_combinations = list(itertools.product(*substituents_lists)) 
        
        
        link.keys()
        
        
        
        # 处理绑定关系，确保绑定位置使用相同的取代基
        new_substituents_lists = []
        processed_indices = set()
        for i, substituents in enumerate(substituents_lists):
            if i in processed_indices:
                continue  # 如果已经处理过这个索引，则跳过
            
            if i in binding_map:
                # 这个位置与其他位置绑定，需要确保使用相同的取代基
                linked_idx = binding_map[i]
                linked_substituents = substituents_lists[linked_idx]
                combined_subs = [(sub,) * 2 for sub in substituents]  # 创建成对的组合
                new_substituents_lists.append(combined_subs)
                processed_indices.add(linked_idx)
            else:
                # 没有绑定关系，正常添加
                new_substituents_lists.append([(sub,) for sub in substituents])
        return new_substituents_lists


    def generate_smiles_combinationss(skeleton, substitution_positions, substituents_lists, bracket_flags,link):
        all_combinations = list(itertools.product(*substituents_lists))  # 将迭代器转换为列表
        
        for num in binding_map.keys():
            if num in substitution_positions:
                
                substitution_positions.remove(num)
        binding_map={}
        
        def process_combination(comb):
            local_new_smile = renumbered_skeleton
            local_offset = 0
            local_next_num = next_num  # 为每个线程创建独立的next_num
            
            print(comb)
            #com：取代基
            for position,substituent, flag in zip(substitution_positions,comb,bracket_flags):
                    #1.主取代基和次取代基完全一样。
                    
                    #2.主取代基和次取代基都不一样，但是相互绑定。出现一个必然出现对应的一个
                    renumbered_sub, local_next_num = Deal.reorder_smiles(substituent, local_next_num)
                    
                    insert_str = f"({renumbered_sub})" if flag == 1 else renumbered_sub# 根据标志决定是否加括号
                    
                    insert_pos = position - 1 + local_offset  # 将化学位置转换为从0开始的索引，并考虑偏移量
                    
                    local_new_smile = local_new_smile[:insert_pos] + insert_str + local_new_smile[insert_pos:]

                    local_offset += len(insert_str)
                    
                    binding_position = binding_map.get(position)
                    if binding_position and substituent != comb[substitution_positions.index(binding_position)]:
                        # 插入绑定的取代基
                        binding_sub = comb[substitution_positions.index(binding_position)]
                        renumbered_binding_sub, _ = Deal.reorder_smiles(binding_sub)
                        binding_insert_str = f"({renumbered_binding_sub})" if bracket_flags[substitution_positions.index(binding_position)] == 1 else renumbered_binding_sub
                        binding_insert_pos = binding_position - 1 + local_offset
                        local_new_smile = local_new_smile[:binding_insert_pos] + binding_insert_str + local_new_smile[binding_insert_pos:]
                        local_offset += len(binding_insert_str)

                    
                    
                    
            return local_new_smile

        # 使用ThreadPoolExecutor来并发处理所有组合
        all_combinations = list(itertools.product(*substituents_lists))  # 将迭代器转换为列表
        
        
        return Executor.ThreadExecutor(process_combination, all_combinations)
    
    def generate_smiles_combinations(skeleton, substitution_positions, substituents_lists, bracket_flags,linked_positions):
        def process_combination(comb):
            local_new_smile = skeleton
            offset = 0  # 因为插入会改变字符串长度，需要一个偏移量
            used_subs = {}  # 用于记录每个位置实际使用的取代基
            for idx, (pos, flag) in enumerate(zip(substitution_positions, bracket_flags)):
                if idx in linked_positions:
                    # 如果当前位置在链接位置列表中，使用链接位置的取代基
                    sub = used_subs[linked_positions[idx]]
                else:
                    # 否则，使用当前组合的取代基
                    sub = comb[idx]
                    used_subs[idx] = sub  # 记录使用的取代基
                
                insert_pos = pos - 1 + offset  # 将化学位置转换为从0开始的索引，并考虑偏移量
                insert_str = f"({sub})" if flag else sub
                local_new_smile = local_new_smile[:insert_pos] + insert_str + local_new_smile[insert_pos:]
                offset += len(insert_str)
            return local_new_smile

        all_combinations = itertools.product(*substituents_lists)
        return [process_combination(comb) for comb in all_combinations]
    
    
    

R1 = ['', 'C', 'O', 'N']
R2 = ['', 'F', 'Cl', 'Br', 'I']
R3 = ['', 'C(=O)O', 'C(=O)OC', 'C(=O)OCC', 'C(=O)OC(C)(C)C']
# 总的分子数为4x5x5=100
#R1=Deal.generate_smiles_combinationss('C1=CC=CC=C1', [2, 4, 5], [R1, R2, R3],[1,1,1],[])

#R1=Deal.generate_smiles_combinationss('C1=CC=CC=C1', [2, 4, 5], [R1, R1, R1],[1,1,1])
a=list(itertools.product(*[R1,R2,R3])) 
print(a)