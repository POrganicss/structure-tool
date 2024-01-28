from typing import List
from itertools import product
import pandas as pd

class CForm:
    def __init__(self):
        self.data = []

    def add(self, values):
        self.data.append(values)

    def Extraadd(self, extra_values):
        self.data.append(extra_values)

    def getParameterNumber(self, parameter_name):
        # Implement this method based on your specific needs
        pass

class Write:
    def WriteToExcel(self, data, filename):
        # Implement this method based on your specific needs
        pass

def MSCompute(elements: List[str], conditions: List[str], M: float) -> int:
    element_weights = {
        "C": C612[0], "H": H11[0], "O": O816[0], "N": N714[0], "S": S1632[0],
        "P": P1531[0], "Si": Si1428[0], "F": F919[0], "Cl": Cl1735[0],
        "Br": Br3579[0], "I": I53127[0]
    }

    def calculate_Hnc(Cn, Nn, On, Fn, Sin, Pn, Sn, Cln, Brn, In):
        return (M - sum(weights * value for element, weights, value in zip(elements, element_weights.values(), [Cn, Nn, On, Fn, Sin, Pn, Sin, Fn, Cln, Brn, In]))) / element_weights["H"]

    Data = []

    for combination in product(*(range(1, int(M / element_weights["C"]) + 1) for _ in elements[1:])):
        combination_dict = dict(zip(elements[1:], combination))
        Hnc = calculate_Hnc(**combination_dict)

        for Hn in range(int(Hnc), int(Hnc) + 1):
            MolMass = sum(weights * value for element, weights, value in zip(elements, element_weights.values(), [combination_dict[element] for element in elements])) + element_weights["H"] * Hn

            D = (combination_dict["C"] + combination_dict["Si"]) * 2 - (Hn + combination_dict["F"] + combination_dict["Cl"] + combination_dict["Br"] + combination_dict["I"]) + combination_dict["N"] + combination_dict["P"] + 2 / 2.0

            if round(D) == D and D >= 0.0:
                CHB = (combination_dict["C"] + Hn) / (sum([combination_dict[element] for element in elements[2:]]) + 1.0)
                CB = combination_dict["C"] / (Hn + 1.0)
                E = abs(MolMass - M)
                ET = abs(MolMass - M) / MolMass * 10000.0

                Mol = CForm()
                Mol.add([combination_dict[element] for element in elements])
                Mol.Extraadd([D, M, MolMass, CHB, CB, E, ET])

                if E < 0.1 and D <= 2 and CHB >= 1.0 and DataDeal.Result(conditions, Mol):
                    RMol = CForm()
                    RMol.add([combination_dict[element] for element in elements])
                    RMol.Extraadd([D, M, MolMass, CHB, CB, E, ET])
                    Data.append(RMol)

    Write().WriteToExcel(Data, f"D:/LCMS-{round(M * 100.0) / 100.0}-Chemkit2.0.xls")
    return len(Data)

if __name__ == "__main__":
    # 调用 MSCompute 函数示例
    elements = ["C", "H", "O", "N", "S", "P", "Si", "F", "Cl", "Br", "I"]
    conditions = ["C=2", "H>4", "O=1"]
    M = 50.0
    result = MSCompute(elements, conditions, M)
    print(f"Number of valid structures: {result}")
