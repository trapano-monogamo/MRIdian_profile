import MyCacher as mc
import MyUtils as utils
import os
import re

def split_table_row(row:str) -> list:
    elements = []
    current = ""
    for c in row:
        if c == "\t":
            print("hey")
            elements.append(current)
            current = ""
        else:
            current += c
    return [e for e in elements if e]


class Reader:
    pkg: str
    run: str
    res_dir: str
    param_num: int
    
    def __init__(self, _pkg: str, _run: str, _res_dir: str, _param_num: int = 3):
        self.pkg = _pkg
        self.run = _run
        self.res_dir = _res_dir
        self.param_num = _param_num

        # process all four centres
        # ...

    def process_centre(self, centre: str):
        # 1. find all profiles of the same centre in pkg
        # 2. read research center table
        # 3. extract fit parameters
        # 4. return parameters

        # MOVE TO ANOTHER FUNCTION
        # 4. run Profile for every one of them with the found params
        #   - Profile(file_path, out_dir, binning, shared_profiles)

        # ~ 1 ~
        filelist = []
        for root, dirs, files in os.walk(self.res_dir):
            for _file in files:
                current_file_path = os.path.join(root, _file).replace("\\", "/")
                if current_file_path.find(centre) != -1:
                    filelist.append(current_file_path)
                else: continue
            break

        # ~ 2 ~
        table_file = ""
        for root, dirs, files in os.walk(f"./out/{self.pkg}/{self.run}/"):
            for _file in files:
                current_file_path = os.path.join(root, _file).replace("\\", "/")
                if current_file_path.find(centre) != -1:
                    table_file = current_file_path
                    break
                else: continue
            break
        table = []
        with open(table_file, "r") as f:
            table = f.read().split('\n')
        
        # ~ 3 ~
        params = []
        for row in table[1:]:
            tokens = [e for e in row.split(' ') if e]
            params.append(tokens[ -self.param_num*2 :])
        params = [e for e in params if e]
        
        # ~ 4 ~