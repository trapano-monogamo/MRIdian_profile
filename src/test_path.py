import os

path = "./res"
filelist = []

for root, dirs, files in os.walk(path):
    for _file in files:
        file_path = os.path.join(root,_file).replace("\\", "/")
        filelist.append(file_path)
        with open(file_path, "r") as f:
            pass

print(filelist)
