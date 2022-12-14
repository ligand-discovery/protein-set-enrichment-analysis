import os

def listdir_util(path):
    for d in os.listdir(path):
        if d.startswith("_"):
            continue
        else:
            yield d