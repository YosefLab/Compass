import json
import os

RESOURCE_DIR = os.path.join("..", "Resources")
MODEL_DIR = os.path.join(RESOURCE_DIR, "Metabolic Models")

def fix(file):
    data = json.load(file)
    if len(data) == 1 and type(data[0]) == list:
        data = data[0]
    
    if type(data[0]) == list and len(data[0]) == 1:
        data = [d[0] for d in data]
    
    for i, d in enumerate(data):
        if type(d) == list and type(d[0]) == int:
            a = "".join(map(chr, d))
            data[i] = a

    file.seek(0)
    json.dump(data, file, indent = 2)
    file.truncate()

model_name = "RECON3_mat"
model_dir = os.path.join(MODEL_DIR, model_name)
model_files = os.path.join(model_dir, "model")
with os.scandir(model_files) as it:
    for entry in it:
        if entry.name == "model.S.json":
            continue
        with open(entry, "r+") as f:
            fix(f)