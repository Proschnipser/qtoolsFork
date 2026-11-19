import glob, sys, os

with open("position.txt","r", encoding="utf-8") as f: 
    positions=f.read()

print(positions)

prefix="trained_"
for first_level in glob.glob(f"{prefix}*"):
    if os.path.isdir(first_level):
        second_level = next((d for d in glob.glob(os.path.join(first_level, "*")) if os.path.isdir(d)),None)
        if second_level:
            print(f"Found: {second_level}")
            abs_path = os.path.abspath(second_level)
            os.system(f"python3 /data/joscha/qtoolsFork/Skripte/histogram.py {abs_path}/weights/ {positions}")