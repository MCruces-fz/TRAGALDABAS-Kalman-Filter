import os
import json
import numpy as np

with open("config/configuration.json", "r") as config_file:
    config = json.load(config_file)

print("-" * 33)
print(r"   _____         _           _ ")
print(r"  |_   _| __ ___| |__   ___ | |")
print(r"    | || '__/ _ \ '_ \ / _ \| |")
print(r"    | || | |  __/ |_) | (_) | |")
print(r"    |_||_|  \___|_.__/ \___/|_|")

print(r"  ____ ___  _   _ _____ ___ ____ ")
print(r" / ___/ _ \| \ | |  ___|_ _/ ___|")
print(r"| |  | | | |  \| | |_   | | |  _ ")
print(r"| |__| |_| | |\  |  _|  | | |_| |")
print(r" \____\___/|_| \_|_|   |___\____|")

print("\n", "-" * 33)
print("\n\n\n")

print('''
P L A N E S   D I S T R I B U T I O N

                                  [mm]   [mm]
T1 # -------------------------- # 1826      0  TOP

T2 # -------------------------- # 1304    522
T3 # -------------------------- #  924    902


T4 # -------------------------- #   87   1739  BOTTOM
                                     0         GROUND
''')

# C H O O S E   V E R S I O N
print("Choose TRAGALDABAS Version:")
print("     1.- Planes -------> T1,        T3, T4")
print("         Heights/mm -> 1826,       924, 87")
print("     2.- Planes -------> T1,   T2,  T3, T4")
print("         Heights/mm -> 1826, 1304, 924, 87")
print("     0.- Write another configuration")
print(f"     c.- {config['vz0']} (CACHED)")
'''
while True:
    version = input("Version (1, 2, 0, c): ")
    print(f"HIT ENTER:{type(version)}")
    if version == "c": break
    try:
        version = int(version)
        if 0 <= version <= 2:
            break
        else:
            print("It must be (2 < version < 1)")
    except ValueError:
        print("It must be an integer")
        pass

if version == 0:
    new_vz0 = input("Write heights in mm separated by commas:\n")
    new_vz0 = np.array(sorted(new_vz0.split(",")), dtype=np.int32).tolist()
    print(new_vz0)
elif version in [1, 2]:
    print(f"Choosen version V{version}")
elif version == "c":
    version = config["vz0"]
print("version!!", version)
'''
while True:
    version = input("Version (1, 2, 0, c): ")
    if version == "1":
        config["vz0"] = [1826, 924, 87]
        break
    elif version == "2":
        config["vz0"] = [1826, 1304, 924, 87]
        break
    elif version == "0":
        new_vz0 = input("Write heights in mm separated by commas:\n")
        config["vz0"] = np.array(sorted(new_vz0.split(",")), dtype=np.int32).tolist()
        break
    elif version in ["c", "C", ""]:
        break
    else:
        pass

print(config["vz0"])
