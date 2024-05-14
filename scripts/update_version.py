import sys

def update_version(file_path, new_version):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            if "version=" in line:
                line = f"version='{new_version}'\n"
            file.write(line)

if __name__ == "__main__":
    update_version('setup.py', sys.argv[1])