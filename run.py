import os

def main():
    stage = "both"

    dataset = "deep"
    datasize = 1

    command = "./main " + stage + " " + dataset + " " + str(datasize)

    os.system("cd build && make -j -B")
    os.system("cd build && " + command)

main()