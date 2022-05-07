import os

def main():
    stage = "both"

    dataset = "sift"
    datasize = 1
    attrsize = 1024
    attrdim = 512

    command = "./main " + stage + " " + dataset + " " + str(datasize) + " " + str(attrsize) + " " + str(attrdim)

    os.system("cd build && make -j -B")
    os.system("cd build && " + command)

main()