
def ReadCSVFile(fileName: str, startIdx: int, endIdx: int) -> list:
    dataSample = []
    with open(fileName) as f:
        for lines in f.readlines():
            strList = lines.split(',')
            if len(strList) > endIdx:
                numberList = []
                for j in range(startIdx, endIdx + 1):
                    numberList = numberList + [float(strList[j])]
                dataSample = dataSample + [numberList]
    return dataSample


def SaveCSVFile(fileName: str, content: list, startIdx: int, endIdx: int):
    with open(fileName, 'w') as f:
        for line in content:
            if len(line) > endIdx:
                for i in range(startIdx, endIdx + 1):
                    if i == endIdx:
                        f.write(str(line[i]) + "\n")
                    else:
                        f.write(str(line[i]) + ", ")
