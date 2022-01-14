from sys import argv
from os.path import isfile
from os import X_OK, access


FLAG = 0

SILVA_ARB = argv[1]
SINA_EXECUTABLE = argv[2]
SORT_ME_RNA_DB1 = argv[3]
SORT_ME_RNA_DB2 = argv[4]
SORT_ME_RNA_TOOL = argv[5]
CLUSTERING_TOOL = argv[6]
RAPID_NJ = argv[7]


if not isfile(SILVA_ARB):
    print('SILVA_ARB: ' + SILVA_ARB + ' not found')
    FLAG = 1

if not isfile(SINA_EXECUTABLE):
    print('SINA_EXECUTABLE: ' + SINA_EXECUTABLE + ' not found')
    FLAG = 1

if not access(SINA_EXECUTABLE, X_OK):
    print('SINA_EXECUTABLE: ' + SINA_EXECUTABLE + ' not executable')
    FLAG = 1

if not isfile(SORT_ME_RNA_DB1):
    print('SORT_ME_RNA_DB1: ' + SORT_ME_RNA_DB1 + ' not found')
    FLAG = 1

if not isfile(SORT_ME_RNA_DB2):
    print('SORT_ME_RNA_DB2: ' + SORT_ME_RNA_DB2 + ' not found')
    FLAG = 1

if not isfile(SORT_ME_RNA_TOOL):
    print('SORT_ME_RNA_TOOL: ' + SORT_ME_RNA_TOOL + ' not found')
    FLAG = 1

if not isfile(CLUSTERING_TOOL):
    print('CLUSTERING_TOOL: ' + CLUSTERING_TOOL + ' not found')
if not access(CLUSTERING_TOOL, X_OK):
    print('CLUSTERING_TOOL: ' + CLUSTERING_TOOL + ' not executable')
    FLAG = 1


if not isfile(RAPID_NJ):
    print('RAPID_NJ: ' + RAPID_NJ + ' not found')
if not access(RAPID_NJ, X_OK):
    print('RAPID_NJ: ' + RAPID_NJ + ' not executable')
    FLAG = 1


if FLAG:
    print('Testing failed, please correct the errors reported.')
else:
    print('Testing successful. All elements of pipeline present.')