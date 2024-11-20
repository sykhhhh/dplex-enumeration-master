# dplex-enumeration-master
## usage
```
g++ main.cpp -O3 -o main
```
### in files BK and DPEnum_without_s
```
./main <file> <k> <l>
```
e.g.
```
./main ./email-Eu-core.bin 2 2
```
### in file DPEnum_with_s
```
./main <file> <k> <l> <s>
```
e.g.
```
./main ./email-Eu-core.bin 2 2 3
```
## input
The input graph should be a binary format. We use ToBin to make a input document become a binary document
usage:
```
g++ ToBin.cpp -o ToBin
./ToBin <graph_file>
```
and it will output a binary document

e.g.
```
./ToBin ./email-Eu-core.txt
```
output: email-Eu-core.bin
