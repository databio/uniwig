# Uniwig

Given a set of bed files, we want to produce 2 wiggle files: one is the track of start coordinates, the other is of end coordinates.

## Compile

Use the included `Makefile` to compile `/bin/uniwig`:

```
make uniwig
```

Test it with:

```
make test
```

Clean (remove `/bin/uniwig`):

```
make clean
```

## usage:

```
uniwig bedfile stepsize smoothSize variableformat > file.wig 

smoothsize 0 - no smoothing
variableFormat 1 - for variable format, 0 - for fixed format

```
Test command to print chromosome to stdout: 

```
./uniwig test.bed 1 5 0
```
