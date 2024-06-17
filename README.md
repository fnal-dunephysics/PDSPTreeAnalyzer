# Compile

```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build . -- install
```

# Run

You may need to setup proxy first
Then,

```
root -l -b -q run_sec_proton.C
```
# Then

You can make your own analyzing classes and corresponding running scripts.
For analyzing classes, please refer to classes in Analyzers directory.
For running script, please refer to files with "run_*.C"
