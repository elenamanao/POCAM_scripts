# data_reader_hdf5

This program reads one or more `.pdd` files and writes a single combined HDF5 output file named `combined_subruns.h5`.

---

## Required files

Before compiling, ensure the following files are present in the same directory:

- `data_reader_hdf5.c`
- `payload_header_drop.h`

The header file is required because it is included directly by the source code. If it is missing or not in the compiler include path, compilation will fail.

---

## HDF5 dependency

The program uses the HDF5 C API (`hdf5.h`). You must have the HDF5 development headers and libraries installed on your system.

On Ubuntu / Debian / WSL systems, install them with:

```bash
sudo apt-get update
sudo apt-get install libhdf5-dev

---

## Compilation

To compile the C source file `data_reader_hdf5.c` into an executable named `data_reader_hdf5`, run the following command from the directory containing the files:

```bash
gcc data_reader_hdf5.c -o data_reader_hdf5 -I/usr/include/hdf5/serial -lhdf5_serial -lm

Depending on the system, the HDF5 library name or installation path may differ, so make sure the path is correct.

---

## Execution

The executable `data_reader_hdf5` takes one or more `.pdd` files as command-line arguments. Each file passed to the program is read and combined internally into a single HDF5 output file.

In other words, the `.pdd` files to be processed must be provided explicitly as arguments when running the program.

To run the program on all `.pdd` files contained in a directory `dir_containing_pdd_files`, execute:

```bash
./data_reader_hdf5 path/to/dir/containing/pdd/files/dir_containing_pdd_files/*.pdd

