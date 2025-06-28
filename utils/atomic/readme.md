# Information
Simple python utility for turning the file assets/atomic/weights.dat into a c++ header which can be included to provide easy access to all atomic weights inside 4DSSE

## Requirments
In order to use this utility you will need

- Python
- Pandas

## Usage
```bash
python format.py <path/to/AME.txt> <path/to/nubase.asc> -o speciesData.h
```
