# Information
Simple python utility for turning the file assets/atomic/weights.dat into a c++ header which can be included to provide easy access to all atomic weights inside 4DSSE

## Requirments
In order to use this utility you will need

- Python
- Pandas

## Usage
```bash
python convertWeightsToHeader.py <path/to/weights.dat> -o atomicWeights.h
```