# BloodBrainClassifier
Classify patients as AD (Alzheimer's) vs healthy controls using extracellular RNA-seq (SILVER-seq)

## setup
Clone this repository and then download the required data.
```
git clone https://github.com/mikecuoco/BloodBrainClassifier && \
cd BloodBrainClassifier && \
scripts/get_data.sh
```

## install
To install the dependencies required to run the code in this repository, execute the following
```
conda env create -f env.yaml
```
