# Mayo Lab Functions

Welcome to the Mayo Lab, we are excited to have you. 
In this repository, members from the lab share their MATLAB scripts/functions for preprocessing, analyzing, and visualizing neuronal and behavioral data. 

### You have a dataset... now what?

#### Let's start by figuring out what is included in the dataset, which is formatted as a MATLAB *structure*.

1. Open MATLAB, make a new script, and load the .mat file you received
```buildoutcfg
# In an empty script, type the following lines.
dataFolder = '/Users/kendranoneman/...' # change to point towards where file is located
fileName = 'combinedMaestroSpkSortFEF.pb10pulsA.mat' # change to your file name

load('-mat',sprintf('%s/%s',dataFolder,fileName)); % raw structure
    
# A variable called "exp" should appear in your workspace
```
<img width="323" alt="Screen Shot 2023-03-20 at 5 57 49 PM" src="https://user-images.githubusercontent.com/37158560/226475890-8e0c124d-3475-4d33-acb3-99e0e322f008.png">
