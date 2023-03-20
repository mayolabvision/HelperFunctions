# Mayo Lab Functions

Welcome to the Mayo Lab, we are excited to have you. 
In this repository, members from the lab share their MATLAB scripts/functions for preprocessing, analyzing, and visualizing neuronal and behavioral data. 

| Term               | Meaning                         |
| :----------------------               | :---                            |  
| `structure`        | MATLAB data type that groups related data using data containers called fields   | 
| `command window`   | enables you to enter individual statements at the MATLAB command line, indicated by the prompt (>>)         | 
| `workspace`   | contains variables that you create or import into MATLAB from data files or other programs        | 
| `trial`   | single instance of the monkey completing a task (of a given condition), is repeated hundreds of times to form an entire session       | 
| `session`   | a day of recording that makes up one .mat file, consists of hundreds to thousands of trials   | 
| `unit`   | isolated data recorded from an electrode contact, can be estimated to a "neuron"  |

--- 
### You have a dataset... now what?

#### Let's start by figuring out what is included in the dataset, which is formatted as a MATLAB *structure*.

##### 1. Open MATLAB, make a new script, and load the .mat file you received
```buildoutcfg
# In an empty script, type the following lines.
dataFolder = '/Users/kendranoneman/...' # change to point towards where file is located
fileName = 'combinedMaestroSpkSortFEF.pb10pulsA.mat' # change to your file name

load('-mat',sprintf('%s/%s',dataFolder,fileName)); % raw structure
    
# A variable called "exp" should appear in your workspace
```
<img width="323" alt="Screen Shot 2023-03-20 at 5 57 49 PM" src="https://user-images.githubusercontent.com/37158560/226475890-8e0c124d-3475-4d33-acb3-99e0e322f008.png">

##### 2. Explore the contents of the *exp* and make sure you understand all of its dimensions.

> In the following example, there is data from **2466 trials**

<img width="371" alt="Screen Shot 2023-03-20 at 6 08 09 PM" src="https://user-images.githubusercontent.com/37158560/226476313-8a26af15-191d-4cf5-a3f4-618eb04c5684.png">

**exp.info**

<img width="349" alt="Screen Shot 2023-03-20 at 6 49 07 PM" src="https://user-images.githubusercontent.com/37158560/226482911-915c79b3-edf4-4d3a-8e69-0e5882f9dfca.png">

> **expName** - name of session, same as file name <br />
> **units** - names of all recorded units 
