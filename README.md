# Mayo Lab Functions

Welcome to the Mayo Lab, we are excited to have you. 
In this repository, members from the lab share their MATLAB scripts/functions for preprocessing, analyzing, and visualizing neuronal and behavioral data. 

| Term               | Meaning                         |
| :----------------------               | :---                            |  
| `structure`        | MATLAB data type that groups related data using data containers called fields   | 
| `command window`   | enables you to enter individual statements at the MATLAB command line, indicated by the prompt (>>)         | 
| `workspace`   | contains variables that you create or import into MATLAB from data files or other programs        | 
| `trial`   | repetitive episode in an experiment, where the same motor task is performed again and again so that a changing behavioral response can be identified with a learning process      | 
| `session`   | a day of recording that makes up one .mat file, consists of hundreds to thousands of trials   | 
| `unit`   | isolated data recorded from an electrode contact, can be estimated to a "neuron"; there can be multiple units on one contact, thus the naming convention ("1a", "1b", "3a", etc...) |
| `signal-to-noise ratio (SNR)`  | main feature that characterizes an ideal extracellular microelectrode for recording brain signals, which is a measure of the fidelity of the received message for the whole frequency band containing useful neural information |

--- 
## You have a dataset... now what?

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
> **units** - names of all recorded units (57 in this case) in alphabetical order <br />
> **rotFactor** - degree to which the target directions are rotated <br />
> **SNRs** - signal-to-noise ratio for each unit, sorted from highest to lowest <br />
> **channels** - names of recorded units, in same order as SNRs <br />

**exp.dataMaestroPlx**

<img width="1317" alt="Screen Shot 2023-03-23 at 2 42 38 PM" src="https://user-images.githubusercontent.com/37158560/227317086-bb9f4f8b-105f-48e6-81ee-d0ac586aa52c.png">

*Each field (row) of the structure is a single trial*

> **trName** - name of trial: [1] = experimenter (e.g. "P" = "Patrick"), [2] = monkey (e.g. "A" = "Aristotle"), [3:4] = session number, [end-3:end] = trial number <br />
> **trType** - condition type, specific to experiment (e.g. p\_d180\_c100\_sp10 = "pursuit"\_"direction"\_"contrast"\_"speed") <br />
> **mstEye** - contains horizontal and vertical eye positions and velocities for entire length of trial (at each ms) <br />
> **tagSection** - contains timing info for trial (stTimeMS = stimulus onset, durMS = length of time after stim onset) <br />
> **units** - spike times for each unit recorded during that trial, relative to trial onset <br />

---
## What are the first steps you should take for preprocessing the data? 

#### 1. Remove empty trials, trials missing important fields, and units that are dropped over the course of a session
```buildoutcfg
# After loading in a exp datafile, run the struct_clean.m function from the \preprocessing folder

exp_clean = struct_clean(exp);
    
# A variable called "exp_clean" should appear in your workspace
```

*The data in exp_clean is structured in the same way as exp, but now there are no empty fields/trials or trials missing important information (like stimulus onset time). Trials are also thrown out if the mean firing rate of the neurons in that trial exceeds 3 standard deviations from the mean and units that drop off over the course of the trial are removed. The units are also now numerically/alphabetically ordered in both exp_clean.dataMaestroPlx(:).units and exp_clean.info.channels.*

#### 2. If microsaccades aren't relevant towards your project, remove trials where they occur around stimulus onset
```buildoutcfg
# Run the detect_msTrials.m function from the \preprocessing folder on the exp_clean structure
# This runs for a single trial and outputs whether a saccade was detected, but you can use cellfun or a for loop to go through all of the trials
# Inputs include stimulus onset time, size of window around stimulus onset to analyze, and thresholds (acceleration and velocity)

exp_clean.dataMaestroPlx(logical(cellfun(@(q) detect_msTrials(struct2cell(q),motionStart,50,100,750,50), {exp_clean.dataMaestroPlx.mstEye}.', 'uni', 1))) = [];
    
# This will reduce the number of trials in "exp_clean"
```

#### 3. Remove trials with conditions you don't want to analyze in this project
```buildoutcfg
# Run the struct_pullConditions.m function from the \preprocessing folder on the exp_clean structure
# Inputs include which columns you want to prune the conditions for, what conditions you want to keep from those columns, and what columns you want to use to label the trial

extract_conditions = {'1fXXX','2fXXX'};
extract_columns = [3 4];
define_columns = 1;

[exp_clean,condition_names] = struct_pullConditions(exp_clean,extract_conditions,extract_columns,define_columns);

# This will reduce the number of trials in "exp_clean"
```


