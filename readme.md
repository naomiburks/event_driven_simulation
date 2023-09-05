# Event-Driven Simulation   

This git repository contains tools to run and analyze time-indexed stochastic processes. In particular, we build tools for event-driven population models while keeping abstractions separate enough to maximally extend to other uses. 

## Setup

### Getting Started

This project is written in python 3 with dependencies managed using pip.

In order to install dependencies, run 
```
pip install -r requirements.txt
```

----

To run scripts, you can either run 
```
py -m src.scripts.[script_name]
```
or run 
```
py setup.py develop
```
and then 
```
py src/scripts/[script_name]
```

The latter method is preferred, since running setup is necessary for tests to work properly. Scripts must be run from the root directory. Running from the src/scripts folder **will not work**.

-----

If you use a text editor, I highly recommend switching to an IDE, especially [VSCode](https://code.visualstudio.com/download). Some benefits of IDEs not typical in text editors:
- test support
- automatic linting (code styling and other coherence checks)
- run support 
- version control
- refactoring

Good IDEs also will not strain your computer's resources significantly more than text editors. There is some cost in getting started with IDEs as it can take a bit of work to configure everything to meet your needs. I've found the quality-of-life improvements well worth it and have made the difference for me in being able to build and maintain scalable projects.  

### Running and Instantiating Models

Models should be instantiated in a parameter-agnostic way. Running the model will depend on parameter inputs. To do this, parameter keys (names) must be built in to a model's instantiation. The state space is implicit. Running a model on states outside its state space will lead to unintended behaviour. Most commonly, this leads to errors.  


To run a model, you must provide:
- An initial state
- A dictionary of parameters
- A duration to run 

### Saving and Loading Data

The module src/tools/io.py provides tools to save and load data. 


### Tests

We use pytest for the testing suite. After running py setup.py develop, run ```pytest``` to run all tests. Or ideally use your IDE and make sure your configuration is such that tests can be discovered and run.

### Linting

We use pylint for linting. I'm not sure the best way to use pylint outside an IDE. 


### What is the point of all those .gitkeep files? 

Git only tracks files, and directories are implicitly determined from them. This means that if a directory contains no files, git will not know about it. So to have git to know about a folder, we put a dummy file in it.  


## Project Structure

### Models

The basic object in this library is a model. A model represents a stochastic process. Models provide run functions. These take in an initial state, parameters, and a duration and output the end state. We implement several subclasses. 

- A PopulationModel is one whose state space is a list of nonnegative integers.  
- An EventModel is one who which can be run via stochastic events (with time-independent rates) using the Gillespie algorithm. 
- A LinearModel is both a PopulationModel, an EventModel, and all its events are linear: each event's rate is proportional to the size of one of the populations. LinearModels are used to describe continuous-time multitiype branching processes.
- An ExponentialPopulationModel one whose run function is deterministic with the reuslt determined via exponentiating an instantaneous generator matrix. It can be used to describe the mean behaviour of many LinearModels.
- methylation.py contains classes specific to modeling methylation. These are all LinearModels.


Model classes contain tools to analyze them accessible to all their subclasses. PopulationModels can do a monte carlo simulations of extinction probabilities. LinearModels can directly calculate the extinction probability for an individual of each type directly. LinearModels also provide a converter function that returns the ExponentialPopulationModel that is the mean behaviour of the LinearModel.

### Events

In order to specify an EventModel, you must provide a list of its events. An event has two important properties:

- An implementation function that says what happens to the state when the event occurs. 
- A rate function that says how often the event occurs given a state and model parameter. 

Events have custom subclasses to make this process cleaner, including: 
- Births (two children with the same type as the parent)
- Deaths
- Transitions (cell type becomes another)


### Data Access

Saving and loading data is done via tools in io.py. 

### Plotting data

Plotting data is done via tools in plot.py.

### Scripts

Files to run experiments are inside the scripts foldre. 



## What next?

I want to implement a model for the M->infinity site limit. This requires learning how to efficiently simulate events with continuously-varying rates, since the methylation dynamics are no longer expressed via random events but instead satisfy a deterministic differential equation (until splitting/death). 