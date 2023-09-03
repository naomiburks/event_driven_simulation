# Event-Driven Simulation   

This git repository contains tools to run and analyze time-indexed stochastic processes. In particular, we build tools for event-driven population models while keeping abstractions separate enough to maximally extend to other uses. 

## Getting Started

With pip installed, run 
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

## Documentation


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

## Miscelaneous

### What is the point of all those .gitkeep files? 

Git only tracks files, and directories are implicitly determined from them. This means that if a directory contains no files, git will not know about it. So to have git to know about a folder, we put a dummy file in it.  
