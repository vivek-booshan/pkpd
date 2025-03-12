# Metabolism PBPK

## Coding Style

In order to make code consistent and easier to debug, here are a few rules (not strict, but highly recommended) to follow when programming (ligand.m is a good example of the following practices)

- How to add to the project
  - everytime you start a new task, create a branch with `git switch -c <branchname>` and make all your edits to that branch. This will avoid any contamination of the main branch and
  make it easier to isolate mistakes and debug. Once you have validated the work, create a pull request to merge the branch
  - create atomic commits, every commit should realistically do one thing and be small. This makes it easier to version control and will help you understand what past you was doing.

- each function should only perform one task, if your function has more than 3 indentations or starts to get unwieldy, it's time to refactor
  - functions should use snake case `function out = this_is_a_function(x1, x2, x3)`
  - it is possible for some functions to literally have just one line of code (Use this carefully as this adds runtime overhead just for a single line)
  - if possible, keep an arguments list in each function, this will save you lots of debugging time
  - minimize nesting (as much as necessary and no more)

- Classes should be labeled with capital letters and ideally a single word `classdef Class`
  - minimize classes used to absolutely necessary groupings as classes also add runtime overhead
  - Methods should be static unless needed otherwise

- For each organ:
  - have a primary (main) function that takes a (time, concentrations, parameter) input and outputs a vector of type double
  - however, all subfunctions should output a double (a single value)

