# Metabolism PBPK

## Coding Style

In order to make code consistent and easier to debug, here are a few rules (not strict, but highly recommended) to follow when programming (ligand.m is a good example of the following practices)

- How to add to the project
  - everytime you start a new task, create a branch with `git switch -c <branchname>` and make all your edits to that branch. This will avoid any contamination of the main branch and
  make it easier to isolate mistakes and debug. Once you have validated the work, create a pull request to merge the branch
  - create atomic commits, every commit should realistically do one thing and be small. This makes it easier to version control and will help you understand what past you was doing.
    - commits should be written as if you're giving commands (Instead of "I added X to Y", write "add X to Y")
    - for style guide, see how previous commits have been made. 

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

## Setting Up Git SSH

First make sure git is installed on your computer (try `git --version` and see if a version number pops up). If not, follow your
specific OS's instructions to setup git. (Automatically setup for MacOS and Linux)

VSCode and/or lazygit are perfect options for managing git

The first step you will need to do is to log into your ugrad account and run the following command:

`ssh-keygen -t rsa -b 2048`

When prompted for the file and passphrase, just hit “Enter” to accept the default option. The output you should see in your terminal should look something like this:

```
[username@hostname ~]$ ssh-keygen -t rsa -b 2048
generating public/private rsa key pair.
enter file in which to save the key (/home/daveho/.ssh/id_rsa): 
created directory '/home/daveho/.ssh'.
enter passphrase (empty for no passphrase): 
enter same passphrase again: 
your identification has been saved in /home/daveho/.ssh/id_rsa
your public key has been saved in /home/daveho/.ssh/id_rsa.pub
the key fingerprint is:
sha256:stviaceiefng/gxlbnjtn89flio+r6gdevq0sqs+hx4 daveho@ugradx
the key's randomart image is:
+---[rsa 2048]----+
|++.....oo ooo.   |
|....+.=.oo..o    |
| + o = o.* +     |
|  b + =.+.=      |
| o b x oso       |
|  e b o o .      |
| . +     o       |
|  .              |
|                 |
+----[sha256]-----+
```
The ssh-keygen command will create a directory in your Linux home directory called `.ssh`. Make sure that this directory is only accessible by you by running the command

`ls -ld ~/.ssh `

You should see output something like the following:

```
[username@hostname ~]$ ls -ld ~/.ssh
drwx------. 2 username users 4 Jan 31 16:55 /home/username/.ssh
```
The permissions `rwx------` mean that only your account can read, write, or execute (search) the contents of the `.ssh` directory, which is what you want to see. If you see different permissions, you can fix them by running the command

`chmod 0700 ~/.ssh `

Copy your public ssh key to the clipboard

Print the contents of your ssh public key by running the command

`cat ~/.ssh/id_rsa.pub `

You should see something like the following (it will probably appear as multiple lines of output in the terminal, but it’s really one long line of text):

`ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCqX2atYK7RtuODlxYZ52TpD1abeA7UxUXk4W39ZKKy3n0bguLOzNOveJNiF7ayGtbirGNBVC/f8snNGpFa8EVjW1Wx+yAVBU0sEAz4h1cYarGUNBhr+SgwGbpFHRDjptkkFpfUu6YoAkY6wv4u4s3396EHR0IttUdOqke9OIKt1nQwr1y30qpyXwLj8nd9s4frmFI4Zo/+Gyux1kYX2kg5C8Iao54HDqTRwSbfww/1KANfF3mjfLI9CI/B5y6C4e+JRa4qoN0dAVJxEeyjo3DztdDm18G1vy2Mo4Od7TvjvA2FirDFnonMknd4QoH0tlwtxk4xzFXjZSW2xEEPWxu9 username@hostname `

In a web browser, go to github.com and log in. Click the menu icon in the upper right-hand corner of the window, and choose “Settings”. Click on “SSH and GPG keys”. Click “New SSH key”.
Paste the copied ssh public key to the “Key” text box, and enter a name for the key in the “Title” text box.


## Cloning This Repository

Optional : create a new folder for github repositories
Run `git clone https://github.com/vivek-booshan/pkpd.git` and this should clone this repository to your local machine
If something goes wrong, ask chatgpt, if that doesn't work, ask me.

You're ready to start coding!!
