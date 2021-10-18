# Final Project Guidelines

Your final project is an opportunity for you to write a “substantial” piece of code that leverages the skills you have learned in the course to date (as well as some you will learn by virtue of doing the project). Your goal is to create a package, which is just a fancy name for code which 
- Is stored as functions and classes  in python files (we’ll call them modules), and can thus
- Be imported into other code for general use, and is
= Installable, which simply means you can import your functions/classes from anywhere on the computer, not just in the folder where the source code is installed. 

You may work **alone** or with **one (1) partner**.  If you work with a partner, we expect to see ~50% contribution from both of you, and will confirm that this is the case via questioning, code reviews, and github commits. We also expect the project to be **more substantive** than a project worked on alone.

## Project Requirments

The following should be considered base requirements, but are subject to modification if you talk to the instructors and provide rationale for why such a change is warranted. 

We want you to pick projects that require a substantial amount of code to pull off. When we say substantial, we mean, roughly, N>2 `classes` (each with multiple methods), and/or roughly N>10 `functions` that work together to solve the problem. Your code must also be organized into `modules`, be fully installable on a system, and be hosted on `Github`. We will help you select projects of a scope such that you shouldn’t need to worry about meeting “code minimum” requirements, and we would never ask you to write longer code when you can achieve something in fewer lines. Thus, you shouldn’t think in terms of “how many lines” or “I need one more function.” Instead, we will do a few things to be sure you’re on the right track. 
- First, we are going to expand the code for the solution to Lab 5 (spectra) into something that would be an appropriate minimum amount of code/complexity for your projects. You can peruse this code at will and qualitatively assess with your own code if you’re in the right ballpark. 
- Second, we are going to have code reviews along the way (see below), during which we will tell you if you are on track to meet our requirements for the project.
- Third, we will suggest (and help) you create project ideas that are modular and extensible. This means if you finish some core functionality and your project is not “long enough,” you can easily add a new feature to it which will then contribute to a more fully featured code and a completed project. 

### Examples of Possible Projects 
- Interactive tools for anything we’ve visualized this semester
- Dataset explorer (for example, using the upcoming pandas stuff)
- Image Zoo -- retrieve images (galaxies, etc), and ask the user to categorize them (store results, maybe train ML?) 
- Retrieve cutouts from Google Mars and perform some feature analysis/detection? (or Google Earth)
- A tool that helps project management, etc., for a club or organization you are part of. 
- Implement a classic game (scrabble, battleship, chess), write an AI to play against. 
- Make a website using [streamlit](http://streamlit.io) to make an element of your research code an interactive tool/widget. ([example by me](https://share.streamlit.io/prappleizer/dslm-tilts/main/calc_tilt.py))
- Write a user-friendly N-Body simulator to set up, e.g., planetary systems, or the Toomre & Toomre galaxy simulations, run them, and visualize the results 
- Use MCMC to write a tool that helps the user fit some astronomical (or other) models to data. Examples: Fit light curves to exoplanet data, fit earth temperature data with a climate model, etc)
- Write a tool that teaches an astronomical concept, that could be used in the classroom (like 255!)

##  Project Stages + Timeline and Deadlines

Because we are not having a midterm, your final project is worth 45% of your grade. We will spread this out as follows: 

1. [**5%, due Oct 29**] Submission of Project Abstract (1 paragraph) describing your intended project 

2. [**5%, due Nov 5**] Submission of Revised Project Abstract after instructor comments + skeleton sketch of the modules/functions/classes you think should be present in your code.  

3. [**10%, due Nov 19**] Code Review I - meet with an instructor and describe what you have so far, issues you’re facing, and next intended steps.  Sections of your project should be functional (and demo-able), and you should have ~50% of the code at least drafted 

4. [**10%, due Dec 3**] Code Review II: - meet with instructor with a draft of the completed project. Elements may still be under construction, but the overall code must meet the requirements by this stage.

5. [**15%, due Dec 10**] Final Presentation - consists of (1) peer review, (2) a demo of your code for the class + pick 1-2 areas of the codebase to explain (ideally the part where you learned something new, or used a tool we haven’t covered in class), and (3) You will also “turn in” your projects by hosting them for the world on Github. 

## Grading Breakdown and Rubrics 

Each section of the final project described above has an associated rubric shown below. As you prepare each phase, consider these as a checklist for what you should have ready by each deadline. We hope that by providing these rubrics, we can alleviate any uncertainty about what we’re expecting to see as you progress through the project, and you should ask any clarifying questions about these rubrics as soon as you have them. 


### Project Abstract (15 possible points, 5% of grade) 

Your project abstract should be a paragraph (4-6 sentences ~250 words) defining the idea for your project. It should indicate what your codebase will do, and why such a code existing would be useful/fun/cool. It should roughly sketch (in 1-2 sentences) the structure of the code. 

**Example:** I propose a codebase `touchdownML` (`tdML`) which uses machine learning to help fantasy football managers to predict the scoring of players on their fantasy team, and thus make informed decisions about who to play and who to bench. This project is motivated by the fact that enormous datasets exist detailing NFL football players every catch, drop, score, etc., with ancillary data (which team was winning, the weather, etc.), and such data is readily available via APIs which make importation to python feasible. I plan to use the `sbi` ML inference package[1] to train a network to learn the connection between the myriad data for players, and the conditions which lead to fantasy scoring. My package will contain several classes, including one to request and parse raw data via APIs, one to “clean” the data and select what will be used in training, and once which connects the ML framework to my data. I will also have a “frontend” (either a class or function), which allows the user to input a player and conditions for an upcoming game, and returns a probability distribution for that player’s fantasy scoring. 

#### Project Abstract Rubric 
|                                                                                            | Not present/unacceptable (0) | Needs Improvement (1) | Acceptable (2)  | Strong (3)  |
|--------------------------------------------------------------------------------------------|------------------------------|-----------------------|-----------------|-------------|
| Complete and on time                                                                       |                              |                       |                 |             |
| Clearly articulates a purpose/problem to be addressed                                      |                              |                       |                 |             |
| Describes the project feasibility and any new tools needed to carry it out                 |                              |                       |                 |             |
| Contains a rough framework of some major “pieces” of code that will be needed.             |                              |                       |                 |             |
| Clearly articulates what the code will do, and how it will be interacted with by the user. |                              |                       |                 |             |

### Revised Project Abstract (9 possible points, 5% of grade)

The instructors will provide detailed feedback on your abstract, including whether we think the project is sufficiently broad (but not too broad), sufficiently detailed, achievable in the time we have, in python, etc. We then want you to do a little research and scratch work starting things up, evaluate our suggestions, and make a revised abstract that addresses any concerns and suggestions made. Since the example above is considered “strong” we don’t provide a revised one here. 

Additionally, by this point we want you to submit a code skeleton that roughly lays out a few major areas of your code. Some experimentation / getting going with your projects will help you complete this portion. An example is provided below. 

Package touchdownML 
- Request module (requests.py)
- - Any functions and classes needed to request data from API web services
- Database module (database.py)
- - Methods to clean incoming data from requests (remove junk data, etc) 
- - Methods to create a singular database for each player containing all known information from past games
- - Methods to create a general database of game conditions
- Fitting module (fitting.py)
- - Methods to organize data in a way used by the ML code 
- - Methods to run (train) the network 
- - Methods to visualize and interpret the results 
- Prediction model (predict.py)
- - Methods to allow users to input players and get predicted points in a clean, OOP way.

#### Rubric for Revised Abstract
|                                                                      | Not present/unacceptable (0) | Needs Improvement (1) | Acceptable (2)  | Strong (3)  |
|----------------------------------------------------------------------|------------------------------|-----------------------|-----------------|-------------|
| Complete and on time.                                                |                              |                       |                 |             |
| Adjustments based on feedback have been implemented                  |                              |                       |                 |             |
| Skeleton sketch is feasible and sensible, and appropriately detailed |                              |                       |                 |             |

### Code Review I (18 possible points, 10% of grade)

*Code review* is a critical part of codebase creation and maintenance, and is carried out at all tech companies and in many scientific collaborations as well. For our code review, you will be meeting with one of the instructors for ~10 minutes. During this time, you are expected to show your progress (which by CR1 should be ~50% code started/constructed, with 25-50% functional). You will take us through the codebase as it stands, indicating the purpose of elements at roughly the class or function level (but not within). You should also be able to demo at least some of the core functionality of the code. Everyone will be required to make a repo for their project on Github and by the CR1 deadline the code you have so far should be pushed there. 

We will also ask you to explain how you did certain things, and you should be comfortable enough with your code to answer. ***This means, primarily, if you find a snippet on StackOverflow, make sure you take the time to figure out (or ask) how it does the thing.***

Finally, code review is also for the coder — we want (and require, see below) you to be able to articulate the challenges/roadblocks/bugs you’re facing, so we can provide advice that helps get you back on track.

In the example I provided above, this might look like showing the database and requests modules, demo-ing the functionality of pulling data and constructing the database. Maybe also (in a notebook/scratch space) showing some initial results from trying out ML training, but incomplete/not yet fully working.

#### Rubric for CR1
|                                                                       | Not present/unacceptable (0) | Needs Improvement (1) | Acceptable (2)  | Strong (3)  |
|-----------------------------------------------------------------------|------------------------------|-----------------------|-----------------|-------------|
| Completed by deadline                                                 |                              |                       |                 |             |
| ~50% of final codebase has been constructed                           |                              |                       |                 |             |
| Several core features are working (demo-d)                            |                              |                       |                 |             |
| Uses modular structure, functions, and classes and is pushed to GH    |                              |                       |                 |             |
| Can answer targeted questions about the codebase                      |                              |                       |                 |             |
| Articulates problem areas or elements not working and currently stuck |                              |                       |                 |             |

### Code Review II (18 possible points, 10% of grade)

Code Review II will mimic CR1. We will require that by CR2, you have met the minimum requirements of the project in terms of code length/complexity/structure (but not all of that code needs to be *working* by CR2). We want to be able to say “Good. If you get all this working, and hosted, you’ve met all the requirements of this project.” You should be able to demo all core features (saving convenience / for-fun features till the end). By this point, all code should be in its modular, near-final organization, and again, pushed to your project repo on Github. Finally, by CR2, we want you to show that your code is installed on your own system, and thus can be imported and used from anywhere.

#### Rubric for CR2 

|                                                                       | Not present/unacceptable (0) | Needs Improvement (1) | Acceptable (2)  | Strong (3)  |
|-----------------------------------------------------------------------|------------------------------|-----------------------|-----------------|-------------|
| Completed by deadline                                                 |                              |                       |                 |             |
| ~100% of final codebase has been constructed                          |                              |                       |                 |             |
| Almost all core features are working (can be demo-d)                  |                              |                       |                 |             |
| Uses modular structure, functions, and classes, and is installable    |                              |                       |                 |             |
| Can answer target questions about codebase                            |                              |                       |                 |             |
| Articulates problem areas or elements not working and currently stuck |                              |                       |                 |             |

### Peer Reivew + Final Presentation (21 possible points, 15% of grade) 

Final Presentation comprises both the finalized version of the project and codebase, as well as your presentation in front of the class and a peer review, in which we will pair you up with another person/group and ask you to install each other’s packages, and follow provided instructions to get them working.
```{note}
You should have both instructions for installation and for the use of your code on Github, either in the Readme file or in a doc/tutorial file).
```
Your project itself has been graded in several capacities in CR1 and CR2. By this point, we will check primarily that your code is ~fully functional (it’s ok if one or two things are still not there), and that it is not just hosted on github, but installable from there (i.e., your peer will download and install your code and attempt to use it). 

We will also check finally that all project requirements are met -- you should already have this by CR2, but if not, this is a chance to get some points back. 

The presentation itself will comprise you (and your partner if you have one) demonstrating your code. These presentations will be ~10 minutes, with the first 5 being a demo of the actual, working use of the code (what a user would see), and the other 5 being a quick tour of the codebase, highlighting 2-3 key areas that make your code tick. We can help you pick which parts of your code to present. 

#### Peer Review + Final Presentation Rubric

|                                                      | Not present/unacceptable (0) | Needs Improvement (1) | Acceptable (2)  | Strong (3)  |
|------------------------------------------------------|------------------------------|-----------------------|-----------------|-------------|
| Present for Peer Review                              |                              |                       |                 |             |
|                     **Project**                      |                              |                       |                 |             |
| Code complete + working+docs/instructions            |                              |                       |                 |             |
| Project hosted on github, installable by instructors |                              |                       |                 |             |
| Project Reqs met                                     |                              |                       |                 |             |
|                   **Presentation**                   |                              |                       |                 |             |
| Clarity of presentation                              |                              |                       |                 |             |
| Familiarity with presented code                      |                              |                       |                 |             |
| Demo works and demonstrates codebase’s function      |                              |                       |                 |             |