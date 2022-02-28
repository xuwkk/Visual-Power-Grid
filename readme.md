# Power System Visualization Demo: On Going

![visualization](visualization.gif)

## Introduction
This repository contains a visualization tool on steady state power system operation, cyber attacks on steady state operation, and defence methods. The visualisation repies on python package [PyPower](https://github.com/rwl/PYPOWER) (installization required)  and [Tkinter](https://docs.python.org/3/library/tkinter.html) (embedded in Python). 

Thre current functions include:

1. Infinitely run the optimal power flow (OPF) and display the active line flow measurement.
2. Do steady state AC state estimation and return the residual for bad data detection (BDD) usage.
3. User's choice to perform two attack types on the grid normal operation, namely 'Random Attack' and 'FDI Attack'.
4. Display both the attacked and normal measurment for comparison.
5. Display the residual plot for past OPF runs.

The functionalities on power system operation (e.g. opf, se, and bdd) are wrapped from our repository  [steady-state-power-system](https://github.com/xuwkk/steady-state-power-system) which maintains different algorithms to solve state estimation and detection algorithms on FDI attacks.

## PyPower

**PYPOWER** is a power flow and Optimal Power Flow (OPF) solver. It is a port of [MATPOWER](http://www.pserc.cornell.edu/matpower/) to the [Python](http://www.python.org/) programming language. Current features include:

- DC and AC (Newton's method & Fast Decoupled) power flow and
- DC and AC optimal power flow (OPF)

## Tkinter

The **Tkinter** package (“Tk interface”) is the standard Python interface to the Tcl/Tk GUI toolkit. To learn Tkinter, please refer to one of the following tutorial:

1. [Tkinter Course - Create Graphic User Interfaces in Python Tutorial](https://www.youtube.com/watch?v=YXPyB4XeYLA&t=29s) from [freeCodeCamp.org](https://www.youtube.com/channel/UC8butISFwT-Wl7EV0hUK0BQ) (beginner level).
2. [Python Tkinter](https://www.youtube.com/playlist?list=PL6lxxT7IdTxGoHfouzEK-dFcwr_QClME_) from [John Philip Jones](https://www.youtube.com/c/johnphilipjones) (beginner level).
3. [GUIs with Tkinter](https://www.youtube.com/playlist?list=PLQVvvaa0QuDclKx-QpC9wntnURXVJqLyk) from [sentdex](https://www.youtube.com/c/sentdex) (intermidiante lavel).
