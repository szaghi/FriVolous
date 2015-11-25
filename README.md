<a name="top"></a>

# FriVolous [![GitHub tag](https://img.shields.io/github/tag/szaghi/FriVolous.svg)]() [![Join the chat at https://gitter.im/szaghi/FriVolous](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/szaghi/FriVolous?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]()
[![Build Status](https://travis-ci.org/szaghi/FriVolous.svg?branch=master)](https://travis-ci.org/szaghi/FriVolous)
[![Coverage Status](https://img.shields.io/codecov/c/github/szaghi/FriVolous.svg)](http://codecov.io/github/szaghi/FriVolous?branch=master)

### FriVolous, Finite Volume block-structured Fortran abstract class

A KISS pure Fortran OOD class providing an abstract container for Finite Volume block-structured numerical computations

- FriVolous is a pure Fortran (KISS) library to aid Finite Volume block-structured computations of numerical spatial operators;
- FriVolous is Fortran 2003+ standard compliant;
- FriVolous is OOP designed;
- FriVolous is a Free, Open Source Project.

#### Table of Contents

- [What is FriVolous?](#what-is-FriVolous)
- [Main features](#main-features)
- [Copyrights](#copyrights)
- [Documentation](#documentation)
  - [A Taste of FriVolous](#a-taste-of-FriVolous)

#### Issues

[![GitHub issues](https://img.shields.io/github/issues/szaghi/FriVolous.svg)]()
[![Ready in backlog](https://badge.waffle.io/szaghi/FriVolous.png?label=ready&title=Ready)](https://waffle.io/szaghi/FriVolous)
[![In Progress](https://badge.waffle.io/szaghi/FriVolous.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/FriVolous)
[![Open bugs](https://badge.waffle.io/szaghi/FriVolous.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/FriVolous)

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-v4.9.2+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/Intel-v12.x+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/g95-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/NAG-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/PGI-not%20tested-yellow.svg)]()

## What is FriVolous?

FriVolous is a user-friendly and Object-Oriented designed class for Finite Volume block-structured numerical computations. In particular, FriVolous allows the easy handling of *metrics* data for the *robust* and efficient computation of numerical **spatial operators** in the framework of Finite Volume Methods (**FVM**). It is based on a simple yet powerful Abstract Data Type (ADT) that is *overloaded* with useful methods for handling the back-end operations necessary to compute a numerical spatial operator independently by the actual *state variable* being numerically integrated.

FriVolous adheres to the [KISS](https://en.wikipedia.org/wiki/KISS_principle) concept.

Go to [Top](#top)

## Main features

+ [x] Pure Fortran implementation;
+ [ ] KISS and user-friendly:
    + [ ] simple API (one main *object*);
    + [ ] easy building and porting on heterogeneous architectures;
+ [ ] comprehensive:
    + [ ] 1D topology;
    + [ ] 2D topology;
    + [ ] 3D topology;
+ [ ] efficient and *non intrusive* (all object methods and operators are *pure* or *elemental*):
    + [ ] threads/processes safe;
+ [x] Tests-Driven Developed ([TDD](https://en.wikipedia.org/wiki/Test-driven_development));
+ [ ] well documented:
    + [ ] complete [API](http://szaghi.github.io/FriVolous/index.html) reference;
    + [ ] comprehensive [wiki](https://github.com/szaghi/FriVolous/wiki):
+ [x] collaborative developed on [GitHub](https://github.com/szaghi/FriVolous);
+ [x] [FOSS licensed](https://github.com/szaghi/FriVolous/wiki/Copyrights);

Any feature request is welcome.

Go to [Top](#top)

## Copyrights

FriVolous is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to FriVolous is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/szaghi/FriVolous/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the FriVolous documentation is contained into its own [wiki](https://github.com/szaghi/FriVolous/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://szaghi.github.io/FriVolous/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

### A Taste of FriVolous

To be written.

Go to [Top](#top)
