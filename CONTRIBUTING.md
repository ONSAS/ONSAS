# Contributing to ONSAS

First off, thank you for considering contributing to ONSAS! It's people like you that make this open-source solver a great tool for the structural engineering and academic community.

This document provides guidelines and instructions for contributing to this repository.

## Table of Contents
1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [How to Contribute](#how-to-contribute)
   - [Reporting Bugs](#reporting-bugs)
   - [Suggesting Enhancements](#suggesting-enhancements)
   - [Submitting Pull Requests](#submitting-pull-requests)
4. [Development Environment Setup](#development-environment-setup)
5. [Coding Style and Guidelines](#coding-style-and-guidelines)
6. [Testing](#testing)

## Code of Conduct
By participating in this project, you are expected to uphold a welcoming and collaborative environment. Please treat all maintainers, contributors, and users with respect.

## Getting Started

Before you begin contributing, please ensure you have a basic understanding of how ONSAS operates. You can familiarize yourself with the solver by checking out the `examples` directory and running some of the predefined structural analysis cases (e.g., the static von mises truss, propeller, or the Euler column model).

## How to Contribute

### Reporting Bugs
If you find a bug in the code, please open an issue in the GitHub repository. When creating an issue, please include:
* A clear and descriptive title.
* The version of ONSAS you are using.
* Your environment details (GNU-Octave or MATLAB version, operating system).
* Steps to reproduce the behavior, ideally with a minimal reproducible script.
* Any relevant error messages or screenshots.

### Suggesting Enhancements
If you have ideas to improve the solver—whether it's a new finite element formulation, a new integration method, or performance optimizations—please open an issue to discuss it before you start coding. This ensures your proposed changes align with the project's roadmap and current academic goals.

### Submitting Pull Requests
1. **Fork the repository** and create your branch from `master`.
2. **Create a descriptive branch name** (e.g., `feature/new-element`, `bugfix/stiffness-matrix`).
3. **Make your changes** following the coding guidelines.
4. **Ensure the test suite passes** completely. 
5. **Update documentation** if your changes affect how users interact with the solver.
6. **Submit a Pull Request** with a clear title and description of the changes you've made. Link any relevant issues.

## Development Environment Setup

To set up a local development environment for ONSAS:

1. Clone your forked repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/ONSAS.git](https://github.com/YOUR_USERNAME/ONSAS.git)
   cd ONSAS
   ```
2. Ensure you have **GNU-Octave** or **MATLAB** installed on your system. 
3. *(Optional but recommended)* Install **Paraview** for visualizing output files and **Gmsh** if you are working on meshing components.
4. Add the `src` folder (and relevant subdirectories) to your MATLAB/Octave path.

## Coding Style and Guidelines

Maintaining a consistent codebase is crucial for a project of this scale. 

* **Language Compatibility**: The code must remain fully compatible with both **GNU-Octave** and **MATLAB**. Avoid using toolboxes or functions that are exclusive to only one of the environments.
* **Formatting and Linting**: The project uses `miss_hit` for code quality and style checking. Please ensure your code complies with the rules defined in the `miss_hit.cfg` file located in the root directory.
* **Documentation**: Document all new functions thoroughly. Include comments explaining the mathematical or structural formulations used, citing relevant literature when applicable. 
* **Variable Naming**: Use descriptive names for variables and matrices (e.g., `stiffnessMatrix`, `nodalDisplacements`) to maintain readability.

## Testing

ONSAS uses an automated testing suite to prevent regressions. Before submitting a pull request, you must ensure that all tests pass.

1. Navigate to the `test` directory within your Octave/MATLAB environment.
2. Run the master test script.
3. If you are adding a new feature or fixing a bug, please add a corresponding test case to the `test` directory to verify your changes.
4. Check that your changes do not break the continuous integration (CI) workflows defined in the `.github/workflows` directory.

***

Thank you for contributing to ONSAS and helping advance open-source structural analysis!


