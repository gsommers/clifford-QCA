# clifford-QCA


Julia code and data used in [Crystalline Quantum Circuits](https://arxiv.org/abs/2210.10808) by Grace M. Sommers, David A. Huse, and Michael J. Gullans, for expressing spacetime translation-invariant unitary Clifford circuits as Clifford quantum cellular automata (CQCA). Functionality includes converting a unit cell of Clifford gates into the Laurent polynomial matrix form, getting recurrence times and operator spreading, and plotting the image of initially local Paulis.

Additional code available upon request. Please contact <gsommers@princeton.edu>.

To use the modules included in this repo, execute the following from the command line:

```export JULIA_LOAD_PATH="~/path/to/repo:"```

or, inside a Julia session, do:

```push!(LOAD_PATH, "~/path/to/repo")```

For example, inside the demo notebook, you would do
```push!(LOAD_PATH, "../")```

Then you can access the functions in these modules through the command `using <module name>`.

Required packages
-----------------
  - `IJulia` (to run demo notebook)
  