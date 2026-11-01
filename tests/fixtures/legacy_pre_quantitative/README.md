# Fixed legacy regression source

These three files are exact `git show` copies from revision
`c85821a833209b45384c0a8835e5d2ead8051daf`, before Phase 8 dispatch integration.
The test executes its headless/phantom implementation against the current SignalModel
with identical explicit inputs, NumPy/Torch seeds and nonzero components. SignalModel
ASTs must match after removing imports; Phase 9 only delays its MatlabRunner import
until actual lipid generation. External anatomy, basis data and MM/lipid providers
are controlled small fixtures; equations, residual water, shim and noise run normally.

All complex total/component arrays, concentrations, axes and relevant metadata are
required to match exactly. This baseline does not require Git history at test time.
It establishes compatibility for the exercised configuration, not every possible
external phantom, lipid runner or platform.
