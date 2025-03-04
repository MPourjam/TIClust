# Version 1.0.5

Minor bugs in CLI arguments resolved.
On failure `uclust working directory | uclust_wd` gets deleted if existed.

# Version 1.0.4

- Arguments passed to CLI gets converts to expected type.

  - *threads*: int
  - *species_thr*: float
  - *genus_thr*: float
  - *family_thr*: float

- F/G/SOTU IDs take 0 as a separator between parent and children IDs
- CLI gets `--version/-v` argument.