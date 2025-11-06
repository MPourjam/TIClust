# Version 1.0.9

- edge taxonomy cases like "tax=;;", "tax=tax=Bacteria;;;Firmicutes;;", and "Zotu111 Bacteria;;;Firmicutes;;" are handled safely


# Version 1.0.8

- `*-All-FullTaxonomy.tab` files will be created only if there are non-bacterial ZOTUs in the input ZOTU table.
- Full-Taxonomy ZOTUs tables file names will have the original ZOTU table name as a prefix.
  

# Version 1.0.7

- Taxonomies in ZOTU table, if given, will be updated and written to the output directories
  - `ZOTU-Table-All-FullTaxonomy.tab`: A version of given ZOTU table with full taxonomy (7 levels) for all ZOTUs.
  - `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab`: A version of given ZOTU table with full taxonomy (7 levels) for only Bacteria ZOTUs.

# Version 1.0.6

- TICAnalysis will take a logger object to replace with default logger.
- Logger format changed to logg file and line of log origin

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
