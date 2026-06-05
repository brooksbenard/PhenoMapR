# Extract ICI PRECOG Column

Resolve a user-supplied `cancer_type` to a single ICI-PRECOG column. The
ICI reference's column names are full cohort labels (e.g.
`"LUAD_PD-L1_Primary_Naive"`, `"KIRC_PD1_Metastatic_Naive"`), and
[`list_cancer_types`](https://brooksbenard.github.io/PhenoMapR/reference/list_cancer_types.md)`("ici_precog")`
returns those labels verbatim, so most callers pass an exact column
name. We accept that directly. A legacy short-form (e.g. `"LUAD"` or
`"LUAD_Metastatic"`) is also still supported via a regex fallback, which
picks the first matching column and warns when there are ties.

## Usage

``` r
extract_ici_column(ici_data, cancer_type)
```
