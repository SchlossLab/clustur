# Validate Count Table

If the count table is already valid nothing will change, otherwise it
will add a new group to the count table file.

## Usage

``` r
validate_count_table(count_table_df)
```

## Arguments

- count_table_df:

  The count table \`data.frame\` object.

## Value

A validated count table \`data.frame\`

## Details

Determines whether user supplied count table is valid

## Examples

``` r
 count_table <- read.delim(example_path("amazon.full.count_table"))
 count_table_valid <- validate_count_table(count_table)
```
