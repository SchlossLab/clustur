# Read count table

This function will read and return your count table. It can take in
sparse and full count tables.

## Usage

``` r
read_count(count_table_path)
```

## Arguments

- count_table_path:

  The file path of your count table.

## Value

a count table \`data.frame\`.

## Examples

``` r
count_table <- read_count(example_path("amazon.full.count_table"))
```
