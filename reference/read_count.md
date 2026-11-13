# Read count table

This function will read and return your count table. It can take in
sparse and full count tables.

## Usage

``` r
read_count(count_table)

# S3 method for class 'strollur'
read_count(count_table)

# S3 method for class 'character'
read_count(count_table)
```

## Arguments

- count_table:

  The file path of your count table or a strollur object.

## Value

a count table \`data.frame\`.

## Examples

``` r
count_table <- read_count(example_path("amazon.full.count_table"))
```
