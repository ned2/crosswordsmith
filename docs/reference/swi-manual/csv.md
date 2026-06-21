
## A.11 library(csv): Process CSV (Comma-Separated Values) data

See also  
RFC 4180

To be done  
\- Implement immediate assert of the data to avoid possible stack overflows.  
- Writing creates an intermediate code-list, possibly overflowing resources. This waits for pure output!

This library parses and generates CSV data. CSV data is represented in Prolog as a list of rows. Each row is a compound term, where all rows have the same name and arity.

\[det\]**csv_read_file**(`+File, -Rows`)  
\[det\]**csv_read_file**(`+File, -Rows, +Options`)  
Read a CSV file into a list of rows. Each row is a Prolog term with the same arity. `Options` is handed to [csv//2](csv.html#csv//2). Remaining options are processed by [phrase_from_file/3](pio.html#phrase_from_file/3). The default separator depends on the file name extension and is `\t` for `.tsv` files and `,` otherwise.

Suppose we want to create a predicate table/6 from a CSV file that we know contains 6 fields per record. This can be done using the code below. Without the option `arity(6)`, this would generate a predicate table/N, where N is the number of fields per record in the data.

``` code
?- csv_read_file(File, Rows, [functor(table), arity(6)]),
   maplist(assert, Rows).
```

\[det\]**csv_read_stream**(`+Stream, -Rows, +Options`)  
Read CSV data from `Stream`. See also [csv_read_row/3](csv.html#csv_read_row/3).

\[det\]**csv**(`?Rows`)`//`  
\[det\]**csv**(`?Rows, +Options`)`//`  
Prolog DCG to‘read/write’CSV data. `Options`:

**separator**(`+Code`)  
The comma-separator. Must be a character code. Default is (of course) the comma. Character codes can be specified using the 0’notation. E.g., using `separator(0';)` parses a semicolon separated file.

**ignore_quotes**(`+Boolean`)  
If `true` (default false), threat double quotes as a normal character.

**strip**(`+Boolean`)  
If `true` (default `false`), strip leading and trailing blank space. RFC4180 says that blank space is part of the data.

**skip_header**(`+CommentLead`)  
Skip leading lines that start with `CommentLead`. There is no standard for comments in CSV files, but some CSV files have a header where each line starts with `#`. After skipping comment lines this option causes [csv//2](csv.html#csv//2) to skip empty lines. Note that an empty line may not contain white space characters (space or tab) as these may provide valid data.

**convert**(`+Boolean`)  
If `true` (default), use [name/2](manipatom.html#name/2) on the field data. This translates the field into a number if possible.

**case**(`+Action`)  
If `down`, downcase atomic values. If `up`, upcase them and if `preserve` (default), do not change the case.

**functor**(`+Atom`)  
Functor to use for creating row terms. Default is `row`.

**arity**(`?Arity`)  
Number of fields in each row. This predicate raises a `domain_error(row_arity(Expected), Found)` if a row is found with different arity.

**match_arity**(`+Boolean`)  
If `false` (default `true`), do not reject CSV files where lines provide a varying number of fields (columns). This can be a work-around to use some incorrect CSV files.

\[nondet\]**csv_read_file_row**(`+File, -Row, +Options`)  
True when `Row` is a row in `File`. First unifies `Row` with the first row in `File`. Backtracking yields the second, ... row. This interface is an alternative to [csv_read_file/3](csv.html#csv_read_file/3) that avoids loading all rows in memory. Note that this interface does not guarantee that all rows in `File` have the same arity.

In addition to the options of [csv_read_file/3](csv.html#csv_read_file/3), this predicate processes the option:

**line**(`-Line`)  
`Line` is unified with the 1-based line-number from which `Row` is read. Note that `Line` is not the physical line, but rather the *logical* record number.

\[det\]**csv_read_row**(`+Stream, -Row, +CompiledOptions`)  
Read the next CSV record from `Stream` and unify the result with `Row`. `CompiledOptions` is created from options defined for [csv//2](csv.html#csv//2) using [csv_options/2](csv.html#csv_options/2). `Row` is unified with `end_of_file` upon reaching the end of the input.

\[det\]**csv_options**(`-Compiled, +Options`)  
`Compiled` is the compiled representation of the CSV processing options as they may be passed into [csv//2](csv.html#csv//2), etc. This predicate is used in combination with [csv_read_row/3](csv.html#csv_read_row/3) to avoid repeated processing of the options.

\[det\]**csv_write_file**(`+File, +Data`)  
\[det\]**csv_write_file**(`+File, +Data, +Options`)  
Write a list of Prolog terms to a CSV file. `Options` are given to [csv//2](csv.html#csv//2). Remaining options are given to [open/4](IO.html#open/4). The default separator depends on the file name extension and is `\t` for `.tsv` files and `,` otherwise.

\[det\]**csv_write_stream**(`+Stream, +Data, +Options`)  
Write the rows in `Data` to `Stream`. This is similar to [csv_write_file/3](csv.html#csv_write_file/3), but can deal with data that is produced incrementally. The example below saves all answers from the predicate data/3 to File.

``` code
save_data(File) :-
   setup_call_cleanup(
       open(File, write, Out),
       forall(data(C1,C2,C3),
              csv_write_stream(Out, [row(C1,C2,C3)], [])),
       close(Out)).
```
