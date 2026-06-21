
## A.46 library(quasi_quotations): Define Quasi Quotation syntax

author  
Jan Wielemaker. Introduction of Quasi Quotation was suggested by Michael Hendricks.

See also  
\- [Why it's nice to be quoted: quasiquoting for haskell](http://www.cs.tufts.edu/comp/150FP/archive/geoff-mainland/quasiquoting.pdf)  
- [Why it's nice to be quoted: quasiquoting for Prolog](https://www.swi-prolog.org/download/publications/quasiquoting.pdf)

Inspired by [Haskell](http://www.haskell.org/haskellwiki/Quasiquotation), SWI-Prolog support *quasi quotation*. Quasi quotation allows for embedding (long) strings using the syntax of an external language (e.g., HTML, SQL) in Prolog text and syntax-aware embedding of Prolog variables in this syntax. At the same time, quasi quotation provides an alternative to represent long strings and atoms in Prolog.

The basic form of a quasi quotation is defined below. Here, `Syntax` is an arbitrary Prolog term that must parse into a *callable* (atom or compound) term and Quotation is an arbitrary sequence of characters, not including the sequence `|}`. If this sequence needs to be embedded, it must be escaped according to the rules of the target language or the‘quoter’must provide an escaping mechanism.

``` code
{|Syntax||Quotation|}
```

While reading a Prolog term, and if the Prolog flag `quasi_quotations` is set to `true` (which is the case if this library is loaded), the parser collects quasi quotations. After reading the final full stop, the parser makes the call below. Here, `SyntaxName` is the functor name of `Syntax` above and `SyntaxArgs` is a list holding the arguments, i.e., `Syntax =.. [SyntaxName|SyntaxArgs]`. Splitting the syntax into its name and arguments is done to make the quasi quotation parser a predicate with a consistent arity 4, regardless of the number of additional arguments.

``` code
call(+SyntaxName, +Content, +SyntaxArgs, +VariableNames, -Result)
```

The arguments are defined as

- `SyntaxName` is the principal functor of the quasi quotation syntax. This must be declared using [quasi_quotation_syntax/1](quasiquotations.html#quasi_quotation_syntax/1) and there must be a predicate SyntaxName/4.

- `Content` is an opaque term that carries the content of the quasi quoted material and position information about the source code. It is passed to with_quasi_quote_input/3.

- `SyntaxArgs` carries the additional arguments of the `Syntax`. These are commonly used to make the parameter passing between the clause and the quasi quotation explicit. For example:

  ``` code
      ...,
      {|html(Name, Address)||
       <tr><td>Name<td>Address</tr>
       |}
  ```

- `VariableNames` is the complete variable dictionary of the clause as it is made available throug [read_term/3](termrw.html#read_term/3) with the option `variable_names`. It is a list of terms `Name = Var`.

- `Result` is a variable that must be unified to resulting term. Typically, this term is structured Prolog tree that carries a (partial) representation of the abstract syntax tree with embedded variables that pass the Prolog parameters. This term is normally either passed to a predicate that serializes the abstract syntax tree, or a predicate that processes the result in Prolog. For example, HTML is commonly embedded for writing HTML documents (see `library(http/html_write)`). Examples of languages that may be embedded for processing in Prolog are SPARQL, RuleML or regular expressions.

The file `library(http/html_quasiquotations)` provides the, suprisingly simple, quasi quotation parser for HTML.

\[det\]**with_quasi_quotation_input**(`+Content, -Stream, :Goal`)  
Process the quasi-quoted `Content` using `Stream` parsed by `Goal`. `Stream` is a temporary stream with the following properties:

- Its initial *position* represents the position of the start of the quoted material.
- It is a text stream, using `utf8` *encoding*.
- It allows for repositioning
- It will be closed after `Goal` completes.

|  |  |
|----|----|
| `Goal` | is executed as `once(Goal)`. `Goal` must succeed. Failure or exceptions from `Goal` are interpreted as syntax errors. |

See also  
[phrase_from_quasi_quotation/2](quasiquotations.html#phrase_from_quasi_quotation/2) can be used to process a quotation using a grammar.

\[det\]**phrase_from_quasi_quotation**(`:Grammar, +Content`)  
Process the quasi quotation using the DCG `Grammar`. Failure of the grammar is interpreted as a syntax error.

See also  
[with_quasi_quotation_input/3](quasiquotations.html#with_quasi_quotation_input/3) for processing quotations from stream.

\[det\]**quasi_quotation_syntax**(`:SyntaxName`)  
Declare the predicate `SyntaxName`/4 to implement the the quasi quote syntax `SyntaxName`. Normally used as a directive.

**quasi_quotation_syntax_error**(`+Error`)  
Report `syntax_error(Error)` using the current location in the quasi quoted input parser.

throws  
`error(syntax_error(Error), Position)`
