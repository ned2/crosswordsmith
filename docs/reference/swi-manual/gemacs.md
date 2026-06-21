
## 2.6 GNU Emacs Interface

SWI-Prolog provides tight integration with GNU Emacs through the `sweep` package. This package embeds SWI-Prolog as a dynamic Emacs module, allowing for Prolog queries to be executed directly from Emacs Lisp. The accompanying Emacs package `sweeprolog.el`, available for installation with the standard Emacs package manager `package.el`, builds on top of this embedding to provide a fully integrated development environment for SWI-Prolog in GNU Emacs.

GNU Emacs ships with by default with a Prolog mode called `prolog.el`. Compared to `sweeprolog.el`, this mode suffers from some problems that arise due to the lack of a proper Prolog parser. The original `prolog.el` by Masanobu Umeda has been included in GNU Emacs since 1989, in 2006 Stefan Monnier added explicit support for SWI-Prolog to `prolog.el`. In 2011, most of the original implementation has been replaced with a new Prolog mode written by initially for the XEmacs port by Stefan Bruda. Bruda's mode was adapted to GNU Emacs by Stefan Monnier, who has been maintaining it along with other GNU Emacs contributor since. Users of this mode may find useful configuration suggestions at [https://www.metalevel.at/pceprolog/](https://www.metalevel.at/pceprolog/).

Other Emacs package that can be useful for working with SWI-Prolog are:

- [https://www.metalevel.at/ediprolog/](https://www.metalevel.at/ediprolog/)  
  Interact with SWI-Prolog directly in Emacs buffers.
- [https://www.metalevel.at/etrace/](https://www.metalevel.at/etrace/)  
  Trace Prolog code with Emacs.
- [https://emacs-lsp.github.io/dap-mode/page/configuration/#swi-prolog](https://emacs-lsp.github.io/dap-mode/page/configuration/#swi-prolog)  
  Debug Adapter Protocol (DAP) support for SWI-Prolog in Emacs via `dap-mode` and the `debug_adapter` pack from [https://github.com/eshelyaron/debug_adapter](https://github.com/eshelyaron/debug_adapter)
- [https://emacs-lsp.github.io/lsp-mode/page/lsp-prolog/](https://emacs-lsp.github.io/lsp-mode/page/lsp-prolog/)  
  Language Server Protocol (LSP) support for SWI-Prolog in Emacs via `lsp-mode` and the `lsp_server` pack from [https://github.com/jamesnvc/lsp_server](https://github.com/jamesnvc/lsp_server)
