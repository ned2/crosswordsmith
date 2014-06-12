crosswordsmith
===============================================================================

A crossword layout generator in Prolog. 

This is a crossword grid generator implemented in SWI Prolog, however
it should work in most Prologs with minimal porting. The solver works
deterministically, laying down one word at a time by placing it such
that it intersects with an existing word, and also such that it
doesn't sit hard up against any letters from another word.

The solver is designed to run online, so it returns the first solution
rather than trying to find them all.  The list of input words is
shuffled first so repeated solutions on different runs are less
likely. The start location for laying down the first word can also be
shuffled or specified in advance. The algorithm does not try guess the
grid length, so this must be specified by the user. It is of course
possible to specify a length for which there are no solutions. Also
note that duplicate solutions exist. These consist of the same layout
of words, but found by laying out words in a different order.

The words and their accompanying solutions are stored in a separate
file, clues.pl, which also allow you to attach a URL to each word, as
this was originally intented for embedding in a web page, where the
clues/solved words can be links.


USAGE
===============================================================================


CLUES FILE
===============================================================================