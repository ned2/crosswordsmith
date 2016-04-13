#!/usr/bin/swipl

:- set_prolog_flag(verbose, silent). 
:- initialization main. 

% crossword.pl - A crossword layout generator in Prolog
% Copyright (C) 2011  Ned Letcher - nedned.net
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
% the Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


% See the file README.mk for background and a general overview..
% 
% This program can be run as an SWI PrologScript like this:
%
% $ ./crossword.pl <grid_length> <start_loc>
% 
% Where grid_length is integer specifying the dimensions of the
% crossword and start_loc specifies where the first word of the
% crossword is placed and can be one of {rand_loc, topleft_down,
% topleft_across, topright, bottomleft}.


% Data structures:
%
% The program uses two simple data structures. The first is a list of
% words that have been placed on the crossword grid with each element
% being a list of attributes of the word. See the placed words utility
% predicate section below for the structure of this list.
%
% The second data structure is as association list used to to store
% the contents of each cell in the grid. The keys are the numbers of
% the cells (1 through to GridLen*GridLen) and the values are the
% contents of the cell.

% Load a file that contains a predicate 'clues' which has a single
% parameter which is a list of clues to be used in the crossword, with
% each clue being a list of the following form [word, clue, link].  It
% is ok to include spaces in the word slot for multi-word answers.
% Link is in case you want to spcify an accompanying URL for the word.
:- include('clues.pl').



% program predicates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main predicate, called when run as a PrologScript.
main :-
    current_prolog_flag(argv, Argv),
    (
     Argv = ['--all', GridLenArg, StartLoc]
    ;
     Argv = ['--all', GridLenArg]
    ),
    atom_number(GridLenArg, GridLen),
    clues(Words),
    all_crossword(GridLen, Words, StartLoc, Num),
    writeln(Num).


main :-
    current_prolog_flag(argv, Argv), 
    clues(ClueWords),
    (
     Argv = ['--shuffle', GridLenArg, StartArg] ->
     % shuffle the input words
     shuffle(ClueWords, UseWords)
     ;
     Argv = [GridLenArg, StartLoc],
     UseWords = ClueWords
    ),
    atom_number(GridLenArg, GridLen),
    ( 
      StartArg == 'rand_loc' ->
      % choose a start location at random
      start_locs(Locs),
      shuffle(Locs, ShuffledLocs),
      member(StartLoc, ShuffledLocs)
    ;    
      % just use supplied start location
      StartLoc = StartArg
    ),
    crossword(GridLen, UseWords, StartLoc),
    halt.


% Top level predicate for solving the crossword with a specified
% starting position.
crossword(GridLen, Words, StartLoc) :-
    find_crossword(GridLen, Words, StartLoc, _Grid, PlacedWords),
    assign_clue_numbers(PlacedWords, NumberedPlacedWords),
    annotate_grid(NumberedPlacedWords, GridLen, Grid),
    print_grid(Grid, GridLen),
    write('\n@\n\n'),
    print_clues(NumberedPlacedWords).


% Top level predicate for finding the number of solutions for the
% crossword for a specific starting position.
all_crossword(GridLen, Words, StartLoc, Num) :-
    length(Sols, Num),
    findall(Grid, find_crossword(GridLen, Words, StartLoc, Grid, _), Sols).


% The driver predicate used to solve the crossword.
find_crossword(GridLen, Words, Loc, Grid, PlacedWords) :-
    init_grid(GridLen, G1),     
    % Get the cell number and direction for start loc
    start_loc(Loc, GridLen, StartNum, StartDir),
    assign_words(Words, [], GridLen, StartNum, StartDir, G1, Grid, PlacedWords).


% Assign all words. The starting location is selected by locating an
% intersecting word from the words already placed.
assign_words([], P, _, _, _, G, G, P).
assign_words(Words, PlacedWords, GridLen, Start, Dir, GIn, GOut, PlacedWordsOut) :-
    member([Word, Clue, Link], Words),
    atom_chars(Word, Letters),
    delete(Letters, ' ', Letters2),
    length(Letters2, WLen),
    % on first pass, Start and Dir will be grounded with the start values
    % then afterwards will be unground, with find_intersecting_word grounding them
    find_intersecting_word(Letters2, WLen, PlacedWords, GridLen, Start, Dir),
    assign_word(Word, Letters2, WLen, Clue, Link, Start, Dir, GridLen, GIn, Placed, G1),
    remove_x([Word, Clue, Link], Words, RemWords),
    assign_words(RemWords, [Placed|PlacedWords], GridLen, _Start, _Dir, G1, GOut, PlacedWordsOut).


% Given a Word and a set of Placed words, locates a candidate 
% start cell and direction that intersects with an existing word.

% No placed words; just use grounded Start and Dir values
find_intersecting_word(_Letters, _WLen, [], _GridLen, _Start, _Dir).

find_intersecting_word(Letters, WLen, PlacedWords, GridLen, Start, Dir) :-
    member([_, _, _, PLetters, _, PDir, _, PStart, _], PlacedWords),
    intersection(Letters, PLetters, Vals),
    list_to_set(Vals, Vals2),
    member(Val, Vals2),
    position(Val, PLetters, PPos),
    position(Val, Letters, Pos),
    calc_num(PDir, GridLen, PPos, PStart, PNum),
    swap_dir(PDir, Dir),
    calc_start(Dir, GridLen, Pos, PNum, Start),
    fits_on_grid(Dir, Start, WLen, GridLen).


assign_word(Word, Letters, WLen, Clue, Link, Start, Dir, GridLen, GIn, Placed, GOut) :-
    % make sure previous cell does not have a letter
    check_prev_cell(Dir, Start, GridLen, GIn),
    assign_letters(Letters, Start, Dir, GridLen, Cells, GIn, GOut),
    Placed = [Word, Clue, Link, Letters, Cells, Dir, WLen, Start, _ClueNum].


% Previous cell before start of word. Make sure it doesn't contain
% anything.
check_prev_cell(Dir, Num, GridLen, G) :-
    (
     % don't check if start letter is start of a row/col
     is_start_cell(Dir, Num, GridLen)
    ;
     % otherwise prev cell must be empty
     prev_cell(Dir, Num, GridLen, Prev),
     get_assoc(Prev, G, empty)
    ), !.


% Next cell after end of word. Make sure it doesn't contain anything.
check_next_cell(Dir, Num, GridLen, G) :-
    prev_cell(Dir, Num, GridLen, Prev),
    (
     % no need to check if prev was end of row/col
     is_end_cell(Dir, Prev, GridLen)
    ;
     % then this cell must be empty
     get_assoc(Num, G, empty)
    ), !.


% Assign each letter of the word, checking that adjacent cells
% are empty to prevent words being placed next to each other.

% Last letter of word, make sure next cell is free
assign_letters([], Num, Dir, GridLen, [], G, G) :- 
    check_next_cell(Dir, Num, GridLen, G).


assign_letters([L|Ls], Num, Dir, GridLen, [Num|RestCells], GIn, GOut) :-
    get_assoc(Num, GIn, X),
    (
     % existing letter in this cell matches letter being placed,
     % nothing needs doing, we can continue to next letter
     X == L,
     G1 = GIn
    ;
     % no letter in this cell, so check adjacent cells are free
     % and then add letter to this cell
     X == empty,
     adj_is_free(Dir, Num, GridLen, GIn),
     put_assoc(Num, GIn, L, G1)
    ), !,
    next_cell(Dir, Num, GridLen, Num2),
    assign_letters(Ls, Num2, Dir, GridLen, RestCells, G1, GOut).


% check that adjacent cells are empty
adj_is_free(down, Num, GridLen, G) :-
    N1 is Num - 1,
    N2 is Num + 1,
    M is Num mod GridLen,
    (
     M == 0 -> % last cell in row
     get_assoc(N1, G, empty)
    ;
     M == 1 -> % first cell in row
     get_assoc(N2, G, empty)
    ;
     get_assoc(N1, G, empty),
     get_assoc(N2, G, empty)
    ), !.

adj_is_free(across, Num, GridLen, G) :-
    N1 is Num - GridLen,
    N2 is Num + GridLen,
    LastCell is (GridLen * GridLen),
    (
     N1 =< 0 -> % before beginning of grid
     get_assoc(N2, G, empty)
    ;
     N2 > LastCell -> % after end of grid
     get_assoc(N1, G, empty)
    ;
     get_assoc(N1, G, empty),
     get_assoc(N2, G, empty)
    ), !.


% Takes the placed words and works out the clue numbers of each
% clue/word. This works by interating through the list of words which
% are sorted by their start numbers.  Clue numbers are then assigned
% in the order that words occur in this list.
assign_clue_numbers(PlacedWords, WordsClues) :-
    % sort placed words by their start number
    map_list_to_pairs(start_is, PlacedWords, Pairs),
    keysort(Pairs, SortedPairs),

    % group clues that have the same start num together (ie a cell
    % that is start of both across and down clue) as these will
    % recieve the same clue number
    group_pairs_by_key(SortedPairs, GroupedPairs),

    % assign the clue numbers, starting at 1
    add_clue_nums(GroupedPairs, 1, WordsClues).


% Updates the placed words list by appending the clue number to the end
% of each word. The list now looks like this:
% placed words -- [word, letters, cells, dir, len, start, clue_num] 

add_clue_nums([], _, []).    

% a word whose start cell only belongs to a down or an across word
add_clue_nums([_-[W]|Rest], ClueNum, [WClue|RestClues]) :-    
    add_clue_word(W, ClueNum, WClue),
    ClueNum2 is ClueNum + 1,
    add_clue_nums(Rest, ClueNum2, RestClues).    
    
% a word whose start cell belongs to both a down and an across word
add_clue_nums([_-[W1,W2]|Rest], [WClue1,WClue2|RestClues]) :-
    add_clue_word(W1, ClueNum, WClue1),
    add_clue_word(W2, ClueNum, WClue2),
    ClueNum2 is ClueNum + 1,
    add_clue_nums(Rest, ClueNum2, RestClues).    


% Create a new grid containing the final values for each cell
annotate_grid(PlacedWords, GridLen, G) :-
    init_grid(GridLen, NewG),
    words_to_grid(PlacedWords, NewG, G).


words_to_grid([], G, G).
words_to_grid([W|Ws], GIn, GOut) :-
    W = [_, _, Link, Letters, Cells, Dir, _, Start, ClueNum],
    letters_to_grid(Letters, Cells, Dir, Start, ClueNum, Link, GIn, G1),
    words_to_grid(Ws, G1, GOut).


% this code could be cleaned up a bit
letters_to_grid([], [], _, _, _, _, Grid, Grid).
letters_to_grid([L|Ls], [C|Cs], Dir, Start, ClueNum, Link, GIn, GOut) :-
    get_assoc(C, GIn, X),
    (
     Start == C ->
     (
      Dir == across,
      Print = 'a'
     ;
      Dir == down,
      Print = 'd'
     ), !
    ;
     Print = 'n'
    ),
    (
     X == empty ->
     (
      Dir == across,
      Val = [ClueNum-x, L, Print, Link]
     ;
      Dir == down,
      Val = [x-ClueNum, L, Print, Link]
     ), !     
    ;
     X = [A-D,L,_P,CurrLink],
     (
      Dir == across,
      Val = [ClueNum-D, L, Print, CurrLink]
     ;
      Dir == down,
      Val = [A-ClueNum, L, Print, Link]
     ), !
    ),
    put_assoc(C, GIn, Val, G1),
    letters_to_grid(Ls, Cs, Dir, Start, ClueNum, Link, G1, GOut).



% Placed Word list utility predicates:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A placed word looks like this: 
% [Word, Clue, Link, Letters, Cells, Dir, Len, Start, ClueNum]


dir_is_across([_,_,_,_,_,across,_,_,_]).

start_is([_,_,_,_,_,_,_,Start,_], Start).

% Assign clue number to it's location in the placed word list  
add_clue_word(Word, ClueNum, WordClue) :-
    Word = [W, C, Lnk, Ls, Cs, D, WLen, S, _],
    WordClue = [W, C, Lnk, Ls, Cs, D, WLen, S, ClueNum].



% grid utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direction(down).
direction(across).


swap_dir(down, across).
swap_dir(across, down).


% The list of start locations.
start_locs([topleft_across, topleft_down, topright, bottomleft]).

% Get the cell number of each possible start location. 
start_loc(topleft_across, _GridLen, 1, across).
start_loc(topleft_down, _GridLen, 1, down).
start_loc(topright, GridLen, GridLen, down).
start_loc(bottomleft, GridLen, StartNum, across) :- 
    StartNum is (GridLen * GridLen) - (GridLen - 1).


% make sure word fits in the row...
fits_on_grid(across, Start, WLen, GridLen) :- 
    M is Start mod GridLen,
    M \== 0,
    Space is GridLen - (M - 1),
    WLen =< Space.

% Make sure word fits in the column...
fits_on_grid(down, Start, WLen, GridLen) :- 
     EndNum is Start + (GridLen * (WLen - 1)),
     EndNum =< GridLen * GridLen.


is_start_cell(across, Num, Length) :- 1 is Num mod Length.
is_start_cell(down, Num, Length) :- Num =< Length.
is_end_cell(across, Num, Length) :- 0 is Num mod Length.
is_end_cell(down, Num, Length) :- Num >= (Length - 1) * Length.


prev_cell(across, Num, _Length, Prev) :- Prev is Num - 1.
prev_cell(down, Num, Length, Prev) :- Prev is Num - Length.


next_cell(across, Num, _Length, Next) :- Next is Num + 1.
next_cell(down, Num, Length, Next) :- Next is Num + Length.


% Calculates the start of a word given the number
% of a letter that occurs at a certain position
calc_start(across, _GridLen, PPos, WNum, Start) :-
    Start is WNum - (PPos - 1).

calc_start(down, GridLen, PPos, WNum, Start) :-
    Start is WNum - (GridLen * (PPos - 1)).


% Calulates the number of a letter in a word given
% the position in the word and the number it starts at
calc_num(across, _GridLen, WPos, WStart, WNum) :-
    WNum is  WStart + (WPos - 1).

calc_num(down, GridLen, WPos, WStart, WNum) :-
    WNum is  WStart + (GridLen * (WPos - 1)).


new_tile(Num, Num-empty).

init_grid(GridLen, Grid) :-
    NumTiles is GridLen * GridLen,
    numlist(1, NumTiles, Tiles),
    maplist(new_tile, Tiles, TupleList),
    list_to_assoc(TupleList, Grid).



% printing predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


print_grid(Grid, Length) :-
    assoc_to_list(Grid, List),
    print_list(List, Length).


print_list([], _).
print_list([K-V|Xs], Length) :-
    (
     V == empty ->
     write('*')
    ;
     V = [A-D,L,P,LK],
     % join the contents of the cell with ',' as a seperator
     atomic_list_concat([A,D,L,P,LK], ',', Str),
     write(Str)
    ),
    write(' '),
    M is K mod Length,
    ( M == 0 -> nl ; true),
    print_list(Xs, Length).


% Partition into down and across lists and then 
% print each separately.
print_clues(Words) :-
    partition(dir_is_across, Words, AcrossWords, DownWords),
    write('Across Clues\n'),
    x_print_clues(AcrossWords),
    nl,
    write('Down Clues\n'),
    x_print_clues(DownWords),
    nl.


% todo: map could definintely be used here, except that
% map needs to return a list, we don't want to.
x_print_clues([]).
x_print_clues([W|Ws]) :-
    W = [Word, Clue, Link, _, _, _, _, _, ClueNum],
    atomic_list_concat([Word, ClueNum, Clue, Link], '|', Str),
    write(Str),
    nl,
    x_print_clues(Ws).



% generic utility predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

head([X|_], X).

length_sort(Words, OrderedWords) :-
    map_list_to_pairs(atom_length, Words, Pairs),
    keysort(Pairs, SortedPairs),
    pairs_values(SortedPairs, SortedWords),
    reverse(SortedWords, OrderedWords).


% finds all positions X occurs in List (over backtracking)
position(X, List, Pos) :- x_position(List, X, 1, Pos).

x_position([], _, _, _) :- false.
x_position([X|_], X, Pos, Pos).
x_position([_|Ys], X, N, Pos) :-
    N2 is N + 1,
    x_position(Ys, X, N2, Pos).


% remove_x(X,L,R) :- R is L with first occurrence of X removed from it.
remove_x(Y,[X|Xs],[X|Tail]) :-
	Y \== X,
	remove_x(Y,Xs,Tail).
remove_x(X,[X|Xs],Xs) :- !.
remove_x(_,[],[]).


%% shuffle(ListIn, ListOut) - randomly shuffles
%% ListIn and unifies it with ListOut
shuffle([], []) :- !.
shuffle(List, [Element|Rest]) :-
    choose(List, Element),
    delete(List, Element, NewList),
    shuffle(NewList, Rest).


%% choose(List, Elt) - chooses a random element
%% in List and unifies it with Elt.
choose(List, Elt) :-
    length(List, Length),
    Index is random(Length),
    nth0(Index, List, Elt).
