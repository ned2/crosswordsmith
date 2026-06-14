% these words can be solved with a grid of length 17
%
% Each entry is [Answer, Metadata]. Answer is the word (spaces allowed for
% multi-word answers; they are stripped before placement). Metadata is an
% optional dict of passthrough data, by convention `clue` and `link`; the
% solver never inspects it and copies it verbatim into the output. The
% metadata element may also be omitted entirely (just [Answer]).

clues([
       ['OMEGA POINT',
        _{clue: 'Transcending entropy',
          link: 'http://en.wikipedia.org/wiki/Omega_Point'}],
       ['FLOW',
        _{clue: 'Autotelic activity',
          link: 'http://en.wikipedia.org/wiki/Flow_(psychology)'}],
       ['GNOSTIC GOSPELS',
        _{clue: 'Some apocrypha',
          link: 'http://en.wikipedia.org/wiki/Gnostic_Gospels'}],
       ['BIAS',
        _{clue: 'Why your brain sucks',
          link: 'http://en.wikipedia.org/wiki/List_of_cognitive_biases'}],
       ['ETERNAL RETURN',
        _{clue: 'Live',
          link: 'http://en.wikipedia.org/wiki/Eternal_return'}],
       ['NARRATIVE FALLACY',
        _{clue: 'Telos everywhere',
          link: 'http://wiki.lesswrong.com/wiki/Narrative_fallacy'}]
      ]).
