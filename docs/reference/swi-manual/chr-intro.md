
## 9.1 Introduction to CHR

Constraint Handling Rules (CHR) is a committed-choice rule-based language embedded in Prolog. It is designed for writing constraint solvers and is particularly useful for providing application-specific constraints. It has been used in many kinds of applications, like scheduling, model checking, abduction, and type checking, among many others.

CHR has previously been implemented in other Prolog systems (SICStus, Eclipse, Yap), Haskell and Java. This CHR system is based on the compilation scheme and runtime environment of CHR in SICStus.

In this documentation we restrict ourselves to giving a short overview of CHR in general and mainly focus on elements specific to this implementation. For a more thorough review of CHR we refer the reader to [Frühwirth, 2009](Bibliography.html#Freuhwirth:2009). More background on CHR can be found at [Frühwirth,](Bibliography.html#chrSite).

In [section 9.2](chr-syntaxandsemantics.html#sec:9.2) we present the syntax of CHR in Prolog and explain informally its operational semantics. Next, [section 9.3](practical.html#sec:9.3) deals with practical issues of writing and compiling Prolog programs containing CHR. [Section 9.4](chr-debugging.html#sec:9.4) explains the (currently primitive) CHR debugging facilities. [Section 9.4.3](chr-debugging.html#sec:9.4.3) provides a few useful predicates to inspect the constraint store, and [section 9.5](chr-examples.html#sec:9.5) illustrates CHR with two example programs. [Section 9.6](chr-compatibility.html#sec:9.6) describes some compatibility issues with older versions of this system and SICStus’CHR system. Finally, [section 9.7](chr-guidelines.html#sec:9.7) concludes with a few practical guidelines for using CHR.
