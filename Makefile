# Convenience targets for crosswordsmith.

.PHONY: test unit golden update-golden

# Full suite: plunit tests + golden-output regression.
test:
	./run_tests.sh

# Just the plunit unit/integration tests.
unit:
	swipl -q tests/run_tests.pl

# Just the golden-output regression check.
golden:
	@diff -u tests/golden/grid_17_topleft_across.txt \
		<(./crossword.pl 17 topleft_across 2>/dev/null) \
		&& echo "golden: OK"

# Regenerate the golden file. Use only after an INTENTIONAL output change,
# and review the diff before committing.
update-golden:
	./crossword.pl 17 topleft_across 2>/dev/null > tests/golden/grid_17_topleft_across.txt
	@echo "Regenerated tests/golden/grid_17_topleft_across.txt"
