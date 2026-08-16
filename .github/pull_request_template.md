<!--
Fill in what applies and delete what does not. A short pull request needs a
short description; nobody is asking for an essay to fix a typo.
-->

## What this changes

<!-- One or two sentences. What is different after this is merged? -->

## Why

<!-- The problem it solves. Link an issue with "Closes #123" if there is one. -->

## How it was verified

<!--
The important part, and the one a reviewer cannot reconstruct.

Not "it should work", but what you ran and what it printed. For example: the
test you added and the fact that it fails without the fix, the alignment or the
VCF you ran through it and the dN/dS or pN/pS you got before and after, or the
comparison run whose numbers changed.

Numbers deserve special care here. Most of this code produces a float, and a
float is a plausible-looking answer whether or not it is the right one. If you
touched site counting, a substitution model, a p-value or a correction, say
what you checked the new number against: a hand computation, the published
example for the model, a run of the previous version, or the golden snapshot.

If it is a change that CI already covers, say which check covers it.
-->

## Checklist

- [ ] `cargo test` passes (the whole suite, not only the test you were working on).
- [ ] `cargo clippy --all-targets -- -D warnings` is clean.
- [ ] New behaviour has a test, or I have said above why it does not.
- [ ] If output changed, the golden snapshot was re-blessed on purpose with
      `BLESS=1 cargo test --test golden`, and I read the diff before committing it.
      A snapshot re-blessed without reading it records the bug instead of catching it.
- [ ] If I added or renamed a CLI flag or an output column, `docs/` says so. The
      `docs_contract` test compares the help text and the column headers against the
      documentation, so this is a gate rather than good manners.
- [ ] No em-dashes anywhere in the diff, in code, comments, docs or commit message.
      The pre-commit hook rejects them; a comma, a colon or a full stop does the job.
