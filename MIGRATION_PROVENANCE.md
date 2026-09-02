# Migration provenance

## Ownership boundary

`orthofinder-results` is a generic, standalone package for interrogating
completed OrthoFinder result sets. It does not contain E3-ligase ranking,
structural, expression, ligandability or chemistry assumptions. The E3 project
may later consume versioned outputs or deep-link to composite group identities,
but it is not a runtime dependency of this package.

## Source authority

- Original repository: `peterthorpe5/E3_project_draft`
- Original branch: `feature/orthofinder-results-package`
- Original commit: `0c3b2fedf4ddcc0e8b2697991115ad7fee085a60`
- Original commit subject: `Release orthofinder-results 0.1.5 report views`
- Original complete tree: `1c0f7d2e97b5182cdab9da02d3d18b020f989b3b`
- Original `orthofinder_results/` subtree:
  `ad7ab32e9a5fe640a98bd84ca91dbefd243caff4`
- Package version: `0.1.5`
- Resource schema version: `2`
- Extraction date: `2026-09-02`

The subdirectory was extracted into the repository root with Git's
history-filtering machinery. Authors, author dates, commit dates and commit
messages were retained. Commit identifiers necessarily changed because each
commit now contains the package subtree at its repository root.

| Original E3 commit | Extracted commit | Subject |
|---|---|---|
| `8d4a030f987f2b0bf9e6af7b7092f74f2032ef30` | `6ecdbf199adcd19d4864391e512498ee5cc0a9a7` | Add generic OrthoFinder results interrogation package |
| `0ecb85eff02ca9ef0474872b946e5e2ce20a665d` | `615ce0cc76e3f35ced676bba4be7754212070f48` | Release orthofinder-results 0.1.2 |
| `809ec67c21d5776a8fb4265a1b069de766b53864` | `c18423d0edae3a830459e72e0928125cddd0e53f` | Fix and streamline OrthoFinder HTML reports |
| `9991b04b5eb096c0f06cffa0f3a2d2291751fd4b` | `3fa26fa14c450b6042dde16c66403794ceb0684f` | Clarify OrthoFinder protein networks |
| `0c3b2fedf4ddcc0e8b2697991115ad7fee085a60` | `c7e9a9708b2045297d76c0c87f36af427ca36776` | Release orthofinder-results 0.1.5 report views |

The extracted v0.1.5 commit has tree
`ad7ab32e9a5fe640a98bd84ca91dbefd243caff4`, exactly matching the source
subtree. All 30 tracked source files were additionally compared byte for byte
before repository-specific metadata was added.

## Extraction validation

The unchanged extracted v0.1.5 tree passed the complete isolated quality gate:

- 59 tests passed;
- 95.62% branch coverage;
- pycodestyle passed;
- pydocstyle passed;
- Ruff passed;
- Python compilation passed; and
- shell syntax checks passed.

Repository-boundary changes added after that equivalence check do not alter
scientific calculations, the resource schema or report rendering.

The original E3 commits remain in the E3 repository as historical provenance.
They must not be rewritten or deleted from that history.
