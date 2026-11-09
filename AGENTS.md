# Repository Instructions

## Code Style

All new code and edits, including wholly new files, must follow the established
WoSX style. Use neighboring original code as the primary reference, especially
`include/wosx/point_estimation/walk_on_spheres.h`,
`include/wosx/point_estimation/walk_on_stars.h`,
`include/wosx/core/geometric_queries.h`, and
`demo_apps/basic_2d/model_problem.h`. Match the language and area being edited;
do not apply C++ formatting rules to CMake, Python, or Slang indiscriminately.

For C++:

- Use four spaces for indentation and no tabs.
- Use `PascalCase` for types, `camelCase` for functions and variables, and
  `snake_case` for filenames, directories, and namespaces where applicable.
- Use descriptive names for program state and settings. Short mathematical names
  are appropriate in formulas when their meaning is clear from the local context.
- Put function-definition opening braces on the next line. Put class, struct,
  namespace, control-flow, and lambda opening braces on the same line.
- Expand control flow and separate statements onto individual lines. Avoid dense
  one-line loops, conditionals, and sequences of statements. Small accessors may
  remain inline when consistent with neighboring code.
- Align wrapped parameter and argument lists with the first parameter or
  argument. Follow local wrapping conventions rather than imposing an unrelated
  formatter preset.
- Attach pointer and reference symbols to the type, such as
  `const PDE<T, DIM>& pde`.
- Write template declarations as `template <typename T, size_t DIM>` and place
  them on their own line. Use the existing uppercase dimension/channel parameter
  convention.
- Use a trailing underscore for constructor parameters that directly initialize
  identically named members, as in `queries_`. Follow the surrounding class's
  member convention: core structs commonly use plain names; demo classes often
  use an `m` prefix.
- Use `#pragma once` in headers. Keep template implementations in headers, with
  declarations first and definitions below an `Implementation` section where
  consistent with the original solver headers. Use `inline` when moving
  non-template definitions into headers requires it.
- Add concise file-level descriptions and comments explaining numerical
  assumptions, algorithm stages, and non-obvious choices. Preserve existing
  license notices; do not invent copyright ownership.
- Follow surrounding include ordering and spacing. Avoid introducing a different
  named style or treating formatter defaults as the project's specification.

For CMake, follow the existing indentation, lowercase commands, `WOSX_` option
names, and section/comment conventions. Preserve the established style of other
languages by consulting their original files.

## Scope and Verification

- Apply these conventions to the entirety of new files, not just modified lines.
- Keep formatting changes to existing files focused on the requested work; avoid
  unrelated reformatting of original code or vendored dependencies.
- Keep style-only changes behavior-preserving: preserve numerical expressions,
  evaluation order, RNG use, geometry decisions, CLI/CSV interfaces, and failure
  handling. A style cleanup is not authorization for an algorithm change.
- Review automated formatting and identifier renames manually, including include
  paths, existing API members, defaults, and string literals.
- Use existing builds and tests appropriate to the changes. For numerical
  reorganizations or renames, compare deterministic outputs when practical;
  exclude timing fields from equality comparisons.
