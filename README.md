# Area and Topology Preserving Polygon Simplification

This project implements a C++17 command line tool that simplifies polygon geometry while preserving topology and reporting area-based metrics as per the assignment.

## Executable

- Program name: `simplify`
- Source entry point: `src/main.cpp`
- Build system: `Makefile`

## Build

From the `data-struct-a2` directory:

```bash
make
```

Clean build artifacts:

```bash
make clean
```

## Usage

```bash
./simplify <input_file> <target_vertices>
```

Example:

```bash
./simplify test_cases/input_wavy_with_three_holes.csv 21
```

## Input CSV format

The input file must contain this header:

```text
ring_id,vertex_id,x,y
```

Rules expected by the parser and validator:

- `ring_id` and `vertex_id` are non-negative integers
- vertices for each ring are ordered by `vertex_id`
- ring 0 is typically the outer ring and additional rings may be holes or islands

## Output format

The program writes to stdout:

1. Simplified polygon CSV rows in the same `ring_id,vertex_id,x,y` schema
2. Assignment metric lines:

```text
Total signed area in input: <scientific>
Total signed area in output: <scientific>
Total areal displacement: <scientific>
```

## Topology behavior

The tool validates topology and rejects invalid edits during simplification.

- Global input topology validation is enabled by default, including very large single-ring inputs
- If needed for runtime experiments, large single-ring startup validation can be skipped via `ATPPS_SKIP_LARGE_SINGLE_RING_INPUT_VALIDATION=1`
- Simplification still applies legality checks during collapse steps

## Test data

Provided datasets and expected outputs are in `test_cases/`.

- Input files: `test_cases/input_*.csv`
- Expected files: `test_cases/output_*.txt`
- Dataset overview: `test_cases/README.md`

## Notes

- The implementation is deterministic for a fixed binary and environment
- The repository currently includes assignment planning notes in `IMPLEMENTATION_PLAN.md`
