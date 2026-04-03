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

## Test case match status

The table below shows whether current output matches the expected output files exactly (byte-for-byte `diff`).

| Case | Target vertices | Status | Notes |
|---|---:|---|---|
| rectangle_with_two_holes | 7 | PASS |  |
| cushion_with_hexagonal_hole | 13 | PASS |  |
| blob_with_two_holes | 17 | FAIL |  |
| wavy_with_three_holes | 21 | FAIL |  |
| lake_with_two_islands | 17 | FAIL |  |
| original_01 | 99 | FAIL |  |
| original_02 | 99 | FAIL |  |
| original_03 | 99 | FAIL |  |
| original_04 | 99 | FAIL |  |
| original_05 | 99 | FAIL |  |
| original_06 | 99 | FAIL |  |
| original_07 | 99 | FAIL |  |
| original_08 | 99 | FAIL |  |
| original_09 | 99 | FAIL | topology error: ring 0: ring self-intersection detected |
| original_10 | 99 | FAIL |  |
