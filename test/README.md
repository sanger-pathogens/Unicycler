# Unicycler tests

This page briefly describes the tests written to check Unicycler. Users of Unicycler are probabaly uninterested in this stuff, but if you are a Unicycler developer, then read on!

All tests are run from the root Unicycler repository (i.e the directory which contains `unicycler-runner.py`).


### Unit tests

These tests check various functions in Unicycler's code and should complete fairly quickly.

To run unit tests:
`python3 -m unittest`


### Random sequence assembly test:

This test:
* generates simple random sequences
* makes fake reads from those sequences
* assembles those reads with Unicycler
* verifies that Unicycler's output looks correct

It runs indefinitely (i.e. will need to be stopped with ctrl-c) and may leave behind a directory called `TEST_TEMP_*` that you'll probably want to delete afterward.

To run random sequence assembly tests:
`python3 test/random_sequence_assembly_test.py`


### Overlap removal test:

This test:
* generates sequences with many repeats
* makes fake reads from those sequences
* assembles those reads with SPAdes
* runs the `AssemblyGraph.remove_all_overlaps` function on the SPAdes graph

It is designed to test two things:
* whether overlap removal is successful
* overlap removal performance

It runs indefinitely (i.e. will need to be stopped with ctrl-c) and may leave behind a directory called `TEST_TEMP_*` that you'll probably want to delete afterward.

To run overlap-removal tests:
`python3 test/overlap_removal_test.py`
