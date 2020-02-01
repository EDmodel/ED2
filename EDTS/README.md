# ED2 test suite

**Author**: Alexey Shiklomanov

## Instructions

_From inside this directory_, run `./run-test.sh TESTNAME [ED2_EXE]`, where:

- `TESTNAME` is the name of a test, corresponding to a file `Templates/ED2IN-TESTNAME`.
`TESTNAME` is of the form `SITENAME[.options]`, where `SITENAME` matches the tags of the driver and initial condition data available as releases at https://github.com/ashiklom/edts-datasets/releases.
The `run-test.sh` script assumes that `SITENAME` is separated from `options` by a period (`.`) -- therefore, `SITENAME` _cannot contain any periods_.
These data are automatically downloaded and unzipped when a test is run (unless they are already downloaded).
The `options` are there to allow different tests to be run with the same site data; 
for instance, you might have two tests from UMBS---`umbs.bg` for a "bare-ground" run, `umbs.ic` for a run with real initial conditions, and `umbs.bg.fcr` for a run from bare ground and also enabling the finite canopy radius submodel.
(Note that unlike `SITENAME`, `options` _can_ contain periods).
- `ED2_EXE` is the optional path to an ED2 executable.
Unless given as an absolute path, it is assumed to be relative to the `EDTS` directory.
By default, it is set to `../ED/build/ed_2.1-dbg`.

An example of a full invocation:

``` sh
./run-test.sh umbs.bg ../ED/build/ed_2.1-opt
```

This will use the `Templates/ED2IN-umbs.bg` input file, and download inputs for UMBS from https://github.com/ashiklom/edts-datasets/releases/download/umbs/umbs.tar.gz.

## Adding additional tests

Assuming your test is called `my-site.opts`:

1. Gather all your site-specific inputs and drivers in a directory called `my-site`.
Make sure any paths referred-to inside this directory are _relative_ to the parent directory of `my-site` (i.e. all paths should start with `my-site/`; usually, this is just the `ED_MET_DRIVER_HEADER` file).
Compress this directory with `tar -czvf my-site.tar.gz my-site` and create a new release in the https://github.com/ashiklom/edts-datasets repository with name _and_ tag both set to `my-site`.

2. Create a `Templates/ED2IN-my-site.opts` file.
To work correctly across systems (including the automated GitHub continuous integration tests), all paths should be _relative_ paths relative to the `EDTS` directory.
Common ED2 inputs should be stored in the directory `common`, and site-specific inputs (downloaded from the GitHub release) should be stored in `my-site`.
Output files should be written to the directory `test-outputs/my-site/`.

3. (Optional) Add `my-site/` to the `EDTS/.gitignore` file to avoid accidentally committing files.

4. Run the test locally to make sure it works: From the `EDTS` directory, run `./run-test.sh my-site.opts /path/to/ed_2.1`

5. (Optional) Add the test to the GitHub Actions continuous integration suite.
In a plain text editor, open the `.github/workflows/ci.yml` file, copy one of the existing test blocks (e.g. `test-umbs-bg`), and change the block name (e.g. `test-umbs-bg: -> test-my-site-opts:`; note that the block name _cannot_ contain periods) and the name of the test to be run at the very end (e.g. `./run-test umbs.bg ... -> ./run-test my-site.opts ...`).
You should be able to keep everything else the same.
Assuming you've done this correctly, the test should run the next time you push your version of ED2 up to GitHub.
