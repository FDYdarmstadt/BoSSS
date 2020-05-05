# Overview on the test architecture in BoSSS

## Test definition

Tests are defined using NUnit 3. 
If the test requires additional files 
(e.g. test-databases, mesh files, BoSSSpad scripts for tutorials, etc.)
they should be marked with the `NUnitFileToCopyHackAttribute`.

![Screenshot of test marked with `NUnitFileToCopyHackAttribute`](TestAttributes.png "Test using additional data")

## Test execution

Tests are executes via the `PublicTestRunner.exe`, resp. the
`InternalTestsRunner.exe` if one has access to the internal parts of 
the repository, too.
The main intention of these runners is to execute a 
large number of tests distributed, on a HPC compute cluster.
This can be e.g. a Microsoft HPC Cluster, or a SLURM-based 
Linux cluster. Therefore, the test runner has two separate execution modes 
(or sub-programs, run `PublicTestRunner.exe help` to obtain the 
respective  syntax):
- local test execution (sub-program `nunit3`)
- deployment of tests at a compute cluster (sub-program `runjobmanager`)


In both cases one is able to continue to work (and compile) 
other parts of the source code, since the 
test runner uses a copy of the respective assembly (exe or dll).

## Examples

- Local execution all tests in `TutorialTests.exe`:
  ```
  PublicTestRunner.exe nunit3 TutorialTests.exe
  ```
- Local execute some specific test:
  ```
  PublicTestRunner.exe nunit3 TutorialTests.exe --test=BoSSS.Application.TutorialTests.AllUpTest.Run__CsharpAndBoSSSpad --result=blabla.xml
  ```
  (i.e. after the filter, there are only 
  standard arguments for nunit3)
- Execute tests on a Batch Processor Queue:
  ```
  PublicTestRunner.exe runjobmanager  TutorialTests.exe queue#1
  ```
  (The 'queue'-argument is optional: 
  after the hash a zero-based index is added to the entries of 
  `~/.BoSSS/etc/BatchProcessorConfig.json`, 
  i.e. to the list of predefined batch processors; default is index 1, i.e. the second entry.)

### Running tests locally using the MiniBatchProcessor
1. If not already done, launch local `MiniBatchProcessor.exe`, preferably not 
   the one compiled in the current repository.
   (A running exe causes a file lock and therefore the repository cannot be compiled as long as the exe is running.)
2. Submit test jobs to the Mini Batch Processor:
   ```
   PublicTestRunner.exe runjobmanager  '*' queue#0
   ```
   (We assume the Mini Batch Processor is the first entry (queue #0) in the 
   local `~/.BoSSS/etc/BatchProcessorConfig.json`).
3. Grab a cup of coffee and wait for results.


## Troubleshooting

Tests fail - thats what they are for.
The question is: if some test fails, what went wrong? 
Therefore, one needs to be able to inspect the respective run.

At first,one has to identify which test, resp. job went wrong.
The `PublicTestRunner.exe' may provide an output like this:
```
...
FinishedSuccessful: Apr23_2350-REL-XNSE_Solver#ScalingViscosityJumpTest_p3 // BoSSS.Application.XNSE_Solver.Tests.UnitTest.ScalingViscosityJumpTest_p3
Failed: Apr23_2350p4-REL-XdgPoisson3#ScalingCircle2D // BoSSS.Application.XdgPoisson3.Tests.ScalingCircle2D at \\terminal03\jenkins\BoSSStstApr23_2350-InternalTestRunner2020Apr24_001752
Only 225 tests/jobs finished successfully -- 1 have other states.
FAILURE.
04/24/2020 01:16:07
 ```
So, it appears that the failed **test** is `BoSSS.Application.XdgPoisson3.Tests.ScalingCircle2D`, the respective **job** is `Apr23_2350p4-REL-XdgPoisson3#ScalingCircle2D`.

Note that ll tests/jobs sent to the batch queue have the same prefix - 
in this case `Apr23_2350`. This allows to find the jobs in the 
Job Manger at first glance. Btw., 2350 stands for 11:50 p.m., which is the time 
at which `PublicTestRunner.exe` started submitting jobs.

For "Failed" tests/jobs, `PublicTestRunner.exe` also specifies the working directory.
One can copy&paste this, e.g. into the file explorer and have a look at the output there.
For "Success" jobs, the directory is deleted, directories of "Failed" jobs remain untouched.

One can now examine the job directory, or also the respective job manager for errors.
Alternatively, you can also view the output in the "allout-*.txt" log file 
produced by `PublicTestRunner.exe`, e.g. `allout-Apr23_2350-REL.txt`, 
and search there (with Ctrl-F of course) for the corresponding job/test name.



