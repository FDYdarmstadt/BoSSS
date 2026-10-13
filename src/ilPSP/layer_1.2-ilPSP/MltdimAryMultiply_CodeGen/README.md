# Tensor multiply code generator

Requires the .NET 10 SDK. This standalone tool is not included in `internal/src/Master.sln`.
It generates `../ilPSP/MultidimensionalArray_CodeGenMultiply.cs`; it does not run tensor products.

From the repository root:

```bash
dotnet build public/src/ilPSP/layer_1.2-ilPSP/MltdimAryMultiply_CodeGen/MltdimAryMultiply_CodeGen.sln -c Release
dotnet run --project public/src/ilPSP/layer_1.2-ilPSP/MltdimAryMultiply_CodeGen -c Release -- /tmp/MultidimensionalArray_CodeGenMultiply.cs
```

An optional positional argument selects the output file. Relative output paths are
relative to the current working directory. Its parent directory must exist.
Without an argument, the tool overwrites the sibling ilPSP generated file, locating
it relative to the standard `bin/<Configuration>/net10.0/` executable directory.
Thus `dotnet run` can be invoked from the repository root or the project directory.
For relocated/published binaries or custom output directories, pass an explicit path.

Output is UTF-8 without BOM and uses CRLF on all platforms. The generator does not
format C# indentation; the existing checked-in file has been formatted separately.
The .NET 10 migration was checked against that file after whitespace normalization:
no generated expressions, methods or dispatch decisions changed.

To regenerate in place after editing the generator:

```bash
dotnet run --project public/src/ilPSP/layer_1.2-ilPSP/MltdimAryMultiply_CodeGen -c Release
```

Review the generated diff, accounting for formatting, then build the master solution
and run the relevant multiply tests in Debug and Release as described in AGENTS.md.
