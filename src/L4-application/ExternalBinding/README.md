# Native hosting for foam-dg

ExternalBinding targets .NET 10. On Linux x64, foam-dg embeds the JIT runtime
through nethost/hostfxr and calls generated `[UnmanagedCallersOnly]` entry points.
The runtime is framework-dependent: install .NET 10, and use the SDK/native
hosting pack when building the C++ host. Mono is no longer used by this bridge.

## Generate and build

`ExternalBinding.CodeGen` discovers `IForeignLanguageProxy` classes in BoSSS's
referenced assemblies. Constructors/methods tagged `[CodeGenExport]`, plus the
foreign-pointer interface methods, define the supported surface. Run:

```bash
dotnet build public/src/L4-application/ExternalBinding.CodeGen -c Release -p:ExcludeNativeExports=true
dotnet public/src/L4-application/ExternalBinding.CodeGen/bin/Release/net10.0/ExternalBinding.CodeGen.dll \
    /absolute/staging-directory public/src/L4-application/ExternalBinding
dotnet publish public/src/L4-application/ExternalBinding -c Release \
    --self-contained false -o /absolute/managed-deployment-directory
```

Review and install the generated C++ headers/sources in foam-dg's
`src/discontinuousGalerkin/BoSSSwrapper/`. Its `createBindings.sh` automates
these steps when given the BoSSS source checkout. It rebuilds both sides after
generation. Commit `NativeExports.g.cs` with export-surface changes so ordinary
solution builds include the adapters. `ExcludeNativeExports=true` omits stale adapters during
the generator bootstrap, allowing removed or renamed exports to regenerate.
Normal builds/publishing include the generated file. The generator overwrites only its named
outputs; it does not delete other files in the destination.

Deploy the entire publish output, including `.runtimeconfig.json`, `.deps.json`
and transitive managed dependencies. Keep `BOSSS_INSTALL` pointed at the separate
native installation (`bin/native/linux/amd64-openmpi`). Hosting the CLR does not
replace BoSSS's native libraries or MPI bootstrap. foam-extend and BoSSS must use
compatible system OpenMPI installations.

## Boundary and ownership

The C++ public class wrappers retain ordinary instance methods. Each call uses a
cached native function pointer to a generated static entry point. Scalars and
pointers cross the boundary directly; managed object references become
pointer-sized opaque tokens backed by normal, unpinned `GCHandle`s. Tokens are
registered and type-checked on the managed side. The pointed-to native buffers
must remain valid for the duration of the synchronous call.

A managed token is distinct from `IForeignLanguageProxy._GetForeignPointer()`,
which stores the C++ wrapper address. A returned object creates an owned root.
If its native wrapper already exists, `_FromHandle` releases only the duplicate
root and reuses that wrapper. Otherwise it adopts the root. Null object returns
produce null native pointers. Native wrappers cannot be copied. Their
destructors clear the foreign pointer, dispose `IDisposable` targets where
applicable, and release the root without throwing across native destruction.
Do not destroy a proxy concurrently with a call using it.

Managed entry points catch exceptions and return a status code. On failure the
native wrapper retrieves per-thread UTF-8 diagnostics into caller-owned storage
and throws `std::runtime_error`; no result is read after a failed call. Never
allow an exception to unwind directly through a UCO/native boundary.

`Globals::Init` discovers hostfxr, starts .NET using ExternalBinding's runtime
configuration and resolves its adapters. Initialization checks a generated ABI
fingerprint before exposing the call surface. Regenerate/rebuild both sides if
it reports a mismatch. The runtime remains alive for the native process
lifetime; closing the hostfxr configuration context does not unload the CLR.
OpenFOAM enables hardware floating-point traps that can also affect .NET's
background JIT. Each transition (including CLR startup) saves the native FP
environment, masks traps for CLR execution and restores the saved environment
on return, including exception paths. This preserves OpenFOAM's native error
traps. The native smoke test enables those traps before startup and checks
restoration. See also [the runtime hosting issue](https://github.com/dotnet/runtime/issues/76257).

`Initializer` retains the existing conditional ownership of MPI initialization
and finalization.

## Checks

Use PublicTestRunner's exact assembly names:

- `ExternalBinding.CodeGen`, test `BoSSS.Application.ExternalBinding.CodeGen.Test.RunGenerator`
- `BoSSS.Application.ExternalBinding`, test `BoSSS.Application.ExternalBinding.GridImportTestSmall.ConvertFOAMGrid`

Run from the respective runner build directory in Debug and Release, as required
by AGENTS.md. For native integration, run the generated `compile.sh`, then
`ExtBindingTest.out /absolute/managed-deployment-directory` with `BOSSS_INSTALL`
set. It checks grid/field exchange, returned wrapper identity and exception
translation. foam-dg's convergence and droplet cases remain the end-to-end
numerical checks. This host implementation currently targets Linux x64.
