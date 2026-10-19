using System;
using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using ilPSP.Connectors;

namespace BoSSS.Application.ExternalBinding {
    /// <summary>
    /// Owns roots crossing the C ABI. Tokens and foreign C++ wrapper addresses have distinct roles.
    /// Each returned root is owned by the caller; disposing a duplicate root must not dispose its target.
    /// </summary>
    internal static class NativeExportRuntime {
        static readonly ConcurrentDictionary<IntPtr, GCHandle> Roots = new ConcurrentDictionary<IntPtr, GCHandle>();
        [ThreadStatic] internal static string LastError;

        internal static IntPtr Root(object target) {
            if (target == null) return IntPtr.Zero;
            var root = GCHandle.Alloc(target, GCHandleType.Normal);
            var token = GCHandle.ToIntPtr(root);
            if (!Roots.TryAdd(token, root)) {
                root.Free();
                throw new InvalidOperationException("Duplicate managed root token.");
            }
            return token;
        }

        internal static T Target<T>(IntPtr token, bool allowNull = false) where T : class {
            if (token == IntPtr.Zero && allowNull) return null;
            if (!Roots.TryGetValue(token, out var root))
                throw new ArgumentException("Unknown or released managed root token.", nameof(token));
            return root.Target as T ?? throw new InvalidCastException("Managed root has an unexpected type.");
        }

        internal static void Release(IntPtr token, bool dispose) {
            if (!Roots.TryRemove(token, out var root))
                throw new ArgumentException("Unknown or released managed root token.", nameof(token));
            try {
                if (dispose) {
                    ((IForeignLanguageProxy)root.Target)._SetForeignPointer(IntPtr.Zero);
                    if (root.Target is IDisposable disposable) disposable.Dispose();
                }
            } finally { root.Free(); }
        }

        internal static int Fail(Exception error) {
            LastError = error.ToString();
            return -1;
        }
    }

    /// <summary>Blittable entry points used by the generated C++ hostfxr wrappers.</summary>
    public static unsafe partial class NativeExports {
        [UnmanagedCallersOnly(CallConvs = new[] { typeof(CallConvCdecl) })]
        public static int Bind(nint handle, nint foreign) {
            try {
                NativeExportRuntime.Target<IForeignLanguageProxy>(handle)._SetForeignPointer(foreign);
                return 0;
            } catch (Exception error) { return NativeExportRuntime.Fail(error); }
        }

        [UnmanagedCallersOnly(CallConvs = new[] { typeof(CallConvCdecl) })]
        public static int Foreign(nint handle, nint* result) {
            try {
                if (result == null) throw new ArgumentNullException(nameof(result));
                *result = NativeExportRuntime.Target<IForeignLanguageProxy>(handle)._GetForeignPointer();
                return 0;
            } catch (Exception error) { return NativeExportRuntime.Fail(error); }
        }

        [UnmanagedCallersOnly(CallConvs = new[] { typeof(CallConvCdecl) })]
        public static int Release(nint handle, int dispose) {
            try {
                NativeExportRuntime.Release(handle, dispose != 0);
                return 0;
            } catch (Exception error) { return NativeExportRuntime.Fail(error); }
        }

        // Copies a per-thread UTF-8 message into caller-owned storage. Returns required bytes, including NUL.
        [UnmanagedCallersOnly(CallConvs = new[] { typeof(CallConvCdecl) })]
        public static int CopyLastError(byte* buffer, int capacity) {
            try {
                var bytes = Encoding.UTF8.GetBytes(NativeExportRuntime.LastError ?? "No managed error details.");
                int required = checked(bytes.Length + 1);
                if (buffer != null && capacity >= required) {
                    Marshal.Copy(bytes, 0, (IntPtr)buffer, bytes.Length);
                    buffer[bytes.Length] = 0;
                }
                return required;
            } catch { return 0; }
        }
    }
}
