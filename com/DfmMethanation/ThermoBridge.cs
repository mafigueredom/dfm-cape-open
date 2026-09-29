using System;
using System.Globalization;
using System.Reflection;
using System.Runtime.InteropServices;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// ICapeThermoMaterial 1.1 plus ICapeThermoCompounds for the compound list.
    /// Those are different COM interfaces on the same TEA stream. SI: K, Pa, mol/s.
    /// </summary>
    internal static class ThermoBridge
    {
        static readonly Guid IidMaterial = new Guid("678C0A9B-7D66-11D2-A67D-00105A42887F");
        static readonly Guid IidCompounds = new Guid("678C0A9D-7D66-11D2-A67D-00105A42887F");
        static readonly Guid IidMaterial10 = new Guid("678C0994-7D66-11D2-A67D-00105A42887F");

        public static double[] ToDoubles(object raw)
        {
            if (raw == null)
                return Array.Empty<double>();
            if (raw is double d)
                return new[] { d };
            if (raw is float f)
                return new[] { (double)f };
            if (raw is int i)
                return new[] { (double)i };
            if (raw is Array arr)
            {
                var o = new double[arr.Length];
                for (int k = 0; k < arr.Length; k++)
                    o[k] = Convert.ToDouble(arr.GetValue(k), CultureInfo.InvariantCulture);
                return o;
            }
            return new[] { Convert.ToDouble(raw, CultureInfo.InvariantCulture) };
        }

        public static string[] ToStrings(object raw)
        {
            if (raw == null)
                return Array.Empty<string>();
            if (raw is string s)
                return new[] { s };
            if (raw is Array arr)
            {
                var o = new string[arr.Length];
                for (int k = 0; k < arr.Length; k++)
                    o[k] = Convert.ToString(arr.GetValue(k), CultureInfo.InvariantCulture);
                return o;
            }
            return new[] { Convert.ToString(raw, CultureInfo.InvariantCulture) };
        }

        public static double GetOverallScalar(object mo, string property, string basis = null)
        {
            var values = GetOverall(mo, property, basis);
            if (values.Length == 0)
                throw new InvalidOperationException("empty " + property);
            return values[0];
        }

        public static double[] GetOverallVector(object mo, string property, string basis)
        {
            return GetOverall(mo, property, basis);
        }

        static double[] GetOverall(object mo, string property, string basis)
        {
            var mat = (ICapeThermoMaterial11)AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
            var b = string.IsNullOrEmpty(basis) ? "undefined" : basis;
            var names = property == "totalFlow"
                ? new[] { "flow", "totalFlow" }
                : new[] { property };
            Exception last = null;
            foreach (var name in names)
            {
                try
                {
                    var values = ToDoubles(mat.GetOverallProp(name, b));
                    CapeTrace.Write("GetOverallProp " + name + " " + b + " n=" + values.Length);
                    return values;
                }
                catch (Exception ex)
                {
                    last = ex;
                    CapeTrace.Write("GetOverallProp " + name + " " + b + " " + ex.Message);
                }
            }
            throw new InvalidOperationException(
                "GetOverallProp " + property + ": " + (last == null ? "failed" : last.Message),
                last);
        }

        public static string[] GetCasNumbers(object mo)
        {
            var compounds = (ICapeThermoCompounds)AsInterface(
                mo, IidCompounds, typeof(ICapeThermoCompounds));
            object compIds, formulae, names, boilTemps, molwts, casnos;
            compounds.GetCompoundList(
                out compIds, out formulae, out names, out boilTemps, out molwts, out casnos);
            var cas = ToStrings(casnos);
            var ids = ToStrings(compIds);
            CapeTrace.Write(
                "GetCompoundList cas=" + string.Join(",", cas) + " ids=" + string.Join(",", ids));
            bool anyCas = false;
            foreach (var c in cas)
            {
                if (!string.IsNullOrWhiteSpace(c))
                    anyCas = true;
            }
            if (!anyCas)
                throw new InvalidOperationException(
                    "TEA returned no CAS numbers. Compound ids: " + string.Join(", ", ids));
            return cas;
        }

        public static void SetOverall(object mo, string property, object value, string basis = null)
        {
            var mat = (ICapeThermoMaterial11)AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
            var b = string.IsNullOrEmpty(basis) ? "undefined" : basis;
            var names = property == "totalFlow"
                ? new[] { "flow", "totalFlow" }
                : new[] { property };
            Exception last = null;
            foreach (var name in names)
            {
                try
                {
                    mat.SetOverallProp(name, b, value);
                    CapeTrace.Write("SetOverallProp " + name + " " + b + " ok");
                    return;
                }
                catch (Exception ex)
                {
                    last = ex;
                    CapeTrace.Write("SetOverallProp " + name + " " + b + " " + ex.Message);
                }
            }
            throw new InvalidOperationException(
                "SetOverallProp " + property + ": " + (last == null ? "failed" : last.Message),
                last);
        }

        public static void FlashTP(object mo)
        {
            object legacy;
            try
            {
                legacy = AsInterface(mo, IidMaterial10, null);
            }
            catch (Exception ex)
            {
                CapeTrace.Write("FlashTP no thermo 1.0 " + ex.Message);
                return;
            }
            try
            {
                legacy.GetType().InvokeMember(
                    "CalcEquilibrium",
                    BindingFlags.InvokeMethod,
                    null,
                    legacy,
                    new object[] { "TP" },
                    CultureInfo.InvariantCulture);
                CapeTrace.Write("FlashTP ok");
            }
            catch (Exception ex)
            {
                CapeTrace.Write("FlashTP " + ex.Message);
            }
        }

        static object AsInterface(object mo, Guid iid, Type typedAs)
        {
            IntPtr unk = IntPtr.Zero;
            IntPtr p = IntPtr.Zero;
            try
            {
                unk = Marshal.GetIUnknownForObject(mo);
                Guid g = iid;
                int hr = Marshal.QueryInterface(unk, ref g, out p);
                if (hr < 0 || p == IntPtr.Zero)
                    throw new InvalidOperationException(
                        "QueryInterface " + iid.ToString("B") + " failed 0x" + hr.ToString("X8"));
                return typedAs == null
                    ? Marshal.GetObjectForIUnknown(p)
                    : Marshal.GetTypedObjectForIUnknown(p, typedAs);
            }
            finally
            {
                if (p != IntPtr.Zero)
                    Marshal.Release(p);
                if (unk != IntPtr.Zero)
                    Marshal.Release(unk);
            }
        }
    }

    [ComImport]
    [Guid("678C0A9D-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsIDispatch)]
    internal interface ICapeThermoCompounds
    {
        void GetCompoundList(
            [MarshalAs(UnmanagedType.Struct)] out object compIds,
            [MarshalAs(UnmanagedType.Struct)] out object formulae,
            [MarshalAs(UnmanagedType.Struct)] out object names,
            [MarshalAs(UnmanagedType.Struct)] out object boilTemps,
            [MarshalAs(UnmanagedType.Struct)] out object molwts,
            [MarshalAs(UnmanagedType.Struct)] out object casnos);
    }

    [ComImport]
    [Guid("678C0A9B-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsIDispatch)]
    internal interface ICapeThermoMaterial11
    {
        [return: MarshalAs(UnmanagedType.Struct)]
        object GetOverallProp(
            [MarshalAs(UnmanagedType.BStr)] string property,
            [MarshalAs(UnmanagedType.BStr)] string basis);

        void SetOverallProp(
            [MarshalAs(UnmanagedType.BStr)] string property,
            [MarshalAs(UnmanagedType.BStr)] string basis,
            [MarshalAs(UnmanagedType.Struct)] object values);
    }
}
