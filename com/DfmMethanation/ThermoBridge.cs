using System;
using System.Globalization;
using System.Reflection;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// Late-bound ICapeThermoMaterial (1.1) with MaterialObject (1.0) fallback.
    /// SI: K, Pa, mol/s, mole fraction. Does not read packed-bed D_z.
    /// </summary>
    internal static class ThermoBridge
    {
        static object Invoke(object mo, string method, params object[] args)
        {
            return mo.GetType().InvokeMember(
                method,
                BindingFlags.InvokeMethod | BindingFlags.GetProperty,
                null,
                mo,
                args,
                CultureInfo.InvariantCulture);
        }

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
            try
            {
                var v = basis == null
                    ? Invoke(mo, "GetOverallProp", property)
                    : Invoke(mo, "GetOverallProp", property, basis);
                var d = ToDoubles(v);
                if (d.Length == 0)
                    throw new InvalidOperationException("empty " + property);
                return d[0];
            }
            catch
            {
                var v = Invoke(mo, "GetProp", property, "overall", Type.Missing, Type.Missing, basis ?? Type.Missing);
                return ToDoubles(v)[0];
            }
        }

        public static double[] GetOverallVector(object mo, string property, string basis)
        {
            try
            {
                return ToDoubles(Invoke(mo, "GetOverallProp", property, basis));
            }
            catch
            {
                return ToDoubles(Invoke(mo, "GetProp", property, "overall", Type.Missing, Type.Missing, basis));
            }
        }

        public static string[] GetCasNumbers(object mo)
        {
            object ids = null, formulae = null, names = null, tb = null, mw = null, cas = null;
            try
            {
                var args = new object[] { ids, formulae, names, tb, mw, cas };
                mo.GetType().InvokeMember(
                    "GetCompoundList",
                    BindingFlags.InvokeMethod,
                    null,
                    mo,
                    args,
                    new ParameterModifier[] { new ParameterModifier(6) },
                    CultureInfo.InvariantCulture,
                    null);
                cas = args[5];
                if (cas != null)
                    return ToStrings(cas);
            }
            catch
            {
                // 1.0: CompIds / GetComponentIds
            }
            try
            {
                return ToStrings(Invoke(mo, "get_ComponentIds"));
            }
            catch
            {
                return ToStrings(Invoke(mo, "GetComponentIds"));
            }
        }

        public static void SetOverall(object mo, string property, object value, string basis = null)
        {
            try
            {
                if (basis == null)
                    Invoke(mo, "SetOverallProp", property, value);
                else
                    Invoke(mo, "SetOverallProp", property, basis, value);
            }
            catch
            {
                Invoke(mo, "SetProp", property, "overall", Type.Missing, Type.Missing, basis ?? Type.Missing, value);
            }
        }

        public static void FlashTP(object mo)
        {
            try
            {
                Invoke(mo, "CalcEquilibrium", "TP", Type.Missing);
            }
            catch
            {
                try { Invoke(mo, "CalcEquilibrium", "TP"); }
                catch { Invoke(mo, "Equilibrium", "TP"); }
            }
        }
    }
}
