using System;
using System.Collections.Generic;
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

        public static void GetTP(object mo, out double temperature, out double pressure)
        {
            object raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
            try
            {
                var mat = (ICapeThermoMaterial11)raw;
                object composition;
                mat.GetOverallTPFraction(out temperature, out pressure, out composition);
                CapeTrace.Write(
                    "GetOverallTPFraction T=" + temperature.ToString(CultureInfo.InvariantCulture)
                    + " P=" + pressure.ToString(CultureInfo.InvariantCulture));
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        public static double GetMolarFlow(object mo)
        {
            var values = GetOverall(mo, "totalFlow", "mole");
            if (values.Length == 0)
                throw new InvalidOperationException("empty flow");
            double sum = 0.0;
            foreach (var v in values)
                sum += v;
            CapeTrace.Write(
                "flow sum " + sum.ToString(CultureInfo.InvariantCulture) + " n=" + values.Length);
            return sum;
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
            object raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
            try
            {
                var mat = (ICapeThermoMaterial11)raw;
                var names = property == "totalFlow"
                    ? new[] { "flow", "totalFlow" }
                    : new[] { property };
                Exception last = null;
                foreach (var name in names)
                {
                    foreach (var b in BasisForms(basis))
                    {
                        try
                        {
                            object results;
                            mat.GetOverallProp(name, b, out results);
                            var values = ToDoubles(results);
                            CapeTrace.Write("GetOverallProp " + name + " " + b + " n=" + values.Length);
                            return values;
                        }
                        catch (Exception ex)
                        {
                            last = ex;
                            ClearError();
                            CapeTrace.Write("GetOverallProp " + name + " " + b + " " + ex.Message);
                        }
                    }
                }
                throw new InvalidOperationException(
                    "GetOverallProp " + property + ": " + (last == null ? "failed" : last.Message),
                    last);
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        public static string[] GetCasNumbers(object mo)
        {
            object raw = AsInterface(mo, IidCompounds, typeof(ICapeThermoCompounds));
            var compounds = (ICapeThermoCompounds)raw;
            object compIds, formulae, names, boilTemps, molwts, casnos;
            CapeTrace.Write("GetCompoundList call");
            try
            {
                compounds.GetCompoundList(
                    out compIds, out formulae, out names, out boilTemps, out molwts, out casnos);
            }
            catch (COMException ex)
            {
                throw new InvalidOperationException(
                    "GetCompoundList 0x" + ex.ErrorCode.ToString("X8") + " " + ex.Message, ex);
            }
            finally
            {
                ReleaseRcw(raw);
            }
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

        [DllImport("oleaut32.dll", PreserveSig = true)]
        static extern int SetErrorInfo(int dwReserved, IntPtr perrinfo);

        public static void ClearError()
        {
            try { SetErrorInfo(0, IntPtr.Zero); }
            catch { }
        }

        public static void ReleaseRcw(object comObject)
        {
            try
            {
                if (comObject != null && Marshal.IsComObject(comObject))
                    Marshal.ReleaseComObject(comObject);
            }
            catch { }
        }

        static Array OneBased(double[] values)
        {
            var arr = Array.CreateInstance(typeof(double), new[] { values.Length }, new[] { 1 });
            for (int i = 0; i < values.Length; i++)
                arr.SetValue(values[i], i + 1);
            return arr;
        }

        public static bool TrySetOverall(object mo, string property, object value, string basis = null)
        {
            try
            {
                SetOverall(mo, property, value, basis);
                return true;
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write("TrySetOverall skip " + property + " " + ex.GetBaseException().Message);
                return false;
            }
        }

        public static void SetOverall(object mo, string property, object value, string basis = null)
        {
            object raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
            try
            {
                var mat = (ICapeThermoMaterial11)raw;
                var names = property == "totalFlow"
                    ? new[] { "totalFlow" }
                    : new[] { property };
                double[] numbers = value is double[] direct
                    ? direct
                    : new[] { Convert.ToDouble(value, CultureInfo.InvariantCulture) };
                // TEA accepted a 0-based array for flow and fraction. A 1-based
                // array came back as 0x80040501, so that form is only a fallback.
                var payloads = new object[] { numbers, OneBased(numbers) };
                Exception last = null;
                foreach (var payload in payloads)
                {
                    var bound = payload is double[] ? "lb=0" : "lb=1";
                    foreach (var name in names)
                    {
                        foreach (var b in BasisForms(basis))
                        {
                            try
                            {
                                mat.SetOverallProp(name, b, payload);
                                CapeTrace.Write(
                                    "SetOverallProp " + name + " " + b + " ok n=" + numbers.Length + " " + bound);
                                return;
                            }
                            catch (Exception ex)
                            {
                                last = ex;
                                ClearError();
                                CapeTrace.Write("SetOverallProp " + name + " " + b + " " + bound + " " + ex.Message);
                            }
                        }
                    }
                }
                throw new InvalidOperationException(
                    "SetOverallProp " + property + ": " + (last == null ? "failed" : last.Message),
                    last);
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        static string[] BasisForms(string basis)
        {
            if (string.IsNullOrEmpty(basis) ||
                basis.Equals("undefined", StringComparison.OrdinalIgnoreCase))
                return new[] { "UNDEFINED", "" };
            if (basis.Equals("mole", StringComparison.OrdinalIgnoreCase))
                return new[] { "Mole", "mole" };
            if (basis.Equals("mass", StringComparison.OrdinalIgnoreCase))
                return new[] { "Mass", "mass" };
            return new[] { basis };
        }

        static readonly Guid IidEquilibrium = new Guid("678C0AA0-7D66-11D2-A67D-00105A42887F");

        public static void LogPhases(object mo, string tag)
        {
            object raw = null;
            try
            {
                raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
                object labels, status;
                ((ICapeThermoMaterial11)raw).GetPresentPhases(out labels, out status);
                CapeTrace.Write(
                    tag + " phases=" + string.Join(",", ToStrings(labels))
                    + " status=" + string.Join(",", ToStrings(status)));
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write(tag + " phases " + ex.GetBaseException().Message);
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        static readonly Guid IidPhases = new Guid("678C0A9E-7D66-11D2-A67D-00105A42887F");

        public static void EnsurePhases(object mo)
        {
            var known = PhaseLabels(mo);
            var candidates = new List<string[]>();
            if (known.Length > 0)
                candidates.Add(known);
            candidates.Add(new[] { "Vapour", "Liquid" });
            candidates.Add(new[] { "Vapor", "Liquid" });
            candidates.Add(new[] { "gas", "liquid" });
            foreach (var labels in candidates)
            {
                if (TrySetPhases(mo, labels))
                    return;
            }
            CapeTrace.Write("SetPresentPhases none accepted");
        }

        static string[] PhaseLabels(object mo)
        {
            object raw = null;
            try
            {
                raw = AsInterface(mo, IidPhases, typeof(ICapeThermoPhases));
                object labels, state, key;
                ((ICapeThermoPhases)raw).GetPhaseList(out labels, out state, out key);
                var names = ToStrings(labels);
                CapeTrace.Write(
                    "GetPhaseList " + string.Join(",", names)
                    + " state=" + string.Join(",", ToStrings(state)));
                return names;
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write("GetPhaseList " + ex.GetBaseException().Message);
                return new string[0];
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        static bool TrySetPhases(object mo, string[] labels)
        {
            var statusText = new string[labels.Length];
            var statusInt = new int[labels.Length];
            for (int i = 0; i < labels.Length; i++)
            {
                statusText[i] = "Cape_UnknownPhaseStatus";
                statusInt[i] = 0;
            }
            var attempts = new[]
            {
                new object[] { labels, statusText },
                new object[] { OneBasedStrings(labels), OneBasedStrings(statusText) },
                new object[] { labels, statusInt },
                new object[] { OneBasedStrings(labels), OneBasedInts(statusInt) }
            };
            object raw = null;
            try
            {
                raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
                var mat = (ICapeThermoMaterial11)raw;
                foreach (var attempt in attempts)
                {
                    try
                    {
                        mat.SetPresentPhases(attempt[0], attempt[1]);
                        CapeTrace.Write("SetPresentPhases " + string.Join(",", labels));
                        return true;
                    }
                    catch (Exception ex)
                    {
                        ClearError();
                        CapeTrace.Write(
                            "SetPresentPhases " + string.Join(",", labels) + " "
                            + ex.GetBaseException().Message);
                    }
                }
                return false;
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        static Array OneBasedStrings(string[] values)
        {
            var arr = Array.CreateInstance(typeof(string), new[] { values.Length }, new[] { 1 });
            for (int i = 0; i < values.Length; i++)
                arr.SetValue(values[i], i + 1);
            return arr;
        }

        static Array OneBasedInts(int[] values)
        {
            var arr = Array.CreateInstance(typeof(int), new[] { values.Length }, new[] { 1 });
            for (int i = 0; i < values.Length; i++)
                arr.SetValue(values[i], i + 1);
            return arr;
        }

        public static bool FlashEquilibrium(object mo)
        {
            object raw = null;
            try
            {
                raw = AsInterface(mo, IidEquilibrium, typeof(ICapeThermoEquilibriumRoutine11));
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write("Flash no equilibrium routine " + ex.GetBaseException().Message);
                return false;
            }
            try
            {
                ((ICapeThermoEquilibriumRoutine11)raw).CalcEquilibrium(
                    "temperature", "pressure", "Unspecified");
                CapeTrace.Write("CalcEquilibrium temperature pressure ok");
                return true;
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write("CalcEquilibrium " + ex.GetBaseException().Message);
                return false;
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        public static void MarkVaporFlashed(object mo, double[] fraction, double temperature, double pressure)
        {
            // CAPE_ATEQUILIBRIUM is 1. COFE rejects an outlet whose phases stay at 0.
            // A single-phase gas makes TEA's own TP flash return an error, so the
            // vapor phase is stored directly after the overall composition.
            object raw = null;
            try
            {
                raw = AsInterface(mo, IidMaterial, typeof(ICapeThermoMaterial11));
                var mat = (ICapeThermoMaterial11)raw;
                TryPhaseProp(mat, "fraction", "Vapor", "Mole", fraction);
                TryPhaseProp(mat, "phaseFraction", "Vapor", "Mole", new[] { 1.0 });
                TryPhaseProp(mat, "temperature", "Vapor", "", new[] { temperature });
                TryPhaseProp(mat, "pressure", "Vapor", "", new[] { pressure });
                mat.SetPresentPhases(new[] { "Vapor", "Liquid" }, new[] { 1, 0 });
                CapeTrace.Write("SetPresentPhases Vapor equilibrium");
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write("MarkVaporFlashed " + ex.GetBaseException().Message);
            }
            finally
            {
                ReleaseRcw(raw);
            }
        }

        static void TryPhaseProp(
            ICapeThermoMaterial11 mat, string property, string phase, string basis, double[] values)
        {
            try
            {
                mat.SetSinglePhaseProp(property, phase, basis, values);
                CapeTrace.Write("SetSinglePhaseProp " + property + " " + phase + " ok");
            }
            catch (Exception ex)
            {
                ClearError();
                CapeTrace.Write(
                    "SetSinglePhaseProp " + property + " " + phase + " "
                    + ex.GetBaseException().Message);
            }
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

    // Dual, so the CLR uses the vtable. IDispatch::Invoke asks TEA for a
    // type library and TEA returns TYPE_E_ELEMENTNOTFOUND. Slots follow
    // the CAPE-OPEN dispids: GetCompoundList is the second method.
    [ComImport]
    [Guid("678C0A9D-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    internal interface ICapeThermoCompounds
    {
        void GetCompoundConstant(object props, object compIds);

        void GetCompoundList(
            [MarshalAs(UnmanagedType.Struct)] out object compIds,
            [MarshalAs(UnmanagedType.Struct)] out object formulae,
            [MarshalAs(UnmanagedType.Struct)] out object names,
            [MarshalAs(UnmanagedType.Struct)] out object boilTemps,
            [MarshalAs(UnmanagedType.Struct)] out object molwts,
            [MarshalAs(UnmanagedType.Struct)] out object casnos);
    }

    // GetOverallProp is dispid 4, SetOverallProp is dispid 10.
    [ComImport]
    [Guid("678C0A9B-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    internal interface ICapeThermoMaterial11
    {
        void ClearAllProps();
        void CopyFromMaterial(object source);
        object CreateMaterial();

        void GetOverallProp(
            [MarshalAs(UnmanagedType.BStr)] string property,
            [MarshalAs(UnmanagedType.BStr)] string basis,
            [MarshalAs(UnmanagedType.Struct)] out object results);

        void GetOverallTPFraction(
            out double temperature,
            out double pressure,
            [MarshalAs(UnmanagedType.Struct)] out object composition);
        void GetPresentPhases(out object phaseLabels, out object phaseStatus);
        void GetSinglePhaseProp(string property, string phaseLabel, string basis, out object results);
        void GetTPFraction(string phaseLabel, out object temperature, out object pressure, out object composition);
        void GetTwoPhaseProp(string property, object phaseLabels, string basis, out object results);

        void SetOverallProp(
            [MarshalAs(UnmanagedType.BStr)] string property,
            [MarshalAs(UnmanagedType.BStr)] string basis,
            [MarshalAs(UnmanagedType.Struct)] object values);

        // Dispid 11. Must stay immediately after SetOverallProp.
        void SetPresentPhases(
            [MarshalAs(UnmanagedType.Struct)] object phaseLabels,
            [MarshalAs(UnmanagedType.Struct)] object phaseStatus);

        // Dispid 12. Must stay immediately after SetPresentPhases.
        void SetSinglePhaseProp(
            [MarshalAs(UnmanagedType.BStr)] string property,
            [MarshalAs(UnmanagedType.BStr)] string phaseLabel,
            [MarshalAs(UnmanagedType.BStr)] string basis,
            [MarshalAs(UnmanagedType.Struct)] object values);
    }

    // GetPhaseList is dispid 3.
    [ComImport]
    [Guid("678C0A9E-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    internal interface ICapeThermoPhases
    {
        int GetNumPhases();
        object GetPhaseInfo(string phaseLabel, string phaseAttribute);
        void GetPhaseList(out object phaseLabels, out object stateOfAggregation, out object keyCompoundId);
    }

    // Dispid 1 on this interface, so it is the first vtable method.
    [ComImport]
    [Guid("678C0AA0-7D66-11D2-A67D-00105A42887F")]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    internal interface ICapeThermoEquilibriumRoutine11
    {
        void CalcEquilibrium(
            [MarshalAs(UnmanagedType.Struct)] object specification1,
            [MarshalAs(UnmanagedType.Struct)] object specification2,
            [MarshalAs(UnmanagedType.BStr)] string solutionType);
    }
}
