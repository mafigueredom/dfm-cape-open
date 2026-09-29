using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using System.Web.Script.Serialization;
using Microsoft.Win32;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// CAPE-OPEN 1.1 Unit Operation. Thin COM wrapper: Material Objects → JSON
    /// → Docker (<c>python -m adr.cape</c>) → product / outputs / reports.
    /// Physics stay out of the PME process. Not CORBA. Not AmsterCHEM Python UO.
    /// </summary>
    [ComVisible(true)]
    [Guid(Guids.Clsid)]
    [ProgId(Guids.ProgId)]
    [ClassInterface(ClassInterfaceType.None)]
    [ComDefaultInterface(typeof(ICapeUnit))]
    public class DfmUnit : ICapeUnit, ICapeUtilities, ICapeIdentification, ICapeUnitReport
    {
        static readonly string[] Species = { "CO2", "H2", "CH4", "H2O", "N2" };
        static readonly Dictionary<string, string> CasToSp = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
        {
            { "124-38-9", "CO2" }, { "1333-74-0", "H2" }, { "74-82-8", "CH4" },
            { "7732-18-5", "H2O" }, { "7727-37-9", "N2" }
        };
        static readonly Dictionary<string, string> SpToCas = new Dictionary<string, string>
        {
            { "CO2", "124-38-9" }, { "H2", "1333-74-0" }, { "CH4", "74-82-8" },
            { "H2O", "7732-18-5" }, { "N2", "7727-37-9" }
        };

        readonly CapeCollection _ports = new CapeCollection();
        readonly CapeCollection _params = new CapeCollection();
        readonly Dictionary<string, CapeParameter> _p = new Dictionary<string, CapeParameter>(StringComparer.OrdinalIgnoreCase);
        readonly JavaScriptSerializer _json = new JavaScriptSerializer { MaxJsonLength = int.MaxValue, RecursionLimit = 256 };
        CapePort _feedAds, _feedPurge, _feedRxn, _feedPurge2, _product;
        Dictionary<string, object> _lastResult;
        static Dictionary<string, object> _sharedResult;
        string _selectedReport = "atom_balance";
        object _simulationContext;
        CapeValidationStatus _val = CapeValidationStatus.CAPE_NOT_VALIDATED;

        public DfmUnit()
        {
            CapeTrace.Guard("DfmUnit", () =>
            {
                ComponentName = Guids.UnitName;
                ComponentDescription =
                    "Isothermal DFM methanation packed bed (cycle-average product).";
                BuildPorts();
                BuildParameters();
            });
        }

        public string ComponentName { get; set; }
        public string ComponentDescription { get; set; }
        public ICapeCollection ports => _ports;
        public ICapeCollection Parameters => _params;
        public CapeValidationStatus ValStatus => _val;
        public object simulationContext
        {
            get => _simulationContext;
            set
            {
                CapeTrace.Write("simulationContext " + (value == null ? "null" : "set"));
                _simulationContext = value;
            }
        }
        public object reports
        {
            get
            {
                var names = new[] { "atom_balance", "outlet_C_raw", "outlet_C_sopdt", "profiles_C_z" };
                var arr = Array.CreateInstance(typeof(string), new[] { names.Length }, new[] { 1 });
                for (int i = 0; i < names.Length; i++)
                    arr.SetValue(names[i], i + 1);
                CapeTrace.Write("reports n=" + names.Length);
                return arr;
            }
        }
        public string selectedReport
        {
            get => _selectedReport;
            set => _selectedReport = value;
        }

        public void Initialize()
        {
            CapeTrace.Write("Initialize");
        }

        public void Terminate()
        {
            CapeTrace.Write("Terminate");
        }

        public void Edit()
        {
            CapeTrace.Guard("Edit", () =>
            {
                var inputs = new List<CapeParameter>();
                foreach (var item in _params.All)
                {
                    var p = (CapeParameter)item;
                    if (p.Mode == CapeParamMode.CAPE_INPUT || p.Mode == CapeParamMode.CAPE_INPUT_OUTPUT)
                        inputs.Add(p);
                }
                using (var dlg = new ParameterEditorForm(inputs))
                {
                    CapeTrace.Write("Edit show");
                    var result = dlg.ShowDialog();
                    CapeTrace.Write("Edit closed " + result);
                }
            });
        }

        public bool Validate(ref string message)
        {
            CapeTrace.Write(">> Validate");
            try
            {
                if (!_feedAds.IsConnected ||
                    !_feedPurge.IsConnected ||
                    !_feedRxn.IsConnected)
                {
                    message = "Feed_ads, Feed_purge and Feed_rxn must be connected.";
                    _val = CapeValidationStatus.CAPE_INVALID;
                    CapeTrace.Write("<< Validate invalid feeds");
                    return false;
                }
                if (!_product.IsConnected)
                {
                    message = "Product must be connected.";
                    _val = CapeValidationStatus.CAPE_INVALID;
                    CapeTrace.Write("<< Validate invalid product");
                    return false;
                }
                var t2 = P("t_purge2").AsDouble();
                if (t2 > 0.0 && !_feedPurge2.IsConnected)
                {
                    message = "t_purge2 > 0 requires Feed_purge2 connected.";
                    _val = CapeValidationStatus.CAPE_INVALID;
                    return false;
                }
                message = "ok";
                _val = CapeValidationStatus.CAPE_VALID;
                CapeTrace.Write("<< Validate ok");
                return true;
            }
            catch (Exception ex)
            {
                CapeTrace.Write("!! Validate" + Environment.NewLine + ex);
                message = ex.Message;
                _val = CapeValidationStatus.CAPE_INVALID;
                return false;
            }
        }

        public void Calculate()
        {
            // A failed HRESULT from this callback is unhandled in COFE and the process exits.
            CapeTrace.Write(">> Calculate");
            try
            {
                string msg = "";
                if (!Validate(ref msg))
                    throw new InvalidOperationException(msg);

                var request = BuildRequest();
                var work = Path.Combine(Path.GetTempPath(), "dfm-cape-" + Guid.NewGuid().ToString("N"));
                Directory.CreateDirectory(work);
                try
                {
                    File.WriteAllText(
                        Path.Combine(work, "cape_request.json"),
                        _json.Serialize(request),
                        new UTF8Encoding(false));
                    var extra = new StringBuilder();
                    var tEnd = P("t_end_override").AsDouble();
                    if (tEnd > 0.0)
                        extra.Append(" --t-end ").Append(tEnd.ToString(CultureInfo.InvariantCulture));
                    var raw = EngineClient.Run(work, P("docker_image").AsString(), extra.ToString());
                    _lastResult = _json.Deserialize<Dictionary<string, object>>(raw);
                    _sharedResult = _lastResult;
                    WriteResultFile();
                    ApplyOutputs();
                    _val = CapeValidationStatus.CAPE_VALID;
                    CapeTrace.Write("<< Calculate");
                }
                finally
                {
                    try { Directory.Delete(work, true); } catch { }
                }
            }
            catch (Exception ex)
            {
                ThermoBridge.ClearError();
                _val = CapeValidationStatus.CAPE_INVALID;
                CapeTrace.Write("!! Calculate" + Environment.NewLine + ex);
            }
        }

        public string ProduceReport()
        {
            CapeTrace.Write("ProduceReport " + _selectedReport);
            if (Result == null)
                return "No Calculate() result yet.\n";
            var fromCatalog = CatalogText(_selectedReport);
            if (!string.IsNullOrEmpty(fromCatalog))
                return fromCatalog;
            if (string.Equals(_selectedReport, "atom_balance", StringComparison.OrdinalIgnoreCase)
                && Result.TryGetValue("report_text", out var t) && t != null)
                return Convert.ToString(t, CultureInfo.InvariantCulture);
            return "Report '" + _selectedReport + "' has no text on this run.\n";
        }

        Dictionary<string, object> Result => _lastResult ?? _sharedResult;

        string CatalogText(string name)
        {
            var entry = Dict(Dict(Dict(Result, "reports"), "catalog"), name);
            if (!entry.TryGetValue("text", out var text) || text == null)
                return null;
            return Convert.ToString(text, CultureInfo.InvariantCulture);
        }

        Dictionary<string, object> BuildRequest()
        {
            var feeds = new Dictionary<string, object>();
            feeds["ads"] = ReadFeed(_feedAds, true);
            feeds["purge"] = ReadFeed(_feedPurge, true);
            feeds["rxn"] = ReadFeed(_feedRxn, true);
            feeds["purge2"] = P("t_purge2").AsDouble() > 0.0
                ? ReadFeed(_feedPurge2, true)
                : null;

            var parameters = new Dictionary<string, object>
            {
                { "geometry", Nested("L", "R", "d_p", "epsilon", "rho_a") },
                { "adsorption", Nested("k_ads", "k_des", "q_e", "K_C", "K_R", "a_csc") },
                { "reaction", Nested("A_rxn", "Ea", "m1", "m2") },
                { "cycle", new Dictionary<string, object>
                    {
                        { "n_cycles", P("n_cycles").AsInt() },
                        { "t_ads", P("t_ads").AsDouble() },
                        { "t_purge", P("t_purge").AsDouble() },
                        { "t_rxn", P("t_rxn").AsDouble() },
                        { "t_purge2", P("t_purge2").AsDouble() }
                    }
                }
            };

            return new Dictionary<string, object>
            {
                { "schema", "cape_request_v1" },
                { "factory_id", "methanation_dfm" },
                { "config_path", "config/dfm_config.pow_v07_10_60_opt-1.json" },
                { "C_t_export", P("C_t_export").AsString() },
                { "balance_tol", P("balance_tol").AsDouble() },
                { "fail_if_unbalanced", P("fail_if_unbalanced").AsBool() },
                { "feeds", feeds },
                { "parameters", parameters }
            };
        }

        Dictionary<string, object> Nested(params string[] names)
        {
            var d = new Dictionary<string, object>();
            foreach (var n in names)
                d[n] = P(n).AsDouble();
            return d;
        }

        Dictionary<string, object> ReadFeed(CapePort port, bool required)
        {
            CapeTrace.Write("ReadFeed " + port.ComponentName);
            var mo = port.MaterialObject;
            if (mo == null)
            {
                if (required)
                    throw new InvalidOperationException(port.ComponentName + " is not connected");
                return null;
            }
            try
            {
                return ReadFeed(mo, port.ComponentName);
            }
            finally
            {
                ThermoBridge.ReleaseRcw(mo);
            }
        }

        Dictionary<string, object> ReadFeed(object mo, string portName)
        {
            var cas = ThermoBridge.GetCasNumbers(mo);
            var yvec = ThermoBridge.GetOverallVector(mo, "fraction", "mole");
            var y = new Dictionary<string, double>();
            int n = Math.Min(cas.Length, yvec.Length);
            for (int i = 0; i < n; i++)
            {
                if (CasToSp.ContainsKey(cas[i]))
                    y[cas[i]] = yvec[i];
            }
            double sum = 0.0;
            foreach (var kv in y)
                sum += kv.Value;
            if (sum <= 0.0)
                throw new InvalidOperationException(portName + ": no mapped CAS compounds");
            if (Math.Abs(sum - 1.0) > 1e-6)
            {
                var keys = new List<string>(y.Keys);
                foreach (var k in keys)
                    y[k] = y[k] / sum;
            }
            var yObj = new Dictionary<string, object>();
            foreach (var kv in y)
                yObj[kv.Key] = kv.Value;
            double pressure;
            ThermoBridge.GetTP(mo, out _, out pressure);
            return new Dictionary<string, object>
            {
                { "F_mol_s", ThermoBridge.GetMolarFlow(mo) },
                { "P_Pa", pressure },
                { "y", yObj }
            };
        }

        void ApplyOutputs()
        {
            var outputs = Dict(_lastResult, "outputs");
            foreach (var name in new[]
                     { "R_C_rel", "R_H_rel", "R_O_rel", "Y_CH4", "n_CH4", "N_CO2_ads" })
            {
                if (outputs.ContainsKey(name) && outputs[name] != null)
                {
                    var number = ToD(outputs, name);
                    SolvedOutputs.Values[name] = number;
                    P(name).value = number;
                }
            }
            if (outputs.ContainsKey("balance_ok"))
            {
                var ok = Convert.ToBoolean(outputs["balance_ok"], CultureInfo.InvariantCulture);
                SolvedOutputs.Values["balance_ok"] = ok;
                P("balance_ok").value = ok;
            }
        }

        void WriteResultFile()
        {
            var path = Path.Combine(Path.GetTempPath(), "dfm-cape-open-last.txt");
            var product = Dict(_lastResult, "product");
            var fss = Dict(product, "F_ss_mol_s");
            var outputs = Dict(_lastResult, "outputs");
            var text = new StringBuilder();
            text.AppendLine("DFM methanation solve finished");
            text.AppendLine("T_K=" + ToD(product, "T_K").ToString(CultureInfo.InvariantCulture));
            text.AppendLine("P_Pa=" + ToD(product, "P_Pa").ToString(CultureInfo.InvariantCulture));
            foreach (var sp in Species)
                text.AppendLine(sp + "_mol_s=" + ToD(fss, sp).ToString("G6", CultureInfo.InvariantCulture));
            foreach (var name in new[] { "R_C_rel", "R_H_rel", "R_O_rel", "Y_CH4", "n_CH4", "N_CO2_ads", "balance_ok" })
            {
                if (outputs.ContainsKey(name) && outputs[name] != null)
                    text.AppendLine(name + "=" + Convert.ToString(outputs[name], CultureInfo.InvariantCulture));
            }
            File.WriteAllText(path, text.ToString());
            CapeTrace.Write("result file " + path);
            CapeTrace.Write(
                "solved CH4_mol_s=" + ToD(fss, "CH4").ToString("G6", CultureInfo.InvariantCulture)
                + " Y_CH4=" + Convert.ToString(outputs.ContainsKey("Y_CH4") ? outputs["Y_CH4"] : "", CultureInfo.InvariantCulture)
                + " balance_ok=" + Convert.ToString(outputs.ContainsKey("balance_ok") ? outputs["balance_ok"] : "", CultureInfo.InvariantCulture));
        }

        static Dictionary<string, object> Dict(Dictionary<string, object> parent, string key)
        {
            if (!parent.TryGetValue(key, out var o) || o == null)
                return new Dictionary<string, object>();
            if (o is Dictionary<string, object> d)
                return d;
            if (o is IDictionary raw)
            {
                var c = new Dictionary<string, object>();
                foreach (DictionaryEntry e in raw)
                    c[Convert.ToString(e.Key, CultureInfo.InvariantCulture)] = e.Value;
                return c;
            }
            return new Dictionary<string, object>();
        }

        static double ToD(Dictionary<string, object> d, string key)
        {
            if (!d.TryGetValue(key, out var o) || o == null)
                return 0.0;
            return Convert.ToDouble(o, CultureInfo.InvariantCulture);
        }

        CapeParameter P(string name) => _p[name];

        void AddP(string name, object value, CapeParamMode mode, string desc = null)
        {
            CapeParameter p;
            if (value is bool b)
                p = new CapeBooleanParameter(name, b, mode, desc);
            else if (value is int n)
                p = new CapeIntegerParameter(name, n, mode, desc);
            else if (value is string s)
                p = new CapeOptionParameter(name, s, Array.Empty<string>(), false, mode, desc);
            else
                p = new CapeRealParameter(name, Convert.ToDouble(value, CultureInfo.InvariantCulture), mode, desc);
            _p[name] = p;
            _params.Add(name, p);
        }

        void AddOption(string name, string value, string[] options, bool restricted, CapeParamMode mode, string desc = null)
        {
            var p = new CapeOptionParameter(name, value, options, restricted, mode, desc);
            _p[name] = p;
            _params.Add(name, p);
        }

        void BuildPorts()
        {
            _feedAds = new CapePort("Feed_ads", CapePortDirection.CAPE_INLET);
            _feedPurge = new CapePort("Feed_purge", CapePortDirection.CAPE_INLET);
            _feedRxn = new CapePort("Feed_rxn", CapePortDirection.CAPE_INLET);
            _feedPurge2 = new CapePort("Feed_purge2", CapePortDirection.CAPE_INLET);
            _product = new CapePort("Product", CapePortDirection.CAPE_OUTLET);
            _ports.Add("Feed_ads", _feedAds);
            _ports.Add("Feed_purge", _feedPurge);
            _ports.Add("Feed_rxn", _feedRxn);
            _ports.Add("Feed_purge2", _feedPurge2);
            _ports.Add("Product", _product);
        }

        void BuildParameters()
        {
            AddP("L", 0.06, CapeParamMode.CAPE_INPUT, "Bed length [m]");
            AddP("R", 0.00225, CapeParamMode.CAPE_INPUT, "Bed radius [m]; D=2R");
            AddP("d_p", 0.0001525, CapeParamMode.CAPE_INPUT, "Particle diameter [m]");
            AddP("epsilon", 0.3865222222, CapeParamMode.CAPE_INPUT, "Bed voidage (typed)");
            AddP("rho_a", 576.3635799, CapeParamMode.CAPE_INPUT, "Apparent density [kg/m3_bed]");
            AddP("k_ads", 3.3333333333333335, CapeParamMode.CAPE_INPUT, "LDF k_ads [1/s]");
            AddP("k_des", 1.7e-6, CapeParamMode.CAPE_INPUT, "H2 desorption k [m3/(kg s)]");
            AddP("q_e", 0.1681, CapeParamMode.CAPE_INPUT, "CSC q_e [mol/kg]");
            AddP("K_C", 245.7369, CapeParamMode.CAPE_INPUT, "CSC K_C (per atm)");
            AddP("K_R", 2.010374, CapeParamMode.CAPE_INPUT);
            AddP("a_csc", 1.33012, CapeParamMode.CAPE_INPUT);
            AddP("A_rxn", 16491.0, CapeParamMode.CAPE_INPUT);
            AddP("Ea", 14560.0, CapeParamMode.CAPE_INPUT, "Ea [J/mol]");
            AddP("m1", 1.7, CapeParamMode.CAPE_INPUT);
            AddP("m2", 5.0, CapeParamMode.CAPE_INPUT);
            AddP("n_cycles", 1, CapeParamMode.CAPE_INPUT, "Number of cycles");
            AddP("t_ads", 400.0, CapeParamMode.CAPE_INPUT, "Adsorption duration [s]");
            AddP("t_purge", 180.0, CapeParamMode.CAPE_INPUT);
            AddP("t_rxn", 300.0, CapeParamMode.CAPE_INPUT);
            AddP("t_purge2", 0.0, CapeParamMode.CAPE_INPUT);
            AddOption("C_t_export", "both", new[] { "raw", "sopdt", "both" }, true, CapeParamMode.CAPE_INPUT,
                "Report series: raw, sopdt, or both");
            AddP("balance_tol", 0.001, CapeParamMode.CAPE_INPUT);
            AddP("fail_if_unbalanced", false, CapeParamMode.CAPE_INPUT);
            AddOption("docker_image", "dfm-methanation-cape:v1", Array.Empty<string>(), false,
                CapeParamMode.CAPE_INPUT, "Docker engine image");
            AddP("t_end_override", 0.0, CapeParamMode.CAPE_INPUT, "0 = full horizon; >0 passed as --t-end");
            AddP("R_C_rel", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("R_H_rel", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("R_O_rel", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("Y_CH4", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("n_CH4", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("N_CO2_ads", 0.0, CapeParamMode.CAPE_OUTPUT);
            AddP("balance_ok", false, CapeParamMode.CAPE_OUTPUT);
        }

        [ComRegisterFunction]
        public static void RegisterFunction(Type t)
        {
            var key = Registry.ClassesRoot.OpenSubKey("CLSID\\{" + Guids.Clsid + "}", true);
            if (key == null)
                return;
            var cats = key.CreateSubKey("Implemented Categories");
            cats.CreateSubKey(Guids.CatidUnit).Close();
            cats.CreateSubKey(Guids.CatidPmc).Close();
            cats.CreateSubKey(Guids.CatidNetCom).Close();
            cats.CreateSubKey(Guids.CatidConsumesThermo).Close();
            cats.CreateSubKey(Guids.CatidThermo10).Close();
            cats.CreateSubKey(Guids.CatidThermo11).Close();
            cats.Close();
            var desc = key.CreateSubKey("CapeDescription");
            desc.SetValue("Name", Guids.UnitName);
            desc.SetValue("Description",
                "DFM methanation packed bed; cycle-average product; FEniCSx via Docker.");
            desc.SetValue("CapeVersion", Guids.CapeVersion);
            desc.SetValue("ComponentVersion", t.Assembly.GetName().Version.ToString());
            desc.SetValue("VendorURL", "https://www.colan.org/");
            desc.SetValue("About", "Doctoral research DFM CAPE-OPEN 1.1 unit.");
            desc.Close();
            key.Close();
        }

        [ComUnregisterFunction]
        public static void UnregisterFunction(Type t)
        {
            try
            {
                var key = Registry.ClassesRoot.OpenSubKey("CLSID\\{" + Guids.Clsid + "}", true);
                if (key == null)
                    return;
                key.DeleteSubKeyTree("Implemented Categories", false);
                key.DeleteSubKeyTree("CapeDescription", false);
                key.Close();
            }
            catch
            {
                // already removed
            }
        }
    }
}
