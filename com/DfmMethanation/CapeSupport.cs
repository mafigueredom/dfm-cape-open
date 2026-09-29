using System;
using System.Collections.Generic;
using System.Globalization;
using System.Runtime.InteropServices;

namespace PhD.DfmMethanation
{
    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
        [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A001")]
    public class CapeCollection : ICapeCollection
    {
        readonly List<object> _items = new List<object>();
        readonly Dictionary<string, object> _byName =
            new Dictionary<string, object>(StringComparer.OrdinalIgnoreCase);

        public void Add(string name, object item)
        {
            _items.Add(item);
            _byName[name] = item;
        }

        public object Item(object index)
        {
            if (index is string s)
                return _byName[s];
            int i = Convert.ToInt32(index);
            if (i < 1 || i > _items.Count)
                throw new IndexOutOfRangeException("CAPE-OPEN collections are 1-based");
            return _items[i - 1];
        }

        public int Count => _items.Count;

        public IEnumerable<object> All => _items;
    }

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
        [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A002")]
    public class CapePort : ICapeUnitPort, ICapeIdentification
    {
        object _connected;

        public CapePort(string name, CapePortDirection dir)
        {
            ComponentName = name;
            ComponentDescription = name;
            direction = dir;
        }

        public string ComponentName { get; set; }
        public string ComponentDescription { get; set; }
        public CapePortType portType => CapePortType.CAPE_MATERIAL;
        public CapePortDirection direction { get; }
        public object connectedObject => _connected;

        public void Connect(object objectToConnect) => _connected = objectToConnect;
        public void Disconnect() => _connected = null;
    }

    /// <summary>
    /// Shared parameter state. Not a COM class; the typed subclasses are.
    /// </summary>
    [ComVisible(false)]
    public abstract class CapeParameter : ICapeIdentification
    {
        object _value;
        readonly object _default;

        protected CapeParameter(string name, object value, CapeParamMode mode, string desc)
        {
            ComponentName = name;
            ComponentDescription = string.IsNullOrEmpty(desc) ? name : desc;
            _default = value;
            _value = value;
            Mode = mode;
            ValStatus = CapeValidationStatus.CAPE_VALID;
        }

        public string ComponentName { get; set; }
        public string ComponentDescription { get; set; }
        public CapeParamMode Mode { get; set; }
        public CapeValidationStatus ValStatus { get; protected set; }
        public object Specification => this;
        // COFE indexes this array while drawing the parameter dialog. An empty
        // vector is an access violation. Nine zeros means dimensionless.
        public object Dimensionality => new double[9];

        public object value
        {
            get => _value;
            set
            {
                if (value == null)
                    return;
                _value = Coerce(value);
                ValStatus = CapeValidationStatus.CAPE_NOT_VALIDATED;
            }
        }

        public bool Validate(ref string message)
        {
            message = "ok";
            ValStatus = CapeValidationStatus.CAPE_VALID;
            return true;
        }

        public void Reset()
        {
            _value = _default;
            ValStatus = CapeValidationStatus.CAPE_NOT_VALIDATED;
        }

        protected object DefaultObject => _default;
        protected abstract object Coerce(object incoming);

        public double AsDouble() => Convert.ToDouble(_value, CultureInfo.InvariantCulture);
        public int AsInt() => Convert.ToInt32(_value, CultureInfo.InvariantCulture);
        public bool AsBool() => Convert.ToBoolean(_value, CultureInfo.InvariantCulture);
        public string AsString() => Convert.ToString(_value, CultureInfo.InvariantCulture);
    }

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
    [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A003")]
    public class CapeRealParameter : CapeParameter, ICapeParameter, ICapeParameterSpec, ICapeRealParameterSpec
    {
        public CapeRealParameter(string name, double value, CapeParamMode mode, string desc = null)
            : base(name, value, mode, desc)
        {
        }

        public CapeParamType Type => CapeParamType.CAPE_REAL;
        public double DefaultValue => Convert.ToDouble(DefaultObject, CultureInfo.InvariantCulture);
        public double LowerBound => -1e30;
        public double UpperBound => 1e30;

        public bool Validate(double value, ref string message)
        {
            message = "ok";
            return true;
        }

        protected override object Coerce(object incoming) =>
            Convert.ToDouble(incoming, CultureInfo.InvariantCulture);
    }

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
    [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A004")]
    public class CapeIntegerParameter : CapeParameter, ICapeParameter, ICapeParameterSpec, ICapeIntegerParameterSpec
    {
        public CapeIntegerParameter(string name, int value, CapeParamMode mode, string desc = null,
            int lower = -1000000000, int upper = 1000000000)
            : base(name, value, mode, desc)
        {
            LowerBound = lower;
            UpperBound = upper;
        }

        public CapeParamType Type => CapeParamType.CAPE_INT;
        public int DefaultValue => Convert.ToInt32(DefaultObject, CultureInfo.InvariantCulture);
        public int LowerBound { get; }
        public int UpperBound { get; }

        public bool Validate(int value, ref string message)
        {
            if (value < LowerBound || value > UpperBound)
            {
                message = ComponentName + " is outside [" + LowerBound + ", " + UpperBound + "].";
                return false;
            }
            message = "ok";
            return true;
        }

        protected override object Coerce(object incoming) =>
            Convert.ToInt32(incoming, CultureInfo.InvariantCulture);
    }

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
    [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A005")]
    public class CapeBooleanParameter : CapeParameter, ICapeParameter, ICapeParameterSpec, ICapeBooleanParameterSpec
    {
        public CapeBooleanParameter(string name, bool value, CapeParamMode mode, string desc = null)
            : base(name, value, mode, desc)
        {
        }

        public CapeParamType Type => CapeParamType.CAPE_BOOLEAN;
        public bool DefaultValue => Convert.ToBoolean(DefaultObject, CultureInfo.InvariantCulture);

        public bool Validate(bool value, ref string message)
        {
            message = "ok";
            return true;
        }

        protected override object Coerce(object incoming)
        {
            if (incoming is bool b)
                return b;
            if (incoming is short s)
                return s != 0;
            if (incoming is int i)
                return i != 0;
            return Convert.ToBoolean(incoming, CultureInfo.InvariantCulture);
        }
    }

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
    [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A006")]
    public class CapeOptionParameter : CapeParameter, ICapeParameter, ICapeParameterSpec, ICapeOptionParameterSpec
    {
        readonly string[] _options;

        public CapeOptionParameter(string name, string value, string[] options, bool restricted,
            CapeParamMode mode, string desc = null)
            : base(name, value ?? "", mode, desc)
        {
            _options = options ?? Array.Empty<string>();
            RestrictedToList = restricted;
        }

        public CapeParamType Type => CapeParamType.CAPE_OPTION;
        public string DefaultValue => Convert.ToString(DefaultObject, CultureInfo.InvariantCulture);
        public object OptionList => _options.Length == 0 ? new[] { DefaultValue } : _options;
        public bool RestrictedToList { get; }

        public bool Validate(string value, ref string message)
        {
            if (!RestrictedToList)
            {
                message = "ok";
                return true;
            }
            foreach (var opt in _options)
            {
                if (string.Equals(opt, value, StringComparison.OrdinalIgnoreCase))
                {
                    message = "ok";
                    return true;
                }
            }
            message = ComponentName + " must be one of: " + string.Join(", ", _options);
            return false;
        }

        protected override object Coerce(object incoming) =>
            Convert.ToString(incoming, CultureInfo.InvariantCulture) ?? "";
    }
}
