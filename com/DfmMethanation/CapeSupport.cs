using System;
using System.Collections.Generic;
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

    [ComVisible(true)]
    [ClassInterface(ClassInterfaceType.None)]
        [Guid("B3E1C0A4-9F2D-4A77-8C11-0DF100C0A003")]
    public class CapeParameter : ICapeParameter, ICapeIdentification
    {
        object _value;

        public CapeParameter(string name, object value, CapeParamMode mode, string desc = null)
        {
            ComponentName = name;
            ComponentDescription = desc ?? name;
            _value = value;
            Mode = mode;
        }

        public string ComponentName { get; set; }
        public string ComponentDescription { get; set; }
        public CapeParamMode Mode { get; }

        public object value
        {
            get => _value;
            set
            {
                if (Mode == CapeParamMode.CAPE_OUTPUT)
                    _value = value;
                else
                    _value = value;
            }
        }

        public double AsDouble() => Convert.ToDouble(_value);
        public int AsInt() => Convert.ToInt32(_value);
        public bool AsBool() => Convert.ToBoolean(_value);
        public string AsString() => Convert.ToString(_value);
    }
}
