using System;
using System.Runtime.InteropServices;

namespace PhD.DfmMethanation
{
    public enum CapePortDirection
    {
        CAPE_INLET = 0,
        CAPE_OUTLET = 1,
        CAPE_INLET_OUTLET = 2
    }

    public enum CapePortType
    {
        CAPE_MATERIAL = 0,
        CAPE_ENERGY = 1,
        CAPE_INFORMATION = 2
    }

    public enum CapeParamMode
    {
        CAPE_INPUT = 0,
        CAPE_OUTPUT = 1,
        CAPE_INPUT_OUTPUT = 2
    }

    public enum CapeParamType
    {
        CAPE_REAL = 0,
        CAPE_INT = 1,
        CAPE_OPTION = 2,
        CAPE_BOOLEAN = 3,
        CAPE_ARRAY = 4
    }

    public enum CapeValidationStatus
    {
        CAPE_VALID = 0,
        CAPE_INVALID = 1,
        CAPE_NOT_VALIDATED = 2
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeIdentification)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeIdentification
    {
        string ComponentName { get; set; }
        string ComponentDescription { get; set; }
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeCollection)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeCollection
    {
        [return: MarshalAs(UnmanagedType.IDispatch)]
        object Item(object index);
        int Count { get; }
    }

    /// <summary>
    /// Vtable order matches the CAPE-OPEN 1.1 type library COFE early-binds.
    /// Specification is first: COFE calls it as soon as the block is dropped.
    /// </summary>
    [ComVisible(true)]
    [Guid(Guids.ICapeParameter)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeParameter
    {
        object Specification { [return: MarshalAs(UnmanagedType.IDispatch)] get; }
        object value { get; set; }
        CapeValidationStatus ValStatus { get; }
        CapeParamMode Mode { get; set; }
        bool Validate(ref string message);
        void Reset();
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeParameterSpec)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeParameterSpec
    {
        CapeParamType Type { get; }
        object Dimensionality { get; }
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeRealParameterSpec)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeRealParameterSpec
    {
        double DefaultValue { get; }
        double LowerBound { get; }
        double UpperBound { get; }
        bool Validate(double value, ref string message);
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeIntegerParameterSpec)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeIntegerParameterSpec
    {
        int DefaultValue { get; }
        int LowerBound { get; }
        int UpperBound { get; }
        bool Validate(int value, ref string message);
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeBooleanParameterSpec)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeBooleanParameterSpec
    {
        bool DefaultValue { get; }
        bool Validate(bool value, ref string message);
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeOptionParameterSpec)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeOptionParameterSpec
    {
        string DefaultValue { get; }
        object OptionList { get; }
        bool RestrictedToList { get; }
        bool Validate(string value, ref string message);
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeUnitPort)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUnitPort
    {
        CapePortType portType { get; }
        CapePortDirection direction { get; }
        // COFE uses this pointer as the material interface. A generic IDispatch
        // vtable is the wrong object and COFE access-violates in its own exe.
        [PreserveSig]
        int get_connectedObject(out IntPtr connectedObject);
        void Connect([MarshalAs(UnmanagedType.IDispatch)] object objectToConnect);
        void Disconnect();
    }

    /// <summary>
    /// Official CAPE-OPEN vtable and dispids. Aspen early-binds this order:
    /// parameters, simulationContext, Initialize, Terminate, Edit.
    /// </summary>
    [ComVisible(true)]
    [Guid(Guids.ICapeUtilities)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUtilities
    {
        [DispId(1)]
        ICapeCollection Parameters { [return: MarshalAs(UnmanagedType.IDispatch)] get; }
        [DispId(2)]
        object simulationContext { [param: MarshalAs(UnmanagedType.IDispatch)] set; }
        [DispId(3)]
        void Initialize();
        [DispId(4)]
        void Terminate();
        [DispId(5)]
        void Edit();
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeUnit)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUnit
    {
        ICapeCollection ports { get; }
        CapeValidationStatus ValStatus { get; }
        void Calculate();
        bool Validate(ref string message);
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeUnitReport)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUnitReport
    {
        object reports { get; }
        string selectedReport { get; set; }
        string ProduceReport();
    }
}
