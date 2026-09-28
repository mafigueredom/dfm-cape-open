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
        object Item(object index);
        int Count { get; }
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeParameter)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeParameter
    {
        object value { get; set; }
        CapeParamMode Mode { get; }
        string ComponentName { get; set; }
        string ComponentDescription { get; set; }
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeUnitPort)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUnitPort
    {
        CapePortType portType { get; }
        CapePortDirection direction { get; }
        object connectedObject { get; }
        void Connect(object objectToConnect);
        void Disconnect();
        string ComponentName { get; set; }
        string ComponentDescription { get; set; }
    }

    [ComVisible(true)]
    [Guid(Guids.ICapeUtilities)]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface ICapeUtilities
    {
        void Initialize();
        void Terminate();
        object simulationContext { set; }
        ICapeCollection Parameters { get; }
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
